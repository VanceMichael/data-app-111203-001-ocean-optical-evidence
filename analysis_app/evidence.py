"""水色反射率与浓度的关联证据流程。

分析单位分三层：
- 观测行：只有 plastic 非空的配对行进入浓度关联；detected=0 表示低于检出限，
  其浓度按删失策略替代，绝不当作零；
- 航次（cruise）：作为可选控制变量，在站点层面以航次指示列参与偏相关；
- 站点（station）：关联与区间估计的独立单元。同一站点的重复行先做站点内
  聚合，区间使用以站点为重采样单元的聚类自助法，重复记录不会被当成独立样本。

本模块只产出关联证据与稳定性对照，不输出污染因果结论或风险阈值。
"""

import hashlib
import random

from .io import fingerprint
from .schema import PRODUCT
from .stats import (benjamini_hochberg, partial_spearman, percentile, spearman,
                    t_two_sided_p)

CENSOR_STRATEGIES = ('as_reported', 'half_lod', 'lod')
CENSOR_LABELS = {
    'as_reported': '按原始记录保留低于检出限的测量值',
    'half_lod': '低于检出限者以检出限的二分之一替代',
    'lod': '低于检出限者以检出限替代',
}
ADJUSTMENTS = ('none', 'chlorophyll', 'chlorophyll_cruise')
ADJUSTMENT_LABELS = {
    'none': '原始斯皮尔曼等级相关',
    'chlorophyll': '控制叶绿素后的偏斯皮尔曼相关',
    'chlorophyll_cruise': '同时控制叶绿素与航次后的偏斯皮尔曼相关',
}
BOOTSTRAP_ITERATIONS_DEFAULT = 1000

REASON_NO_PAIRING = '截点后没有配对浓度记录'
REASON_SINGLE_STATION = '仅有一个有效站点，无法估计站点间关联'
REASON_CONSTANT_BAND = '该波段反射率在站点间为常数，关联不可估计'
REASON_CONSTANT_CONCENTRATION = '替代后浓度在站点间为常数，关联不可估计'
REASON_CONSTANT_CHLOROPHYLL = '叶绿素在站点间为常数，控制叶绿素的偏相关不可估计'
REASON_CONSTANT_CRUISE = '航次指示在站点层面无变化，控制航次的偏相关不可估计'
REASON_DEGREES = '控制变量后自由度不足，关联不可估计'
REASON_EMPTY_INPUT = '输入数据为空'
REASON_NO_WINDOW_ROWS = '截点后没有任何观测行'
REASON_BOOTSTRAP_EMPTY = '自助重采样均未产生有效估计，区间不可用'


def _censor_value(row, strategy):
    if row['detected'] == 0:
        if strategy == 'half_lod':
            return row['detection_limit'] / 2.0
        if strategy == 'lod':
            return row['detection_limit']
    return row['plastic']


def _sub_seed(seed, *parts):
    payload = '|'.join([str(seed), *parts]).encode('utf-8')
    digest = hashlib.sha256(payload).hexdigest()
    return int(digest[:16], 16)


def _control_columns(stations_order, station_payload, cruise_order):
    chlorophyll = [station_payload[station]['chlorophyll_mean'] for station in stations_order]
    # 航次用 n-1 个指示列（省略首个航次作为参考），避免份额列和为 1 的虚拟变量陷阱。
    cruise_columns = []
    for cruise in cruise_order[1:]:
        cruise_columns.append(
            [station_payload[station]['cruise_share'].get(cruise, 0.0)
             for station in stations_order])
    return chlorophyll, cruise_columns


def _estimate_adjustment(x, y, chlorophyll, cruise_columns, adjustment):
    """返回 (相关系数, 控制变量个数, 不可估计原因或 None)。"""
    if adjustment == 'none':
        rho = spearman(x, y)
        return rho, 0, None if rho is not None else REASON_CONSTANT_BAND
    if len(set(chlorophyll)) < 2:
        return None, None, REASON_CONSTANT_CHLOROPHYLL
    if adjustment == 'chlorophyll':
        rho = partial_spearman(x, y, [chlorophyll])
        control_count = 1
    else:
        active_cruise = [column for column in cruise_columns if len(set(column)) >= 2]
        if not active_cruise:
            return None, None, REASON_CONSTANT_CRUISE
        rho = partial_spearman(x, y, [chlorophyll] + active_cruise)
        control_count = 1 + len(active_cruise)
    if rho is None:
        return None, control_count, REASON_DEGREES
    return rho, control_count, None


def _bootstrap_interval(seed, sensor, band, strategy, adjustment, stations_order,
                        station_payload, cruise_order, iterations):
    rng = random.Random(_sub_seed(seed, sensor, band, strategy, adjustment))
    size = len(stations_order)
    replicates = []
    failed = 0
    for _ in range(iterations):
        picked = rng.choices(stations_order, k=size)
        x = [station_payload[station]['reflectance_mean'] for station in picked]
        y = [station_payload[station]['y_mean'][strategy] for station in picked]
        chlorophyll = [station_payload[station]['chlorophyll_mean'] for station in picked]
        cruise_columns = [[station_payload[station]['cruise_share'].get(cruise, 0.0)
                           for station in picked] for cruise in cruise_order]
        rho, _, reason = _estimate_adjustment(x, y, chlorophyll, cruise_columns, adjustment)
        if rho is None:
            failed += 1
        else:
            replicates.append(rho)
    if not replicates:
        return {'iterations': iterations, 'successful': 0, 'failed': failed,
                'ci_low': None, 'ci_high': None, 'reason': REASON_BOOTSTRAP_EMPTY,
                'method': 'percentile_5_95', 'resampling_unit': 'station'}
    ordered = sorted(replicates)
    return {'iterations': iterations, 'successful': len(ordered), 'failed': failed,
            'ci_low': percentile(ordered, 0.05), 'ci_high': percentile(ordered, 0.95),
            'reason': None, 'method': 'percentile_5_95', 'resampling_unit': 'station'}


def _build_cell(sensor, band, strategy, adjustment, window_rows, seed, iterations):
    cell_rows = [row for row in window_rows
                 if row['sensor'] == sensor and row['band'] == band]
    paired_rows = [row for row in cell_rows if row['plastic'] is not None]
    audit_rows_per_station = {}
    grouped = {}
    for row in paired_rows:
        station = row['station']
        audit_rows_per_station.setdefault(station, 0)
        audit_rows_per_station[station] += 1
        grouped.setdefault(station, []).append(row)
    stations_order = sorted(grouped)
    base = {'sensor': sensor, 'band': band, 'censor_strategy': strategy,
            'adjustment': adjustment, 'n_rows': len(paired_rows),
            'n_stations': len(stations_order),
            'stations': stations_order,
            'rows_per_station': {station: audit_rows_per_station[station]
                                 for station in stations_order},
            'point_estimate': None, 'p_value': None, 'adjusted_p_value': None,
            'bootstrap': None, 'estimable': False, 'not_estimable_reason': None}
    if not paired_rows:
        base['not_estimable_reason'] = REASON_NO_PAIRING
        return base, None
    if len(stations_order) < 2:
        base['not_estimable_reason'] = REASON_SINGLE_STATION
        return base, None

    cruise_order = sorted({row['cruise'] for row in paired_rows})
    station_payload = {}
    for station in stations_order:
        members = grouped[station]
        count = len(members)
        share = {}
        for cruise in cruise_order:
            share[cruise] = sum(1 for row in members if row['cruise'] == cruise) / count
        station_payload[station] = {
            'reflectance_mean': sum(row['reflectance'] for row in members) / count,
            'chlorophyll_mean': sum(row['chlorophyll'] for row in members) / count,
            'y_mean': {strategy: sum(_censor_value(row, strategy) for row in members) / count},
            'cruise_share': share,
        }
    x = [station_payload[station]['reflectance_mean'] for station in stations_order]
    y = [station_payload[station]['y_mean'][strategy] for station in stations_order]
    if len(set(x)) < 2:
        base['not_estimable_reason'] = REASON_CONSTANT_BAND
        return base, None
    if len(set(y)) < 2:
        base['not_estimable_reason'] = REASON_CONSTANT_CONCENTRATION
        return base, None
    chlorophyll, cruise_columns = _control_columns(stations_order, station_payload, cruise_order)
    rho, control_count, reason = _estimate_adjustment(x, y, chlorophyll,
                                                      cruise_columns, adjustment)
    if rho is None:
        base['not_estimable_reason'] = reason
        return base, None
    base['point_estimate'] = rho
    base['p_value'] = t_two_sided_p(rho, len(stations_order), control_count)
    base['bootstrap'] = _bootstrap_interval(
        seed, sensor, band, strategy, adjustment, stations_order,
        station_payload, cruise_order, iterations)
    base['n_cruises'] = len(cruise_order)
    base['control_count'] = control_count
    base['estimable'] = True
    return base, base['p_value']


def _cross_sensor_view(cells):
    views = []
    keys = {(cell['band'], cell['censor_strategy'], cell['adjustment']) for cell in cells}
    for band, strategy, adjustment in sorted(keys):
        members = [cell for cell in cells
                   if cell['band'] == band and cell['censor_strategy'] == strategy
                   and cell['adjustment'] == adjustment]
        directions = {}
        positive, negative, zero, missing = [], [], [], []
        for cell in sorted(members, key=lambda item: item['sensor']):
            if not cell['estimable']:
                directions[cell['sensor']] = 'not_estimable'
                missing.append(cell['sensor'])
                continue
            if cell['point_estimate'] > 0:
                sign = 'positive'
                positive.append(cell['sensor'])
            elif cell['point_estimate'] < 0:
                sign = 'negative'
                negative.append(cell['sensor'])
            else:
                sign = 'zero'
                zero.append(cell['sensor'])
            directions[cell['sensor']] = sign
        views.append({'band': band, 'censor_strategy': strategy,
                      'adjustment': adjustment, 'directions': directions,
                      'positive_sensors': positive, 'negative_sensors': negative,
                      'zero_sensors': zero, 'not_estimable_sensors': missing,
                      'direction_conflict': bool(positive and negative)})
    return views


def _family_adjustment(cells):
    families = []
    keys = {(cell['censor_strategy'], cell['adjustment']) for cell in cells}
    for strategy, adjustment in sorted(keys):
        members = [cell for cell in cells
                   if cell['censor_strategy'] == strategy
                   and cell['adjustment'] == adjustment]
        tested = [cell for cell in members if cell['estimable']
                  and cell['p_value'] is not None]
        adjusted = benjamini_hochberg([cell['p_value'] for cell in tested])
        member_report = []
        for cell, q_value in zip(tested, adjusted):
            cell['adjusted_p_value'] = q_value
            member_report.append({'sensor': cell['sensor'], 'band': cell['band'],
                                  'p_value': cell['p_value'],
                                  'adjusted_p_value': q_value})
        families.append({'censor_strategy': strategy, 'adjustment': adjustment,
                         'family_size_total': len(members),
                         'family_size_tested': len(tested),
                         'method': 'benjamini_hochberg',
                         'members': sorted(member_report,
                                           key=lambda item: (item['sensor'], item['band']))})
    return families


def _global_counts(records, window_rows):
    window_paired = [row for row in window_rows if row['plastic'] is not None]
    stations = [row['station'] for row in window_rows]
    rows_by_cell = {}
    for row in window_rows:
        key = f"{row['sensor']}|{row['band']}"
        rows_by_cell[key] = rows_by_cell.get(key, 0) + 1
    return {
        'input_rows': len(records),
        'window_rows': len(window_rows),
        'excluded_after_cutoff_rows': len(records) - len(window_rows),
        'paired_rows': len(window_paired),
        'unpaired_rows': len(window_rows) - len(window_paired),
        'detected_paired_rows': sum(1 for row in window_paired if row['detected'] == 1),
        'censored_paired_rows': sum(1 for row in window_paired if row['detected'] == 0),
        'stations_total': len(set(stations)),
        'stations_with_pairing': len({row['station'] for row in window_paired}),
        'repeated_station_rows': len(stations) - len(set(stations)),
        'cruises_total': len({row['cruise'] for row in window_rows}),
        'rows_by_sensor_band': dict(sorted(rows_by_cell.items())),
    }


def _narrative(result):
    lines = []
    params = result['parameters']
    counts = result['counts']
    lines.append('一、数据口径与样本去向')
    if not counts:
        lines.append('输入数据为空，未进入任何统计步骤，无样本去向可报告。')
        lines.append('二、关联结论')
        lines.append('本次数据不可估计：' + result['not_estimable_reason'] + '。')
        lines.append('三、稳定性与限制')
        lines.append('未产生关联估计，因而无稳定性对照；本证据包不构成污染因果结论，也不提供风险阈值推断。')
        return '\n'.join(lines)
    lines.append(
        f"截点 day<={params['cutoff_day']} 共保留 {counts['window_rows']} 行观测"
        f"（截点后剔除 {counts['excluded_after_cutoff_rows']} 行），其中配对浓度行 "
        f"{counts['paired_rows']} 行、未配对 {counts['unpaired_rows']} 行；"
        f"配对行中低于检出限 {counts['censored_paired_rows']} 行，按删失策略替代而非记零。"
        f"有效站点 {counts['stations_with_pairing']} 个、航次 {counts['cruises_total']} 个；"
        f"同站重复行 {counts['repeated_station_rows']} 行已先在站点内聚合，"
        "区间以站点为唯一重采样单元，重复记录不视作独立样本。")
    if not result['estimable']:
        lines.append('二、关联结论')
        lines.append('本次数据不可估计：' + result['not_estimable_reason'] + '。')
        lines.append('三、稳定性与限制')
        lines.append('未产生关联估计，因而无稳定性对照；本证据包不构成污染因果结论，也不提供风险阈值推断。')
        return '\n'.join(lines)

    cells = result['cells']
    estimable = [cell for cell in cells if cell['estimable']]
    lines.append('二、关联方向与跨传感器冲突')
    raw = [cell for cell in estimable if cell['adjustment'] == 'none'
           and cell['censor_strategy'] == 'as_reported']
    for cell in sorted(raw, key=lambda item: (item['sensor'], item['band'])):
        boot = cell['bootstrap']
        if boot is not None and boot['ci_low'] is not None:
            interval = f"90% 自助区间 [{boot['ci_low']:.3f}, {boot['ci_high']:.3f}]，"
        else:
            interval = '站点重采样区间不可用（原因见 JSON bootstrap.reason），'
        q_value = cell['adjusted_p_value']
        q_text = f"{q_value:.3f}" if q_value is not None else '不可用'
        lines.append(
            f"{cell['sensor']} / {cell['band']}：原始斯皮尔曼相关 "
            f"{cell['point_estimate']:.3f}，{interval}"
            f"BH 校正后 p={q_text}（基于 {cell['n_stations']} 个站点）。")
    conflicts = [view for view in result['cross_sensor'] if view['direction_conflict']]
    if conflicts:
        first = conflicts[0]
        lines.append(
            f"存在跨传感器方向冲突（如 {first['band']} 波段、{first['censor_strategy']} 策略、"
            f"{first['adjustment']} 口径：正向 {','.join(first['positive_sensors']) or '无'}；"
            f"负向 {','.join(first['negative_sensors']) or '无'}），采样点位决策不应只采纳单一传感器方向。")
    else:
        lines.append('各可估计组合的跨传感器方向一致，未出现正负冲突；部分组合不可估计时以 JSON 中原因字段为准。')

    lines.append('三、稳定性对照')
    lines.append('删失敏感性覆盖三种口径：' + '；'.join(
        f"{key}（{label}）" for key, label in CENSOR_LABELS.items()) + '。')
    spread = []
    for sensor in sorted({cell['sensor'] for cell in estimable}):
        for band in sorted({cell['band'] for cell in estimable}):
            values = [cell['point_estimate'] for cell in estimable
                      if cell['sensor'] == sensor and cell['band'] == band]
            if values:
                spread.append(f"{sensor}/{band} 全口径范围 [{min(values):.3f}, {max(values):.3f}]")
    lines.append('估计稳定性：' + '；'.join(spread) + '。')
    controlled = [cell for cell in estimable if cell['adjustment'] == 'chlorophyll']
    if controlled:
        lines.append(
            f"控制叶绿素后共 {len(controlled)} 个组合仍可估计，方向是否保持请对照 JSON 中 "
            "adjustment=chlorophyll 的 point_estimate；同时控制航次的组合若标记不可估计，"
            "原因是站点层面航次指示无变化或自由度不足。")

    lines.append('四、限制')
    lines.append(
        '本证据包仅描述反射率与浓度的统计关联及其对删失处理、叶绿素控制、站点聚类的敏感性，'
        '属于观测性关联证据，不构成污染因果结论，也不提供任何风险阈值或采样点位的直接建议；'
        f"全部结果由固定种子 {params['seed']} 与 {params['bootstrap_iterations']} 次站点聚类自助复算，"
        '改变输入行顺序不改变同种子结论。')
    return '\n'.join(lines)


def build_evidence(records, cutoff_day, seed, bootstrap_iterations=BOOTSTRAP_ITERATIONS_DEFAULT):
    if not isinstance(cutoff_day, int) or isinstance(cutoff_day, bool):
        raise ValueError('截点必须是整数 day')
    if not isinstance(seed, int) or isinstance(seed, bool):
        raise ValueError('种子必须是整数')
    if not isinstance(bootstrap_iterations, int) or bootstrap_iterations <= 0:
        raise ValueError('自助次数必须为正整数')
    ordered = sorted(records, key=lambda row: row['observation_id'])
    result = {
        'product': PRODUCT,
        'schema_version': 1,
        'analysis': 'evidence',
        'input_fingerprint': fingerprint(ordered),
        'parameters': {
            'cutoff_day': cutoff_day,
            'seed': seed,
            'bootstrap_iterations': bootstrap_iterations,
            'bootstrap_unit': 'station',
            'bootstrap_method': 'cluster_percentile_5_95',
            'multiple_testing': 'benjamini_hochberg_within_strategy_and_adjustment',
        },
        'censor_strategies': [{'key': key, 'description': label}
                              for key, label in CENSOR_LABELS.items()],
        'adjustments': [{'key': key, 'description': label}
                        for key, label in ADJUSTMENT_LABELS.items()],
        'counts': {}, 'estimable': False, 'not_estimable_reason': None,
        'cells': [], 'bh_families': [], 'cross_sensor': [], 'narrative_zh': ''}

    if not ordered:
        result['not_estimable_reason'] = REASON_EMPTY_INPUT
        result['narrative_zh'] = _narrative(result)
        return result
    window_rows = [row for row in ordered if row['day'] <= cutoff_day]
    result['counts'] = _global_counts(ordered, window_rows)
    if not window_rows:
        result['not_estimable_reason'] = REASON_NO_WINDOW_ROWS
        result['narrative_zh'] = _narrative(result)
        return result
    if result['counts']['paired_rows'] == 0:
        result['not_estimable_reason'] = REASON_NO_PAIRING
        result['narrative_zh'] = _narrative(result)
        return result

    sensors = sorted({row['sensor'] for row in window_rows})
    bands = sorted({row['band'] for row in window_rows})
    cells = []
    for sensor in sensors:
        for band in bands:
            for strategy in CENSOR_STRATEGIES:
                for adjustment in ADJUSTMENTS:
                    cell, _ = _build_cell(sensor, band, strategy, adjustment,
                                          window_rows, seed, bootstrap_iterations)
                    cells.append(cell)
    result['cells'] = cells
    result['bh_families'] = _family_adjustment(cells)
    result['cross_sensor'] = _cross_sensor_view(cells)
    result['estimable'] = any(cell['estimable'] for cell in cells)
    if not result['estimable']:
        reasons = sorted({cell['not_estimable_reason'] for cell in cells
                          if cell['not_estimable_reason']})
        result['not_estimable_reason'] = reasons[0]
    result['narrative_zh'] = _narrative(result)
    return result
