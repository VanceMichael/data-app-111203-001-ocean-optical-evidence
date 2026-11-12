import copy
import csv
import json
import random
import subprocess
import sys
import tempfile
import unittest
from itertools import count
from pathlib import Path

from analysis_app.cli import save_text
from analysis_app.evidence import (ADJUSTMENTS, CENSOR_STRATEGIES, REASON_CONSTANT_BAND,
                                   REASON_NO_PAIRING, REASON_SINGLE_STATION, build_evidence)
from analysis_app.io import load_records
from analysis_app.stats import (average_ranks, benjamini_hochberg, partial_spearman,
                                regularized_incomplete_beta, spearman, t_two_sided_p)

FIELDNAMES = ['observation_id', 'station', 'sensor', 'cruise', 'day', 'reflectance',
              'plastic', 'chlorophyll', 'latitude', 'detected', 'band', 'detection_limit']
SEED = 20261112
ITERATIONS = 200


def make_row(identifier, station, sensor, cruise, day, reflectance, plastic,
             chlorophyll, band, detected=1, detection_limit=2.0, latitude=20.0):
    return {'observation_id': identifier, 'station': station, 'sensor': sensor,
            'cruise': cruise, 'day': day, 'reflectance': reflectance,
            'plastic': '' if plastic is None else plastic, 'chlorophyll': chlorophyll,
            'latitude': latitude, 'detected': detected, 'band': band,
            'detection_limit': detection_limit}


def write_csv(directory, rows):
    path = Path(directory) / 'evidence_input.csv'
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)
    return path


def monotone_rows(station_count=8, sensors=('sensor-0',), bands=('blue',),
                  cruises=('C0', 'C1'), conflict_sensor=None, noise_seed=31):
    """构造站点层面反射率与浓度单调同向（或某传感器反向）的合成观测。"""
    # 叶绿素站点排序与反射率刻意不同，保证取秩后仍可做偏相关。
    shuffled = list(range(station_count))
    random.Random(noise_seed + 7).shuffle(shuffled)
    station_noise = {index: 0.35 * shuffled[index]
                     for index in range(station_count)}
    # 浓度在站点均值层面带固定交替扰动，刻意打乱相邻站点秩次，
    # 避免完美单调导致 rho=±1 时 t 近似无定义。
    plastic_jitter = {index: 0.32 if index % 2 == 0 else -0.32
                      for index in range(station_count)}
    rows = []
    seq = count()
    for station_index in range(station_count):
        station = f'S{station_index:02d}'
        for cruise_index, cruise in enumerate(cruises):
            day = cruise_index * 2
            for sensor in sensors:
                for band in bands:
                    base = station_index + 1
                    direction = -1 if sensor == conflict_sensor else 1
                    reflectance = 0.01 * base + (0.0001 * hash_band(band))
                    plastic = (0.5 * base * direction + 5.0
                               + plastic_jitter[station_index])
                    chlorophyll = 0.3 * base + 0.1 * cruise_index + station_noise[station_index]
                    rows.append(make_row(f'obs-{next(seq):04d}', station, sensor,
                                         cruise, day, reflectance, plastic,
                                         chlorophyll, band))
    return rows


def hash_band(band):
    return 1 if band == 'green' else 0


class StatsTest(unittest.TestCase):
    def test_average_ranks_with_ties(self):
        self.assertEqual(average_ranks([3, 1, 1, 5, 3]), [3.5, 1.5, 1.5, 5.0, 3.5])

    def test_spearman_monotone_and_constant(self):
        self.assertEqual(spearman(list(range(10)), list(range(10))), 1.0)
        self.assertEqual(spearman(list(range(10)), list(range(10))[::-1]), -1.0)
        self.assertIsNone(spearman([1, 1, 1], [1, 2, 3]))

    def test_partial_spearman_removes_confounder(self):
        rng = random.Random(7)
        z = [rng.gauss(0, 1) for _ in range(300)]
        x = [value + rng.gauss(0, 0.01) for value in z]
        y = [value + rng.gauss(0, 0.01) for value in z]
        self.assertGreater(spearman(x, y), 0.99)
        self.assertLess(abs(partial_spearman(x, y, [z])), 0.3)

    def test_t_p_value_known_points(self):
        self.assertAlmostEqual(t_two_sided_p(0.0, 20, 0), 1.0)
        r = 2.22813885 / (10 + 2.22813885 ** 2) ** 0.5
        self.assertAlmostEqual(t_two_sided_p(r, 12, 0), 0.05, places=4)
        self.assertIsNone(t_two_sided_p(0.5, 3, 2))

    def test_beta_function_exact_value(self):
        self.assertAlmostEqual(regularized_incomplete_beta(2, 3, 0.5), 0.6875)

    def test_bh_adjustment_keeps_none_and_caps_at_one(self):
        adjusted = benjamini_hochberg([0.01, 0.02, 0.03, 0.04, 0.05])
        for value in adjusted:
            self.assertAlmostEqual(value, 0.05)
        with_none = benjamini_hochberg([0.04, None, 0.01, 0.03, 0.02])
        self.assertIsNone(with_none[1])
        for index, value in enumerate(with_none):
            if index != 1:
                self.assertAlmostEqual(value, 0.04)


class EvidenceTest(unittest.TestCase):
    def run_evidence(self, rows, cutoff=99, seed=SEED, iterations=ITERATIONS):
        with tempfile.TemporaryDirectory() as directory:
            path = write_csv(directory, rows)
            records = load_records(path)
        return build_evidence(records, cutoff, seed, iterations)

    def cell(self, result, sensor, band, strategy='as_reported', adjustment='none'):
        return next(cell for cell in result['cells']
                    if cell['sensor'] == sensor and cell['band'] == band
                    and cell['censor_strategy'] == strategy
                    and cell['adjustment'] == adjustment)

    def test_full_pipeline_structure_on_fixture(self):
        records = load_records('data/observations.csv')
        result = build_evidence(records, 7, SEED, ITERATIONS)
        self.assertTrue(result['estimable'])
        self.assertEqual(len(result['cells']), 3 * 2 * len(CENSOR_STRATEGIES) * len(ADJUSTMENTS))
        self.assertEqual(len(result['bh_families']), len(CENSOR_STRATEGIES) * len(ADJUSTMENTS))
        self.assertEqual(result['counts']['window_rows'], 192)
        self.assertEqual(result['counts']['excluded_after_cutoff_rows'], 96)
        self.assertEqual(result['counts']['paired_rows'] + result['counts']['unpaired_rows'],
                         result['counts']['window_rows'])
        cell = self.cell(result, 'sensor-0', 'blue')
        self.assertEqual(cell['n_stations'], 24)
        self.assertIsNotNone(cell['point_estimate'])
        self.assertIsNotNone(cell['adjusted_p_value'])
        self.assertEqual(cell['bootstrap']['resampling_unit'], 'station')
        self.assertEqual(cell['bootstrap']['iterations'], ITERATIONS)
        self.assertLess(cell['bootstrap']['ci_low'], cell['point_estimate'])
        self.assertGreater(cell['bootstrap']['ci_high'], cell['point_estimate'])

    def test_row_order_does_not_change_conclusion(self):
        records = load_records('data/observations.csv')
        baseline = build_evidence(records, 7, SEED, ITERATIONS)
        shuffled = copy.deepcopy(records)
        random.Random(99).shuffle(shuffled)
        shuffled_again = copy.deepcopy(records)
        random.Random(12345).shuffle(shuffled_again)
        self.assertEqual(build_evidence(shuffled, 7, SEED, ITERATIONS), baseline)
        self.assertEqual(build_evidence(shuffled_again, 7, SEED, ITERATIONS), baseline)

    def test_seed_changes_bootstrap_but_not_point_estimate(self):
        records = load_records('data/observations.csv')
        first = build_evidence(records, 7, 1, ITERATIONS)
        second = build_evidence(records, 7, 2, ITERATIONS)
        c1 = self.cell(first, 'sensor-0', 'blue')
        c2 = self.cell(second, 'sensor-0', 'blue')
        self.assertEqual(c1['point_estimate'], c2['point_estimate'])
        self.assertNotEqual((c1['bootstrap']['ci_low'], c1['bootstrap']['ci_high']),
                            (c2['bootstrap']['ci_low'], c2['bootstrap']['ci_high']))

    def test_unpaired_rows_never_enter_concentration_association(self):
        rows = monotone_rows()
        rows[0]['plastic'] = ''
        result = self.run_evidence(rows)
        cell = self.cell(result, 'sensor-0', 'blue')
        self.assertEqual(cell['n_rows'], result['counts']['paired_rows'])
        self.assertEqual(result['counts']['unpaired_rows'], 1)
        self.assertEqual(sum(cell['rows_per_station'].values()), cell['n_rows'])

    def test_censored_is_not_zero_and_strategies_differ(self):
        # 四个站点各一行配对；两个低于检出限站点的检出限不同，
        # 三种替代口径给出不同秩次，证明 detected=0 未被当作零。
        rows = [
            make_row('o1', 'S00', 'sensor-0', 'C0', 0, 0.01, 0.4, 1.0,
                     'blue', detected=0, detection_limit=2.0),
            make_row('o2', 'S01', 'sensor-0', 'C0', 0, 0.02, 0.8, 1.1,
                     'blue', detected=0, detection_limit=10.0),
            make_row('o3', 'S02', 'sensor-0', 'C0', 0, 0.03, 3.0, 1.2, 'blue'),
            make_row('o4', 'S03', 'sensor-0', 'C0', 0, 0.04, 9.0, 1.3, 'blue'),
        ]
        result = self.run_evidence(rows)
        as_reported = self.cell(result, 'sensor-0', 'blue', 'as_reported')
        half = self.cell(result, 'sensor-0', 'blue', 'half_lod')
        lod = self.cell(result, 'sensor-0', 'blue', 'lod')
        self.assertEqual(result['counts']['censored_paired_rows'], 2)
        for cell in (as_reported, half, lod):
            self.assertTrue(cell['estimable'])
        estimates = {as_reported['point_estimate'], half['point_estimate'],
                     lod['point_estimate']}
        self.assertEqual(len(estimates), 3)
        # 没有任何口径把删失值处理成零：half_lod 最小替代值为 1.0，lod 为 2.0；
        # 三个口径的秩次结构各不相同（0.4 < 1.0 < 3.0 < 9.0 等），相关系数可手算复核。
        self.assertAlmostEqual(half['point_estimate'], 0.8)
        self.assertNotEqual(lod['point_estimate'], half['point_estimate'])
        self.assertNotEqual(as_reported['point_estimate'], half['point_estimate'])

    def test_cutoff_filters_rows(self):
        rows = monotone_rows()
        result = self.run_evidence(rows, cutoff=1)
        self.assertLess(result['counts']['window_rows'], len(rows))
        self.assertEqual(result['counts']['excluded_after_cutoff_rows'],
                         len(rows) - result['counts']['window_rows'])

    def test_station_repeats_are_aggregated_not_independent(self):
        base = monotone_rows(station_count=6)
        duplicated = base + [dict(row, observation_id=f'dup-{index}')
                             for index, row in enumerate(base[:10])]
        result = self.run_evidence(duplicated)
        cell = self.cell(result, 'sensor-0', 'blue')
        self.assertEqual(cell['n_stations'], 6)
        self.assertGreaterEqual(min(cell['rows_per_station'].values()), 2)
        single = self.run_evidence(base)
        single_cell = self.cell(single, 'sensor-0', 'blue')
        self.assertEqual(cell['point_estimate'], single_cell['point_estimate'])
        self.assertEqual(cell['bootstrap'], single_cell['bootstrap'])

    def test_constant_band_is_not_estimable_with_reason(self):
        rows = monotone_rows(station_count=6, bands=('blue', 'green'))
        for row in rows:
            if row['band'] == 'blue':
                row['reflectance'] = 0.05
        result = self.run_evidence(rows)
        cell = self.cell(result, 'sensor-0', 'blue')
        self.assertFalse(cell['estimable'])
        self.assertEqual(cell['not_estimable_reason'], REASON_CONSTANT_BAND)
        green = self.cell(result, 'sensor-0', 'green')
        self.assertTrue(green['estimable'])

    def test_single_effective_station_reports_reason(self):
        rows = monotone_rows(station_count=3)
        for row in rows:
            row['station'] = 'S00'
            row['observation_id'] = f"only-{row['observation_id']}"
        result = self.run_evidence(rows)
        cell = self.cell(result, 'sensor-0', 'blue')
        self.assertFalse(cell['estimable'])
        self.assertEqual(cell['not_estimable_reason'], REASON_SINGLE_STATION)

    def test_no_paired_concentration_reports_reason(self):
        rows = monotone_rows(station_count=4)
        for row in rows:
            row['plastic'] = ''
        result = self.run_evidence(rows)
        self.assertFalse(result['estimable'])
        self.assertEqual(result['not_estimable_reason'], REASON_NO_PAIRING)
        self.assertTrue(all(not cell['estimable'] for cell in result['cells']))
        self.assertIn(REASON_NO_PAIRING, result['narrative_zh'])

    def test_empty_input_reports_reason(self):
        result = build_evidence([], 9, SEED, ITERATIONS)
        self.assertFalse(result['estimable'])
        self.assertEqual(result['counts'], {})
        self.assertTrue(result['narrative_zh'])

    def test_cutoff_before_all_rows(self):
        result = self.run_evidence(monotone_rows(), cutoff=-1)
        self.assertFalse(result['estimable'])
        self.assertEqual(result['counts']['window_rows'], 0)

    def test_cross_sensor_direction_conflict(self):
        rows = monotone_rows(station_count=8, sensors=('sensor-0', 'sensor-1', 'sensor-2'),
                             conflict_sensor='sensor-2')
        result = self.run_evidence(rows)
        view = next(view for view in result['cross_sensor']
                    if view['censor_strategy'] == 'as_reported'
                    and view['adjustment'] == 'none')
        self.assertTrue(view['direction_conflict'])
        self.assertEqual(view['positive_sensors'], ['sensor-0', 'sensor-1'])
        self.assertEqual(view['negative_sensors'], ['sensor-2'])
        self.assertIn('冲突', result['narrative_zh'])

    def test_bh_family_groups_strategy_and_adjustment(self):
        result = self.run_evidence(monotone_rows(sensors=('sensor-0', 'sensor-1'),
                                                 bands=('blue', 'green')))
        family = next(family for family in result['bh_families']
                      if family['censor_strategy'] == 'as_reported'
                      and family['adjustment'] == 'none')
        self.assertEqual(family['family_size_total'], 4)
        self.assertEqual(family['family_size_tested'], 4)
        self.assertEqual(len(family['members']), 4)
        for member in family['members']:
            self.assertGreaterEqual(member['adjusted_p_value'], member['p_value'])
        for cell in result['cells']:
            if cell['estimable'] and cell['p_value'] is not None:
                self.assertIsNotNone(cell['adjusted_p_value'])

    def test_cruise_control_with_single_cruise_is_not_estimable(self):
        rows = monotone_rows(station_count=6, cruises=('C0',))
        result = self.run_evidence(rows)
        cell = self.cell(result, 'sensor-0', 'blue', adjustment='chlorophyll_cruise')
        self.assertFalse(cell['estimable'])
        from analysis_app.evidence import REASON_CONSTANT_CRUISE
        self.assertEqual(cell['not_estimable_reason'], REASON_CONSTANT_CRUISE)
        controlled = self.cell(result, 'sensor-0', 'blue', adjustment='chlorophyll')
        self.assertTrue(controlled['estimable'])

    def test_invariant_cruise_dummy_is_dropped_but_other_cruise_kept(self):
        # C0 行覆盖全部站点（其站点层面指示列恒为正），C1 只覆盖部分站点；
        # 偏相关应自动剔除常数列，控制列数少于 1 + 航次数。
        rows = monotone_rows(station_count=8, cruises=('C0',))
        extra = monotone_rows(station_count=4, cruises=('C1',), noise_seed=73)
        for index, row in enumerate(extra):
            row['observation_id'] = f'extra-{index:04d}'
            row['station'] = f'S{index:02d}'
        result = self.run_evidence(rows + extra)
        cell = self.cell(result, 'sensor-0', 'blue', adjustment='chlorophyll_cruise')
        self.assertTrue(cell['estimable'])
        self.assertEqual(cell['n_cruises'], 2)
        self.assertLess(cell['control_count'], cell['n_cruises'] + 1)
        self.assertEqual(cell['control_count'], 2)

    def test_chlorophyll_control_runs_and_audit_trail_conserved(self):
        records = load_records('data/observations.csv')
        result = build_evidence(records, 11, SEED, ITERATIONS)
        controlled = self.cell(result, 'sensor-1', 'green', adjustment='chlorophyll')
        self.assertTrue(controlled['estimable'])
        self.assertEqual(controlled['control_count'], 1)
        raw = self.cell(result, 'sensor-1', 'green', adjustment='none')
        self.assertEqual(controlled['n_rows'], raw['n_rows'])
        self.assertEqual(controlled['n_stations'], raw['n_stations'])
        total_audit = sum(
            cell['n_rows'] for cell in result['cells']
            if cell['censor_strategy'] == 'as_reported' and cell['adjustment'] == 'none')
        self.assertEqual(total_audit, result['counts']['paired_rows'])

    def test_invalid_arguments(self):
        records = load_records('data/observations.csv')
        with self.assertRaisesRegex(ValueError, '截点'):
            build_evidence(records, 1.5, SEED, ITERATIONS)
        with self.assertRaisesRegex(ValueError, '种子'):
            build_evidence(records, 7, 'seed')
        with self.assertRaisesRegex(ValueError, '自助次数'):
            build_evidence(records, 7, SEED, 0)

    def test_narrative_contains_required_sections_and_no_threshold_claim(self):
        records = load_records('data/observations.csv')
        narrative = build_evidence(records, 7, SEED, ITERATIONS)['narrative_zh']
        for section in ('数据口径与样本去向', '关联方向与跨传感器冲突', '稳定性对照', '限制'):
            self.assertIn(section, narrative)
        self.assertIn('不构成污染因果结论', narrative)
        self.assertNotIn('阈值', narrative.replace('风险阈值', ''))

    def test_json_is_finite_and_serializable(self):
        result = build_evidence(load_records('data/observations.csv'), 7, SEED, ITERATIONS)
        payload = json.dumps(result, ensure_ascii=False, sort_keys=True, allow_nan=False)
        decoded = json.loads(payload)
        self.assertEqual(decoded['parameters']['seed'], SEED)

    def test_cli_evidence_end_to_end(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / 'nested' / 'evidence.json'
            narrative = Path(directory) / 'narrative.txt'
            run = subprocess.run(
                [sys.executable, '-m', 'analysis_app', 'evidence',
                 '--input', 'data/observations.csv', '--cutoff', '7',
                 '--seed', str(SEED), '--iterations', '80',
                 '--output', str(output), '--narrative', str(narrative)],
                text=True, capture_output=True)
            self.assertEqual(run.returncode, 0, run.stderr)
            on_disk = json.loads(output.read_text(encoding='utf-8'))
            on_stdout = json.loads(run.stdout)
            self.assertEqual(on_disk, on_stdout)
            self.assertEqual(narrative.read_text(encoding='utf-8').strip(),
                             on_disk['narrative_zh'])
            self.assertEqual([item.name for item in output.parent.iterdir()],
                             ['evidence.json'])

    def test_cli_evidence_failure_leaves_no_output(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / 'evidence.json'
            run = subprocess.run(
                [sys.executable, '-m', 'analysis_app', 'evidence',
                 '--input', str(Path(directory) / 'missing.csv'), '--cutoff', '7',
                 '--seed', str(SEED), '--output', str(output)],
                text=True, capture_output=True)
            self.assertEqual(run.returncode, 2)
            self.assertFalse(output.exists())

    def test_atomic_text_output(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'note.txt'
            save_text(path, '第一段')
            save_text(path, '第二段')
            self.assertEqual(path.read_text(encoding='utf-8').strip(), '第二段')


if __name__ == '__main__':
    unittest.main()
