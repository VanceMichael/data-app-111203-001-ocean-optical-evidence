"""无第三方依赖的关联统计基础函数。

所有函数只接受有限实数；输入长度不足或变量无变化时返回 None，由上层
转换为明确的不可估计原因，而不是抛出异常或静默给出伪数值。
"""

import math


def average_ranks(values):
    """返回平均秩（结取平均秩次），最小取值为 1。"""
    order = sorted(range(len(values)), key=lambda index: values[index])
    ranks = [0.0] * len(values)
    start = 0
    while start < len(values):
        end = start + 1
        while end < len(values) and values[order[end]] == values[order[start]]:
            end += 1
        tied_rank = (start + 1 + end) / 2.0
        for index in range(start, end):
            ranks[order[index]] = tied_rank
        start = end
    return ranks


def _mean(values):
    return sum(values) / len(values)


def pearson(x, y):
    """皮尔逊相关系数；任一序列为常数或长度不足时返回 None。"""
    n = len(x)
    if n != len(y) or n < 2:
        return None
    mean_x = _mean(x)
    mean_y = _mean(y)
    dx = [value - mean_x for value in x]
    dy = [value - mean_y for value in y]
    sxx = sum(value * value for value in dx)
    syy = sum(value * value for value in dy)
    if sxx <= 0.0 or syy <= 0.0:
        return None
    sxy = sum(a * b for a, b in zip(dx, dy))
    return sxy / math.sqrt(sxx * syy)


def spearman(x, y):
    """斯皮尔曼等级相关系数；常数变量返回 None。"""
    if len(x) != len(y) or len(x) < 2:
        return None
    return pearson(average_ranks(x), average_ranks(y))


def _solve_linear(matrix, vector):
    """高斯消元（部分选主元）解 A x = b；奇异时返回 None。"""
    n = len(vector)
    work = [list(row) + [vector[index]] for index, row in enumerate(matrix)]
    for col in range(n):
        pivot = max(range(col, n), key=lambda row: abs(work[row][col]))
        if abs(work[pivot][col]) < 1e-12:
            return None
        work[col], work[pivot] = work[pivot], work[col]
        pivot_value = work[col][col]
        for row in range(col + 1, n):
            factor = work[row][col] / pivot_value
            if factor == 0.0:
                continue
            for j in range(col, n + 1):
                work[row][j] -= factor * work[col][j]
    solution = [0.0] * n
    for row in range(n - 1, -1, -1):
        remainder = work[row][n] - sum(work[row][j] * solution[j] for j in range(row + 1, n))
        solution[row] = remainder / work[row][row]
    return solution


def _residualize(target, controls):
    """对控制变量做线性投影，返回残差；设计矩阵共线时返回 None。"""
    n = len(target)
    k = len(controls)
    matrix = [[0.0] * k for _ in range(k)]
    vector = [0.0] * k
    for i in range(k):
        vector[i] = sum(controls[i][row] * target[row] for row in range(n))
        for j in range(i, k):
            value = sum(controls[i][row] * controls[j][row] for row in range(n))
            matrix[i][j] = value
            matrix[j][i] = value
    coefficients = _solve_linear(matrix, vector)
    if coefficients is None:
        ridge = 1e-10 * sum(matrix[i][i] for i in range(k)) / k
        if ridge > 0.0:
            for i in range(k):
                matrix[i][i] += ridge
            coefficients = _solve_linear(matrix, vector)
    if coefficients is None:
        return None
    return [target[row] - sum(coefficients[j] * controls[j][row] for j in range(k))
            for row in range(n)]


def partial_correlation(x, y, controls):
    """控制 controls（列的列表）后的偏相关；无可控制列时退化为皮尔逊。"""
    if not controls:
        return pearson(x, y)
    residual_x = _residualize(x, controls)
    residual_y = _residualize(y, controls)
    if residual_x is None or residual_y is None:
        return None
    return pearson(residual_x, residual_y)


def partial_spearman(x, y, raw_controls):
    """先对全部变量取秩再做偏相关；控制列为常数时返回 None。"""
    ranked_x = average_ranks(x)
    ranked_y = average_ranks(y)
    ranked_controls = []
    for column in raw_controls:
        if len(set(column)) < 2:
            return None
        ranked_controls.append(average_ranks(column))
    return partial_correlation(ranked_x, ranked_y, ranked_controls)


def _beta_continued_fraction(a, b, x):
    qab = a + b
    qap = a + 1.0
    qam = a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < 1e-30:
        d = 1e-30
    d = 1.0 / d
    result = d
    for m in range(1, 201):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < 1e-30:
            d = 1e-30
        c = 1.0 + aa / c
        if abs(c) < 1e-30:
            c = 1e-30
        d = 1.0 / d
        result *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < 1e-30:
            d = 1e-30
        c = 1.0 + aa / c
        if abs(c) < 1e-30:
            c = 1e-30
        d = 1.0 / d
        delta = d * c
        result *= delta
        if abs(delta - 1.0) < 1e-14:
            break
    return result


def regularized_incomplete_beta(a, b, x):
    """正则化不完全贝塔函数 I_x(a, b)，用于 t 分布尾概率。"""
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    log_prefix = (math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
                  + a * math.log(x) + b * math.log1p(-x))
    prefix = math.exp(log_prefix)
    if x < (a + 1.0) / (a + b + 2.0):
        return prefix * _beta_continued_fraction(a, b, x) / a
    return 1.0 - prefix * _beta_continued_fraction(b, a, 1.0 - x) / b


def t_two_sided_p(correlation, sample_size, control_count):
    """相关系数的 t 近似双侧 p 值（斯皮尔曼沿用同一近似）。"""
    if correlation is None or abs(correlation) >= 1.0:
        return None
    degrees = sample_size - 2 - control_count
    if degrees <= 0:
        return None
    t_stat = abs(correlation) * math.sqrt(degrees / (1.0 - correlation * correlation))
    x = degrees / (degrees + t_stat * t_stat)
    return regularized_incomplete_beta(degrees / 2.0, 0.5, x)


def benjamini_hochberg(p_values):
    """BH 校正；None 位置不参与校正并原样返回。"""
    valid = [(index, value) for index, value in enumerate(p_values) if value is not None]
    adjusted = [None] * len(p_values)
    if not valid:
        return adjusted
    ordered = sorted(valid, key=lambda item: item[1])
    count = len(ordered)
    running = 1.0
    for reverse_position in range(count, 0, -1):
        index, value = ordered[reverse_position - 1]
        running = min(running, value * count / reverse_position)
        adjusted[index] = min(1.0, running)
    return adjusted


def percentile(sorted_values, quantile):
    """线性插值分位数（quantile 取 0 到 1），输入需已升序。"""
    if not sorted_values:
        return None
    if len(sorted_values) == 1:
        return sorted_values[0]
    position = quantile * (len(sorted_values) - 1)
    lower = int(math.floor(position))
    upper = min(lower + 1, len(sorted_values) - 1)
    fraction = position - lower
    return sorted_values[lower] * (1.0 - fraction) + sorted_values[upper] * fraction
