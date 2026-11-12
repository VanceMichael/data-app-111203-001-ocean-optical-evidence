# ocean-optical-evidence

海洋观测团队按传感器汇总水色反射率与配对采样覆盖率，原始浓度单位为毫克每立方米。

## 使用

运行需要 Python 3.11 或更新版本，没有第三方运行依赖。数据在本地处理，不请求外部接口。

```bash
python3 -m analysis_app summary
python3 -m analysis_app summary --input data/observations.csv --output output/summary.json
```

标准输出是完整 JSON，输出文件使用临时文件原子替换；失败退出码为 2，不覆盖既有结果。空数据集保留零计数与空分组，缺失数值不进入数值分母。数据指纹由按 observation_id 排序后的解析记录计算。

## 关联证据流程（evidence）

`evidence` 子流程回答“反射率与浓度是否存在统计关联、该关联对删失处理和控制变量是否稳定”，不输出污染因果结论或风险阈值。

```bash
python3 -m analysis_app evidence --cutoff 7 --seed 20261112 \
  --output output/evidence.json --narrative output/narrative.md
```

参数：`--cutoff` 为包含 `day<=cutoff` 的截点（整数）；`--seed` 为站点聚类自助法的固定种子（整数）；`--iterations` 为自助次数（默认 1000）；`--output` 写机器可读 JSON，`--narrative` 写采样会阅读的中文说明，二者均为临时文件原子替换，标准输出始终打印完整 JSON。

分析单位与口径：

- 只有 `plastic` 非空的配对观测行进入浓度关联；`detected=0` 表示低于检出限而非浓度为零，提供三种删失替代策略的敏感性对照：`as_reported`（保留原始测量值）、`half_lod`（以检出限二分之一替代）、`lod`（以检出限替代）。
- 同一站点的重复行先在站点内取均值，站点是关联与区间估计的独立单元；90% 区间使用以站点为唯一重采样单元的聚类自助法（5/95 百分位），同站重复记录不会被当作独立样本而缩窄区间。
- 每个传感器×波段×删失策略分别给出原始斯皮尔曼等级相关、控制叶绿素后的偏斯皮尔曼相关，以及同时控制叶绿素与航次（n-1 个航次指示列）的偏斯皮尔曼相关；p 值采用 t 近似，并在每个（删失策略, 调整口径）检验族内对同批传感器×波段检验统一做 BH 校正。
- 跨传感器方向冲突单独列出：同一口径下同时出现正向与负向点估计时标记 `direction_conflict=true`。
- 每一步保留样本去向（截点剔除、未配对、低于检出限、有效站点数、每站行数）；站点重复、全常数波段、仅一个有效站点、无配对浓度等情形输出明确的中文不可估计原因，而不是静默给出数值。
- 全部随机性来自按（种子, 传感器, 波段, 策略, 调整口径）派生的独立随机流；记录在进入流程前按 `observation_id` 排序，改变输入行顺序不改变相同种子的结论。

中文说明固定包含四节：数据口径与样本去向、关联方向与跨传感器冲突、稳定性对照、限制，并明确声明不构成污染因果结论、不提供风险阈值推断。

## 验证与容器

```bash
python3 -m unittest discover -s tests -v
python3 -m compileall -q analysis_app
docker build -t ocean-optical-evidence .
docker run --rm ocean-optical-evidence
```

容器内也可执行 `python -m unittest discover -s tests -v`。镜像构建阶段取得 Python 镜像，产品启动、测试和演示均无需下载。小规模样例的完整单元测试和 summary 可以连续执行，不需要等待定时任务。

## 代码与数据

- `analysis_app/io.py`：有类型的 CSV 读取、完整性校验及内容指纹。
- `analysis_app/schema.py`：字段类型、主体与聚合维度。
- `analysis_app/domain.py`：当前业务核算函数。
- `analysis_app/summary.py`：确定性的分组与数值口径。
- `analysis_app/stats.py`：无第三方依赖的等级相关、偏相关、t 近似 p 值、BH 校正与分位数。
- `analysis_app/evidence.py`：站点聚类的关联证据流程、删失敏感性与中文说明。
- `analysis_app/cli.py`：本地命令及原子输出。
- `docs/data.md`：数据契约。

演示数据全部为本地合成，不对应现实机构、设备或人员；不能用来推断真实业务效果。
