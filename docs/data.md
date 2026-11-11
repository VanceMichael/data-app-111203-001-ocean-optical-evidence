# 数据契约

海洋观测团队按传感器汇总水色反射率与配对采样覆盖率，原始浓度单位为毫克每立方米。

主测量 `reflectance` 的单位为无量纲，主体字段是 `station`，核算维度是 `sensor`。`observation_id` 是文件内唯一观测编号，主体允许重复观测。数值只接受有限值；空单元格解析为 None，不能直接当成零。额外列不进入当前核算。

| 字段 | 类型 | 约束 |
| --- | --- | --- |
| `band` | str | 必填，blue 或 green |
| `detection_limit` | float | 必填，浓度检出限 |
| `observation_id` | str | 必填 |
| `station` | str | 必填 |
| `sensor` | str | 必填 |
| `cruise` | str | 必填 |
| `day` | int | 必填 |
| `reflectance` | float | 必填 |
| `plastic` | float | 允许空值，空值不是零 |
| `chlorophyll` | float | 必填 |
| `latitude` | float | 必填 |
| `detected` | int | 必填 |

时间字段 day、week、tick、hour、minute、cycle、start_day、reported_day、known_at、visible_at、opened、observed 均以各字段名所述的相对时间单位表示整数；不同字段并不必然代表同一种时间。缺少的登记码与浓度等信息保持为空。示例只表达产品数据结构，不能代替实际调查资料。

## 业务字段语义

plastic 为毫克每立方米，空值表示未配对。detected=0 表示低于检出限，detected=1 表示可检出。反射率为无量纲；chlorophyll 为合成浓度代理，latitude 为度。station 可跨传感器和航次重复。

band 标识波段，每个传感器均包含多个波段；detection_limit 与 plastic 单位相同，是该记录的检出限。
