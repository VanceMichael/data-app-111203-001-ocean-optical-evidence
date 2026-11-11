# 数据契约

海洋观测团队按传感器汇总水色反射率与配对采样覆盖率，原始浓度单位为毫克每立方米。

主测量 `reflectance` 的单位为无量纲，主体字段是 `station`，核算维度是 `sensor`。`observation_id` 是文件内唯一观测编号，主体允许重复观测。数值只接受有限值；空单元格解析为 None，不能直接当成零。额外列不进入当前核算。

| 字段 | 类型 | 约束 |
| --- | --- | --- |
| `observation_id` | str | 必填 |
| `station` | str | 必填 |
| `sensor` | str | 必填 |
| `cruise` | str | 允许空值，空值不是零 |
| `day` | int | 允许空值，空值不是零 |
| `reflectance` | float | 必填 |
| `plastic` | float | 允许空值，空值不是零 |
| `chlorophyll` | float | 允许空值，空值不是零 |
| `latitude` | float | 允许空值，空值不是零 |
| `detected` | int | 允许空值，空值不是零 |

时间字段 day、week、tick、hour、minute、cycle、start_day、reported_day、known_at、visible_at、opened、observed 均以各字段名所述的相对时间单位表示整数；不同字段并不必然代表同一种时间。缺少的登记码与浓度等信息保持为空。示例只表达产品数据结构，不能代替实际调查资料。
