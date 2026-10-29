# ACCEPTANCE — 三维数字孪生胃

逐项实测验收。每一行的"实测"列都是从生成的报告里读出来的，不是从叙述里抄的。
未达成的条目**不删不藏**，集中列在最后一节。

## A. 几何

| # | 要求 | 实测 | 判定 | 说明 |
|---|---|---|---|---|
| A1 | 真实分割网格 | TotalSegmentator v2.0.1 `stomach`，受试者 s0236 | **达成** | CC-BY-4.0，仅取掩膜，未下载影像 |
| A2 | 面数 13 万–16 万 | 156,492 | **达成** |  |
| A3 | glb ≤ 7 MB | 5.09 MB | **达成** |  |
| A4 | 分区尽量从本受试者形态划出 | 胃底=实测水平面；角切迹=实测曲率峰 κ=0.37788 rad/mm（13.89× 中位）；贲门领与幽门管=声明值 | **部分达成** | 贲门领 20 mm 与幽门管 20 mm 为构造，已声明 |
| A5 | 报实测容积并与空腹/餐后对照 | 433.3 mL | **达成（照实报）** | 落在两参考区间之间，CT 未声明禁食状态 |
| A6 | 与小肠页在幽门处相接 | 接触带 6 体素，幽门点到质心 0.35 mm | **达成** | 接触带很薄，照实报 |
| A7 | 与肝脏页在门静脉处相接 | 门静脉汇合开口锚在 `portal_vein_and_splenic_vein` 标签表面 | **达成** |  |

## B. 管道树

| # | 要求 | 实测 | 判定 | 说明 |
|---|---|---|---|---|
| B1 | CCO 生长 | 6 棵树 / 2438 段 / 847 分叉 / 853 终端 / 总长 10242.2 mm | **达成** |  |
| B2 | Murray 每个分叉最大相对残差为 0 | 六棵树全部 `0.0`（共核验 847 个分叉） | **达成** | 结构性立方律，0 是构造属性而非拟合结果，已声明 |
| B3 | 多树不相交 | 15 对树全部非负；全局最小表面间隙 **0.7475 mm**（right_gastric ↔ right_gastroepiploic） | **达成** | 不等于求解器目标，是真实测量 |
| B4 | 不静默删除段 | 删除 0 段；推出可通行空间的节点 0 个 | **达成** | 未接入目标 653 个，已逐树报出 |
| B5 | 节点位移有上界并报最大值 | 生长后松弛**未运行**（0 次扫掠）；骨架求解位移上限 20.0 mm，实际最大 0.0 mm | **达成** | 生长时即以保守半径施加净空，故无需事后松弛 |

### 逐树

| 树 | 段 | 分叉 | 终端 | 直径 mm | 总长 mm | Murray 残差 | 流量 mL/min | 主干口径参考已核实 |
|---|---|---|---|---|---|---|---|---|
| left_gastric | 453 | 176 | 177 | 0.527–2.9588 | 1793.6 | `0.0` | 11.073 | 是 |
| right_gastric | 332 | 118 | 119 | 0.357–1.7561 | 1390.8 | `0.0` | 6.376 | 是 |
| left_gastroepiploic | 450 | 176 | 177 | 0.3948–2.2167 | 1605.9 | `0.0` | 5.552 | **否（假设值）** |
| right_gastroepiploic | 472 | 180 | 181 | 0.4517–2.5551 | 1850.9 | `0.0` | 12.736 | **否（假设值）** |
| short_gastric | 381 | 117 | 118 | 0.3016–1.4794 | 1778.8 | `0.0` | 16.034 | **否（假设值）** |
| gastric_veins | 350 | 80 | 81 | 0.7084–3.0652 | 1822.2 | `0.0` | 51.771 | **否（假设值）** |

### 跨树最小表面间隙（全部 15 对）

| 树对 | 实测间隙 mm | 发生级别 | 等于求解器目标 |
|---|---|---|---|
| right_gastric ↔ right_gastroepiploic | 0.7475 | right_gastric_trunk ↔ right_gastroepiploic_trunk | 否 |
| left_gastroepiploic ↔ short_gastric | 1.1267 | left_gastroepiploic_course ↔ short_gastric_terminal | 否 |
| right_gastric ↔ gastric_veins | 1.4082 | right_gastric_course ↔ gastric_veins_trunk | 否 |
| short_gastric ↔ gastric_veins | 1.662 | short_gastric_terminal ↔ gastric_veins_terminal | 否 |
| left_gastric ↔ gastric_veins | 1.6736 | left_gastric_terminal ↔ gastric_veins_terminal | 否 |
| left_gastric ↔ left_gastroepiploic | 2.1063 | left_gastric_terminal ↔ left_gastroepiploic_terminal | 否 |
| left_gastric ↔ short_gastric | 2.1821 | left_gastric_terminal ↔ short_gastric_terminal | 否 |
| left_gastric ↔ right_gastric | 2.2161 | left_gastric_course ↔ right_gastric_terminal | 否 |
| right_gastric ↔ short_gastric | 2.2854 | right_gastric_terminal ↔ short_gastric_terminal | 否 |
| left_gastroepiploic ↔ right_gastroepiploic | 2.6573 | left_gastroepiploic_terminal ↔ right_gastroepiploic_terminal | 否 |
| right_gastroepiploic ↔ short_gastric | 3.0366 | right_gastroepiploic_terminal ↔ short_gastric_terminal | 否 |
| right_gastroepiploic ↔ gastric_veins | 3.3799 | right_gastroepiploic_terminal ↔ gastric_veins_terminal | 否 |
| left_gastroepiploic ↔ gastric_veins | 4.5123 | left_gastroepiploic_trunk ↔ gastric_veins_terminal | 否 |
| right_gastric ↔ left_gastroepiploic | 6.0773 | right_gastric_terminal ↔ left_gastroepiploic_terminal | 否 |
| left_gastric ↔ right_gastroepiploic | 8.1275 | left_gastric_terminal ↔ right_gastroepiploic_terminal | 否 |

**标注 `at_solver_target`：** 骨架求解阶段的最小间隙 `0.25` mm **恰好等于**求解器打磨目标 0.25 mm。
这是求解器的属性，**不是对几何的独立测量**。最终报告的 15 对间隙都在该阶段之后、
在生长完成的真实半径上重新测量，没有一个等于求解器目标。

## C. 流动

| # | 要求 | 实测 | 判定 | 说明 |
|---|---|---|---|---|
| C1 | v(s)=Q(s)/A(s) 烘进路径 | 118 条通路 × 64 采样点 | **达成** |  |
| C2 | 逐采样点零逆流审计 | 7,552 个采样点，逆流 **0** 个、零速 0 个、不前进 0 个 | **达成** | 是逐点计数，不是断言 |
| C3 | 峰值流速落在已发表范围并引文献 | 腹腔干 PSV 90–190 cm/s 已引；**命名胃动脉分支未找到已发表流速范围** | **部分达成** | 无参考值者按实测报，不与臆造范围比较 |
| C4 | 显示增益声明 | 血管层 **1×**（无增益） | **达成** | 运动层单独声明，微观视图只有空间放大 |

| 树 | 时均主干 cm/s | 峰值上界 cm/s | 已发表参考 | 判定 |
|---|---|---|---|---|
| left_gastric | 2.775 | 18.448 | — | 无公开参考值 |
| right_gastric | 6.368 | 42.338 | — | 无公开参考值 |
| left_gastroepiploic | 2.667 | 17.732 | — | 无公开参考值 |
| right_gastroepiploic | 7.118 | 47.325 | — | 无公开参考值 |
| short_gastric | 19.85 | 131.976 | — | 无公开参考值 |
| gastric_veins | 11.945 | 79.418 | 10.0–25.0 | **落在范围内** |

## D. 动画（用户点名要）

| # | 要求 | 声明值 | 浏览器内实测 | 判定 |
|---|---|---|---|---|
| D1 | 胃窦收缩波频率 | 3.0 cpm | 3 cpm | **达成** |
| D2 | 传播速度 | 胃体 3.0 → 终末 7.5 mm/s | 中位 7.248 (p10 7.244 · p90 7.248) mm/s | **达成** |
| D3 | 是行进波不是整体缩放 | — | t=0 未位移采样环 44/112（0.3929）| **达成** |
| D4 | 壁顶点真的在动 | — | 探针 4 个，相位间移动 2 个，最大位移 2.5076 mm | **达成** |
| D5 | 回推 retropulsion | 筛分尺寸 2.0 mm | 浏览器计数：通过率 1（事件 17）| **部分达成** |

**D5 必须说清楚。** 上表那个"通过率 1.0"**不能**当作生理通过率读。
交付包里的 `mobile_report.json` 是用 **8 秒**的回推审计窗口测的，
而一个蠕动波周期是 **20 秒**——窗口比一个周期还短，所以大颗粒还没来得及被推回胃体，
审计只数到了小颗粒的"通过"事件（17 / 12 / 12 个），回推事件计数为 0。
**机制本身是实现了的**（剖面视图里五档颗粒随壁运动平流、按筛分尺寸分流，
计数由浏览器完成），但这一项的**定量证据不完整**。
`tools/benchmark_stomach3d.py` 里的窗口已改为 **1800 步 / 60 秒 = 三个完整波周期**，
但那次重测作业在 Modal 冷启动阶段被 15 分钟空闲看门狗杀掉，**没有跑完**，
所以随包交付的仍是 8 秒窗口的那份报告。照实列入未达成条目。
| D6 | 容受性舒张 | 胃底容积变化、压力近似不变 | 定性，**未建压力–容积曲线** | **部分达成** |

### 一个未通过的独立核验

the mean of the piecewise profile over the antrum is 3.88 mm/s against an independently published antral average of 5.9 mm/s - a DISAGREEMENT of -34%. The profile was built only from the corpus value and the terminal-antral acceleration, so this is an independent check and it does not pass.  The likely cause is geometric: this subject's antrum measures 72.5 mm along the centreline against a published antral length of 30-40 mm, so the slow proximal portion dominates the average, whereas in a 35 mm antrum almost the whole length would lie inside the accelerating terminal 28 mm.  The profile was NOT retuned to hit 5.9 mm/s.

`mean_antral_within_20pct = false`。

## E. 解剖剖面与微观单元

| # | 要求 | 实测/声明 | 判定 |
|---|---|---|---|
| E1 | 真实比例微观结构 | 泌酸腺总深 800 µm（胃小凹 200 + 腺体 600，比 1:4） | **达成** |
| E2 | 放大倍数必须声明 | 腺体柱 **60×**；黏液凝胶嵌套 **8×**，合计 **480×** | **达成** |
| E3 | pH 梯度真实比例着色 | 屏障厚 106.0 µm，场景内 6.36 mm（嵌套 50.88 mm），48 层切片 | **达成** |
| E4 | 质子浓度梯度 | 由声明端点算出：pH 1.5→7.0 为 **3.16e+05** 倍，pH 3.5→7.0 为 **3.16e+03** 倍 | **达成（照算出来报）** |

### 关于"约 200 µm"这个数字

the micro view is built at the HUMAN antral thickness (106 um) because it is the only human in-vivo measurement in that study.  The frequently quoted ~200 um is the rat corpus value and is shown alongside it, labelled as rat, rather than being used as if human.

- 人体**胃窦**在体实测：**106.0 µm**（最小 50）—— 本页微观视图采用此值
- 大鼠**胃体**在体实测：189.0 µm（牢固黏附层约 80 µm）—— 这才是"~200 µm"的出处，**不是人体胃体值**
- 大鼠**胃窦**：274.0 µm

## F. 分子层

| # | 要求 | 实测 | 判定 |
|---|---|---|---|
| F1 | 接真实单细胞数据 | CELLxGENE Census `2025-11-08`，`tissue_general == 'stomach'`，`is_primary_data == True` | **达成** |
| F2 | 规模 | 131 细胞类型 × 61 基因，557,791 细胞 / 141 donor / 8 数据集 | **达成** |
| F3 | 读取路线 | 连续 soma_joinid 区块 4 块 × 150,000 行，**失败 0 块**，耗时 9890.1 s | **达成** |
| F4 | 基因面板与其它器官一致 | 61/61 解析，顺序与 kidney/small_intestine 一致，缺失 [] | **达成** |
| F5 | ATP4A/ATP4B 作为附加基因单独报 | `ATP4A`, `ATP4B`，单独文件 `stomach_proton_pump_expression.json` | **达成** |

**面板外附加基因的声明：** NOT part of the 61-gene ADME panel used for heart, liver, kidney and small intestine.  Added for the stomach only, because the gastric proton pump is the organ's defining pharmacological target (PPI target). Read in the same Census X pass, reported separately so the ADME panel is not silently widened.

### 关键细胞类型是否在 Census 胃队列中存在

| 细胞类型 | Census 标签 | 细胞数 | 状态 |
|---|---|---|---|
| 壁细胞 | `parietal cell` | 1,802 | **存在** |
| 主细胞 (Census: peptic cell) | `peptic cell` | 8,840 | **存在** |
| G 细胞 | `type G enteroendocrine cell` | 1,892 | **存在** |
| 表面黏液细胞 | `foveolar cell of stomach` | 41,796 | **存在** |
| 颈黏液细胞 | `mucous neck cell of gastric gland` | 26,574 | **存在** |
| ECL 细胞 | `enterochromaffin-like cell` | 0 | **缺失** |
| D 细胞 | `type D enteroendocrine cell` | 0 | **缺失** |

**缺失的关键细胞类型：ECL 细胞、D 细胞。** 没有用别的细胞类型或组织平均值顶替。

特别说明：本队列里**确实存在** `type EC enteroendocrine cell`（肠嗜铬细胞，5-羟色胺）、`type X enteroendocrine cell` 与 `P/D1 enteroendocrine cell`（胃饥饿素细胞），但它们**都不是**泌酸腺里那个分泌组胺的**肠嗜铬样（ECL）细胞**。把 EC 当 ECL 填进表里属于捏造，所以 ECL 行按缺失报告。同理，泛化的 `enteroendocrine cell` 标签混合了多个谱系，**不**用作 D 细胞的替身。

### 关键基因实测值与膜面

| 基因 | 蛋白 | 膜面 | 细胞类型 | 平均表达 log1p CP10K | 表达比例 % |
|---|---|---|---|---|---|
| `ATP4A` | H+/K+-ATPase α | **顶侧 apical** | parietal cell | 2.6147 | 77.52 |
| `ATP4A` | H+/K+-ATPase α | **顶侧 apical** | peptic cell | 0.0648 | 3.81 |
| `ATP4A` | H+/K+-ATPase α | **顶侧 apical** | type G enteroendocrine cell | 0.0111 | 1.16 |
| `ATP4A` | H+/K+-ATPase α | **顶侧 apical** | foveolar cell of stomach | 0.0177 | 1.30 |
| `ATP4B` | H+/K+-ATPase β | **顶侧 apical** | parietal cell | 2.6134 | 84.07 |
| `ATP4B` | H+/K+-ATPase β | **顶侧 apical** | peptic cell | 0.0881 | 4.33 |
| `ATP4B` | H+/K+-ATPase β | **顶侧 apical** | type G enteroendocrine cell | 0.0055 | 0.63 |
| `ATP4B` | H+/K+-ATPase β | **顶侧 apical** | foveolar cell of stomach | 0.0117 | 0.93 |
| `CYP3A4` | CYP3A4 | 胞内（ER） | parietal cell | 0.0016 | 0.06 |
| `CYP3A4` | CYP3A4 | 胞内（ER） | peptic cell | 0.0006 | 0.08 |
| `CYP3A4` | CYP3A4 | 胞内（ER） | type G enteroendocrine cell | 0.0017 | 0.16 |
| `CYP3A4` | CYP3A4 | 胞内（ER） | foveolar cell of stomach | 0.0052 | 0.60 |
| `CYP3A5` | CYP3A5 | 胞内（ER） | parietal cell | 0.4411 | 31.63 |
| `CYP3A5` | CYP3A5 | 胞内（ER） | peptic cell | 0.0392 | 2.79 |
| `CYP3A5` | CYP3A5 | 胞内（ER） | type G enteroendocrine cell | 0.1627 | 15.59 |
| `CYP3A5` | CYP3A5 | 胞内（ER） | foveolar cell of stomach | 0.9131 | 48.97 |
| `ABCB1` | P-gp / MDR1 | 顶侧 apical | parietal cell | 0.0000 | 0.00 |
| `ABCB1` | P-gp / MDR1 | 顶侧 apical | peptic cell | 0.0010 | 0.05 |
| `ABCB1` | P-gp / MDR1 | 顶侧 apical | type G enteroendocrine cell | 0.0011 | 0.16 |
| `ABCB1` | P-gp / MDR1 | 顶侧 apical | foveolar cell of stomach | 0.0010 | 0.08 |
| `ABCG2` | BCRP | 顶侧 apical | parietal cell | 0.0012 | 0.22 |
| `ABCG2` | BCRP | 顶侧 apical | peptic cell | 0.0011 | 0.11 |
| `ABCG2` | BCRP | 顶侧 apical | type G enteroendocrine cell | 0.0309 | 3.33 |
| `ABCG2` | BCRP | 顶侧 apical | foveolar cell of stomach | 0.0060 | 0.48 |
| `SLCO2B1` | OATP2B1 | 顶侧 apical | parietal cell | 0.0006 | 0.06 |
| `SLCO2B1` | OATP2B1 | 顶侧 apical | peptic cell | 0.0013 | 0.12 |
| `SLCO2B1` | OATP2B1 | 顶侧 apical | type G enteroendocrine cell | 0.0003 | 0.05 |
| `SLCO2B1` | OATP2B1 | 顶侧 apical | foveolar cell of stomach | 0.0076 | 0.77 |
| `SLC15A1` | PEPT1 | 顶侧 apical | parietal cell | 0.0009 | 0.11 |
| `SLC15A1` | PEPT1 | 顶侧 apical | peptic cell | 0.0021 | 0.27 |
| `SLC15A1` | PEPT1 | 顶侧 apical | type G enteroendocrine cell | 0.0466 | 4.86 |
| `SLC15A1` | PEPT1 | 顶侧 apical | foveolar cell of stomach | 0.0086 | 0.64 |
| `SLC22A1` | OCT1 | 基底侧 basolateral | parietal cell | 0.0017 | 0.17 |
| `SLC22A1` | OCT1 | 基底侧 basolateral | peptic cell | 0.0001 | 0.01 |
| `SLC22A1` | OCT1 | 基底侧 basolateral | type G enteroendocrine cell | 0.0008 | 0.05 |
| `SLC22A1` | OCT1 | 基底侧 basolateral | foveolar cell of stomach | 0.0007 | 0.06 |
| `SLC22A3` | OCT3 | 基底侧 basolateral | parietal cell | 0.0062 | 0.67 |
| `SLC22A3` | OCT3 | 基底侧 basolateral | peptic cell | 0.0033 | 0.46 |
| `SLC22A3` | OCT3 | 基底侧 basolateral | type G enteroendocrine cell | 0.0030 | 0.53 |
| `SLC22A3` | OCT3 | 基底侧 basolateral | foveolar cell of stomach | 0.0164 | 1.92 |
| `ABCC3` | MRP3 | 基底侧 basolateral | parietal cell | 0.1340 | 11.76 |
| `ABCC3` | MRP3 | 基底侧 basolateral | peptic cell | 0.0409 | 3.13 |
| `ABCC3` | MRP3 | 基底侧 basolateral | type G enteroendocrine cell | 0.0587 | 6.87 |
| `ABCC3` | MRP3 | 基底侧 basolateral | foveolar cell of stomach | 0.6295 | 38.36 |


## G. 性能与移动端

测试硬件：**measured on a data-centre NVIDIA A10G, not a laptop and not a real phone**。

Two numbers per view.  frame_stats comes from the animation loop and in headless Chromium reflects the compositor cadence as much as the scene - headless has NO VSYNC, so requestAnimationFrame is uncapped and must not be read as a user-visible frame rate.  render_cost times back-to-back renderer.render() calls each followed by gl.finish(), which is the real per-frame cost, and compares it with the 16.67 ms 60 fps budget.  Read render_cost.

| 视口 | unmasked renderer | 视图 | ms/帧 | 对 16.67 ms 余量 | 帧时 p95 ms | DOM 溢出 | 水平溢出 px |
|---|---|---|---|---|---|---|---|
| desktop | `ANGLE (NVIDIA, Vulkan 1.4.341 (NVIDIA NVIDIA A10G (0x00002237)), NVIDIA)` | organ | 0.2525 | 66.02× | 4.6 | 0 | 0 |
| desktop | `ANGLE (NVIDIA, Vulkan 1.4.341 (NVIDIA NVIDIA A10G (0x00002237)), NVIDIA)` | vessels | 0.255 | 65.37× | 4.3 | 0 | 0 |
| desktop | `ANGLE (NVIDIA, Vulkan 1.4.341 (NVIDIA NVIDIA A10G (0x00002237)), NVIDIA)` | section | 0.1825 | 91.34× | 0.4 | 0 | 0 |
| desktop | `ANGLE (NVIDIA, Vulkan 1.4.341 (NVIDIA NVIDIA A10G (0x00002237)), NVIDIA)` | gland | 3.635 | 4.59× | 4.9 | 0 | 0 |
| desktop | `ANGLE (NVIDIA, Vulkan 1.4.341 (NVIDIA NVIDIA A10G (0x00002237)), NVIDIA)` | molecular | 0.1375 | 121.24× | 0.3 | 0 | 0 |
| phone390 | `ANGLE (NVIDIA, Vulkan 1.4.341 (NVIDIA NVIDIA A10G (0x00002237)), NVIDIA)` | organ | 0.2425 | 68.74× | 2.5 | 0 | 0 |
| phone390 | `ANGLE (NVIDIA, Vulkan 1.4.341 (NVIDIA NVIDIA A10G (0x00002237)), NVIDIA)` | vessels | 0.24 | 69.46× | 2.2 | 0 | 0 |
| phone390 | `ANGLE (NVIDIA, Vulkan 1.4.341 (NVIDIA NVIDIA A10G (0x00002237)), NVIDIA)` | section | 0.1775 | 93.92× | 0.5 | 0 | 0 |
| phone390 | `ANGLE (NVIDIA, Vulkan 1.4.341 (NVIDIA NVIDIA A10G (0x00002237)), NVIDIA)` | gland | 3.3775 | 4.94× | 4.4 | 0 | 0 |
| phone390 | `ANGLE (NVIDIA, Vulkan 1.4.341 (NVIDIA NVIDIA A10G (0x00002237)), NVIDIA)` | molecular | 0.1425 | 116.98× | 0.3 | 0 | 0 |
| phone360 | `ANGLE (NVIDIA, Vulkan 1.4.341 (NVIDIA NVIDIA A10G (0x00002237)), NVIDIA)` | organ | 0.245 | 68.04× | 2.4 | 0 | 0 |
| phone360 | `ANGLE (NVIDIA, Vulkan 1.4.341 (NVIDIA NVIDIA A10G (0x00002237)), NVIDIA)` | vessels | 0.2525 | 66.02× | 2.2 | 0 | 0 |
| phone360 | `ANGLE (NVIDIA, Vulkan 1.4.341 (NVIDIA NVIDIA A10G (0x00002237)), NVIDIA)` | section | 0.1775 | 93.92× | 0.5 | 0 | 0 |
| phone360 | `ANGLE (NVIDIA, Vulkan 1.4.341 (NVIDIA NVIDIA A10G (0x00002237)), NVIDIA)` | gland | 3.1725 | 5.25× | 4.4 | 0 | 0 |
| phone360 | `ANGLE (NVIDIA, Vulkan 1.4.341 (NVIDIA NVIDIA A10G (0x00002237)), NVIDIA)` | molecular | 0.1375 | 121.24× | 0.4 | 0 | 0 |

`gpu_measured = true` —— 页面拿到的是真实 GPU renderer 字符串，不是 SwiftShader。

审计门槛全部通过：`problems: []`。逐项：
- unmasked renderer 全部为真实 GPU 字符串，`gpu_measured = true`
- 着色器链接：每个视口 12 个程序，**链接失败 0 个**
- 粒子逐步审计：**逐单步**分类，三个视口合计 1380600 步，**逆流 0 步**、非有限值 0 个；路径回绕（s -= L）单独计数 7342 次
- 水平溢出 0 px、DOM 越界元素 0 个（三个视口 × 五个视图全部为 0）

## H. 未达成 / 未测量 / 合成 / 超范围未调参

这一节是本文件的重点。以下条目**没有**在别处被淡化处理。

1. **胃窦长度超出已发表范围且未调参**：实测 72.5 mm，参考 30.0–40.0 mm。角切迹按形态实测定位，没有为了把胃窦压进 30–40 mm 而移动它。
2. **蠕动波速的独立核验未通过**：分段剖面在胃窦上的平均速度 3.882 mm/s，独立发表的胃窦平均值 5.9 mm/s，相差 −34%。未重调剖面。
3. **胃动脉总入流是假设而非实测**：取心脏模型 CO 的 1.0% = 51.771 mL/min。`published_whole_stomach_flow_located = false` —— 本次构建没有找到可用于校验的人体全胃血流实测值。
4. **命名胃动脉分支没有已发表流速参考范围**：胃左/右、网膜左/右、胃短动脉的峰值或平均流速范围本次构建未找到，相应行按实测报出而**不与任何臆造范围比较**。
5. **峰值流速只是上界，不是预测**：峰/均比 6.6488 测于主动脉根部并原样应用到远端内脏分支，真实分支是被衰减过的。
6. **反算终末口径全部超出已发表小动脉范围**：6 棵树的反算终末直径均大于 50.0–150.0 µm。10² 个终末的立方律级联在锚定公开主干口径时到不了小动脉口径，**没有发明任何系数去凑**。
7. **653 个目标点未能接入**（其中 640 个属胃静脉树）。胃静脉树最终 81 个终端，低于请求的 180 个。这是保守净空约束的代价，**不是被隐藏的失败**，也没有靠删段来掩盖。
8. **三条主干口径没有找到已发表直径**：胃网膜左、胃网膜右、胃短动脉，以及胃静脉根部。这四个值标为假设值（`reference_located = false`）。腹腔干 7.0 mm、胃左 2.8 mm、胃右 1.6 mm 有文献。
9. **骨架求解阶段的最小间隙 0.25 mm 恰好等于求解器目标 0.25 mm**，已标 `at_solver_target`。这是求解器的属性，不是对几何的独立测量。
10. **生长后松弛被刻意关闭**（0 次扫掠）。每一种配置都会把一棵树的主干拖进另一棵树的终端里，重新制造相互穿插。净空改为在生长时以保守（最终）半径施加。
11. **`subdivide()` 的横向弓形偏置被关闭**（bow = 0）。它是装饰性曲率，但施加在所有净空检查之后，是最后一个造成相互穿插的来源。
12. **壁/腔分界是构造的**：源数据只有一个实心 `stomach` 标签，壁厚 4.0 mm 是公开充盈 CT 值，**非本受试者实测**。腔与壁的体积都继承这个假设。
13. **腹腔干与每一条胃动脉、胃静脉都是构造的**。只有开口（腹腔干起始水平、脾门上下极、门静脉汇合、肝门侧、胃十二指肠侧）是 contour_anchored。
14. **胃黏膜皱襞是构造的**：12 条、幅度 2.5 mm，源数据无此结构。
15. **腔内 pH 的区域值是文献值，不是本受试者实测**，泌酸腺密度分布是教科书区域组织学。两者都只能定性读。
16. **黏液屏障的 pH 剖面不是输运方程的解**：两个已发表端点之间在 pH 上线性插值。真实的近黏膜剖面是 S 形的，由 H⁺ 反向扩散与 HCO₃⁻ 分泌的平衡决定。`pH_profile_is_a_transport_solution = false`。
17. **"约 200 µm 屏障"是大鼠胃体值**，不是人体胃体值。本页按人体胃窦实测 106.0 µm 建，并把大鼠值标注为大鼠值并列显示。
18. **腺体口径与间距是假设值**（50.0 µm / 70.0 µm），由此算出的腺体密度 204.1/mm² 是**结果而不是实测**；本次构建没有找到人体胃腺密度的实测值。
19. **固体与液体半排空时间没有核实到单一来源**（`reference_located = false`），只能当量级用。共识标准用的是滞留率阈值（2 h > 60%、4 h > 10%），那两个值有来源。
20. **容受性舒张只做了定性**：有胃底容积变化与"压力近似不变"的说明，**没有任何压力–容积曲线**。
21. **没有任何标定系数**：本项目不把表达量或细胞丰度换算成酸分泌通量、电导或清除率。分子层与 pH 场的映射只是定性的。
22. **没有药物浓度、没有 PBPK 房室、没有 pH 依赖吸收的计算**。PPI 与幽门螺杆菌改变胃内 pH 进而影响其它药物吸收，只作为**文献陈述**放在界面上（`computed_by_this_project = false`）。
23. **大弯:小弯的逐站迹线定义结构上无法复现教科书比值**（实测 1.12），因为两条迹线都按同一条中心线参数化。轮廓弧定义给出 2.31，但投影会压扁器官，两个数都是真实三维弧长的**下界**。
24. **心率、心输出量、平均动脉压不是用户的个人数据**：来自本项目心脏模型，其驱动数据为 PhysioNet nsrdb/16265 公开记录。所有文件中 `is_user_personal_data = false`。用户的 Apple Health 导出始终未找到。
25. **幽门–十二指肠接触带只有 6 个体素**，很薄；这是幽门锚点的全部实测依据，照实报出。
26. **胃标签有第二个 4 体素的连通分量**（最大分量占比 0.999969），未做清理，照实报出。
27. **未建模**：胃淋巴、跨黏膜水分运动、幽门括约肌的独立肌层、迷走神经与肠神经系统的解剖、黏膜肌层、固有肌层的三层走向（斜/环/纵）、胃小凹的真实开口形态。
28. **没有做口腔、食管、结肠**——本任务只是胃。食管与十二指肠只作为**邻接锚定物**出现，按 26% 不透明度的上下文渲染，不是本页建模的器官。
29. 分子层：关键细胞类型 ECL 细胞、D 细胞 在 Census 胃队列中不存在，按缺失报告，未用替代值

## I. 引文核实等级

本次构建的文献值分两级标注，**区分"论文身份可核"与"数字可核"**：

| 等级 | 含义 |
|---|---|
| **A** | 论文标题与 DOI/URL 出现在**本会话留存的 web_search 结果文件**里，可直接回查 |
| **B** | 论文身份可核（同上），但**具体数值是从检索摘要片段读来的**，那些片段的正文**没有**留在本会话记录里，本次构建**也没有**回到全文逐一核对 |

**A 级的证据在哪里。** 检索结果的正文片段没有留在对话记录里（只留下 URL 与标题），
所以"论文身份可核"这一说法本身也必须可核。下表每一行的标识符都能在本会话的
`web_search` 结果文件里 grep 到，文件名一并列出：

| 标识符 | 留存的检索结果文件 |
|---|---|
| `doi:10.1152/ajpgi.00125.2010` | `srvtoolu_019PSRVCGDQznXboxndaJFNt.json` |
| `PMC5130547 / doi:10.1152/ajpgi.00255.2016` | `srvtoolu_019PSRVCGDQznXboxndaJFNt.json` |
| `doi:10.1152/ajpgi.2001.280.5.G922` | `srvtoolu_014x4StMWGUFy45G5SUTkfxP.json` |
| `doi:10.1152/ajpgi.2001.280.1.G121` | `srvtoolu_014x4StMWGUFy45G5SUTkfxP.json` |
| `Gastroenterology 0016-5085(93)90339-E` | `srvtoolu_01Pr6D1mSVHz8pGsd6L2draj.json` |
| `PMC12528509` | `srvtoolu_013YUbkoDfrhcVf4bD7vhjmM.json` |
| `PMC13440139` | `srvtoolu_013YUbkoDfrhcVf4bD7vhjmM.json` |
| `PMC10519692` | `srvtoolu_013YUbkoDfrhcVf4bD7vhjmM.json` |
| `PMID 14717069` | `srvtoolu_013YUbkoDfrhcVf4bD7vhjmM.json` |
| `Appl Sci 2024;14:12033` | `srvtoolu_013YUbkoDfrhcVf4bD7vhjmM.json` |
| `PMC12129437` | `srvtoolu_01HY2L1T3ewyiwWxc6DyRucK.json` |
| `PathologyOutlines stomachnormalhistology` | `srvtoolu_01Pr6D1mSVHz8pGsd6L2draj.json` |
| `Gastric Emptying Scintigraphy (snmjournals 47/2/111)` | `srvtoolu_01EEvTzLqwA8p1Y7MgwktWwD.json` |

这些文件位于本会话的 tool-results 目录下。**没有一个标识符是凭记忆写出来的**——
每一个都是从上述文件里的 `url` / `title` 字段取回的。

| 数值 | 来源 | 等级 |
|---|---|---|
| 慢波 3 cpm、胃体 3.0 mm/s、胃窦 5.9 mm/s、起搏区位置 | "Origin and propagation of human gastric slow-wave activity defined by high-resolution mapping", AJP-GI 2010, doi:10.1152/ajpgi.00125.2010 | 论文身份 **A** / 数值 **B** |
| 终末胃窦加速 3.3→7.5 mm/s、距幽门 6 mm 停止激动 | "Functional physiology of the human terminal antrum defined by high-resolution electrical mapping and computational model" (PMC5130547 / doi:10.1152/ajpgi.00255.2016) | 论文身份 **A** / 数值 **B** |
| 人胃窦黏液凝胶 106 µm、大鼠胃体 189 µm、大鼠胃窦 274 µm | "The adherent gastrointestinal mucus gel layer: thickness and physical state in vivo", AJP-GI 2001;280:G922 | 论文身份 **A** / 数值 **B** |
| 上皮表面 pH ≈ 7（腔内 pH ≈ 2 时） | "Intraluminal acid and gastric mucosal integrity", AJP-GI 2001;280:G121 | 论文身份 **A** / 数值 **B** |
| 胃体黏膜 0.6–0.8 mm、胃小凹 0.2 mm、腺体 0.6 mm、胃窦 0.8/0.4 mm | Gastroenterology 1993;104:1310（Ménétrier 系列中的正常值表） | 论文身份 **A** / 数值 **B** |
| 腹腔干 ~7 mm、胃左 ~2.7–3.6 mm、胃右 ~1.6–1.9 mm | PMC12528509 / PMC13440139 / PMC10519692 / Appl Sci 2024;14:12033 | 论文身份 **A** / 数值 **B** |
| 腹腔干发自 T12 的占比、胃左动脉起点位于 T12–L1 的占比 | PMID 14717069 / PMC13440139 | 论文身份 **A** / 具体百分比 **B** |
| 腹腔干容积流量与 PSV 正常范围 | "Gastrointestinal Ultrasound: Measurements and Normal Findings" (PMC12129437) | 论文身份 **A** / 数值 **B** |
| 排空滞留阈值（2 h >60%、4 h >10%） | 胃排空闪烁显像多学会共识 | 论文身份 **A** / 数值 **B** |
| 胃窦/幽门为远端 3–4 cm、贲门黏膜 1–15 mm、胃小凹:腺体比 | PathologyOutlines 胃正常组织学条目（转引 Am J Surg Pathol 2002;26:1207） | 论文身份 **A** / 数值 **B** |

**没有一个数字是 A 级数值。** 也就是说：本页引用的每一个文献数值都值得在全文里复核一遍；
论文本身都是可以点开的，但"7.0 ± 0.5"这类精确到小数位与样本量的写法，
在本次构建里**只到摘要片段为止**。先前版本的脚本里写过"7.06 mm（SD 1.33，范围 3.00–10.30）、260 例 CT 血管成像"、"385 例 CT 系列"、"T12 占 79.8%"这类精确数字，**已经按本节的标准收窄**，
因为那些精确值无法从本会话记录里回查。

## J. 补测与更正

构建过程中被发现并更正的事项，列此以免只在代码注释里留痕：

1. **中心线原先从贲门起算**，导致整个胃底穹窿落在器官自身弧长坐标之外，
   近端大弯完全没有被覆盖。改为从**胃底顶点**起算，贲门成为其上的一站
   （s = 34.1 mm）。
2. **角切迹第一次"检出"是数值假象**：κ = 1.636 rad/mm，来自平滑前迹线上
   近乎重合的相邻点使分母趋零。改为先平滑并按弧长重采样，再在固定 8 mm
   物理窗口上求转折，并加显著度判据（≥ 2× 该迹线中位转折）。
   更正后 κ = 0.37788 rad/mm，显著度 13.89×。
3. **大弯:小弯比值第一次报成 1.1**，看上去像个解剖错误。实际是定义问题：
   逐站迹线必然同长。补做了教科书定义（C 平面轮廓弧），两套都报，
   并写明逐站定义**结构上**不可能复现教科书比值。
4. **胃静脉与胃右动脉最初相互穿插 −2.52 mm**：两者被给了同一条小弯走行。
   改为声明的伴行偏置，并把两条脾动脉来源（胃短、胃网膜左）的开口
   分到脾门上下极——真实解剖里它们本来就不从同一点发出。
5. **生长后松弛越修越坏**：把一棵树的主干拖进另一棵树的终端。
   最终改为四阶段架构（铺骨架 → 以主干半径求解骨架互净空 → 冻结骨架、
   以保守最终半径生长 → 不做事后松弛）。
6. **最后一个穿插源是装饰性的弓形偏置**：`subdivide(bow=0.08)` 在所有净空
   检查之后把一段 30 mm 的血管横移 2.4 mm。关掉后 15 对树全部非负。
7. **"700 个采样环"这个说法被更正了**。运动剖面确实按 700 站采样，
   但几何是绑到中心线的 116 个点（2 mm 间距）上的，所以实际只有 **112** 个
   不同的站位承载黏膜顶点，未位移环数是在这 112 个里数的。
   声称测了 700 个环会夸大分辨率：700 个环摊在本器官上是 0.3 mm 间距，
   比 1.5 mm 的源体素还细。小肠页的"700 环"框架不能直接搬过来。
8. **逆流计数第一次报出 1124/445/451 步，是计数口径错误而不是真的逆流**。
   跨 30 步做首尾快照无法把"路径回绕（s -= L，一次大的负跳）"和真逆流分开，
   因为动画循环同时在推进同一批粒子，回绕之后再前进会落在 (-L, 0) 区间内任意位置。
   改为**逐单步**分类后，回绕被单独计数，逆流为 **0**。
9. **回推审计第一次在错误的视图里取数**：`retro.step()` 只在剖面视图运行，
   在器官视图读计数器只会得到全 0。改为由审计直接驱动 `api.retroStep()`。
10. **心脏面板是 51 个离子通道基因，不是 61 个 ADME 基因**：
   `stomach_proton_pump_expression.json` 的面板声明最初把心脏也算进
   "61 基因"那一组，已更正为"肝、肾、小肠"。
11. **过于精确的文献数字被收窄**。脚本里原先写着"腹腔干 7.06 mm（SD 1.33，范围 3.00–10.30）、260 例"、"385 例 CT 系列"、"T12 占 79.8%"。
   这些精确到小数位与样本量的值来自检索摘要片段，而片段正文没有留在会话记录里，
   本次构建也没有回到全文核对，**无法回查**。已按 I 节的 A/B 分级收窄为可核的论文身份加量级表述。几何与血管的任何实测数字都不受影响。
12. **渲染重新配色这一步当时没有被验证**。那个 cell 无论 `str.replace` 是否命中都会打印 "restyled"，而同一 cell 里的 grep 只检查了 `main.js`。
   事后已在 `materials.js` 与**交付的 tar 包内**同时 grep 确认：新的分区配色（胃底 0.80/0.33/0.30、胃窦 0.52/0.34/0.56、幽门 0.30/0.40/0.62）与 `roughness: 0.78, metalness: 0.0` 都在包里，剖面截图上也能看到紫色胃窦与琥珀色贲门，确认生效。
13. **Vulkan 后端的注释有歧义**。原注释写"vulkan 优先、gl-egl 作为回退"，读起来像有自动故障转移，但代码只是一个可用 `ST_ANGLE` 覆盖的单一 flag，**没有**运行时检测 SwiftShader 再重试的逻辑。注释已改为如实描述；判定始终以页面实际拿到的 unmasked renderer 字符串为准。
