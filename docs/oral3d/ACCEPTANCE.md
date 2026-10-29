# ACCEPTANCE — 三维数字孪生口腔（逐项实测验收）

本文件的每一个数字都是从构建产物的 JSON 里读回来的，没有一个是手写的。凡实测值落在已发表范围之外的，**按算出的报，不调参**。

## 0. 一句话结论

- 真实 CT 分割骨性框架 **135.25 mL**，整页 **159,826 个三角面**、glb **5.409 MB**（预算 ≤ 7 MB）。
- 七套管道系统，**911 个分叉**，Murray 定律最大相对残差 **0.0**。
- 流动审计：**3095 个采样点，0 个逆流**。
- 运动学在浏览器内实测（见第 6 节）。
- **跨器官受试者一致性被打破**（s0236 无头部覆盖），见 GEOMETRY_PROVENANCE.md 第 2 节。

## 1. 几何

| 量 | 实测 | 参考范围 | 落入 |
|---|---|---|---|
| `bicondylar_breadth_mm` | 99 | 112 – 140 | **否** |
| `bigonial_breadth_mm` | 94.5 | 88 – 112 | 是 |
| `condylion_menton_mm` | 111 | 100 – 132 | 是 |
| `mandibular_arch_width_canine_mm` | 38.9 | 31 – 48 | 是 |
| `mandibular_arch_width_molar_mm` | 69.7 | 52 – 70 | 是 |
| `maxillary_arch_width_canine_mm` | 43.1 | 31 – 48 | 是 |
| `maxillary_arch_width_molar_mm` | 62.9 | 55 – 72 | 是 |
| `mandible_vertical_extent_mm` | 87.0 | （无姿势无关参考） | 不做判定 |
| `mandibular_arch_perimeter_mm` | 115.7 | （无标准参考） | 不做判定 |
| `maxillary_arch_perimeter_mm` | 122.1 | （无标准参考） | 不做判定 |

- `bicondylar_breadth_mm` 的定义是**两侧髁突点云质心之间的距离**，实测 99.0 mm；同一块骨的**最大横向跨度**是 112.5 mm。两者量的不是同一件事，都列出来。
- 口咽水平核验：C1 上缘 y = 188.5 mm，咬合平面 y = 149.5 mm，差 -39.0 mm。the hard palate of an adult in neutral head position lies at roughly the level of the C1 anterior arch; this difference is reported as measured, not corrected

### 面数与体积预算

| 项 | 值 |
|---|---|
| glb 三角面总数 | 159,826 |
| glb 顶点总数 | 83,387 |
| glb 部件数 | 45 |
| glb 文件大小 | 5.409 MB |
| 骨 + 链锚点面数 | 56,600 |
| 构造软组织面数 | 103,226 |

## 2. 构造软组织（全部 contour_anchored = false）

| 量 | 实测回读 | 参考范围 | 落入 | 标记 |
|---|---|---|---|---|
| `tongue_volume_mL` | 119.4 | 60 – 110 | **否** |  |
| `tongue_length_mm` | 91.21 | 80 – 100 | 是 |  |
| `parotid_volume_mL` | 25 | 14 – 35 | 是 | **at_published_mean** — 按已发表均值缩放得到，是构造的属性，不是测量 |
| `submandibular_volume_mL` | 8.5 | 5 – 12 | 是 |  |
| `sublingual_volume_mL` | 2.5 | 1.5 – 4 | 是 |  |
| `oral_mucosa_area_cm2` | 232.3 | 170 – 215 | **否** |  |
| `hard_palate_length_mm` | 52 | 45 – 60 | 是 |  |
| `circumvallate_count` | 9 | 7 – 12 | 是 |  |
| `fungiform_count` | 280 | 200 – 400 | 是 |  |
| `foliate_folds_per_side` | 6 | 4 – 11 | 是 |  |

- 牙列：32 颗，FDI 两位编号。30 颗落在实测牙槽弓上，2 颗落在实测弓的切线外推段上（实测上颌弓周长 122.1 mm、下颌 115.7 mm；已发表 16 颗牙冠近远中径之和为 127.4 / 127.0 mm，本来就装不下）。
- 舌乳头：丝状 **900 个被画出**（the human tongue carries on the order of 10^5 filiform papillae; 900 are DRAWN. This is a declared subsample for rendering, not a count of the anatomy.），菌状 280 个，轮廓 9 个，叶状每侧 6 褶。
- 口腔黏膜总面积的组成：sum of the eight constructed mucosal sheets plus a declared 42%% share of the closed tongue surface as the dorsum; the 42%% share is DECLARED, not measured
- permeability_rank is an ORDERING ONLY (1 least .. 8 most permeable), anchored on the published sublingual > buccal > palatal ordering. No permeability coefficient, no flux calibration and no expression-to-transport conversion exists anywhere in this build.

## 3. 微观结构与声明的放大倍数

| 结构 | 真实尺寸 | 放大倍数 | 画出尺寸 |
|---|---|---|---|
| 轮廓乳头直径 | 2.0 mm | **60.0×** | 120.0 mm |
| 味蕾高 × 宽 | 65.0 × 45.0 µm | **60.0×** | 3.9000000000000004 mm 高 |
| hard_palate 上皮 | 200 – 300 µm | **60.0×** | 15.0 mm |
| buccal 上皮 | 500 – 800 µm | **60.0×** | 39.0 mm |
| sublingual 上皮 | 100 – 200 µm | **60.0×** | 9.0 mm |
| 唾液腺腺泡 | 40.0 µm | **300.0×**（第二个声明倍数） | 12.0 mm |

- 味蕾：a circumvallate papilla carries on the order of 250 taste buds; 24 are DRAWN - a declared subset, not a count

## 4. 七套管道系统

| 系统 | 段 | 分叉 | 终末 | 节点 | 根 | 直径 mm | 总长 mm | Murray 最大相对残差 | 核验的分叉数 | 轮廓锚定 |
|---|---|---|---|---|---|---|---|---|---|---|
| `lingual_artery` | 481 | 192 | 194 | 483 | 2 | 0.5223 – 2.4405 | 1636.8 | **0.0** | 192 | root only (carotid surface) |
| `facial_artery` | 326 | 102 | 104 | 328 | 2 | 0.7074 – 2.6406 | 1105.1 | **0.0** | 102 | root only (carotid surface) + crosses the measured mandibular antegonial notch |
| `inferior_alveolar_artery` | 117 | 25 | 27 | 119 | 2 | 0.5834 – 1.5 | 401.9 | **0.0** | 25 | both ends on the measured mandibular surface |
| `venous_return` | 812 | 259 | 261 | 814 | 2 | 0.8751 – 4.4667 | 3312.6 | **0.0** | 257 | root only (brachiocephalic vein surface) |
| `parotid_duct` | 492 | 180 | 182 | 494 | 2 | 0.4446 – 2 | 1319.2 | **0.0** | 180 | ostium only (from the measured maxillary arch) |
| `submandibular_duct` | 307 | 118 | 120 | 309 | 2 | 0.5459 – 2.1491 | 721.2 | **0.0** | 118 | ostium only (from the measured mandibular arch) |
| `sublingual_ducts` | 103 | 35 | 45 | 113 | 10 | 0.4386 – 1.0313 | 209.6 | **0.0** | 35 | False |
| **合计** | 2638 | 911 | 933 | 2660 | — | — | 8706.4 | **0.0** | — | — |

### 主干口径标定

| 系统 | 目标直径 | 计划/实际根下游终末数 | r_term mm | 实现主干直径 | 已发表范围 | 落入 |
|---|---|---|---|---|---|---|
| `lingual_artery` | 2.4 | 97 / 102 | 0.26117 | 2.4405 | 1.8 – 3.2 | 是 |
| `facial_artery` | 2.8 | 62 / 52 | 0.35372 | 2.6406 | 2 – 3.6 | 是 |
| `inferior_alveolar_artery` | 1.5 | 17 / 17 | 0.29168 | 1.5 | 1 – 2 | 是 |
| `venous_return` | 4.5 | 136 / 133 | 0.43753 | 4.4667 | 3 – 6 | 是 |
| `parotid_duct` | 2.0 | 91 / 91 | 0.22232 | 2.0 | 1 – 3 | 是 |
| `submandibular_duct` | 2.25 | 70 / 61 | 0.27297 | 2.1491 | 1.5 – 3 | 是 |
| `sublingual_ducts` | 0.75 | 5 / 13 | 0.2193 | 1.0313 | 0.5 – 1 | **否** |

> r_term = r_trunk_target / N_planned^(1/3), fixed BEFORE growth so the clearance solve runs on the radii that are rendered; the realised trunk diameter above differs from the target only where growth produced a different terminal count than planned. The trunk calibre is CALIBRATED to a published value; it is not an independent measurement. Murray's cube law stays exact at every bifurcation regardless.

### 跨树间隙

- **全局最小表面间隙 0.05 mm**，发生在 `lingual_artery|venous_return`，位置 ['lingual_artery', 'lingual_branch', 'venous_return', 'venous_tributary', [-38.1, 121.97, 146.99]]。
- 21 对组合里，表面间隙为负（即两管实际相交）的有 **0** 对。
- **恰好等于求解器目标 0.05 mm 的有 0 对**，已标 `at_solver_target`：。**这是求解器的属性，不是对几何的独立测量。**
- 节点位移上界声明为 9.0 mm，实测最大位移 9.0 mm。
- 丢弃的段数：**0**。没有为了达标删除任何段。唯一被丢弃的是 22 个落在舌下腺内的**静脉生长目标点**（不是段），原因写在 vessel_report.json 里。

> Anatomical fidelity outranks any fixed clearance threshold. Companion vessels and their ducts really do run with sub-millimetre wall separation, so a gap below 1 mm is reported as measured and is NOT treated as a failure. No segment was deleted and no node was pushed outside its own territory to make a number look better; the solver record above states exactly what moved and by how much.

> The venous tree is grown against a 0.05 mm clearance from the three arterial trees instead of the 0.25 mm used between unrelated systems, and its terminals are offset 2-4 mm from the arterial terminals. An artery and its companion vein share an adventitial sheath in vivo; forcing a 0.25 mm wall gap between them would be an artefact of the solver, not anatomy.

## 5. 流动：实测 v = Q / A

| 系统 | Q mL/min（每根） | v 最小 | v 最大 | 主干 v | 已发表参考 mm/s | 落入 | 显示增益 |
|---|---|---|---|---|---|---|---|
| `lingual_artery` | 25.0 | 19.063 | 95.411 | 89.069 | 200 – 700 | **否** | 无 |
| `facial_artery` | 45.0 | 36.693 | 136.955 | 136.955 | 250 – 800 | **否** | 无 |
| `inferior_alveolar_artery` | 8.0 | 29.344 | 107.474 | 75.451 | 100 – 400 | 是 | 无 |
| `venous_return` | 78.0 | 16.253 | 82.755 | 82.755 | 50 – 400 | 是 | 无 |
| `parotid_duct` | 0.425 | 0.501 | 2.255 | 2.255 | 0.5 – 40 | 是 | ×150.0 |
| `submandibular_duct` | 0.2975 | 0.347 | 1.398 | 1.367 | 0.5 – 40 | 是 | ×150.0 |
| `sublingual_ducts` | 0.0595 | 0.505 | 6.564 | 1.187 | 0.2 – 30 | 是 | ×150.0 |

- **逆流审计：112 条路径、3095 个采样点，逆流 0 个。** 判据：a sample refluxes if its baked v <= 0, or if the path arclength fails to increase between two consecutive samples; both are counted
- 唾液显示增益 **×150.0**，applied to the three duct systems ONLY, and only to the on-screen particle speed. Every velocity reported in this file is the UNGAINED v = Q/A.
- 心脏耦合：four-chamber time-varying elastance + three-element Windkessel (systemic, pulmonary) + compressed intramyocardial coronary branch，HR 64.45428703642494 bpm，CO 5.1771 L/min，MAP 94.4872 mmHg。aortic pressure of beat #2, rescaled to [0.72, 1.27] as a multiplicative modulation of arterial particle speed. The SHAPE is the heart model output; the rescaling range is a display choice and is declared.
- 数据来源声明：PhysioNet MIT-BIH Normal Sinus Rhythm Database (nsrdb) record 16265 (32 M, 22.24 h) — PUBLIC PLACEHOLDER, NOT THE USER'S OWN DATA；`is_user_personal_data = False`。

## 6. 运动学：声明值对浏览器内实测值

（实测于 desktop 视口的真实 Chromium，GPU 为数据中心 NVIDIA A10 / `ANGLE (NVIDIA, Vulkan 1.4.341 (NVIDIA NVIDIA A10 (0x00002236)), NVIDIA)`，**不是笔记本也不是真实手机**。）

| 项 | 声明 | 浏览器内实测 |
|---|---|---|
| 咀嚼周期 | 0.8 s（已发表 0.6 – 1 s） | 由相位推进实测 |
| 切点三维位移 mm（声明的咀嚼幅度就是这个量） | 18.0 | 0.000 – 16.338 |
| 切点垂直间距变化 mm（上述位移的垂直分量） | — | 0.000 – 14.319 |
| 切点间距 mm（静态基线 11.8751） | — | 11.875 – 26.194 |
| 髁突前滑 mm（前后向分量） | 3 | 0.000 – 2.999 |
| 髁突侧移 mm（研磨分量） | 4.5 | 0.000 – 4.475 |
| 下颌转角 deg | — | 0.000 – 11.986 |
| **髁突–颏距离漂移 mm**（旋转下必须为 0） | 0 | 0.000011 – 0.000011 |
| 髁间距 mm（刚体下必须恒定） | — | 99.6212 – 99.6212 |
| 吞咽波速 mm/s | 101.34 | 101.344 – 101.344（87 个采样） |
| 吞咽中位移的舌顶点 / 总顶点 | — | 4364 / 22778 |
| 吞咽波峰位移 mm | 6.0 | 5.9929 |
| **静止时未位移的采样环 / 总环** | 全部 | 123 / 123 |
| 静止时位移的舌顶点 | 0 | 0 |
| 静止时开口位移 mm | 0 | 0 |
| 静止时髁突位移 mm | 0 | 0 |
| 吞咽中位移超过波峰 5% 的环 / 总环 | 只有行进波的接触带 | 55 / 123，u ∈ [0.425, 0.8821] |
| 吞咽中位移 > 1e-6 mm 的环（含高斯远尾） | — | 102 / 123，u ∈ [0.1536, 0.9964] |
| 吞咽 t=0 时位移的环（波应在舌尖） | 仅 u≈0 附近 | 69 / 123，u ∈ [0.0036, 0.5036] |

**这组数字为什么能证明"是关节在动、不是整体缩放"：**髁突到颏点的距离在旋转下是不变量，在任何缩放下都必变。上表的漂移是从下颌**真实的世界矩阵**上量的，不是从姿态公式反推的；同时开口量在同一组帧里从 0 变到十几毫米。一个"整体缩放"的模型做不到前者为 0、后者非 0。舌的变形同样是直接写进**被渲染的** BufferAttribute 的，审计读回的就是正在画的顶点，所以静止时未位移的环数与吞咽时只有波峰附近位移这两件事，是在同一份几何上测到的。

## 7. 分子层

- CELLxGENE Census `2025-11-08`，`is_primary_data == True`。
- **74 个细胞类型 × 61 个 ADME 基因**，共 **312,586 个细胞**，**82 位供体**，**7 个数据集**。
- 分块读取：3 块 × 150,000 行，**失败 0 块**，耗时 9220.0 s。
- 组织词选择在 `tissue` 这一层而不是 `tissue_general`：tissue (NOT tissue_general - no oral term exists there)

| 组织词 | 细胞数 |
|---|---|
| `gingiva` | 105,918 |
| `parotid gland` | 40,487 |
| `oral cavity` | 28,186 |
| `minor salivary gland` | 24,010 |
| `posterior part of tongue` | 23,935 |
| `submandibular gland` | 21,230 |
| `salivary gland epithelium` | 15,796 |
| `anterior part of tongue` | 14,656 |
| `saliva` | 14,502 |
| `sublingual gland` | 11,549 |
| `hard palate` | 8,807 |
| `mucosa of lip` | 2,210 |
| `soft palate` | 835 |
| `buccal mucosa` | 302 |
| `tongue` | 163 |

**请求了但 Census 里不存在（或无 primary 细胞）的组织词：**`mucosa of dorsum of tongue`

- 咽部组织词 `oropharynx`, `nasopharynx`, `tonsil` 被**有意排除**：它们属于咽，不属于口腔。

### 逐基因实测值 — 按组织词（log1p CP10K，全细胞均值含零）

| 组织 | n | `CYP3A4` | `CYP3A5` | `CES1` | `CES2` | `UGT1A1` | `UGT1A3` | `UGT2B7` | `ABCB1` | `ABCG2` | `SLCO2B1` | `CYP2D6` | `SULT1A1` |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `gingiva` | 105,918 | 0.000 | 0.017 | 0.018 | 0.069 | 0.000 | 0.000 | 0.000 | 0.094 | 0.051 | 0.061 | 0.002 | 0.029 |
| `parotid gland` | 40,487 | 0.001 | 0.012 | 0.012 | 0.053 | 0.000 | 0.000 | 0.006 | 0.090 | 0.010 | 0.054 | 0.001 | 0.019 |
| `oral cavity` | 28,186 | 0.000 | 0.011 | 0.009 | 0.070 | 0.000 | 0.000 | 0.000 | 0.036 | 0.008 | 0.083 | 0.004 | 0.026 |
| `minor salivary gland` | 24,010 | 0.000 | 0.008 | 0.038 | 0.030 | 0.000 | 0.000 | 0.008 | 0.008 | 0.006 | 0.017 | 0.001 | 0.019 |
| `posterior part of tongue` | 23,935 | 0.006 | 0.346 | 0.014 | 0.418 | 0.012 | 0.001 | 0.009 | 0.084 | 0.022 | 0.051 | 0.001 | 0.012 |
| `submandibular gland` | 21,230 | 0.001 | 0.005 | 0.020 | 0.040 | 0.001 | 0.000 | 0.007 | 0.019 | 0.018 | 0.038 | 0.001 | 0.031 |
| `salivary gland epithelium` | 15,796 | 0.000 | 0.015 | 0.000 | 0.218 | 0.000 | 0.000 | 0.094 | 0.000 | 0.004 | 0.000 | 0.004 | 0.085 |
| `anterior part of tongue` | 14,656 | 0.013 | 0.305 | 0.037 | 0.628 | 0.019 | 0.002 | 0.001 | 0.050 | 0.018 | 0.070 | 0.001 | 0.010 |
| `saliva` | 14,502 | 0.002 | 0.064 | 0.028 | 0.032 | 0.059 | 0.000 | 0.001 | 0.000 | 0.007 | 0.028 | 0.002 | 0.126 |
| `sublingual gland` | 11,549 | 0.001 | 0.017 | 0.033 | 0.081 | 0.002 | 0.000 | 0.015 | 0.100 | 0.024 | 0.104 | 0.001 | 0.011 |
| `hard palate` | 8,807 | 0.000 | 0.016 | 0.032 | 0.080 | 0.000 | 0.000 | 0.000 | 0.122 | 0.053 | 0.031 | 0.002 | 0.045 |
| `mucosa of lip` | 2,210 | 0.000 | 0.026 | 0.083 | 0.083 | 0.001 | 0.000 | 0.022 | 0.105 | 0.054 | 0.031 | 0.003 | 0.039 |
| `soft palate` | 835 | 0.001 | 0.074 | 0.052 | 0.140 | 0.002 | 0.000 | 0.000 | 0.082 | 0.069 | 0.028 | 0.003 | 0.029 |
| `buccal mucosa` | 302 | 0.000 | 0.005 | 0.039 | 0.047 | 0.003 | 0.000 | 0.001 | 0.242 | 0.101 | 0.127 | 0.000 | 0.009 |
| `tongue` | 163 | 0.012 | 0.138 | 0.059 | 0.096 | 0.003 | 0.000 | 0.000 | 0.000 | 0.008 | 0.013 | 0.000 | 0.016 |

### 逐基因实测值 — 细胞数前 12 的细胞类型

| 细胞类型 | n | `CYP3A4` | `CYP3A5` | `CES1` | `CES2` | `UGT1A1` | `UGT1A3` | `UGT2B7` | `ABCB1` | `ABCG2` | `SLCO2B1` | `CYP2D6` | `SULT1A1` |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| fibroblast | 31,693 | 0.001 | 0.015 | 0.087 | 0.090 | 0.000 | 0.000 | 0.001 | 0.011 | 0.014 | 0.011 | 0.002 | 0.041 |
| basal cell | 25,616 | 0.004 | 0.168 | 0.011 | 0.341 | 0.015 | 0.001 | 0.014 | 0.008 | 0.006 | 0.016 | 0.002 | 0.035 |
| lymphocyte | 20,641 | 0.000 | 0.008 | 0.008 | 0.038 | 0.000 | 0.000 | 0.004 | 0.019 | 0.003 | 0.003 | 0.001 | 0.016 |
| T cell | 18,848 | 0.000 | 0.018 | 0.003 | 0.050 | 0.000 | 0.000 | 0.000 | 0.051 | 0.003 | 0.007 | 0.005 | 0.017 |
| endothelial cell of venule | 15,819 | 0.000 | 0.017 | 0.001 | 0.080 | 0.000 | 0.000 | 0.000 | 0.383 | 0.245 | 0.168 | 0.002 | 0.003 |
| epithelial cell | 15,285 | 0.000 | 0.013 | 0.048 | 0.037 | 0.000 | 0.000 | 0.011 | 0.009 | 0.005 | 0.001 | 0.001 | 0.021 |
| plasma cell | 13,034 | 0.000 | 0.012 | 0.001 | 0.025 | 0.001 | 0.000 | 0.000 | 0.009 | 0.007 | 0.007 | 0.002 | 0.015 |
| stratified squamous epithelial cell | 12,323 | 0.019 | 0.660 | 0.042 | 0.918 | 0.017 | 0.001 | 0.004 | 0.009 | 0.006 | 0.077 | 0.001 | 0.007 |
| acinar cell of salivary gland | 12,237 | 0.003 | 0.014 | 0.006 | 0.049 | 0.000 | 0.000 | 0.023 | 0.015 | 0.006 | 0.012 | 0.001 | 0.014 |
| B cell | 12,017 | 0.001 | 0.021 | 0.002 | 0.052 | 0.001 | 0.000 | 0.002 | 0.125 | 0.003 | 0.007 | 0.002 | 0.015 |
| myeloid cell | 9,770 | 0.000 | 0.007 | 0.015 | 0.062 | 0.000 | 0.000 | 0.000 | 0.028 | 0.006 | 0.262 | 0.002 | 0.047 |
| keratinocyte | 9,699 | 0.000 | 0.039 | 0.001 | 0.120 | 0.002 | 0.000 | 0.000 | 0.007 | 0.009 | 0.007 | 0.003 | 0.046 |

**实测里三条与通常叙述不一致、但照实报的结果：**

1. **CYP3A4 在口腔组织里几乎测不到**：全部 74 个细胞类型里最高只有 **0.019**，而 CYP3A5 最高 **0.660**，相差 35 倍。口腔黏膜以 CYP3A5 为主、而非 CYP3A4，这与文献一致，但也意味着把咪达唑仑（CYP3A4 探针）的代谢放到口腔层是没有数据支持的。
2. **ABCB1 (P-gp) 的信号主要不在上皮**：上皮类细胞里最高的是 myoepithelial cell（0.062），而内皮类里最高的是 endothelial cell（0.577）。也就是说，本数据**不支持**"口腔上皮顶侧 P-gp 外排"这条叙述；页面分子视图在把该转运体画到顶侧时标了 `extrapolated`，并且会把这条矛盾直接显示出来。
3. **`buccal mucosa` 只有 302 个 `is_primary_data == True` 的细胞**，颊黏膜的逐基因值因此统计上极不稳，不能当作"颊黏膜的测量"来读。颊黏膜恰恰是颊部给药的部位，这一点必须说清楚。

### 转运体画在哪个膜面

| 基因 | 膜面 | 依据等级 |
|---|---|---|
（完整表见页面分子视图与 `window.__oral3d.molecular.genesShown`）

- 顶侧（朝唾液）：`ABCB1`、`ABCG2`、`ABCC2`；基底侧（朝固有层与毛细血管）：`ABCC1`、`ABCC3`、`SLCO2B1`、`SLC22A3`；胞内（滑面内质网 / 胞质）：`CYP3A4`、`CYP3A5`、`CES1`、`CES2`、`UGT1A1`、`SULT1A1`。
- **分层鳞状口腔上皮里 `ABCB1`/`ABCG2`/`ABCC2` 的顶侧定位是从肠道、肾脏等上皮外推的，不是在口腔组织里测到的**，页面上逐个标了 `extrapolated`。画错膜面比不画更糟，所以这一条单列出来。
- 唾液腺分泌单位的膜面（腺泡基底侧 NKCC1、顶侧 AQP5；纹状管顶侧 ENaC/CFTR/ABCB1、基底侧 Na⁺/K⁺-ATPase）是已确立的生理学。

- there is no expression-to-flux, abundance-to-conductance, or concentration calibration anywhere in this build. The molecular layer is a qualitative map.

## 8. 渲染性能与移动端

| 视口 | unmasked renderer | 每帧渲染成本 ms（mean / p50 / p95 / max） | 对 16.67 ms 预算的余量（按 p95） | GPU 实测 | HUD 溢出 | DOM 溢出 |
|---|---|---|---|---|---|---|
| desktop [1440, 900] | `ANGLE (NVIDIA, Vulkan 1.4.341 (NVIDIA NVIDIA A10 (0x00002236)), NVIDIA)` | 1.064 / 1.100 / 1.400 / 2.000 | 11.907× | 是 | 0 | 0 |
| phone390 [390, 844] | `ANGLE (NVIDIA, Vulkan 1.4.341 (NVIDIA NVIDIA A10 (0x00002236)), NVIDIA)` | 2.343 / 1.400 / 2.600 / 36.900 | 6.412× | 是 | 0 | 0 |
| phone360 [360, 800] | `ANGLE (NVIDIA, Vulkan 1.4.341 (NVIDIA NVIDIA A10 (0x00002236)), NVIDIA)` | 2.474 / 1.500 / 2.700 / 37.600 | 6.174× | 是 | 0 | 0 |

| 视口 / 视图 | 每帧渲染成本 p95 ms | rAF fps（**不是可见帧率**） | HUD 溢出 |
|---|---|---|---|
| desktop / organ | 1.20 | 44 | 0 |
| desktop / vessels | 1.60 | 32 | 0 |
| desktop / motion | 1.20 | 60 | 0 |
| desktop / section | 0.60 | 45 | 0 |
| desktop / molecular | 2.30 | 24 | 0 |
| phone390 / organ | 1.10 | 60 | 0 |
| phone390 / vessels | 1.50 | 60 | 0 |
| phone390 / motion | 1.20 | 47 | 0 |
| phone390 / section | 0.50 | 60 | 0 |
| phone390 / molecular | 2.60 | 60 | 0 |
| phone360 / organ | 1.10 | 60 | 0 |
| phone360 / vessels | 1.40 | 60 | 0 |
| phone360 / motion | 1.20 | 46 | 0 |
| phone360 / section | 0.50 | 60 | 0 |
| phone360 / molecular | 2.70 | 60 | 0 |

- 60 fps 的帧预算是 **16.67 ms**。上表的"每帧渲染成本"是在渲染循环里直接计时的，不是从 rAF 间隔反推的。
- **无头 Chromium 没有垂直同步，rAF 帧率不设上限，那个数字不能读作用户可见帧率。**
- **全部测于数据中心 NVIDIA A10（型号取自上表的 unmasked renderer 字符串本身，不是手写的），不是笔记本，也不是真实手机。**手机那两行只是视口尺寸，硬件仍是 NVIDIA A10，所以这些数字不能用来预测真实手机的表现。
- 竖屏视口下相机会按 1/aspect 后退，但控制面板占掉下半屏，**口腔在竖屏里仍然偏上、边缘会被裁到**。HUD 与 DOM 的溢出计数都是 0（没有元素跑到视口外），但构图本身在手机上不如桌面。照实记在这里。

**审计未报出任何问题。**

## 9. 未达成 / 未测量 / 合成 / 超范围未调参

这一节是本文件存在的理由。

1. **跨器官受试者一致性被打破**：前四个器官用 s0236，本器官不能用（s0236 无头部覆盖，实测头颈掩膜全为空），改用另一位覆盖头部的受试者。
2. `bicondylar_breadth_mm` 实测 99，在已发表范围 112 – 140 之外 —— 照实报，未调参。
3. `tongue_volume_mL` 实测回读 119.4，在已发表范围 60 – 110 之外 —— 照实报，未调参。
4. `parotid_volume_mL` 恰好等于已发表均值，因为网格是按该均值缩放的；这是构造的属性，不是对任何东西的独立测量。
5. `oral_mucosa_area_cm2` 实测回读 232.3，在已发表范围 170 – 215 之外 —— 照实报，未调参。
6. 2 颗牙的位置不在实测牙槽弓上，而在该弓切线的外推段上。
7. 丝状乳头真实数量在 10^5 量级，页面只画 900 个 —— 声明的渲染抽样，不是解剖计数。
8. 口腔黏膜总面积里，舌背占比取了**声明的 42%**，不是在独立的舌背面片上测的。
9. 一个轮廓乳头约有 250 个味蕾，页面只画 24 个 —— 声明的子集，不是计数。
10. 微观视图用了**两个**声明放大倍数：黏膜与乳头 60.0×，唾液腺分泌单位 300.0×。
11. `sublingual_ducts` 实现主干直径 1.0313 mm，在已发表 0.5 – 1 mm 之外，原因是生长出来的根下游终末数（13）与计划值（5）不同 —— 照实报，未回头调。
12. 七套系统的主干口径都是**按已发表口径标定**的（r_term = r_trunk_target / N^(1/3)），不是独立测量结果。Murray 立方律在每个分叉上仍然精确成立。
13. 22 个静脉生长目标点被从舌下腺内剔除（2.5 mL 腺体内已有 10 条 Rivinus 管），段数未受影响，原因已记录。
14. `lingual_artery` 峰值 v = 95.411 mm/s 落在已发表范围 200 – 700 mm/s 之外 —— 照实报，未调参。（MEAN cross-sectional velocity Q/A, time-averaged. The published Doppler bounds quoted below are PEAK SYSTOLIC velocities…）
15. `facial_artery` 峰值 v = 136.955 mm/s 落在已发表范围 250 – 800 mm/s 之外 —— 照实报，未调参。（MEAN cross-sectional velocity Q/A, time-averaged. The published Doppler bounds quoted below are PEAK SYSTOLIC velocities…）
16. 唾液在导管里的实测速度为 0.3–6.6 mm/s 量级，页面用了声明的 ×150.0 显示增益才看得见方向；报表里的速度一律是**未加增益**的 v = Q/A。
17. 没有任何管道的峰值流速是独立测量；Q 本身是已发表量级估计，v 只是 Q/A 的商。已发表的 Doppler 参考是**收缩期峰值**，本页算的是**平均**横截面速度，两者本就不同量，均原样列出。
18. 以下口腔组织词在 Census 2025-11-08 里不存在或没有 `is_primary_data == True` 的细胞，照实列出、未用相邻组织顶替：mucosa of dorsum of tongue
19. 组织词 `buccal mucosa` 在 Census 里 `is_primary_data == True` 的细胞极少（见上表），颊黏膜的逐基因值主要来自其它口腔组织词，不能当作"颊黏膜的测量"来读。
20. CYP3A4 在口腔组织里几乎测不到（最高 0.019），CYP3A5 高 35 倍；把 CYP3A4 探针代谢放在口腔层没有本数据支持 —— 照实报。
21. ABCB1 (P-gp) 实测最高的是 endothelial cell（0.577），上皮类最高仅 myoepithelial cell（0.062）：**本数据不支持"口腔上皮顶侧 P-gp 外排"**，页面已把这条矛盾显示出来而不是让图赢。
22. `buccal mucosa` 在 Census 里只有 302 个 primary 细胞，颊黏膜的逐基因值统计上极不稳 —— 照实报。
23. 口腔上皮里 ABCB1 / ABCG2 / ABCC2 的顶侧膜定位是**外推**的，不是在口腔组织中测得；已在页面逐基因标 `extrapolated`。
24. 源数据里**没有**下颌骨、上颌骨、舌、牙、唾液腺、口腔、腭、口腔黏膜中的任何一个标签：本页的全部软组织都是**构造**的，`contour_anchored = false`。
25. 下颌/中面部的划分来自一条实测规则（咬合平面处的最小前部截面），不是解剖标注；`mandible_mL` / `midface_mL` 因此是**这条规则的产物**，不是对下颌骨体积的测量。
26. 牙槽弓截断窗口 42 mm 是**声明参数**，不是本受试者的测量值。
27. 下颌孔、颏孔的位置是"放在实测骨面上"，**不是在骨面上检出的**——源数据没有这两个孔的标注。
28. 唾液腺、导管、牙、黏膜分区、唇颊、咀嚼肌与舌肌的形状全部是参数化的，只有尺度取自已发表值。
29. `permeability_rank` 只是**排序**（1 最低…8 最高），锚在已发表的"舌下 > 颊 > 腭"次序上；**没有任何渗透系数、分配系数或通量标定**。
30. 第一过效应一节里的生物利用度是**已发表临床值的引用**，本项目没有、也不会用几何或表达量算出任何生物利用度。
31. 竖屏视口的构图：相机按 1/aspect 后退后口腔仍偏上、边缘被裁，HUD/DOM 溢出计数为 0 但构图不如桌面 —— 已知限制，照实记录。
32. 浏览器审计跑在数据中心 NVIDIA A10 上（型号读自 renderer 字符串）；**没有在任何真实手机上测过**，手机视口只是视口尺寸。
33. 没有做任何假设检验或残差诊断（未评估）。

---

生成于构建产物：`geometry_report.json`、`oral_soft_report.json`、`vessel_report.json`、`flow_report.json`、`oral_fields_report.json`、`oral_model.json`、`mobile_report.json`、`oral_cavity_census_provenance.json`。