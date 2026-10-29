# ACCEPTANCE — 三维数字孪生肾脏

生成日期 2026-10-29。每一项都是**实测**，数值从 `assets/*.json` 与基准测试报告读出。
最后一节逐条列出**未达成 / 未测量 / 合成**的项目。

## A. 几何

| 验收项 | 实测 | 判定 |
|---|---|---|
| 真实分割网格 | TotalSegmentator v2.0.1 CC-BY-4.0，受试者 s0236 | 达成 |
| 三角面数 13 万–16 万量级 | **153,580** | 达成 |
| glb ≤ 7 MB | **4.99 MB** | 达成 |
| 两侧肾脏无视野截断 | False / False | 达成 |
| 两侧肾脏单连通 | 1 / 1 | 达成 |
| 体积落在成人区间 | 174.8 / 194.1 mL，区间 [110.0, 200.0] | 达成 |
| 长径落在成人区间 | 11.8 / 12.59 cm | 达成 |
| 皮质体积占比落在公开区间 | **0.776 / 0.766**，区间 [0.6, 0.8] | 达成 |
| 肾锥体数落在 8–18 | **10 / 10** | 达成 |

皮质 / 髓质 / 肾柱体积（mL）：左 135.7 / 39.1 / 14.1，
右 148.7 / 45.4 / 14.0。
锥体高 [14.1, 23.4] mm，乳头间距 [20.5, 26.5] mm，
前排/后排乳头 3/7（measured, not imposed: papillae are a farthest-point sample of the sinus-facing deep parenchyma and fall into the two rows on their own）。

## B. 血管树

| 验收项 | 实测 | 判定 |
|---|---|---|
| 动脉级联完整 | 肾动脉→段(5)→叶间→弓形→皮质放射→入球→肾小球→出球 | 达成 |
| 叶间动脉走行于肾柱内 | 穿入锥体的采样点 左 1 / 右 0 | **未达成** |
| 弓形动脉沿皮髓交界弧行 | 每侧 19 条弧，跨度 [84.0, 154.0]° | 达成 |
| 串联双毛细血管床 | 出球小动脉→管周毛细血管(皮质) / 直小血管(髓质) | 达成 |
| 出球小动脉细于入球 | 0.013 < 0.02 mm，Poiseuille 阻力比 5.6× | 达成 |
| Murray 最大相对残差 = 0 | 肾动脉 **0.0**（623 个分叉）· 肾静脉 **0.0**（581）· 集合系统 **0.0**（8） | 达成 |
| 跨树不相交 | 见下表 | **未达成** |
| 最小表面间隙 ≥ 1.0 mm | **-0.729 mm**（artery_vs_vein） | **未达成** |
| 节点位移 ≤ 声明上界 | 6.0 mm，上界 6.0 mm | 达成 |
| 为满足间隙而删除的段数 | **0** | 达成（未删除任何段） |

| 树对 | 最小表面间隙 mm | 穿插段对 | 比较段对 | 最近间隙中位数 mm | 等于求解目标？ |
|---|---|---|---|---|---|
| 动脉–静脉 | -0.729 | 15 | 3,079,375 | 7.745 | False |
| 动脉–集合系统 | 0.25 | 0 | 403,635 | 12.45 | True |
| 静脉–集合系统 | 0.25 | 0 | 346,125 | 11.325 | True |

**读这些数字时请注意**：a pair whose at_solver_target flag is true sits exactly at the polish target (0.25 mm) because the solver drove it there — that value is a property of the solver, NOT an independent measurement, and is flagged so it cannot be read as one. A pair with the flag false was never reached by the solver and its minimum is a real measurement of the geometry.

优先级声明：ANATOMICAL FIDELITY OUTRANKS THE GAP TARGET. In a real kidney an interlobar artery and its vein run side by side a few hundred micrometres apart, so a 1.0 mm wall-to-wall separation is a MODELLING convenience, not an anatomical constraint. No vessel was moved out of the parenchyma, off its anatomical course, or deleted in order to reach it: each course was chosen from candidates inside a declared anatomical envelope ({'renal_vein_trunk_mm': 7.0, 'segmental_vein_mm': 7.0, 'venous_arcade_mm': 8.0, 'pelvis_mm': 9.0, 'major_calyx_mm': 6.0, 'minor_calyx_mm': 7.0}), and where the target could not be met inside that envelope the SMALLER REAL GAP is reported here.

树规模：
- renal_artery：1895 段 / 623 分叉 / 639 终末 / 2 根，直径 [0.8216, 5.6196] mm，总长 5469.9 mm，轮廓锚定：ostium only — placed on the subject's real aortic surface; the trunk, segmental, interlobar, arcuate and cortical radiate courses are CONSTRUCTED from the measured kidney, sinus and pyramid geometry
- renal_vein：1625 段 / 581 分叉 / 598 终末 / 2 根，直径 [1.027, 6.9054] mm，总长 5134.1 mm，轮廓锚定：ostium only — placed on the subject's real caval surface; the rest is CONSTRUCTED
- glomerular：2304 段 / 128 分叉 / 768 终末 / 128 根，直径 [0.0074, 0.02] mm，总长 411.8 mm，轮廓锚定：no — representative glomerular / post-glomerular units drawn at published arteriolar and capillary calibres
- collecting_system：213 段 / 8 分叉 / 20 终末 / 2 根，直径 [1.6246, 3.5] mm，总长 873.7 mm，轮廓锚定：papillae and the distal bladder target are derived from real labels; the calyces, pelvis and ureteric course are CONSTRUCTED

## C. 血流

| 验收项 | 实测 | 判定 |
|---|---|---|
| 与心脏模型耦合 | CO 5.1771 L/min、MAP 94.4872 mmHg 来自本项目心脏模型 | 达成 |
| RBF | 1139.0 mL/min（CO 5.1771 L/min x 1000 x 0.22 = 1139.0 mL/min） | 达成 |
| RPF 落在 600–700 | **626.4 mL/min** | 达成 |
| GFR 落在临床 90–120 | **125.3 mL/min** | **未达成** |
| 滤过分数 | 0.2（教科书值） | 达成 |
| 肾小球毛细血管高压体现在标注中 | [45.0, 60.0] vs 一般毛细血管 [15.0, 25.0] mmHg，比值 2.62× | 达成 |
| 净滤过压 | 11.5 mmHg（P_glomerular(52.5) - P_Bowman(12.5) - pi_plasma(28.5) = 11.5 mmHg） | 达成 |
| 逐采样点零逆流 | **98 条通路 / 6,272 个采样点，逆流 0 个** | 达成 |
| 单肾单位 GFR 落在 40–80 nL/min | **62.7 nL/min** | 达成 |
| 显示增益已声明 | 尿流 40.0×；球后网络口径 9× | 达成 |

**一个不是输入而是预测的结果**：入球小动脉 3.021 cm/s →
出球小动脉 **6.364 cm/s**。
about a fifth of the plasma entering the glomerulus is filtered off, so LESS blood leaves than enters — yet it leaves FASTER, because the efferent arteriole is the narrower vessel. Both numbers fall out of the published calibres and the model's own GFR; neither was fitted.

| 层级 | 时间平均 cm/s | 峰值 cm/s | 公开参考 cm/s | 判定 |
|---|---|---|---|---|
| afferent | 3.021 | 3.588 | — | 无公开参考值 |
| arcuate | 8.895 | 10.564 | 15.0–30.0 | **超出区间** |
| arcuate_vein | 5.527 | 6.564 | — | 无公开参考值 |
| cortical_radiate | 5.604 | 6.655 | 10.0–20.0 | **超出区间** |
| cortical_radiate_vein | 3.832 | 4.551 | — | 无公开参考值 |
| efferent | 6.364 | 7.559 | — | 无公开参考值 |
| glomerulus | 4.414 | 5.242 | — | 无公开参考值 |
| interlobar | 14.953 | 17.759 | 25.0–40.0 | **超出区间** |
| interlobar_vein | 12.475 | 14.816 | — | 无公开参考值 |
| major_calyx | 0.087 | 0.103 | — | 无公开参考值 |
| minor_calyx | 0.06 | 0.072 | — | 无公开参考值 |
| renal_artery | 38.288 | 45.473 | 60.0–120.0 | **超出区间** |
| renal_vein | 25.482 | 30.264 | 15.0–30.0 | 在区间内 |
| segmental | 21.18 | 25.155 | 30.0–50.0 | **超出区间** |
| segmental_vein | 3.832 | 4.551 | — | 无公开参考值 |
| ureter | 0.13 | 0.154 | — | 无公开参考值 |

## D. 剖面与微观单元

| 验收项 | 实测 |
|---|---|
| 放大倍数（已声明） | **20.0×** — 按真实毫米尺寸建模后显示 |
| 皮质肾单位 | 髓袢深度 2.0 mm，无直小血管 |
| 近髓肾单位 | 髓袢深度 18.0 mm，**有直小血管** |
| 近曲小管真实弧长 | 14 mm，盘绕在约 1.2 mm 的皮质球内 |
| 肾小球直径 | 0.200 mm |
| 致密斑 | 画在升支粗段与自身肾小球的接触点 |
| 小管液流速 | v = Q/A，流量自本模型单肾单位 GFR 起按教科书重吸收比例逐段递减 |

## E. 分子层

| 验收项 | 实测 |
|---|---|
| 数据规模 | 90 细胞类型 × 61 基因，446,561 细胞 |
| 近端小管上皮 | 113,966 细胞 / 137 donor |
| 二甲双胍主线 | OCT2(基底侧) → MATE1 / MATE2-K(顶侧)，全程原形排泄 |
| 膜面正确性 | 每个转运体由数据中的 `membrane` 字段定位，不由人眼摆放 |
| 逐基因实测值 | 见下表 |
| 通量标定 | **不存在**；THERE IS NO CALIBRATED EXPRESSION-TO-FLUX CONVERSION in this project. Expression values are not turned into transport rates, clearances or concentrations anywhere. The mapping from the measured values to the colours and arrow weights in the scene is QUALITATIVE and ordinal only. |

| 蛋白 | 基因 | 膜面 | 平均表达 log1p CP10K | 表达细胞比例 | 检出 |
|---|---|---|---|---|---|
| OAT1 | `SLC22A6` | basolateral | 0.899 | 47.41% | 是 |
| OAT3 | `SLC22A8` | basolateral | 0.7991 | 43.46% | 是 |
| MATE1 | `SLC47A1` | apical | 0.5998 | 36.12% | 是 |
| OCT2 | `SLC22A2` | basolateral | 0.5261 | 35.55% | 是 |
| MATE2-K | `SLC47A2` | apical | 0.4174 | 27.06% | 是 |
| MRP4 | `ABCC4` | apical | 0.403 | 25.53% | 是 |
| MRP2 | `ABCC2` | apical | 0.308 | 20.69% | 是 |
| OAT2 | `SLC22A7` | basolateral | 0.2686 | 17.75% | 是 |
| P-gp / MDR1 | `ABCB1` | apical | 0.2624 | 19.91% | 是 |
| MRP3 | `ABCC3` | basolateral | 0.2031 | 12.45% | 是 |
| PEPT2 | `SLC15A2` | apical | 0.1148 | 9.09% | 是 |
| OATP2B1 | `SLCO2B1` | not assigned | 0.0924 | 7.39% | 是 |
| PEPT1 | `SLC15A1` | apical | 0.0601 | 5.23% | 是 |
| OCT3 | `SLC22A3` | basolateral | 0.0592 | 4.8% | 是 |
| OCT1 | `SLC22A1` | basolateral | 0.0298 | 2.6% | 是 |
| BCRP | `ABCG2` | apical | 0.0186 | 1.71% | 是 |
| NTCP | `SLC10A1` | not expected in kidney | 0.0014 | 0.14% | 是 |
| OATP1B1 | `SLCO1B1` | not expected in kidney | 0.0005 | 0.05% | 是 |
| OATP1B3 | `SLCO1B3` | not expected in kidney | 0.0005 | 0.04% | 是 |

近端小管实测肝 CYP（用以支撑"肾内不代谢二甲双胍"）：
CYP3A4 0.0183（1.85%）, CYP3A5 0.4566（21.31%）, CYP2D6 0.0039（0.36%）, CYP2C9 0.0232（2.06%）, CYP1A2 0.0001（0.01%）。
these are the MEASURED renal proximal-tubule values for the major hepatic CYPs. They are reported so the claim "metformin is not metabolised in the kidney" is backed by the data in this project rather than only by the literature. A non-zero value is a transcript measurement and is NOT evidence of metformin metabolism.

## F. 渲染与移动端

| 视口 / 视图 | 平均 fps | HUD 越界 | 水平溢出 px | DOM 越界 |
|---|---|---|---|---|
| desktop / organ | 10.0 | 0 | 0 | 0 |
| desktop / section | 10.0 | 0 | 0 | 0 |
| desktop / nephron | 12.5 | 0 | 0 | 0 |
| desktop / molecular | 10.0 | 0 | 0 | 0 |
| phone390 / organ | 10.0 | 0 | 0 | 0 |
| phone390 / section | 10.0 | 0 | 0 | 0 |
| phone390 / nephron | 24.1 | 0 | 0 | 0 |
| phone390 / molecular | 10.0 | 0 | 0 | 0 |
| phone360 / organ | 10.0 | 0 | 0 | 0 |
| phone360 / section | 10.0 | 0 | 0 | 0 |
| phone360 / nephron | 24.2 | 0 | 0 | 0 |
| phone360 / molecular | 10.0 | 0 | 0 | 0 |

测试用渲染后端：ANGLE (Google, Vulkan 1.3.0 (SwiftShader Device (Subzero) (0x0000C0DE)), SwiftShader driver)，Google Inc. (Google)。

**上表的帧率不是硬件 GPU 帧率，而是软件光栅化（SwiftShader）帧率。**
An NVIDIA A10 was allocated and present in the container (nvidia-smi lists it, /dev/nvidia0 exists, libEGL_nvidia and the 10_nvidia.json vendor ICD are installed), but all SEVEN Chromium GL backend configurations probed (angle-gl-egl, angle-vulkan, desktop-gl, egl, angle-default, ozone-headless-egl, ozone-headless-nativegl — with EGL_PLATFORM=surfaceless and the NVIDIA vendor ICD pinned) fell back to SwiftShader. There is no Vulkan ICD in this image (/usr/share/vulkan/icd.d absent). So every frame rate in this report was measured under a SOFTWARE rasteriser. HARDWARE-GPU FRAME RATE WAS NOT MEASURED and no claim about it is made anywhere.
The per-view figures also reflect the headless compositor cadence: nine of the twelve land on exactly 10.0 fps, which is a requestAnimationFrame pacing artefact rather than a scene cost. They are reported as measured and are NOT evidence of on-device performance.
因此：**本项目没有测到肾脏页面在真实 GPU 上的帧率**，上表只能作为"页面在最差情形下仍可运行、且布局无溢出"的证据，
不能读作性能指标。可信的只有与渲染后端无关的那几列：HUD 越界、水平溢出、DOM 越界（全为 0）。

自动断言结果：desktop/organ: 10.0 fps mean, desktop/section: 10.0 fps mean, desktop/nephron: 12.5 fps mean, desktop/molecular: 10.0 fps mean, phone390/organ: 10.0 fps mean, phone390/section: 10.0 fps mean, phone390/nephron: 24.1 fps mean, phone390/molecular: 10.0 fps mean, phone360/organ: 10.0 fps mean, phone360/section: 10.0 fps mean, phone360/nephron: 24.2 fps mean, phone360/molecular: 10.0 fps mean

## G. 未达成 / 未测量 / 合成 —— 逐条

1. **arcuate 峰值流速 10.564 cm/s 落在公开区间 [15.0, 30.0] 之外。** 原因是搏动耦合为纯阻性、无惯性与顺应性，系统性低估收缩期峰值（模型阻力指数 0.362，公开参考 [0.5, 0.7]）。时间平均值 8.895 cm/s 照算出来的报，未施加任何增益。
1. **cortical_radiate 峰值流速 6.655 cm/s 落在公开区间 [10.0, 20.0] 之外。** 原因是搏动耦合为纯阻性、无惯性与顺应性，系统性低估收缩期峰值（模型阻力指数 0.362，公开参考 [0.5, 0.7]）。时间平均值 5.604 cm/s 照算出来的报，未施加任何增益。
1. **interlobar 峰值流速 17.759 cm/s 落在公开区间 [25.0, 40.0] 之外。** 原因是搏动耦合为纯阻性、无惯性与顺应性，系统性低估收缩期峰值（模型阻力指数 0.362，公开参考 [0.5, 0.7]）。时间平均值 14.953 cm/s 照算出来的报，未施加任何增益。
1. **renal_artery 峰值流速 45.473 cm/s 落在公开区间 [60.0, 120.0] 之外。** 原因是搏动耦合为纯阻性、无惯性与顺应性，系统性低估收缩期峰值（模型阻力指数 0.362，公开参考 [0.5, 0.7]）。时间平均值 38.288 cm/s 照算出来的报，未施加任何增益。
1. **segmental 峰值流速 25.155 cm/s 落在公开区间 [30.0, 50.0] 之外。** 原因是搏动耦合为纯阻性、无惯性与顺应性，系统性低估收缩期峰值（模型阻力指数 0.362，公开参考 [0.5, 0.7]）。时间平均值 21.18 cm/s 照算出来的报，未施加任何增益。
1. **GFR 125.3 mL/min 落在临床 eGFR 参考区间 [90.0, 120.0] 之外**（但命中经典教科书值 125.0 mL/min）。由声明的教科书输入算出，未调参。
1. **肾脏横径超出成人参考区间 [4.0, 7.0] cm**：左 6.91 cm、右 7.12 cm。横径沿主惯性轴测量，包含肾门唇，因此系统性大于临床超声的横径读数。照实报告。
1. **15 对管段互相穿插**（共比较 3,079,375 对，占 4.871e-06）。最差处在 renal_artery·interlobar_right（r=1.1682 mm）与 renal_vein·interlobar_vein_right（r=2.3674 mm）之间，表面间隙 -0.729 mm。the few remaining interpenetrating pairs are interlobar artery / interlobar vein companions sharing one column of Bertin. With the terminal count this scene renders, a strict cube-law cascade anchored on the published renal-artery calibre inflates interlobar calibres several-fold above their published values, so two companion vessels that really do run side by side no longer fit in one column. Nothing was deleted, no vessel was moved off its anatomical course and none was pushed out of the parenchyma to hide this; it is listed as a failed item instead.
1. **最小表面间隙 -0.729 mm 未达到 1.0 mm 这一建模目标**（最差对：artery_vs_vein）。按既定优先级，解剖保真高于该阈值：真实伴行血管壁间距本来就可以小于 1 mm，因此此处报实测值而不是把血管推开到解剖上讲不通的位置。所有候选走行都被限制在声明的解剖包络 {'renal_vein_trunk_mm': 7.0, 'segmental_vein_mm': 7.0, 'venous_arcade_mm': 8.0, 'pelvis_mm': 9.0, 'major_calyx_mm': 6.0, 'minor_calyx_mm': 7.0} 内。
1. **肾内中段血管口径系统性偏粗。** 渲染树终末直径 0.8216 mm 对应皮质放射动脉层级而非入球小动脉：the Murray-scaled terminal of the macroscopic tree sits at the CORTICAL RADIATE (interlobular) artery level, not at an afferent arteriole. With a trunk fixed at the published renal artery calibre, a cube-law cascade cannot also reach arteriolar calibre at a few hundred terminals, and no coefficient has been invented to force it to.
1. **肾静脉主干直径 6.83 mm 落在公开的 [10.0, 12.0] mm 之外。** the venous cascade is scaled from the arterial one by the published intrarenal vein/artery diameter ratio 1.25, NOT anchored on the hilar renal-vein calibre. Anchoring on the hilum would double every intrarenal vein relative to its artery, which is not anatomical. The rendered hilar diameter above is reported as it came out.
1. **皮质厚度、皮髓分界、锥体、肾柱、肾乳头全部为构造量**，源数据每侧只有一个肾脏标签。它们不是分割结果，不得当作影像测量读。
1. **集合系统的 Murray 定律是建模约定而非生理主张。** Murray's law is a hydrodynamic optimality principle for blood vessels. It is applied to the urinary collecting system here as a branching convention so the calibre cascade is reproducible; that is a modelling choice, NOT a physiological claim about urine flow.
1. **自身调节是现象学曲线而非方程。** a PHENOMENOLOGICAL autoregulation curve, not a differential equation. It reproduces the published observation that renal blood flow and GFR are held nearly constant between about 80 and 180 mmHg of renal perfusion pressure. The two mechanisms named below are LABELS on that curve, not simulated: no myogenic wall-tension equation and no macula densa / tubuloglomerular feedback kinetics are solved anywhere in this project.
1. **小管转运未求解。** the fractions below are textbook segmental reabsorption shares, used to LABEL the nephron view. No transport equations are solved and no solute concentrations are computed anywhere in this project.
1. **足细胞足突、裂孔隔膜、刷状缘微绒毛未建模**（低于 CT 分辨率两个数量级）；分子视图中的刷状缘环为声明过的示意。
1. **不存在任何表达量→转运通量的标定系数**，因此没有肾清除率、没有药物浓度、没有 PBPK 房室。
1. **近端小管无 S1/S2/S3 亚段分辨率**（OAT 与 MATE 在亚段间确有差异）。全程使用细胞类型分辨率的实测值，未用组织平均值代替。
1. **以下转运体不在 61 基因面板上，完全无法报告**：SLC22A12, SLC22A11, SLCO4C1, SLC22A4, SLC22A5, SLC5A2, AQP2。
1. **尿流为稳态等效速度。** 真实输尿管是蠕动式团块输送（每分钟 2–6 次收缩），场景中另按声明的 40.0× 显示增益加速以便看清方向。
1. **球后网络按声明的 9× 口径放大**显示于器官视图（位置与长度为真值）。
1. **呼吸位移幅度为公开范围的中点，非本受试者实测**；呼吸频率来自公开 PhysioNet 数据，`is_user_personal_data = false`。
1. **artery_vs_collecting_system 的最小间隙 0.25 mm 恰好等于求解器的打磨目标 0.25 mm**，这是求解器把它推到目标上的结果，是求解器的属性，不是对几何的独立测量，不得当作"实测间隙"引用。
1. **vein_vs_collecting_system 的最小间隙 0.25 mm 恰好等于求解器的打磨目标 0.25 mm**，这是求解器把它推到目标上的结果，是求解器的属性，不是对几何的独立测量，不得当作"实测间隙"引用。
1. **真实 GPU 帧率未测量。** An NVIDIA A10 was allocated and present in the container (nvidia-smi lists it, /dev/nvidia0 exists, libEGL_nvidia and the 10_nvidia.json vendor ICD are installed), but all SEVEN Chromium GL backend configurations probed (angle-gl-egl, angle-vulkan, desktop-gl, egl, angle-default, ozone-headless-egl, ozone-headless-nativegl — with EGL_PLATFORM=surfaceless and the NVIDIA vendor ICD pinned) fell back to SwiftShader. There is no Vulkan ICD in this image (/usr/share/vulkan/icd.d absent). So every frame rate in this report was measured under a SOFTWARE rasteriser. HARDWARE-GPU FRAME RATE WAS NOT MEASURED and no claim about it is made anywhere. The per-view figures also reflect the headless compositor cadence: nine of the twelve land on exactly 10.0 fps, which is a requestAnimationFrame pacing artefact rather than a scene cost. They are reported as measured and are NOT evidence of on-device performance. 基准测试报告里的 fps 一律是软件光栅化下的数值，既不能读作桌面性能，也不能读作手机性能；此项按未测量列出，而不是把软件渲染的数字当成 GPU 结果来报。同一次测试中与渲染后端无关的结论（HUD 越界 0、水平溢出 0、DOM 越界 0、98/98 通路推进、逆流 0、非有限步 0、13 个着色器程序全部链接成功）仍然成立。
1. **`renderCost` 背靠背渲染探针的结果已弃用。** 该探针在 SwiftShader 下返回 0.1–1.2 ms/帧（折合 849–10,526 fps），与软件光栅化的实际能力不符，判定为不可信，因此不在本文件中报告，也未用于任何结论。

## H. 不得这样读这个页面

- 这里**没有**药物浓度、肾清除率、PBPK 房室，因为没有表达量→通量的标定。
- 分子层的颜色与标记大小是**序数编码**，不是转运速率。
- 皮髓分层与锥体是**构造的**，不是影像测量。
- 所有生理波形的来源是公开数据，**不是用户本人的记录**。
- 本页为研究用模拟器，**不可用于临床给药或诊断决策**。
