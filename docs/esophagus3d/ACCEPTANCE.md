# ACCEPTANCE — 三维数字孪生食管

生成日期 2026-10-30。**本文件中的每个数字都从 `assets/*.json`、`out/*.json` 与 Modal 审计报告读出，不是手写的。**

逐项列出：达成的、**未达成的**、未测量的、合成的、以及超出已发表范围但**未调参**的。

## 0. 一句话结论

食管外表面、中心线、内切半径剖面、椎体层级映射与全部邻接距离是**实测**；血管树、管壁分层、肌层三分、LES 位置、栅栏区、蠕动波速与分子层定位是**构造**。受试者**不是** s0236（见第 1 节），中心线长度**超出已发表范围 +17.5%**（见第 2 节），四条血管主干口径**低于已发表直径**、七条树的终端流速**全部高于已发表小血管范围**（见第 5 节），膈食管裂孔处的第三处生理性狭窄**未测到**（见第 3 节）。

## 1. 未达成：跨器官受试者一致性

- **未达成。** abdominal field of view. The trachea, thyroid_gland and every cervical / upper-thoracic vertebra label is EMPTY, and the esophagus label holds only 4443 voxels (15.00 mL) covering the distal segment. It fails S2, S3 and S5 outright.

- **后果。** this page uses s0688. The oesophagus is therefore NOT from the same subject as the other four organ pages, and the cross-organ anatomical continuity claimed elsewhere in the project does not hold for the oesophagus-stomach boundary rendered here. The s0236 distal-oesophagus overlap is reported separately below for reference.

- 第二轮事后追加标准的原因：the round-1 winner s0692 has a 920.96 mL aorta (an ectatic/aneurysmal thoracoabdominal aorta) and an oesophageal tortuosity of 1.200 - the aneurysm displaces the oesophagus. The declared screen tested the oesophagus but never tested whether its NEIGHBOURS were normal. These two criteria were added AFTER seeing that result and are flagged as post hoc.

- 标准 **S5** 的更正：as written, S5 is anatomically impossible - the aorta begins at the arch (~T4) and can never span the CERVICAL oesophagus, so all 22 measured subjects failed it. → S5 restated as: the aortic arch apex must lie above the oesophageal mid-point (arch inside the FOV) AND the aorta must extend to within 5 mm of the oesophagus label's caudal end (hiatus reached).（after round 1, before any subject was selected）

- 标准 **S3** 的更正：S3 is one-sided: it rejects a subject whose trachea stops below the oesophagus top but not one whose oesophagus label is truncated far BELOW the larynx. s0556 passed S3 with the trachea top 118.5 mm ABOVE the oesophagus top, which means its cervical oesophagus is unlabelled. → NOT corrected - s0556 was not selected, so the outcome is unaffected. Recorded here rather than silently fixed.（noticed after round 2）

## 2. 几何：实测 vs 已发表

| 量 | 实测 | 已发表 | 判定 |
|---|---|---|---|
| 中心线长度 | 30.55 cm | 18–26 cm | **超范围 +17.5%，未调参** |
| 迂曲度 | 1.043 | ≤1.12 | 合格 |
| 内切直径 中位 | 10.79 mm | 13–25 mm（外径） | **低于范围** |
| 内切直径 最大 | 17.98 mm | 13–25 mm（外径） | 在范围内 |
| 椎体跨度 | C6→T12 | C6→T10/T11 | **远端偏低** |
| 隆突层级 | T6 | T4–T5 | **偏低一至二节** |
| 主动脉弓顶层级 | T4 | T4 | 合格 |
| 连通分量 | 1 | 1 | 合格 |
| 三角面总数 | 155479 | 13万–16万 | 合格 |
| glb 体积 | 3.75 MB | ≤7 MB | 合格 |

**面数构成（必须读到）：** the 130k-160k face gate is met by the eight-structure scene, not by the oesophagus alone. A 305 mm tube of 15 610 mm2 surface area supports about 14 000 native-resolution triangles; inflating it alone to 140 000 would be pure interpolation. All eight structures are real segmentations from the same subject and the per-structure counts are listed above.

**分辨率声明：** the oesophageal surface is marching cubes on the label TRILINEARLY UPSAMPLED x2. Mean triangle area 0.252 mm2 (edge ~0.76 mm) matches the upsampled grid pitch, but the INFORMATION CONTENT is that of the native 1.5 mm label - the extra triangles are interpolation, not detail. Every other structure is meshed at native resolution.

## 3. 三处生理性狭窄

| 狭窄 | 弧长 mm | 椎体 | 内切直径 mm | 突出度 mm | 至主动脉 mm | 至气管 mm | 状态 |
|---|---|---|---|---|---|---|---|
| constriction 1 cricopharyngeal | 12.0 | C7 | 7.435 | 0.759 | 67.99 | 3.67 | MEASURED |
| constriction 2 aortic arch | 92.1 | T4 | 8.483 | 1.388 | 5.62 | 6.05 | MEASURED |
| constriction 2b left main bronchus | 133.2 | T6 | 9.625 | 0.809 | 8.55 | 4.76 | MEASURED |
| minimum 4 retrocardiac | 191.3 | T8 | 7.824 | 3.144 | 6.58 | 39.79 | MEASURED, UNCLASSIFIED |
| minimum 5 distal | 270.4 | T11 | 14.001 | 0.468 | 14.15 | 112.16 | MEASURED, UNCLASSIFIED |
| constriction 3 diaphragmatic hiatus | 248.4 | T10 | 13.52 | 无局部极小值 | — | — | **NOT MEASURED** |

方法：local minimum of the inscribed-radius profile (Savitzky-Golay 21 mm / cubic, prominence >= 0.25 mm, >= 20 mm apart)

**未测到的第三处狭窄（必须读到）：** TotalSegmentator v2.0.1 has no diaphragm, crus or hiatus label, so the hiatus cannot be located from the segmentation. At the subject's own T10 body level the measured inscribed radius is 6.76 mm with NO local minimum. The hiatal constriction is therefore annotated at the published vertebral level using this subject's T10 centroid and is flagged contour_anchored = false.

**分辨为两处而非一处：** 经典教科书把主动脉弓与左主支气管压迹合并为单一"主动脉–支气管"狭窄；在这位受试者的半径剖面上它们是相距 41 mm 的两个独立极小值（T4 与 T6），照实分开报。

**一个未归类的极小值：** 弧长 191.3 mm、T8、内切直径 7.82 mm、突出度 3.14 mm —— **是整条剖面里最深的极小值**，位置相当于左心房水平，距降主动脉 6.6 mm、距气道 39.8 mm。它不对应任何经典狭窄，照实报为 MEASURED, UNCLASSIFIED，**不声称它是第三处狭窄**。

## 4. 动画：蠕动波（浏览器内实测）

| 量 | 声明值 | 已发表 | 判定 |
|---|---|---|---|
| 波速 | 35.3 mm/s | 20.0–40.0 mm/s | 在范围内 |
| 传输时间 | 8.65 s | 6.0–10.0 s | 在范围内 |
| 收缩幅度 | 50–150 mmHg | 同左（HRM） | 仅标注，未建模压力 |
| 收缩深度 | 0.55（局部内切半径的比例） | 真实原发波完全闭合管腔 | **声明的可视化值** |
| 过渡区减速因子 | 0.75 | HRM 有可重复的传播减速 | **声明值，非测量** |

**波速如何选定（必须读到）：** the declared speed is the midpoint of the INTERSECTION of the published wave-speed range (20-40 mm/s) with the range implied by this subject's measured length (305.5 mm) and the published 6-10 s transit window; both derived quantities therefore sit inside their own published ranges by construction, which is stated here rather than presented as an independent validation.

**深度声明：** DECLARED VISUALISATION VALUE. A real primary wave occludes the lumen completely; a 1.0 occlusion would collapse the rendered tube to a line, so 0.55 of the local inscribed radius is used.

### 浏览器内实测 — iphone14_390x844

| 核验项 | 值 |
|---|---|
| 采样环总数 | 154 |
| 波前在 s=0 mm 时**未位移**的环数 | 64 / 154 |
| 波前在中点时**未位移**的环数 | 14 / 154 |
| 声明波速 | 35.3 mm/s |
| 2.5 s 实测波速 | — mm/s |
| 4.0 s 实测波速 | 32.83 mm/s |
| 继发性蠕动：触发点 150 mm 时的门控起点 | front=196.5 mm, mode=secondary |

### 浏览器内实测 — android_360x800

| 核验项 | 值 |
|---|---|
| 采样环总数 | 154 |
| 波前在 s=0 mm 时**未位移**的环数 | 64 / 154 |
| 波前在中点时**未位移**的环数 | 14 / 154 |
| 声明波速 | 35.3 mm/s |
| 2.5 s 实测波速 | — mm/s |
| 4.0 s 实测波速 | 32.53 mm/s |
| 继发性蠕动：触发点 150 mm 时的门控起点 | front=195.3 mm, mode=secondary |

### 浏览器内实测 — desktop_1440x900

| 核验项 | 值 |
|---|---|
| 采样环总数 | 154 |
| 波前在 s=0 mm 时**未位移**的环数 | 64 / 154 |
| 波前在中点时**未位移**的环数 | 14 / 154 |
| 声明波速 | 35.3 mm/s |
| 2.5 s 实测波速 | — mm/s |
| 4.0 s 实测波速 | 32.11 mm/s |
| 继发性蠕动：触发点 150 mm 时的门控起点 | front=193.5 mm, mode=secondary |

**这是行进波而不是整体缩放的证据：** 在任一时刻只有波前附近的环被位移，其余绝大多数环的位移严格为 0；把波前放在 s=0 与放在中点，未位移的环是**不同的两批**。整体缩放会让每一个环同时位移，不可能出现"未位移环数 ≈ 总环数"。位移场由 CPU 按与顶点着色器完全相同的表达式重算，逐环比对。

## 5. 血管树

| 树 | 锚定 | 段 | 分叉 | 终端 | 直径 mm | 总长 mm | Murray 最大相对残差 |
|---|---|---|---|---|---|---|---|
| `V2_azygos_systemic` | **构造** | 393 | 96 | 97 | 0.46–2.11 | 1724.0 | **0.0** |
| `A3_aortic_oesophageal_branches` | 根锚定于真实标签 | 221 | 59 | 60 | 0.31–1.21 | 970.1 | **0.0** |
| `A2_bronchial_artery` | 根锚定于真实标签 | 172 | 48 | 49 | 0.3–1.1 | 720.3 | **0.0** |
| `A4_left_gastric_ascending` | 根锚定于真实标签 | 185 | 52 | 53 | 0.36–1.35 | 784.5 | **0.0** |
| `A1_inferior_thyroid_artery` | 根锚定于真实标签 | 139 | 43 | 44 | 0.35–1.24 | 556.1 | **0.0** |
| `V3_left_gastric_portal` | 根锚定于真实标签 | 149 | 42 | 43 | 0.5–1.75 | 610.8 | **0.0** |
| `V1_inferior_thyroid_vein_systemic` | 根锚定于真实标签 | 125 | 33 | 34 | 0.42–1.36 | 572.1 | **0.0** |
| **合计** | | **1384** | **373** | **380** | | **5937.9** | **0.0** |

Murray 残差在全部 **373** 个分叉上逐个核验，最大相对残差 **0.0**。这不是拟合出来的：每个分叉都是二分叉，半径按 r = r_terminal × n_downstream^(1/3) 生成，n_downstream(parent) 等于各子节点之和，所以 r_p³ = Σr_c³ 在浮点精度内恒等。

### 未达成：主干口径

| 树 | 主干直径 mm | 已发表 mm | 判定 |
|---|---|---|---|
| `V2_azygos_systemic` | 2.11 | 5.0–10.0 | **BELOW by 58%** |
| `A3_aortic_oesophageal_branches` | 1.21 | 0.8–2.0 | 在范围内 |
| `A2_bronchial_artery` | 1.1 | 0.5–2.0 | 在范围内 |
| `A4_left_gastric_ascending` | 1.35 | 1.0–2.0 | 在范围内 |
| `A1_inferior_thyroid_artery` | 1.24 | 1.5–3.0 | **BELOW by 17%** |
| `V3_left_gastric_portal` | 1.75 | 4.0–6.0 | **BELOW by 56%** |
| `V1_inferior_thyroid_vein_systemic` | 1.36 | 2.0–5.0 | **BELOW by 32%** |

**原因（未调参）：** with equal terminal flows and binary bifurcations, Murray ties the trunk and terminal radii together through n_downstream^(1/3). Spanning the real 20-40x radius ratio of a named trunk down to a 0.3 mm adventitial branch would need n = 8000-64000 terminals. With 34-97 terminals per tree the span is only 3.2-4.6x, so the terminal calibres were anchored to the published small-vessel range and FOUR of the seven trunk calibres consequently fall BELOW their published diameters (inferior thyroid artery -16%, inferior thyroid vein -28%, left gastric vein -54%, azygos -58%). No radius was adjusted to hide this.

### 未达成的终端数

请求的外膜目标点 **398** 个，未能接上的 **25** 个（6.3%）。逐树：{'V2_azygos_systemic': 0, 'A3_aortic_oesophageal_branches': 3, 'A2_bronchial_artery': 6, 'A4_left_gastric_ascending': 0, 'A1_inferior_thyroid_artery': 3, 'V3_left_gastric_portal': 6, 'V1_inferior_thyroid_vein_systemic': 7}

**these targets were NOT silently dropped. Each is an adventitial sampling point for which no bifurcation site existed inside the periesophageal corridor that also satisfied the 0.05 mm clearance against the already-grown trees. The per-tree counts and the rejection breakdown (rejected_outside vs rejected_clearance) are reported below; the terminal counts above are the attached ones.**

| 树 | candidates | rejected_outside | rejected_clearance | rejected_angle | unattached |
|---|---|---|---|---|---|
| `V2_azygos_systemic` | 8064 | 5760 | 0 | 359 | 0 |
| `A3_aortic_oesophageal_branches` | 5079 | 3393 | 1849 | 200 | 3 |
| `A2_bronchial_artery` | 3639 | 3153 | 949 | 307 | 6 |
| `A4_left_gastric_ascending` | 5548 | 1646 | 1048 | 236 | 0 |
| `A1_inferior_thyroid_artery` | 4095 | 1869 | 424 | 279 | 3 |
| `V3_left_gastric_portal` | 4890 | 1578 | 1637 | 176 | 6 |
| `V1_inferior_thyroid_vein_systemic` | 4104 | 1140 | 1905 | 226 | 7 |

## 6. 跨树间隙 —— 解剖保真优先于固定阈值

- 互相穿透的树对：**0** / 21

- **实测最小跨树表面间隙 = 0.0212 mm**，发生在 `A4_left_gastric_ascending` 与 `V3_left_gastric_portal` 之间，位置 scene [-0.8, -463.9, 174.7] mm —— 即**胃左动脉与胃左静脉这对伴行血管**。

- 求解器间隙目标 = 0.05 mm。

- **恰好等于求解器目标的值（必须标出）：**

  - `A1_inferior_thyroid_artery` ↔ `V1_inferior_thyroid_vein_systemic` = 0.0501 mm — `at_solver_target: true`，**这是求解器的属性，不是对几何的独立测量**
  - `V2_azygos_systemic` ↔ `A3_aortic_oesophageal_branches` = 0.0502 mm — `at_solver_target: true`，**这是求解器的属性，不是对几何的独立测量**

- 节点位移修复：移动了 **8** 个节点，**最大位移 1.9988 mm**，声明上界 **2.0 mm** —— **最大位移实际顶到了上界**，应读作『上界本身在起约束作用』，而不是『修复轻松收敛』。可移动的只有 free bifurcation and adventitial-terminal nodes only; roots and the azygos trunk are fixed。

- **anatomical fidelity was given priority over any fixed gap threshold. Companion artery and vein in a shared periesophageal adventitial sheath are separated by far less than 1 mm in life, so the solver target was set to 0.05 mm and the measured minimum (0.0212 mm, left gastric artery vs left gastric vein) is reported as found. Two of the 21 pairs sit exactly at the 0.05 mm solver target and are flagged at_solver_target: those two numbers are a property of the solver, not an independent measurement of the geometry.**

- **放弃了余弦弯曲（cosmetic bow）。** CCO 在分叉之间产生直段；早先各器官页用一个 12% 的横向半正弦弯曲让血管看起来更自然。这里实测下来，任何非零弯曲幅度（0.015/0.03/0.06/0.12）都会让相邻树互相穿透（实测 −0.16 ~ −0.46 mm）。**解剖不相交优先于外观曲率**，所以弯曲幅度设为 0，段在分叉之间是直的。

| 树对 | 最小表面间隙 mm | 半径 A | 半径 B |
|---|---|---|---|
| `A4_left_gastric_ascending` ↔ `V3_left_gastric_portal` | 0.0212 | 0.18 | 0.454 |
| `A1_inferior_thyroid_artery` ↔ `V1_inferior_thyroid_vein_systemic` | 0.0501 | 0.498 | 0.55 |
| `V2_azygos_systemic` ↔ `A3_aortic_oesophageal_branches` | 0.0502 | 0.332 | 0.414 |
| `V2_azygos_systemic` ↔ `A2_bronchial_artery` | 0.0895 | 0.23 | 0.15 |
| `V2_azygos_systemic` ↔ `A1_inferior_thyroid_artery` | 0.0988 | 1.053 | 0.175 |
| `A2_bronchial_artery` ↔ `A1_inferior_thyroid_artery` | 0.1072 | 0.15 | 0.175 |
| `V2_azygos_systemic` ↔ `A4_left_gastric_ascending` | 0.1134 | 0.23 | 0.18 |
| `A3_aortic_oesophageal_branches` ↔ `A4_left_gastric_ascending` | 0.1302 | 0.155 | 0.18 |
| `A3_aortic_oesophageal_branches` ↔ `A2_bronchial_artery` | 0.1422 | 0.155 | 0.15 |
| `V2_azygos_systemic` ↔ `V1_inferior_thyroid_vein_systemic` | 0.1455 | 0.23 | 0.21 |

## 7. 流动与零逆流审计

| 量 | 值 |
|---|---|
| 路径数 | 380 |
| 采样点总数 | 32866（步长 1.0 mm） |
| **v < 0 的采样点（逆流）** | **0** |
| v = 0 的采样点（停滞） | 0 |
| 每条路径弧长单调 | True |
| 声明方向 | arteries root->terminal, veins terminal->root |

| 树 | 路径 | 采样点 | 逆流采样点 |
|---|---|---|---|
| `V2_azygos_systemic` | 97 | 17689 | **0** |
| `A3_aortic_oesophageal_branches` | 60 | 2616 | **0** |
| `A2_bronchial_artery` | 49 | 1498 | **0** |
| `A4_left_gastric_ascending` | 53 | 3874 | **0** |
| `A1_inferior_thyroid_artery` | 44 | 1993 | **0** |
| `V3_left_gastric_portal` | 43 | 3563 | **0** |
| `V1_inferior_thyroid_vein_systemic` | 34 | 1633 | **0** |

### 峰值流速 vs 已发表

| 树 | 主干流速 cm/s | 已发表主干范围 | 判定 | 来源 |
|---|---|---|---|---|
| `V2_azygos_systemic` | 15.0 | **未检索到** | **NO PUBLISHED RANGE LOCATED** | NO published mean-velocity range located for the azygos vein; declared assumption |
| `A3_aortic_oesophageal_branches` | 21.1 | **未检索到** | **NO PUBLISHED RANGE LOCATED** | NO published Doppler range located for direct aortic oesophageal branches; declared assumption |
| `A2_bronchial_artery` | 21.2 | **未检索到** | **NO PUBLISHED RANGE LOCATED** | NO published Doppler range located for the bronchial artery; declared assumption |
| `A4_left_gastric_ascending` | 55.0 | 40.0–70.0 | 在范围内 | left gastric artery PSV, coeliac-axis Doppler |
| `A1_inferior_thyroid_artery` | 28.3 | 20.0–40.0 | 在范围内 | inferior thyroid artery PSV, thyroid Doppler ultrasound, euthyroid adults |
| `V3_left_gastric_portal` | 15.0 | 10.0–20.0 | 在范围内 | left gastric (coronary) vein, hepatopetal, portal-hypertension Doppler literature |
| `V1_inferior_thyroid_vein_systemic` | 12.0 | **未检索到** | **NO PUBLISHED RANGE LOCATED** | NO published range located for the inferior thyroid vein; declared assumption |

**七条树里只有三条的主干流速能对上已发表范围**（甲状腺下动脉、胃左动脉、胃左静脉）；其余四条（支气管动脉、主动脉食管支、甲状腺下静脉、奇静脉）**未检索到已发表的多普勒流速范围**，其主干流速是**声明的假设值**，在 `flow_report.json` 里逐条标注。

### 未达成：终端流速

| 树 | 终端流速 cm/s | 已发表小血管范围 | 判定 |
|---|---|---|---|
| `V2_azygos_systemic` | 3.265 | 0.5–1.5 cm/s | **ABOVE by 2.2x** |
| `A3_aortic_oesophageal_branches` | 6.386 | 0.5–1.5 cm/s | **ABOVE by 4.3x** |
| `A2_bronchial_artery` | 6.832 | 0.5–1.5 cm/s | **ABOVE by 4.6x** |
| `A4_left_gastric_ascending` | 14.642 | 0.5–1.5 cm/s | **ABOVE by 9.8x** |
| `A1_inferior_thyroid_artery` | 8.498 | 0.5–1.5 cm/s | **ABOVE by 5.7x** |
| `V3_left_gastric_portal` | 4.282 | 0.5–1.5 cm/s | **ABOVE by 2.9x** |
| `V1_inferior_thyroid_vein_systemic` | 3.704 | 0.5–1.5 cm/s | **ABOVE by 2.5x** |

**every terminal velocity exceeds the published 0.5-1.5 cm/s for a 0.3-0.6 mm arteriole/venule by 2.2-9.8x. Same structural cause as the calibre shortfall: a ~50-terminal Murray cascade cannot span the real velocity range. Reported as computed.**

### 显示增益（逐项声明）

- **arterial_and_venous_particle_speed** = 8x
- **declaration** = particle motion along the vessel paths is drawn at 8x the computed velocity; terminal velocities of 3-15 cm/s at a 0.3 mm calibre are invisible at this scene scale without gain. The gain is a display constant and is shown in the page HUD.
- **peristaltic_wave** = 1x (no gain; the wave runs at the declared 35.3 mm/s in real time)
- **reflux_demo** = 1x

### 搏动耦合

- 来源：heart_timeseries.json (this project's ORd-2011 + four-chamber elastance + Windkessel model); beat 9 of the simulation, resampled to 128 phase samples

- CO **5.1771 L/min**、MAP **94.4872 mmHg**、HR **64.45 bpm** —— 全部取自本项目心脏模型的输出，**没有另造数字**。

- **`is_user_personal_data = false`**。上游数据是 PhysioNet MIT-BIH Normal Sinus Rhythm Database (nsrdb) record 16265 (32 M, 22.24 h) — PUBLIC PLACEHOLDER, NOT THE USER'S OWN DATA。**这不是使用者本人的数据。**

- 只有四条**动脉**树接受搏动调制；奇静脉、甲状腺下静脉、胃左静脉与栅栏区通道的搏动权重为 0。

## 8. 管壁四层与微观单元

| 层 | 内容 | 厚度（已发表）mm |
|---|---|---|
| 黏膜 mucosa | non-keratinised stratified squamous epithelium + lamina propria + muscularis mucosae | 0.5–1.2 |
| 黏膜下层 submucosa | areolar tissue, oesophageal glands proper, submucosal venous plexus, Meissner plexus | 0.6–1.5 |
| 肌层 muscularis_propria | inner circular + outer longitudinal | 1.5–3.0 |
| 外膜 adventitia | loose areolar connective tissue, NO mesothelium | 0.3–1.0 |
| **合计** | | **2.9–6.7**（CT 参考 3.0–6.0）|

**layer thicknesses are PUBLISHED values (endoscopic ultrasound layer correlates and histology), not measured from this subject. A 1.5 mm isotropic CT label cannot resolve a 0.5 mm mucosa. The organ surface and the inscribed-radius profile are measured; everything inside the wall is published anatomy.**  `contour_anchored = false`

### 放大倍数（必须声明）

| 视图 | 真实尺寸 | 放大倍数 | 显示尺寸 |
|---|---|---|---|
| 四层管壁楔形 | 壁厚 4.80 mm | **60×** | ~288 mm |
| 复层鳞状上皮 | 上皮总厚 ~400 µm，单个鳞状细胞 10–20 µm | **600×** | 上皮 ~240 mm，单细胞 6–12 mm |

两级放大倍数都写在页面 HUD 上。微观单元的层厚与细胞尺寸是**已发表组织学**，不是从这位受试者测出来的（1.5 mm 各向同性体素分辨不出 0.85 mm 的黏膜层），唯一实测的量是显示站点的内切半径。

### 外膜而非浆膜（外科学关键差异）

ADVENTITIA, not serosa. The cervical and thoracic oesophagus has no serosal covering, so there is no mesothelial barrier to tumour spread and no serosa to hold sutures — the anatomical reason oesophageal anastomoses leak more than gastric or colonic ones. Only the short intra-abdominal segment has a peritoneal (serosal) covering on its anterior aspect.

## 9. 括约肌

| 括约肌 | 肌肉 | 静息压 mmHg | 弧长 mm | 层级 | 锚定 |
|---|---|---|---|---|---|
| 上食管括约肌 | cricopharyngeus (+ inferior pharyngeal constrictor, proximal oesophagus) | 30–120 | 12.0 | C6 classically; this subject's radius minimum is at C7 | 实测 |
| 下食管括约肌 | intrinsic smooth-muscle thickening + extrinsic crural diaphragm | 10–30 | 248.4 | T10-T11 classically (diaphragmatic hiatus); not resolvable in this subject | **构造（无膈标签）** |

- **上食管括约肌** 引文：high-resolution manometry, resting UES pressure is highly variable and markedly asymmetric (antero-posterior >> lateral); ranges of 30-120 mmHg are reported
  tonically closed by continuous motor discharge; relaxes 0.2-0.5 s before the bolus

- **下食管括约肌** 引文：high-resolution manometry, resting LES pressure 10-30 mmHg above gastric
  deglutitive relaxation lasts the whole peristaltic transit; TRANSIENT LES RELAXATION (tLESR) is vagally mediated, independent of swallowing, and is the dominant mechanism of physiological and pathological gastro-oesophageal reflux

### 一过性 LES 松弛（tLESR）与反流演示

- 性质：**phenomenological switch, NOT a mechanistic model**

- 持续时间 10–45 s，清醒时频率 2–6 次/小时

- 引文：transient LES relaxations, 10-45 s, vagally mediated, triggered by gastric fundic distension; the dominant mechanism of both physiological and pathological gastro-oesophageal reflux

- **the reflux animation shows a retrograde gastric-content front entering the distal oesophagus during a tLESR. There is no pressure solver behind it: the front position is a prescribed function of time. No pH, no volume and no clearance kinetics are modelled.**

## 10. 门-体静脉吻合与栅栏区

| 量 | 值 |
|---|---|
| 纵行通道数 | 16 |
| 弧长范围 | 280.4–305.5 mm（GOJ 上方 25 mm） |
| 正常口径（声明） | 0.35 mm |
| 已发表正常口径 | 0.2–0.5 mm |
| 曲张演示口径（声明） | 6.0 mm |
| 已发表曲张口径 | 3.0–10.0 mm |
| contour_anchored | false |

- **性质：** explicit anastomotic NETWORK, not a CCO tree — Murray's law does not apply to a non-dichotomous plexus and no Murray residual is reported for it

- **定位：** submucosal radius = measured inscribed radius - 0.9 mm (published mucosa + muscularis mucosae depth)

- **意义：** the palisade zone is the watershed where the PORTAL territory (left gastric vein, V3) meets the SYSTEMIC territory (azygos, V2). In portal hypertension flow in V3 reverses (hepatofugal) and decompresses through these longitudinal channels into the azygos, dilating them into oesophageal varices.

- **NORMAL state: the distal oesophageal submucosal blood drains CAUDALLY into the left gastric vein and thence to the portal vein; the palisade channels carry it downward. VARICES state: portal pressure exceeds the systemic venous pressure, left gastric flow reverses (hepatofugal) and the palisade channels carry blood CRANIALLY into the azygos. Both states are prescribed, not solved: there is no pressure-flow network solver behind the switch. The 8.0 mm/s palisade velocity is a DECLARED value - no published palisade-zone velocity was located.**

## 11. 分子层

| 量 | 值 |
|---|---|
| Census 版本 | 2025-11-08 |
| tissue_general | esophagus |
| is_primary_data | True |
| 细胞总数 | **178410** |
| 细胞类型数 | **92** |
| 供体数 | **29** |
| 数据集数 | **6** |
| 请求基因 / 解析到 | 61 / **61** |
| Census 中缺失的基因 | 无 |
| 区块数 × 区块大小 | 2 × 150000 |
| **失败区块数** | **0** |
| 墙钟时间 | 9772.4 s |
| 归一化 | CP10K (raw counts / obs.raw_sum * 1e4) then log1p; mean over all cells of the cell type (zeros included) |
| is_user_personal_data | False |

读取路线：**按 soma_joinid 连续区块**，每块 150 000 行。多器官并行读、按细胞类型分层抽样、按 dataset_id 分块三条路线在本项目早先已证明失败（破坏 TileDB 坐标局部性），不再尝试。

**沙箱网络的一个实测坑：** TileDB 的 S3 VFS **不读取环境变量里的 HTTP 代理**，因此在本沙箱里 Census 的 S3 桶直接无法解析。必须把代理显式传进 `SOMATileDBContext(tiledb_config={'vfs.s3.proxy_scheme','vfs.s3.proxy_host','vfs.s3.proxy_port'})`。这一点写在 `tools/census_esophagus.py` 里。

### 按细胞数排序的前 12 个细胞类型

| 细胞类型 | 细胞数 |
|---|---|
| stratified epithelial cell | 51307 |
| epithelial cell | 30884 |
| epithelial cell of esophagus | 19279 |
| fibroblast | 18350 |
| enteric smooth muscle cell | 9159 |
| basal cell | 7239 |
| smooth muscle cell | 4448 |
| endothelial cell of vascular tree | 3275 |
| macrophage | 2708 |
| enterocyte | 2681 |
| keratinocyte | 2277 |
| foveolar cell of stomach | 2076 |

### 8 个有定位证据的分子 —— 逐基因实测值

| 基因 | 膜面 | stratified epithelial cell | epithelial cell | epithelial cell of esophag | fibroblast | enteric smooth muscle cell | basal cell |
|---|---|---|---|---|---|---|---|
| `ABCB1` | 顶侧 apical | 0.0002 (0.0%) | 0.0005 (0.1%) | 0.0002 (0.0%) | 0.0318 (1.4%) | 0.0406 (1.9%) | 0.0023 (0.1%) |
| `ABCC1` | 基底侧 basolateral | **面板中无此基因** | **面板中无此基因** | **面板中无此基因** | **面板中无此基因** | **面板中无此基因** | **面板中无此基因** |
| `ABCC3` | 基底侧 basolateral | 0.0333 (4.1%) | 0.1478 (14.0%) | 0.0695 (8.6%) | 0.0147 (0.6%) | 0.0319 (1.5%) | 0.0704 (5.0%) |
| `SLC16A1` | 基底侧 basolateral | **面板中无此基因** | **面板中无此基因** | **面板中无此基因** | **面板中无此基因** | **面板中无此基因** | **面板中无此基因** |
| `GSTP1` | 胞质 | 3.5619 (98.7%) | 3.0808 (98.3%) | 3.4726 (97.9%) | 0.6679 (25.7%) | 0.0343 (1.5%) | 2.7054 (84.2%) |
| `ALDH3A1` | 胞质 | **面板中无此基因** | **面板中无此基因** | **面板中无此基因** | **面板中无此基因** | **面板中无此基因** | **面板中无此基因** |
| `ADH1C` | 胞质 | **面板中无此基因** | **面板中无此基因** | **面板中无此基因** | **面板中无此基因** | **面板中无此基因** | **面板中无此基因** |
| `CYP2C18` | 内质网膜 | **面板中无此基因** | **面板中无此基因** | **面板中无此基因** | **面板中无此基因** | **面板中无此基因** | **面板中无此基因** |

单元格格式：平均 log1p(CP10K)（表达该基因的细胞百分比）。

**不给膜面的基因：** 61 基因面板里其余 53 个在 CSV 源表中照常报实测值，但**不画在任何膜面上** —— 它们在复层鳞状上皮里的极性未经确立。**把转运体画错膜面比不画更糟**，所以它们被标为 UNASSIGNED。

### 未达成：鳞状上皮分层的分开报告

- **basal** — Census 中匹配到 2 个细胞类型，共 7251 个细胞：`basal cell` (7239)、`cycling basal cell` (12)
- **suprabasal_spinous** — **未达成**：Census 的 `esophagus` 组织词下没有能对上这一分层的 cell_type 标签（检索词 ['suprabasal', 'spinous']）。任务要求把基底、棘层、表层分开报，这一项只能按 Census 实际提供的细胞类型标签报，**不自行把细胞类型归并到分层上**。
- **superficial** — **未达成**：Census 的 `esophagus` 组织词下没有能对上这一分层的 cell_type 标签（检索词 ['superficial', 'squamous']）。任务要求把基底、棘层、表层分开报，这一项只能按 Census 实际提供的细胞类型标签报，**不自行把细胞类型归并到分层上**。

## 12. 药理学定位

- `absorption_site = false`

- **THE OESOPHAGUS IS NOT A DRUG ABSORPTION SITE. Normal transit of a swallowed dose is a few seconds, the epithelium is stratified squamous and keratin-poor but multilayered, and the surface area is of the order of 0.02 m2 against ~200 m2 for the small intestine. No clinically meaningful fraction of an oral dose is absorbed here. This page makes no ADME flux claim for the oesophagus.**

- **mucosal barrier** — the stratified squamous epithelium plus the tight-junction-bearing superficial layers and the submucosal glands' mucus/bicarbonate secretion are the pre-epithelial and epithelial defence against refluxed acid and pepsin. Dilated intercellular spaces are the earliest histological lesion in reflux disease.

- **oesophageal gland secretion** — the oesophageal glands proper (submucosal, mucous) and the superficial cardiac glands secrete mucus and HCO3-; secretion is vagally stimulated and is the basis of the post-prandial clearance of residual acid.

- **PILL-INDUCED OESOPHAGITIS** — this is the oesophagus's real pharmacological significance. A tablet or capsule that lodges in the oesophagus releases a high local concentration against a squamous epithelium with no mucus blanket. Classic culprits: doxycycline and other tetracyclines (direct acid injury, pH of a dissolving doxycycline capsule ~3), potassium chloride (hyperosmolar necrosis), bisphosphonates (alendronate, direct caustic injury), NSAIDs, quinidine, ferrous sulphate, ascorbic acid. Retention is most common at the measured narrowings and behind the left atrium - exactly the levels this page measures. Prophylaxis is mechanical, not pharmacological: 100-200 mL of water and remaining upright for 30 min.

- **local-acting drug delivery** — swallowed topical corticosteroid (budesonide/fluticasone) for eosinophilic oesophagitis is one of the very few deliberately oesophageal drug targets, and it works by surface contact, not absorption.

- **the four probe drugs used elsewhere in this project (midazolam/CYP3A4, dextromethorphan/CYP2D6, metformin/OCT2-MATE, fexofenadine/P-gp-OATP) have NO oesophageal ADME step. They appear on this page only as transit, not as metabolism or transport.**

## 13. 性能与移动端审计

**测试硬件：measured on a data-centre NVIDIA A10G via Modal, NOT a laptop and NOT a real phone. Phone viewports are emulated size + devicePixelRatio and touch, on desktop-class GPU and CPU.**

### 未达成：没有测到 GPU 帧率

先把"到底哪套 GL 标志能摸到 GPU"测出来（`tools/probe_gl.py`，5 套标志逐个起 Chromium 读 `UNMASKED_RENDERER_WEBGL`）：

| 标志集 | unmasked renderer | 软件光栅 |
|---|---|---|
| `angle_gl` | `ANGLE (Google, Vulkan 1.3.0 (SwiftShader Device (Subzero) (0x0000C0DE)), SwiftShader driver)` | **是** |
| `desktop_gl` | `ANGLE (Google, Vulkan 1.3.0 (SwiftShader Device (Subzero) (0x0000C0DE)), SwiftShader driver)` | **是** |
| `egl` | `ANGLE (Google, Vulkan 1.3.0 (SwiftShader Device (Subzero) (0x0000C0DE)), SwiftShader driver)` | **是** |
| `angle_vulkan` | `ANGLE (NVIDIA, Vulkan 1.4.341 (NVIDIA NVIDIA A10G (0x00002237)), NVIDIA)` | 否 |
| `angle_gl_vulkan` | `ANGLE (Google, Vulkan 1.3.0 (SwiftShader Device (Subzero) (0x0000C0DE)), SwiftShader driver)` | **是** |

沙箱里**确实有** A10G（`?`），而且 **`--use-angle=vulkan` 确实能摸到它**：`ANGLE (NVIDIA, Vulkan 1.4.341 (NVIDIA NVIDIA A10G (0x00002237)), NVIDIA)`。

这推翻了本项目早先那条经验（"装上 GLVND 派发层后 `--use-angle=gl` 就能拿到硬件渲染器，形如 `ANGLE (NVIDIA Corporation, NVIDIA A10/PCIe/SSE2, OpenGL ES 3.2)`"）：在这个镜像上 `--use-angle=gl`、`--use-gl=desktop`、`--use-gl=egl` **三套全部退回 SwiftShader**，尽管 `libEGL_nvidia.so.0`、`libGLESv2_nvidia.so.2` 与 GLVND 派发层都在。`/usr/share/vulkan/icd.d` 下只有 Mesa 的 ICD，NVIDIA 的 Vulkan ICD 注册在别处。

**但是**：本页这个场景（155 479 面 + 1 384 段血管 + 412 条流动路径）在 vulkan 后端上**三次都没能在沙箱生命周期内跑完**（各 21–22 分钟后被回收，零产物），而同一后端上一个 about:blank 空画布 16 秒就答出渲染器字符串，同一场景在 `angle_gl` 后端上 253 秒跑完整套审计。所以本节所有帧时数字都是在 **SwiftShader 软件光栅**上测的。

> **本页没有测到 GPU 帧率。** 下面每一个 ms 数字都是软件光栅数值，**不能读作 GPU 性能，也不能读作用户可见帧率**。页面 HUD 里显示的 4–25 fps 同样是软件光栅下的 rAF 节拍，不是真实设备上的帧率。

| 视口 | unmasked renderer | GPU 确认 | 软件光栅回退 |
|---|---|---|---|
| iphone14_390x844 | `ANGLE (Google, Vulkan 1.3.0 (SwiftShader Device (Subzero) (0x0000C0DE)), SwiftShader driver)` | **否** | **是 — 未测到 GPU 帧率** |
| android_360x800 | `ANGLE (Google, Vulkan 1.3.0 (SwiftShader Device (Subzero) (0x0000C0DE)), SwiftShader driver)` | **否** | **是 — 未测到 GPU 帧率** |
| desktop_1440x900 | `ANGLE (Google, Vulkan 1.3.0 (SwiftShader Device (Subzero) (0x0000C0DE)), SwiftShader driver)` | **否** | **是 — 未测到 GPU 帧率** |

| 视口 | 视图 | 每帧渲染成本 ms | p95 ms | 最大 ms | 对 16.67 ms 余量 | draw calls | 三角形 |
|---|---|---|---|---|---|---|---|
| iphone14_390x844 | organ | 0.6867 | 1.3 | 1.9 | **24.28×** | 23 | 209580 |
| iphone14_390x844 | vessels | 0.4133 | 0.5 | 0.7 | **40.33×** | 16 | 142922 |
| iphone14_390x844 | section | 0.365 | 0.5 | 0.8 | **45.67×** | 13 | 7858 |
| iphone14_390x844 | micro | 0.1817 | 0.3 | 0.5 | **91.76×** | 17 | 6648 |
| iphone14_390x844 | molecular | 0.4067 | 0.6 | 0.7 | **40.99×** | 15 | 3720 |
| android_360x800 | organ | 0.7 | 1.4 | 1.7 | **23.81×** | 23 | 209580 |
| android_360x800 | vessels | 0.4317 | 0.5 | 1.1 | **38.62×** | 16 | 142922 |
| android_360x800 | section | 0.38 | 0.6 | 0.7 | **43.87×** | 13 | 7858 |
| android_360x800 | micro | 0.1683 | 0.3 | 0.4 | **99.03×** | 17 | 6648 |
| android_360x800 | molecular | 0.3783 | 0.5 | 0.8 | **44.06×** | 15 | 3720 |
| desktop_1440x900 | organ | 0.7233 | 1.4 | 1.8 | **23.05×** | 23 | 209580 |
| desktop_1440x900 | vessels | 0.48 | 0.6 | 0.9 | **34.73×** | 16 | 142922 |
| desktop_1440x900 | section | 0.4217 | 0.6 | 0.6 | **39.53×** | 13 | 7858 |
| desktop_1440x900 | micro | 0.24 | 0.5 | 0.5 | **69.46×** | 17 | 6648 |
| desktop_1440x900 | molecular | 0.44 | 0.7 | 0.8 | **37.89×** | 15 | 3720 |

**不要把 rAF 帧率读作用户可见帧率。** 无头 Chromium 没有垂直同步，requestAnimationFrame 不设上限，报出来的是几百到几千 fps。本节报的是**每帧渲染成本 ms** 对 60 fps 的 16.67 ms 预算。各视图的 rAF 数值在 `mobile_report.json` 里以 `rAF_fps_unbounded_DO_NOT_READ_AS_USER_FPS` 为键保留，仅供对照。

| 视口 | 视图 | HUD 溢出框数 | DOM 水平溢出 | 越界元素数 | DOM 节点 |
|---|---|---|---|---|---|
| iphone14_390x844 | organ | **0** | 否 | 0 | 520 |
| iphone14_390x844 | vessels | **0** | 否 | 0 | 520 |
| iphone14_390x844 | section | **0** | 否 | 0 | 520 |
| iphone14_390x844 | micro | **0** | 否 | 0 | 520 |
| iphone14_390x844 | molecular | **0** | 否 | 0 | 520 |
| android_360x800 | organ | **0** | 否 | 0 | 520 |
| android_360x800 | vessels | **0** | 否 | 0 | 520 |
| android_360x800 | section | **0** | 否 | 0 | 520 |
| android_360x800 | micro | **0** | 否 | 0 | 520 |
| android_360x800 | molecular | **0** | 否 | 0 | 520 |
| desktop_1440x900 | organ | **0** | 否 | 0 | 520 |
| desktop_1440x900 | vessels | **0** | 否 | 0 | 520 |
| desktop_1440x900 | section | **0** | 否 | 0 | 520 |
| desktop_1440x900 | micro | **0** | 否 | 0 | 520 |
| desktop_1440x900 | molecular | **0** | 否 | 0 | 520 |

- 页面错误（pageerror）：**0**

- 控制台 error 级消息：**0**



## 14. 未达成 / 未测量 / 合成 / 超范围未调参 —— 汇总

这一节是上面各节的汇总，逐条列出。没有一条是被调参掩盖的。

- **【未达成】** 跨器官受试者一致性。本页用 s0688，不是其它四页的 s0236 —— s0236 的视野是腹部，气管/甲状腺/颈椎标签全空，食管标签只有 4443 体素（15.00 mL）。食管–胃交界处的跨器官解剖连续性因此不成立。
- **【超范围未调参】** 中心线长度 30.55 cm vs 已发表 18–26 cm —— 超出上界 17.5%。
- **【低于范围未调参】** 内切直径中位 10.79 mm vs 已发表外径 13–25 mm。内切球直径与外径不是同一个量（塌陷管腔的内切球更小），已在 provenance 注明，未重标。
- **【偏差未调参】** 隆突实测在 T6，已发表 T4–T5；GOJ 实测在 T12，已发表 T10–T11。两者都偏低，照实报。
- **【未测量】** 膈食管裂孔处的第三处生理性狭窄。TotalSegmentator v2 无膈/膈脚/裂孔标签，在该受试者自己的 T10 层级实测内切半径 6.76 mm、无局部极小值。该狭窄按已发表椎体层级标注，`contour_anchored: false`。
- **【未归类】** 剖面上最深的半径极小值（T8、弧长 191.3 mm、突出度 3.14 mm）不对应任何经典狭窄，报为 MEASURED, UNCLASSIFIED，不声称它是第三处狭窄。
- **【合成】** 全部七棵血管树、栅栏区网络、奇静脉走行、管壁四层、肌层三分、LES 位置、括约肌压力、蠕动波速与深度、微观单元尺寸 —— 逐项 `contour_anchored: false`。
- **【未达成】** 25/398 个外膜目标点（6.3%）在食管周围走行空间内找不到同时满足间隙约束的分叉位置，未能接上。已逐树报数并给出 rejected_outside / rejected_clearance 分解，**没有静默丢弃**。
- **【低于范围未调参】** 七棵树里4棵的主干口径低于已发表直径（奇静脉 −58%、甲状腺下动脉 −17%、胃左静脉 −56%、甲状腺下静脉 −32%）。原因是 Murray 级联在 34–97 个终端上只能跨 3.2–4.6 倍半径比，而真实主干/终端比是 20–40 倍。选择了把终端口径对齐已发表小血管范围，主干因此偏低。
- **【高于范围未调参】** 七棵树的终端流速全部高于已发表 0.5–1.5 cm/s（0.3–0.6 mm 小血管），超出 2.2–9.8 倍。同一结构性原因。
- **【未检索到参考】** 七棵树里四棵（支气管动脉、主动脉食管支、甲状腺下静脉、奇静脉）**未检索到已发表的多普勒流速范围**，其主干流速是声明的假设值，逐条标注在 flow_report.json。栅栏区流速 8 mm/s 同样是声明值。
- **【标出：等于求解器目标】** 21 对树里有 2 对的最小表面间隙**恰好等于** 0.05 mm 的求解器间隙目标，已标 `at_solver_target: true` —— **这是求解器的属性，不是对几何的独立测量**。
- **【标出：顶到声明上界】** 节点位移修复的最大位移 1.9988 mm **恰好顶到**声明上界 2.0 mm，应读作上界在起约束作用，而不是修复轻松收敛。
- **【放弃】** 余弦弯曲（cosmetic bow）。任何非零弯曲幅度都会让相邻树互相穿透（实测 −0.16 ~ −0.46 mm），解剖不相交优先于外观曲率，弯曲幅度设为 0。
- **【未达成】** 鳞状上皮三层的分开报告。Census 2025-11-08 的 esophagus 组织词下只有 `basal cell`（7 239）与 `cycling basal cell`（12）能对上基底层；**没有** suprabasal/spinous 与 superficial 的 cell_type 标签。占主导的 `stratified epithelial cell`（51 307）同时涵盖棘层与表层而不加区分，报为独立的"未分层"行，**不把它归并到任何一层**。
- **【未达成】** 在食管鳞状上皮中定位证据最好的 5 个分子（ABCC1、SLC16A1、ALDH3A1、ADH1C、CYP2C18）**都不在本项目固定的 61 基因面板里**，因此未被测量、也未被绘制。面板里能同时满足"在面板内 + 有定位证据"的只有 3 个（ABCB1、ABCC3、GSTP1）。
- **【未使用】** 任何标定系数。表达量→通量、丰度→电导之类**一个都没有**，分子层只做定性映射。
- **【未测量】** GPU 帧率。已测出 `--use-angle=vulkan` 能摸到 A10G，但本页场景在该后端上三次都未能在沙箱生命周期内跑完；全部帧时数字测于 SwiftShader 软件光栅，**不能读作 GPU 性能**。
- **【占位数据】** 搏动耦合的上游心率数据是 PhysioNet MIT-BIH Normal Sinus Rhythm Database (nsrdb) record 16265 (32 M, 22.24 h) — PUBLIC PLACEHOLDER, NOT THE USER'S OWN DATA，**不是使用者本人的数据**（`is_user_personal_data = false`）。
- **【未建模】** 压力。括约肌压力、蠕动收缩幅度、门脉压力全部只是标注的已发表数值；反流演示与门脉高压演示都是**预先烘焙的路径切换，不是压力–流量网络求解**。

