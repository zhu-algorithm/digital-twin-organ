# GEOMETRY_PROVENANCE — 三维数字孪生肾脏

生成日期 2026-10-29。本文件中的每个数字都从 `assets/*.json` 读出，不是手写的。

## 1. 唯一几何来源

| 项 | 值 |
|---|---|
| 数据集 | TotalSegmentator dataset v2.0.1 (Wasserthal et al.) |
| DOI | `10.5281/zenodo.10047292` |
| 许可 | **CC-BY-4.0** — 署名类，**可在公开仓再分发** |
| 受试者 | `s0236` |
| 体素 | [1.5, 1.4999998807907104, 1.5] mm |
| 取数方式 | HTTP range extraction of 7 mask members from the Zenodo zip; no CT image volume downloaded, nothing derived from greyscale |
| 使用的结构 | `kidney_left`, `kidney_right`, `aorta`, `inferior_vena_cava`, `urinary_bladder`, `kidney_cyst_left`, `kidney_cyst_right` |

**不下载 CT 影像本体。** 本项目不含任何基于影像灰度的内容：所有几何都来自二值掩膜。
对 23.5 GB 的 Zenodo 压缩包只做 HTTP Range 读取，按需取出上述掩膜成员。

### 已被否决、未使用的来源
- Gibson 等多器官腹部 CT 分割（Zenodo `10.5281/zenodo.1169361`）：其自带条款叠加 TCIA
  Pancreas-CT 与 BTCV 挑战赛条款，**再分发权利不明**，不用。
- 任何带 NC（非商业）条款的网格：与本仓 MIT 许可不兼容，不用。
- BodyParts3D、3D-IRCADb-01：不可达且非纯署名许可，不用。

## 2. 受试者筛选（按肾脏自己的标准重新核验）

与本项目肝脏使用同一受试者 `s0236` 以保证跨器官解剖一致性，但**并未因此免检**。
筛选标准（全部为实测，无一为假设）：

1. 两侧 `kidney_left` / `kidney_right` 标签均存在且**完整落在视野内**
2. 每侧肾脏为**单连通分量**
3. 主轴（上下极）长径落在成人区间 [9.0, 13.0] cm
4. 实质体积落在成人区间 [110.0, 200.0] mL
5. `aorta`、`inferior_vena_cava`、`urinary_bladder` 三个标签均存在
   （血管开口与输尿管远端锚点取自它们）
6. `kidney_cyst_*` 体积可忽略（< 2.0 mL）——
   囊肿会挤占它所在的锥体

筛查 **18 例，2 例通过**。
`s0236` 通过；脚本见 `tools/screen_kidney_subjects.py`，逐例结果在 `out/kidney_subject_screen.json`。

| 受试者 | 左肾 | 右肾 | 囊肿 mL | 结果 |
|---|---|---|---|---|
| `s0236` | 174.8 mL / 11.80 cm / 1 comp / FOV-edge 0 | 194.1 mL / 12.59 cm / 1 comp / FOV-edge 0 | 0.0 | **pass** |
| `s0000` | label absent | label absent | 0.0 | fail |
| `s0084` | label absent | label absent | 0.0 | fail |
| `s0169` | label absent | label absent | 0.0 | fail |
| `s0252` | 184.7 mL / 12.35 cm / 1 comp / FOV-edge 0 | 182.0 mL / 12.23 cm / 1 comp / FOV-edge 0 | 0.0 | **pass** |
| `s0335` | label absent | label absent | 0.0 | fail |
| `s0411` | label absent | label absent | 0.0 | fail |
| `s0487` | 1.4 mL / 2.83 cm / 1 comp / FOV-edge 1 | label absent | 0.0 | fail |
| `s0567` | label absent | label absent | 0.0 | fail |
| `s0643` | label absent | label absent | 0.0 | fail |
| `s0714` | 27.0 mL / 6.29 cm / 1 comp / FOV-edge 1 | 37.4 mL / 6.34 cm / 1 comp / FOV-edge 1 | 0.0 | fail |
| `s0791` | 1.1 mL / 3.29 cm / 1 comp / FOV-edge 1 | label absent | 0.0 | fail |
| `s0868` | label absent | label absent | 0.0 | fail |
| `s0953` | label absent | label absent | 0.0 | fail |
| `s1029` | 213.7 mL / 12.46 cm / 1 comp / FOV-edge 0 | 205.0 mL / 12.53 cm / 1 comp / FOV-edge 0 | 0.0 | fail |
| `s1110` | 112.3 mL / 8.78 cm / 1 comp / FOV-edge 1 | 216.6 mL / 11.52 cm / 1 comp / FOV-edge 1 | 0.0 | fail |
| `s1189` | 90.5 mL / 9.77 cm / 1 comp / FOV-edge 0 | 106.0 mL / 10.18 cm / 1 comp / FOV-edge 0 | 0.0 | fail |
| `s1273` | 266.4 mL / 14.64 cm / 1 comp / FOV-edge 0 | 222.1 mL / 12.86 cm / 1 comp / FOV-edge 0 | 0.0 | fail |

## 3. 实测几何量

| 量 | 左 | 右 | 参考区间 | 判定 |
|---|---|---|---|---|
| 实质体积 mL | 174.8 | 194.1 | [110.0, 200.0] | 均在区间内 |
| 长径 cm | 11.8 | 12.59 | [9.0, 13.0] | 均在区间内 |
| 横径 cm | 6.91 | 7.12 | [4.0, 7.0] | 左 在 / 右 **超出** |
| 前后径 cm | 5.39 | 5.74 | — | 无参考区间 |
| 连通分量 | 1 | 1 | 1 | 通过 |
| 触及视野边缘 | False | False | false | 通过 |
| 肾门凹陷深度 mm | 49.86 | 45.79 | — | 实测，无参考 |
| 肾窦体积 mL | 37.89 | 43.07 | — | 实测，无参考 |
| 肾门→主动脉 mm | 51.2 | 34.21 | — | 实测 |
| 肾门→下腔静脉 mm | 76.84 | 19.71 | — | 实测 |

**本例右肾（194.1 mL）大于左肾（174.8 mL）**，
与"左肾通常略大于右肾"的人群趋势相反。照实报告，未做任何处理。

肾窦与肾门凹陷的测量方式：per-slice convex-hull deficit of the kidney label (measured on the label, not assumed)。

## 4. 源数据中不存在、因而必须构造的部分

源发布每侧肾脏只有**一个标签**。以下标注**完全不存在**：
- `renal_artery`
- `renal_vein`
- `renal_cortex`
- `renal_medulla`
- `renal_pelvis`
- `calyx`
- `ureter`

后果（the vascular cascade, the collecting system and the corticomedullary architecture are CONSTRUCTED, not contour-anchored; only the capsule, the hilar notch, the renal sinus and the aortic/caval/bladder anchors are measured）：

| 结构 | 轮廓锚定？ | 依据 |
|---|---|---|
| 肾包膜（外表面） | **是** | `kidney_*` 标签本身 |
| 肾门凹陷、肾窦 | **是** | 逐冠状层凸包亏损，在标签上量出 |
| 肾动脉开口 | **是（仅开口）** | 该受试者真实主动脉表面上距肾门最近点 |
| 肾静脉开口 | **是（仅开口）** | 该受试者真实下腔静脉表面上距肾门最近点 |
| 输尿管远端锚点 | **是** | 真实 `urinary_bladder` 标签上极 |
| 肾乳头位置 | 由实测几何导出 | 深部实质中最靠近实测肾窦的面，最远点采样 |
| 皮质 / 髓质 / 皮髓交界 | **否，构造** | 外包膜距离场 + 公开皮质厚度 8.0 mm |
| 髓质锥体 | **否，构造** | 肾乳头 Voronoi 划分 |
| Bertin 肾柱 | **否，构造** | 相邻锥体之间声明半宽 1.6 mm 的隔 |
| 肾动脉主干 / 段 / 叶间 / 弓形 / 皮质放射 | **否，构造** | 实测肾窦、肾柱、锥体底几何 |
| 肾静脉全树 | **否，构造** | 同上 |
| 肾小球与球后网络 | **否，构造** | 公开小动脉与毛细血管口径 |
| 肾盏 / 肾盂 / 输尿管走行 | **否，构造** | 实测肾乳头与膀胱锚点之间 |

**皮质厚度的来源必须讲清楚**：midpoint of the published adult 6-10 mm range; NOT measurable from a single-label segmentation and NOT measured on this subject。

## 5. 面数与文件

| 项 | 值 |
|---|---|
| 三角面合计 | 153,580 |
| 顶点合计 | 76,870 |
| `kidney.glb` | 4.99 MB |
| 逐部件面数 | capsule_left 31,608, capsule_right 34,972, medulla_left 11,000, medulla_right 11,000, cut_left 20,000, cut_right 20,000, cut_medulla_left 7,000, cut_medulla_right 7,000, aorta 6,000, ivc 5,000 |

烘进网格的逐顶点场：
`_FIELD0 = ['capsule_depth_01', 'pyramid_id', 'medulla_01', 'concavity']`，
`_FIELD1 = ['column_of_bertin', 'cmj_distance_norm', 'corticomedullary_zone', 'side']`。
其中**实测的只有** ['capsule_depth_01', 'concavity']，
其余 ['pyramid_id', 'medulla_01', 'column_of_bertin', 'cmj_distance_norm', 'corticomedullary_zone'] 均为构造量。

## 6. 与心脏模型的耦合

心输出量、平均动脉压、心率与主动脉压波形全部取自本项目已验证的心脏模型，
**不另造数字**：CO 5.1771 L/min，MAP 94.4872 mmHg，
HR 64.45 bpm。呼吸频率 0.18323 Hz
（10.99 次/分）沿用同一来源。

该来源是 **PhysioNet nsrdb 记录 16265 的公开数据**，
**不是用户本人的任何记录**：所有文件中 `is_user_personal_data = false`。
