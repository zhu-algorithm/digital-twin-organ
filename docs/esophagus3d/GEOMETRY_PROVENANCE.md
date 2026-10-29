# GEOMETRY_PROVENANCE — 三维数字孪生食管

生成日期 2026-10-30。本文件中的每个数字都从 `assets/*.json` 与 `out/*.json` 读出，不是手写的。

## 1. 唯一几何来源

| 项 | 值 |
|---|---|
| 数据集 | TotalSegmentator dataset v2.0.1 (Wasserthal J. et al.) |
| DOI | `10.5281/zenodo.10047292` |
| 许可 | **CC-BY-4.0** — 署名类，**可在公开仓再分发** |
| 受试者 | `s0688` |
| 体素 | [1.5, 1.5, 1.5] mm |
| 网格 | [241, 241, 234] |
| 该版本结构数 | 117 |
| 取数方式 | HTTP range extraction of 24 mask members from the Zenodo zip; no CT image volume downloaded, nothing derived from greyscale |
| 使用的结构 | `esophagus`, `trachea`, `aorta`, `stomach`, `thyroid_gland`, `spinal_cord`, `vertebrae`, `heart`, `portal_vein_and_splenic_vein`, `brachiocephalic_vein_left`, `liver`, `lung_lower_lobe_left`, `lung_upper_lobe_right`, `subclavian_artery_left`, `common_carotid_artery_left`, `brachiocephalic_trunk`, `atrial_appendage_left`, `spleen` |

**不下载 CT 影像本体。** 本项目不含任何基于影像灰度的内容：所有几何都来自二值掩膜。对 23.5 GB 的 Zenodo 压缩包只做 HTTP Range 读取，按需取出上述掩膜成员。

### 该版本标签清单的实测核实

从 zip 中央目录枚举，该版本每位受试者的 `segmentations/` 下共 **117** 个结构标签文件。**没有** head_glands_cavities 等附加任务的输出，因此没有腮腺、下颌下腺、咽腔、喉等头颈腺体与腔隙标签 —— 这一点是查出来的，不是假设的。本器官需要而该版本**不提供**的标签：

- `diaphragm`
- `diaphragmatic crus`
- `oesophageal hiatus`
- `azygos vein`
- `hemiazygos vein`
- `bronchial artery`
- `inferior thyroid artery`
- `left gastric artery`
- `left gastric vein`
- `coeliac trunk`
- `submucosal venous plexus`
- `muscularis propria`
- `mucosa`
- `oesophageal glands`

**后果：** the oesophageal outer surface, its centreline, the inscribed-radius profile, the vertebral levels and every neighbour distance are MEASURED. Every vessel, the wall layering, the muscle-type transition, both sphincter locations/pressures and the palisade zone are CONSTRUCTED from published anatomy and are individually flagged contour_anchored = false.

### 已被否决、未使用的来源
- Gibson 等多器官腹部 CT 分割（Zenodo `10.5281/zenodo.1169361`）：条款叠加 TCIA Pancreas-CT 与 BTCV，再分发权利不明，不用。
- 任何带 NC（非商业）条款的网格：与本仓 MIT 许可不兼容，不用。
- BodyParts3D、3D-IRCADb-01：不可达且非纯署名许可，不用。

## 2. 受试者筛选 —— 为什么**不是** s0236

本项目的心脏、肝脏、肾脏、小肠四页都用受试者 `s0236`。跨器官解剖一致性是既定偏好，所以 s0236 **被第一个拿来按食管自己的标准复核** —— 它不及格。

| s0236 实测 | 值 |
|---|---|
| `esophagus` 标签体素 | 4443 |
| `esophagus` 标签体积 | 15.00 mL |
| `trachea` 标签体素 | 0 |
| `thyroid_gland` 标签体素 | 0 |
| 颈椎至 T9 椎体标签体素 | 0 |

**不及格原因：** abdominal field of view. The trachea, thyroid_gland and every cervical / upper-thoracic vertebra label is EMPTY, and the esophagus label holds only 4443 voxels (15.00 mL) covering the distal segment. It fails S2, S3 and S5 outright.

**后果（必须读到）：** this page uses s0688. The oesophagus is therefore NOT from the same subject as the other four organ pages, and the cross-organ anatomical continuity claimed elsewhere in the project does not hold for the oesophagus-stomach boundary rendered here. The s0236 distal-oesophagus overlap is reported separately below for reference.

### 第一轮筛选标准（在任何测量之前声明）

- **S1** — largest connected component (26-connectivity) holds >= 95% of labelled voxels
- **S2** — per-slice centre-of-mass centreline length >= 200 mm
- **S3** — trachea label top no more than 5 mm below the oesophagus label top (cervical oesophagus inside the FOV)
- **S4** — minimum surface gap from the distal 25 mm of the oesophagus label to the stomach label <= 10 mm
- **S5** — aorta label spans the whole oesophageal cranio-caudal range
- **S6** — maximum inscribed diameter <= 25 mm (exclude disease-driven dilatation)

预筛：compressed size of esophagus.nii.gz > 3x the s0236 value AND trachea.nii.gz > 2000 B AND stomach.nii.gz > 2000 B (free: read from the zip central directory, no download)

预筛命中 **62** 例，实测 **22** 例，第一轮通过 **14** 例。

### 标准的更正（照实记录，不静默修补）

**S5** — 问题：as written, S5 is anatomically impossible - the aorta begins at the arch (~T4) and can never span the CERVICAL oesophagus, so all 22 measured subjects failed it.

  更正：S5 restated as: the aortic arch apex must lie above the oesophageal mid-point (arch inside the FOV) AND the aorta must extend to within 5 mm of the oesophagus label's caudal end (hiatus reached).（after round 1, before any subject was selected）

**S3** — 问题：S3 is one-sided: it rejects a subject whose trachea stops below the oesophagus top but not one whose oesophagus label is truncated far BELOW the larynx. s0556 passed S3 with the trachea top 118.5 mm ABOVE the oesophagus top, which means its cervical oesophagus is unlabelled.

  更正：NOT corrected - s0556 was not selected, so the outcome is unaffected. Recorded here rather than silently fixed.（noticed after round 2）

### 第二轮：事后追加的标准

**追加原因（必须读到）：** the round-1 winner s0692 has a 920.96 mL aorta (an ectatic/aneurysmal thoracoabdominal aorta) and an oesophageal tortuosity of 1.200 - the aneurysm displaces the oesophagus. The declared screen tested the oesophagus but never tested whether its NEIGHBOURS were normal. These two criteria were added AFTER seeing that result and are flagged as post hoc.

- **S7** — aorta label volume between 90 and 350 mL (declared band for a non-ectatic thoracoabdominal aorta; not a single published figure)
- **S8** — oesophageal centreline tortuosity (path / straight-line end-to-end) <= 1.12

第二轮通过 **9** 例；决胜规则：largest labelled voxel count among passing subjects (most surface support for meshing) → 选中 **`s0688`**。

| 受试者 | 食管体素 | 中心线 mm | 迂曲度 | 主动脉 mL | S7 | S8 | 结果 |
|---|---|---|---|---|---|---|---|
| `s0692` | 22015 | 306.5 | 1.2 | 921.0 | ✗ | ✗ | fail |
| `s1361` | 18075 | 290.4 | 1.185 | 822.8 | ✗ | ✗ | fail |
| `s0688` | 17497 | 305.5 | 1.043 | 309.8 | ✓ | ✓ | **pass** |
| `s0751` | 17447 | 286.8 | 1.05 | 497.8 | ✗ | ✓ | fail |
| `s1135` | 17298 | 256.8 | 1.056 | 239.4 | ✓ | ✓ | **pass** |
| `s0918` | 17164 | 282.8 | 1.115 | 321.8 | ✓ | ✓ | **pass** |
| `s0550` | 15665 | 282.7 | 1.1 | 312.6 | ✓ | ✓ | **pass** |
| `s1152` | 15630 | 259.5 | 1.128 | 192.1 | ✓ | ✗ | fail |
| `s0657` | 15088 | 289.9 | 1.042 | 229.8 | ✓ | ✓ | **pass** |
| `s0556` | 14645 | 275.5 | 1.067 | 269.6 | ✓ | ✓ | **pass** |
| `s0771` | 13486 | 286.4 | 1.06 | 203.2 | ✓ | ✓ | **pass** |
| `s1304` | 13428 | 265.2 | 1.101 | 333.9 | ✓ | ✓ | **pass** |
| `s1082` | 13192 | 327.3 | 1.132 | 362.1 | ✗ | ✗ | fail |
| `s1088` | 10980 | 294.4 | 1.113 | 308.0 | ✓ | ✓ | **pass** |

### 选中受试者的实测值

| 量 | 值 |
|---|---|
| esophagus_voxels | 17497 |
| esophagus_mL | 59.05 |
| centreline_mm | 305.5 |
| tortuosity | 1.043 |
| aorta_mL | 309.8 |
| max_inscribed_diameter_mm | 18.5 |
| gap_to_stomach_mm | 1.5 |
| cranial_level | C6 |
| caudal_level | T12 |
| connected_components | 1 |

## 3. 实测 vs 已发表参考

| 量 | 实测 | 已发表参考 | 判定 |
|---|---|---|---|
| 中心线长度 | 30.55 cm | 18–26 cm | **超出上界 +17.5%，照实报，未调参** |
| 迂曲度 | 1.043 | ≤1.12（本页筛选带） | 合格 |
| 椎体跨度 | C6 → T12 | C6 → T10/T11 | 远端偏低一至两节 |
| 内切直径 中位 | 10.79 mm | 13–25 mm（外径） | 低于范围 |
| 内切直径 最大 | 17.98 mm | 13–25 mm（外径） | 在范围内 |
| 标签体积 | 59.05 mL | — | 实测 |
| 连通分量 | 1 | 1 | 合格 |
| 至胃标签最小表面间隙 | 1.50 mm | — | 贲门锚定成立 |
| 至主动脉最小距离 | 5.27 mm | — | 实测 |
| 至气管最小距离 | 3.63 mm | — | 实测 |
| 至甲状腺最小距离 | 10.91 mm | — | 实测 |
| 隆突 | 弧长 128.2 mm / T6 | T4–T5 | 偏低 |
| 主动脉弓顶 | 弧长 83.1 mm / T4 | T4 | 合格 |
| 食管周围可用走行空间 | 368.3 mL | — | 实测 |

**长度超范围的说明（不调参）：** MEASURED 30.55 cm vs published 18-26 cm -> ABOVE the published range by 17.5%. Reported as computed; no parameter was adjusted. The label spans C6 to T12 (cricopharyngeal region to the gastro-oesophageal junction, 1.50 mm minimum surface gap to the stomach label), so it is not obviously over-extended at either end. Candidate contributors, none of which this project can separate: in-vivo supine CT length exceeds cadaveric and endoscopic estimates; this subject's GOJ label sits at T12, lower than the usual T10-T11; and a centreline through a curved tube is longer than a cranio-caudal measurement (measured tortuosity 1.043).

**内切直径的说明：** median 10.80 mm and max 17.98 mm, both BELOW the 13-25 mm published OUTER diameter range at the median. The label follows the outer wall of a COLLAPSED resting oesophagus, and the maximum inscribed sphere of a flattened slit is smaller than the outer transverse diameter, so these two quantities are not directly comparable. Reported as measured with that caveat rather than rescaled.

## 4. 网格

| 结构 | 面 | 顶点 | 表面积 mm² | 网格体积 mL | 标签体积 mL | 上采样 | 平均三角形 mm² |
|---|---|---|---|---|---|---|---|
| `esophagus` | 61999 | 31036 | 15610.2 | 56.39 | 59.05 | 2× | 0.2518 |
| `trachea` | 15999 | 8029 | 11942.4 | 55.24 | 60.54 | 1× | 0.7464 |
| `aorta` | 22000 | 11022 | 41154.0 | 297.7 | 309.81 | 1× | 1.8706 |
| `stomach` | 14999 | 7553 | 22792.1 | 231.97 | 242.44 | 1× | 1.5196 |
| `thyroid_gland` | 3482 | 1793 | 2194.1 | 6.28 | 8.47 | 1× | 0.6301 |
| `spinal_cord` | 6000 | 3021 | 17833.0 | 66.89 | 75.77 | 1× | 2.9722 |
| `vertebrae` | 22000 | 10960 | 124118.6 | 556.12 | 589.66 | 1× | 5.6418 |
| `heart` | 9000 | 4502 | 51364.9 | 538.96 | 552.97 | 1× | 5.7072 |
| **合计** | **155479** | | | | | | |

门槛 13万–16万：**合格**；glb 3.75 MB（门槛 ≤ 7.0 MB）。食管本体占 61999 面（39.9%）。

**分辨率声明：** the oesophageal surface is marching cubes on the label TRILINEARLY UPSAMPLED x2. Mean triangle area 0.252 mm2 (edge ~0.76 mm) matches the upsampled grid pitch, but the INFORMATION CONTENT is that of the native 1.5 mm label - the extra triangles are interpolation, not detail. Every other structure is meshed at native resolution.

**面数构成声明：** the 130k-160k face gate is met by the eight-structure scene, not by the oesophagus alone. A 305 mm tube of 15 610 mm2 surface area supports about 14 000 native-resolution triangles; inflating it alone to 140 000 would be pure interpolation. All eight structures are real segmentations from the same subject and the per-structure counts are listed above.

## 5. 逐项 contour_anchored

| 结构 / 量 | contour_anchored |
|---|---|
| esophagus_outer_surface | **true（实测）** |
| centreline | **true（实测）** |
| inscribed_radius_profile | **true（实测）** |
| vertebral_level_map | **true（实测）** |
| neighbour_distances | **true（实测）** |
| gastro_oesophageal_junction | **true（实测）** |
| constriction_1_cricopharyngeal | **true（实测）** |
| constriction_2_aortic_arch | **true（实测）** |
| constriction_2b_left_main_bronchus | **true（实测）** |
| constriction_3_diaphragmatic_hiatus | false（构造） |
| upper_oesophageal_sphincter_location | **true（实测）** |
| lower_oesophageal_sphincter_location | false（构造） |
| wall_four_layer_stack | false（构造） |
| muscle_type_thirds | false（构造） |
| sphincter_pressures | false（构造） |
| arterial_trees | false（构造） |
| venous_trees | false（构造） |
| azygos_course | false（构造） |
| palisade_zone | false（构造） |
| peristaltic_wave_speed | false（构造） |
| micro_cross_section_unit | false（构造） |
| molecular_layer | false（构造） |

