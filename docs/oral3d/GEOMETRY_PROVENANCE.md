# GEOMETRY_PROVENANCE — 三维数字孪生口腔

## 1. 唯一几何来源

- **数据集** TotalSegmentator dataset v2.0.1 (Wasserthal et al.)
- **DOI** `10.5281/zenodo.10047292`
- **许可** **CC-BY-4.0** — attribution only, redistribution permitted。可再分发，本仓库仅需署名。
- **受试者** `s1304`
- **体素** 1.5 × 1.5 × 1.5 mm
- **取用方式** HTTP range extraction of 13 mask members from the Zenodo zip; no CT image volume downloaded, nothing derived from greyscale
- **用到的掩膜** `skull`, `vertebrae_C1`, `vertebrae_C2`, `vertebrae_C3`, `esophagus`, `trachea`, `thyroid_gland`, `brachiocephalic_vein_left`, `brachiocephalic_vein_right`, `superior_vena_cava`, `heart`, `common_carotid_artery_left`, `common_carotid_artery_right`
- **方位核验** RAS +y = anterior (trachea centroid 13.9 mm anterior of oesophagus centroid)

没有下载任何 CT 影像本体。本项目**不含任何基于影像灰度的内容**：所有锚定量都来自二值掩膜的几何，而不是 HU 值。

### 被否决、未使用的来源

- Gibson et al., Zenodo `10.5281/zenodo.1169361` — 条款叠加 TCIA Pancreas-CT 与 BTCV，不用。
- 任何含 NC（非商业）条款的网格 — 与本仓库 MIT 许可冲突，不用。
- BodyParts3D、3D-IRCADb-01 — 不用。
- TotalSegmentator 的附加任务（含头颈腺体与腔隙的 `head_glands_cavities` 等）是**模型权重，不是数据**，不属于本 CC-BY-4.0 数据集，未使用。

## 2. 受试者选择：为什么不是 s0236

本项目的心脏、肝脏、肾脏、小肠都用受试者 `s0236`。本器官**无法**沿用。
`s0236` 的 `meta.csv` 记录为 `ct abdomen-pelvis`，实测其头颈部掩膜：

| 掩膜 | s0236 体素数 |
|---|---|
| `skull` | 0（空） |
| `brain` | 0（空） |
| `vertebrae_C1` | 0（空） |
| `vertebrae_C2` | 0（空） |
| `vertebrae_C3` | 0（空） |
| `vertebrae_C4` | 0（空） |
| `trachea` | 0（空） |
| `thyroid_gland` | 0（空） |
| `common_carotid_artery_left` | 0（空） |
| `common_carotid_artery_right` | 0（空） |

头部完全在扫描野之外，口腔**不可能**锚定到 s0236。**跨器官受试者一致性在这一个器官上被打破，这一点在页面正文里也写明了。**

### 本器官自己的筛选标准（先定标准，再筛）

1. `skull` 掩膜非空且最大连通体积 ≥ 600 mL（排除只切到颅顶的扫描）
2. 下颌弓在扫描野内：`skull` 掩膜在其下部轴位层面上形成明确的弓形
3. `vertebrae_C1`、`C2`、`C3` 均非空（确认口咽水平在野内）
4. 吞咽链锚点齐全：`trachea`、`esophagus`、`thyroid_gland` 非空
5. 静脉链锚点齐全：`brachiocephalic_vein_left/right`、`superior_vena_cava`、`heart` 非空
6. `pathology == no_pathology` 且 `pathology_location == no_location`
7. 各向同性 1.5 mm 体素

### 实测值（被选中的受试者）

| 标准 | 实测 |
|---|---|
| `skull` 标签体积 | 998.2 mL |
| `skull` 连通分量数 / 丢弃的碎点 | 4 / 8 体素 |
| 颅盖汇合层面 z | 194.5 mm |
| 口腔骨性框架体积 | 135.25 mL |
| 框架连通分量 / 丢弃体素 | 1 / 0 |
| 咬合平面 z | 149.5 mm |
| 下颌部分体积 | 85.96 mL |
| 中面部部分体积 | 49.29 mL |
| `trachea` | 47.64 mL |
| `esophagus` | 46.99 mL |
| `thyroid_gland` | 35.82 mL |
| `vertebrae_C1` | 18.8 mL |
| `vertebrae_C2` | 19.38 mL |
| `vertebrae_C3` | 13.75 mL |
| `brachiocephalic_vein_left` | 21.49 mL |
| `brachiocephalic_vein_right` | 10.6 mL |
| `superior_vena_cava` | 33.66 mL |
| `heart` | 1012.17 mL |
| `common_carotid_artery_left` | 9.46 mL |
| `common_carotid_artery_right` | 7.12 mL |

## 3. 提取规则（每一条都是可复现的实测规则，不是手工描点）

- **颅盖汇合层面**：first axial level at which the upward-grown retained area jumps >=4x and the retained posterior extent jumps >=30 mm; this is the temporomandibular/zygomatic arch level at which the cranial vault joins。实测 z = 194.5 mm。
- **咬合平面**：local minimum of the anterior (y > 173.5 mm RAS, the incisor region) cross-sectional area of the frame, searched over its middle 60%。实测 z = 149.5 mm。
- **下颌 / 中面部分割**：以咬合平面上下的体素为两个种子，在骨性框架内做竞争性膨胀（测地分水），不是 z 平面硬切——硬切会把升支和髁突从下颌上削掉。
- **牙槽弓曲线**：在上述两个层面带内按极角分箱取外缘，再截断到声明的 42.0 mm 弓深窗口内。**这个窗口是声明参数，不是本受试者的测量值**；窗口内测到的宽度是。

## 4. 哪些是轮廓锚定的，哪些不是

| 结构 | contour_anchored |
|---|---|
| `mandible` | **true** |
| `midface` | **true** |
| `vertebrae_C1` | **true** |
| `vertebrae_C2` | **true** |
| `vertebrae_C3` | **true** |
| `trachea` | **true** |
| `esophagus` | **true** |
| `thyroid_gland` | **true** |
| `brachiocephalic_vein_left` | **true** |
| `brachiocephalic_vein_right` | **true** |
| `superior_vena_cava` | **true** |
| `heart` | **true** |
| `common_carotid_artery_left` | **true** |
| `common_carotid_artery_right` | **true** |
| 上下牙槽弓曲线、髁突轴、颏点、咬合平面 | **true**（在上述骨面上实测） |
| `buccal_left` | **false** |
| `buccal_right` | **false** |
| `floor_of_mouth` | **false** |
| `gingiva_lower` | **false** |
| `gingiva_upper` | **false** |
| `hard_palate` | **false** |
| `labial` | **false** |
| `masseter_left` | **false** |
| `masseter_right` | **false** |
| `micro_acinus` | **false** |
| `micro_mucosa_buccal` | **false** |
| `micro_mucosa_hard_palate` | **false** |
| `micro_mucosa_sublingual` | **false** |
| `micro_papilla` | **false** |
| `papillae_circumvallate` | **false** |
| `papillae_filiform` | **false** |
| `papillae_foliate` | **false** |
| `papillae_fungiform` | **false** |
| `parotid_left` | **false** |
| `parotid_right` | **false** |
| `soft_palate` | **false** |
| `sublingual_left` | **false** |
| `sublingual_right` | **false** |
| `submandibular_left` | **false** |
| `submandibular_right` | **false** |
| `teeth_lower` | **false** |
| `teeth_upper` | **false** |
| `temporalis_left` | **false** |
| `temporalis_right` | **false** |
| `tongue` | **false** |
| `uvula` | **false** |
| 32 颗牙的**位置** | 其中 30 颗在实测弓上，2 颗在实测弓的切线外推段上 |
| 32 颗牙的**尺寸** | **false** — Ash MM & Nelson SJ, Wheeler's Dental Anatomy, Physiology and Occlusion, 9th ed. - table of average crown dimensions |
| `lingual_artery` | root only (carotid surface) |
| `facial_artery` | root only (carotid surface) + crosses the measured mandibular antegonial notch |
| `inferior_alveolar_artery` | both ends on the measured mandibular surface |
| `venous_return` | root only (brachiocephalic vein surface) |
| `parotid_duct` | ostium only (from the measured maxillary arch) |
| `submandibular_duct` | ostium only (from the measured mandibular arch) |
| `sublingual_ducts` | False |

## 5. 管道树锚点的来源

- `carotid_top_*` — measured: the superior-most vertex of the common_carotid_artery mask surface (the carotid bifurcation level), contour-anchored
- `bcv_*` — measured: the superior-most vertex of the brachiocephalic vein mask surface, contour-anchored
- `mandibular_foramen_*` — measured: the most medial mandibular surface vertex in a 8 mm band midway between the occlusal plane and the condyle. Anatomically correct level, but the foramen itself is NOT labelled in the source - it is placed on the measured surface, not detected on it.
- `mental_foramen_*` — measured surface vertex under the second premolar; same caveat as the mandibular foramen
- `parotid_ostium_*` — buccal mucosa opposite the maxillary second molar (tooth 17/27), whose position comes from the measured maxillary arch
- `caruncle_*` — sublingual caruncle beside the mandibular central incisor (tooth 41/31), from the measured mandibular arch
