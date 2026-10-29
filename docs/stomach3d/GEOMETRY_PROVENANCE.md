# GEOMETRY_PROVENANCE — 三维数字孪生胃

## 1. 唯一几何来源

| 项 | 值 |
|---|---|
| 数据集 | TotalSegmentator dataset v2.0.1 (Wasserthal et al. 2023) |
| DOI | 10.5281/zenodo.10047292 |
| 许可 | **CC-BY-4.0** — attribution only; redistribution of derivatives permitted |
| 受试者 | `s0236` |
| 体素 | [1.5, 1.4999998807907104, 1.5] mm |
| 体数据尺寸 | [252, 198, 308] |
| 取用方式 | HTTP range extraction of 13 mask members from the Zenodo zip; no CT image volume downloaded, nothing derived from greyscale |
| 用到的结构 | `stomach`, `esophagus`, `duodenum`, `liver`, `spleen`, `pancreas`, `aorta`, `portal_vein_and_splenic_vein`, `small_bowel`, `colon`, `inferior_vena_cava`, `vertebrae_T12`, `vertebrae_L1` |

**没有下载 CT 影像本体。** 本项目不含任何基于影像灰度的内容：所有几何都来自二值掩膜。

### 被否决、未使用的来源
- Gibson et al., Zenodo `10.5281/zenodo.1169361` —— 条款叠加 TCIA Pancreas-CT 与 BTCV，未用。
- 任何含 NC（非商业）条款的网格 —— 与本仓 MIT 许可不兼容，未用。
- BodyParts3D、3D-IRCADb-01 —— 未用。

### v2.0.1 标签清单的核实结果
该发布每位受试者有 117 个分割文件（已逐一列出核实）。**确认存在** `stomach`、
`esophagus`、`duodenum`、`liver`、`spleen`、`pancreas`、`aorta`、
`portal_vein_and_splenic_vein`、`small_bowel`、`colon`、`inferior_vena_cava`、
`vertebrae_T12`、`vertebrae_L1`。**确认不存在**：腹腔干、任何一条胃动脉或胃静脉、
胃壁/胃腔分层、胃黏膜、胃腺、贲门/胃底/胃体/胃窦/幽门亚标签。
这份发布的压缩包只含 total 任务的标签；头颈腺体与腔隙等附加任务（如 `head_glands_cavities`）
**不在**这个 zip 里，因此本器官无法从本数据集取得任何口腔或咽部结构。

## 2. 本器官自己的受试者筛查

判据 6 项，结果 **通过**：

1. stomach label present and is one dominant connected component (largest component >= 95% of label volume)
2. label does not touch the field-of-view edge
3. label volume between the fasted lower bound and the published capacity upper bound (25-1500 mL)
4. a measurable oesophago-gastric contact band exists (cardia anchor)
5. a measurable pyloro-duodenal contact band exists (pylorus anchor)
6. both great-vessel anchors present: aorta and portal_vein_and_splenic_vein labels non-empty

| 实测量 | 值 |
|---|---|
| `connected_components` | 2 |
| `component_voxel_sizes` | [128386, 4] |
| `largest_component_fraction` | 0.999969 |
| `touches_fov_edge` | False |
| `label_volume_mL` | 433.3 |
| `egj_contact_voxels` | 84 |
| `pyloroduodenal_contact_voxels` | 6 |
| `aorta_voxels` | 35167 |
| `portal_splenic_vein_voxels` | 8485 |
| `vertebra_T12_voxels` | 20790 |

第二个连通分量只有 4 个体素（最大分量占比 0.999969），照实报出而不是四舍五入成 1。
幽门–十二指肠接触带只有 6 个体素，**很薄**，照实报出——这是本页幽门锚点的全部实测依据。

## 3. 逐项 contour_anchored

| 结构 | contour_anchored | 依据 |
|---|---|---|
| 胃浆膜轮廓 | **true** | `stomach` 标签本身 |
| 食管–胃连接（贲门口） | **true** | `stomach` ∩ 膨胀后的 `esophagus`，取接触带质心 |
| 幽门–十二指肠连接 | **true** | `stomach` ∩ 膨胀后的 `duodenum` |
| 胃底边界 | **true** | 经实测贲门口的水平面（标准解剖学定义，施加在本受试者自己的标签上） |
| 角切迹（胃体/胃窦界） | **true** | C 平面内小弯迹线曲率峰 κ=0.37788 rad/mm，显著度 13.89× 中位值，越过阈值 |
| 中心线 | **true** | 胃标签内的测地中轴路径（胃底顶点→幽门） |
| 大弯/小弯迹线 | **true** | 逐站极值表面体素 + C 平面轮廓弧（两套定义都报） |
| 腹腔干**起始水平** | **true** | 本受试者自己的 `vertebrae_T12` 标签 |
| 腹腔干**血管本身** | false | 源数据无此标签，构造 |
| 脾门（上/下极） | **true** | `spleen` 标签表面距胃最近的带内上下极体素 |
| 门静脉汇合开口 | **true** | `portal_vein_and_splenic_vein` 标签表面 |
| 肝门/幽门侧肝动脉开口 | **true** | `liver` 标签距实测幽门最近的体素 |
| 胃十二指肠动脉开口 | **true** | `duodenum` 标签距实测幽门最近的体素 |
| 胃左/右动脉 | false | 源数据无标签，CCO 构造 |
| 胃网膜左/右动脉 | false | 源数据无标签，CCO 构造 |
| 胃短动脉 | false | 源数据无标签，CCO 构造 |
| 胃静脉 | false | 源数据无标签，CCO 构造（逆行） |
| 胃壁/胃腔分界 | false | 单一实心标签，按公开充盈 CT 壁厚 4.0 mm 腐蚀 |
| 贲门领宽度 | false | 声明 20.0 mm 带宽 |
| 幽门管长度 | false | 声明 20.0 mm |
| 胃黏膜皱襞 | false | 12 条 / 2.5 mm，构造 |
| 胃腺与黏液屏障 | false | 显微尺度，远低于 1.5 mm 体素；按公开形态计量构造 |
| 腔内 pH 区域值 | false | 公开文献值，非本受试者实测 |
| 泌酸腺密度分布 | false | 教科书区域组织学，构造 |

## 4. 实测尺寸对参考范围

| 量 | 实测 | 参考 | 判定 |
|---|---|---|---|
| 胃标签体积 | 433.3 mL | 空腹 25.0–50.0 mL；容量 1000.0–1500.0 mL | **两区间之间**，CT 未声明禁食状态，两者都不作为判据 |
| 胃窦中心线长 | 72.5 mm | 30.0–40.0 mm | **超出范围（未调参）** |
| 大弯:小弯（轮廓弧定义） | 2.31 | 教科书定性：大弯显著长于小弯 | 方向一致 |
| 大弯:小弯（逐站迹线定义） | 1.12 | — | 结构上必然接近 1，**不能**复现教科书比值，已说明 |
| 三角面 | 156,492 | 130,000–160,000 | **落在预算内** |
| glb 大小 | 5.09 MB | ≤ 7 MB | **通过** |
| 胃–脾最小表面间隙 | 3.35 mm | — | 实测 |

## 5. 坐标系与跨器官一致性

场景坐标系 = （患者左、上、前），单位 mm，与本项目心脏、肝脏、肾脏、小肠完全一致；
`ras_to_scene` = [[-1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]]。
五个器官共用受试者 `s0236`，所以幽门→十二指肠、
胃静脉→门静脉这两处衔接在同一坐标系下直接对得上。
