# 数据署名 / Dataset attribution

上游单细胞数据来自 **CZ CELLxGENE Discover**。其数据提交政策规定公开发布的数据一律以 **CC-BY 4.0** 提供——
可自由下载、分享与使用，唯一义务是**向原始数据贡献者署名**。本文件就是为履行该义务而生成的。

三维几何来自 **TotalSegmentator v2.0.1**（Wasserthal et al.，DOI 10.5281/zenodo.10047292，CC-BY-4.0）。

Census 版本 `2025-11-08`，共 **130 个不同数据集**，全部已对上其 collection 的引文与 DOI。

## 一处必须讲明的口径差异

本表按各器官的 `tissue_general` 过滤重新枚举数据集。其中 **7 个器官的细胞数与已交付的表达层完全一致**；
**肾脏不一致**（表达层 446,561 对本表 1,463,942）——肾脏的表达层是本项目最早一批、
走的是与后续器官不同的采集路线，覆盖的是该过滤条件下的一个子集。因此肾脏一行应读作
**「该器官过滤条件下的全部贡献数据集」**，而非「表达层确切用到的那些」。其余器官两者相同。

## 逐器官数据集

### oral_cavity

（本表未覆盖：该器官按组织词枚举时无对应 `tissue_general` 条目）

### esophagus — 6 个数据集 / 178,410 细胞

| 数据集 | 细胞数 | 供体 | 研究 | DOI |
|---|---:|---:|---|---|
| `b07fb54c` | 87,947 | 6 | scRNA-seq assessment of the human lung, spleen, and esophagus tissue s | [10.1186/s13059-019-1906-x](https://doi.org/10.1186/s13059-019-1906-x) |
| `4ed927e9` | 34,173 | 3 | Single-nucleus cross-tissue molecular reference maps to decipher disea | [10.1126/science.abl4290](https://doi.org/10.1126/science.abl4290) |
| `7bb64315` | 32,700 | 15 | Single-cell RNA sequencing unifies developmental programs of Esophagea | [10.1158/2159-8290.cd-22-0824](https://doi.org/10.1158/2159-8290.cd-22-0824) |
| `2adb1f8a` | 11,364 | 2 | Construction of a human cell landscape at single-cell level | [10.1038/s41586-020-2157-4](https://doi.org/10.1038/s41586-020-2157-4) |
| `19053a82` | 7,510 | 1 | Single-cell integration reveals metaplasia in inflammatory gut disease | [10.1038/s41586-024-07571-1](https://doi.org/10.1038/s41586-024-07571-1) |
| `9968be68` | 4,716 | 2 | Charting human development using a multi-endodermal organ atlas and or | [10.1016/j.cell.2021.04.028](https://doi.org/10.1016/j.cell.2021.04.028) |

### stomach — 8 个数据集 / 557,791 细胞

| 数据集 | 细胞数 | 供体 | 研究 | DOI |
|---|---:|---:|---|---|
| `19053a82` | 237,365 | 61 | Single-cell integration reveals metaplasia in inflammatory gut disease | [10.1038/s41586-024-07571-1](https://doi.org/10.1038/s41586-024-07571-1) |
| `7bb64315` | 229,962 | 59 | Single-cell RNA sequencing unifies developmental programs of Esophagea | [10.1158/2159-8290.cd-22-0824](https://doi.org/10.1158/2159-8290.cd-22-0824) |
| `53d208b0` | 33,064 | 3 | Tabula Sapiens | [10.1126/science.abl4896](https://doi.org/10.1126/science.abl4896) |
| `2adb1f8a` | 22,577 | 5 | Construction of a human cell landscape at single-cell level | [10.1038/s41586-020-2157-4](https://doi.org/10.1038/s41586-020-2157-4) |
| `9968be68` | 13,632 | 3 | Charting human development using a multi-endodermal organ atlas and or | [10.1016/j.cell.2021.04.028](https://doi.org/10.1016/j.cell.2021.04.028) |
| `f7c1c579` | 12,106 | 4 | A human cell atlas of fetal gene expression | [10.1126/science.aba7721](https://doi.org/10.1126/science.aba7721) |
| `6725ee8e` | 7,600 | 5 | Integrated human endoderm-derived organoids cell atlas (HEOCA) | [10.1038/s41588-025-02182-6](https://doi.org/10.1038/s41588-025-02182-6) |
| `a77176ec` | 1,485 | 1 | Transcriptome Analysis of Archived Tumor Tissues by Visium, GeoMx DSP, | [10.1101/2024.11.01.621259](https://doi.org/10.1101/2024.11.01.621259) |

### small_intestine — 22 个数据集 / 1,434,955 细胞

| 数据集 | 细胞数 | 供体 | 研究 | DOI |
|---|---:|---:|---|---|
| `a37f857c` | 201,072 | 42 | The landscape of immune dysregulation in Crohn’s disease revealed thro | [10.1016/j.immuni.2023.01.002](https://doi.org/10.1016/j.immuni.2023.01.002) |
| `1b350d0a` | 195,115 | 21 | Multimodal profiling reveals tissue-directed signatures of human immun | [10.1038/s41590-025-02241-4](https://doi.org/10.1038/s41590-025-02241-4) |
| `fe4b89d5` | 154,136 | 42 | The landscape of immune dysregulation in Crohn’s disease revealed thro | [10.1016/j.immuni.2023.01.002](https://doi.org/10.1016/j.immuni.2023.01.002) |
| `80a2c5b6` | 144,292 | 14 | Cells of the human intestinal tract mapped across space and time | [10.1038/s41586-021-03852-1](https://doi.org/10.1038/s41586-021-03852-1) |
| `19053a82` | 127,362 | 29 | Single-cell integration reveals metaplasia in inflammatory gut disease | [10.1038/s41586-024-07571-1](https://doi.org/10.1038/s41586-024-07571-1) |
| `6725ee8e` | 125,887 | 48 | Integrated human endoderm-derived organoids cell atlas (HEOCA) | [10.1038/s41588-025-02182-6](https://doi.org/10.1038/s41588-025-02182-6) |
| `6ceeaa86` | 89,849 | 13 | Pediatric FGID lamina propria of ileum biopsies. | [10.7554/eLife.91792](https://doi.org/10.7554/eLife.91792) |
| `0f4865d5` | 75,695 | 40 | The landscape of immune dysregulation in Crohn’s disease revealed thro | [10.1016/j.immuni.2023.01.002](https://doi.org/10.1016/j.immuni.2023.01.002) |
| `9968be68` | 64,619 | 9 | Charting human development using a multi-endodermal organ atlas and or | [10.1016/j.cell.2021.04.028](https://doi.org/10.1016/j.cell.2021.04.028) |
| `53d208b0` | 42,036 | 5 | Tabula Sapiens | [10.1126/science.abl4896](https://doi.org/10.1126/science.abl4896) |
| `b46237d1` | 41,702 | 9 | Single-Cell Sequencing of Developing Human Gut Reveals Transcriptional | [10.1016/j.devcel.2020.11.010](https://doi.org/10.1016/j.devcel.2020.11.010) |
| `3f32121d` | 36,359 | 6 | Insulin is expressed by enteroendocrine cells during human fetal devel | [10.1038/s41591-021-01586-1](https://doi.org/10.1038/s41591-021-01586-1) |
| `1b9d8702` | 28,741 | 10 | Cross-tissue immune cell analysis reveals tissue-specific features in  | [10.1126/science.abl5197](https://doi.org/10.1126/science.abl5197) |
| `9d5df009` | 24,787 | 13 | Spatiotemporal analysis of human intestinal development at single-cell | [10.1016/j.cell.2020.12.016](https://doi.org/10.1016/j.cell.2020.12.016) |
| `8e47ed12` | 22,502 | 15 | Single-Cell Sequencing of Developing Human Gut Reveals Transcriptional | [10.1016/j.devcel.2020.11.010](https://doi.org/10.1016/j.devcel.2020.11.010) |
| `b6ee0939` | 19,200 | 1 | Epithelial zonation along the mammalian small intestine defines five d | [10.1038/s41556-023-01337-z](https://doi.org/10.1038/s41556-023-01337-z) |
| `2adb1f8a` | 13,597 | 3 | Construction of a human cell landscape at single-cell level | [10.1038/s41586-020-2157-4](https://doi.org/10.1038/s41586-020-2157-4) |
| `019c7af2` | 7,336 | 3 | A proximal-to-distal survey of healthy adult human small intestine and | [10.1016/j.jcmgh.2022.02.007](https://doi.org/10.1016/j.jcmgh.2022.02.007) |
| `248cf7bf` | 5,980 | 2 | Single-cell transcriptome analysis reveals differential nutrient absor | [10.1084/jem.20191130](https://doi.org/10.1084/jem.20191130) |
| `a77176ec` | 5,335 | 1 | Transcriptome Analysis of Archived Tumor Tissues by Visium, GeoMx DSP, | [10.1101/2024.11.01.621259](https://doi.org/10.1101/2024.11.01.621259) |
| `ee195b7d` | 5,200 | 2 | Charting human development using a multi-endodermal organ atlas and or | [10.1016/j.cell.2021.04.028](https://doi.org/10.1016/j.cell.2021.04.028) |
| `7bb64315` | 4,153 | 4 | Single-cell RNA sequencing unifies developmental programs of Esophagea | [10.1158/2159-8290.cd-22-0824](https://doi.org/10.1158/2159-8290.cd-22-0824) |

### colon — 25 个数据集 / 1,365,592 细胞

| 数据集 | 细胞数 | 供体 | 研究 | DOI |
|---|---:|---:|---|---|
| `16023185` | 370,115 | 62 | Spatially organized multicellular immune hubs in human colorectal canc | [10.1016/j.cell.2021.08.003](https://doi.org/10.1016/j.cell.2021.08.003) |
| `19053a82` | 203,644 | 47 | Single-cell integration reveals metaplasia in inflammatory gut disease | [10.1038/s41586-024-07571-1](https://doi.org/10.1038/s41586-024-07571-1) |
| `80a2c5b6` | 111,645 | 13 | Cells of the human intestinal tract mapped across space and time | [10.1038/s41586-021-03852-1](https://doi.org/10.1038/s41586-021-03852-1) |
| `518d9049` | 101,164 | 22 | The landscape of immune dysregulation in Crohn’s disease revealed thro | [10.1016/j.immuni.2023.01.002](https://doi.org/10.1016/j.immuni.2023.01.002) |
| `5ce42b38` | 67,815 | 22 | The landscape of immune dysregulation in Crohn’s disease revealed thro | [10.1016/j.immuni.2023.01.002](https://doi.org/10.1016/j.immuni.2023.01.002) |
| `e40c6272` | 58,644 | 26 | HTAN VUMC - Differential pre-malignant programs and microenvironment c | [10.1016/j.cell.2021.11.031](https://doi.org/10.1016/j.cell.2021.11.031) |
| `6725ee8e` | 53,217 | 12 | Integrated human endoderm-derived organoids cell atlas (HEOCA) | [10.1038/s41588-025-02182-6](https://doi.org/10.1038/s41588-025-02182-6) |
| `d6dfdef1` | 50,280 | 23 | HTAN VUMC - Differential pre-malignant programs and microenvironment c | [10.1016/j.cell.2021.11.031](https://doi.org/10.1016/j.cell.2021.11.031) |
| `9d5df009` | 49,095 | 12 | Spatiotemporal analysis of human intestinal development at single-cell | [10.1016/j.cell.2020.12.016](https://doi.org/10.1016/j.cell.2020.12.016) |
| `829a3cd1` | 34,955 | 28 | Progressive plasticity during colorectal cancer metastasis | [10.1038/s41586-024-08150-0](https://doi.org/10.1038/s41586-024-08150-0) |
| `4dd00779` | 34,772 | 12 | Intra- and Inter-cellular Rewiring of the Human Colon during Ulcerativ | [10.1016/j.cell.2019.06.029](https://doi.org/10.1016/j.cell.2019.06.029) |
| `2adb1f8a` | 32,505 | 3 | Construction of a human cell landscape at single-cell level | [10.1038/s41586-020-2157-4](https://doi.org/10.1038/s41586-020-2157-4) |
| `6cf3634d` | 30,814 | 22 | The landscape of immune dysregulation in Crohn’s disease revealed thro | [10.1016/j.immuni.2023.01.002](https://doi.org/10.1016/j.immuni.2023.01.002) |
| `53d208b0` | 30,084 | 5 | Tabula Sapiens | [10.1126/science.abl4896](https://doi.org/10.1126/science.abl4896) |
| `b252b015` | 29,443 | 4 | MSK SPECTRUM – Ovarian cancer mutational processes drive site-specific | [10.1038/s41586-022-05496-1](https://doi.org/10.1038/s41586-022-05496-1) |
| `2872f4b0` | 27,692 | 5 | Distinct microbial and immune niches of the human colon | [10.1038/s41590-020-0602-z](https://doi.org/10.1038/s41590-020-0602-z) |
| `b46237d1` | 21,147 | 9 | Single-Cell Sequencing of Developing Human Gut Reveals Transcriptional | [10.1016/j.devcel.2020.11.010](https://doi.org/10.1016/j.devcel.2020.11.010) |
| `b617ee1b` | 18,837 | 7 | Dissecting novel myeloid-derived cell states through single-cell RNA-S | [10.1038/s41467-024-49916-4](https://doi.org/10.1038/s41467-024-49916-4) |
| `9968be68` | 13,053 | 3 | Charting human development using a multi-endodermal organ atlas and or | [10.1016/j.cell.2021.04.028](https://doi.org/10.1016/j.cell.2021.04.028) |
| `6a270451` | 7,976 | 52 | HTAN VUMC - Differential pre-malignant programs and microenvironment c | [10.1016/j.cell.2021.11.031](https://doi.org/10.1016/j.cell.2021.11.031) |
| `019c7af2` | 5,254 | 3 | A proximal-to-distal survey of healthy adult human small intestine and | [10.1016/j.jcmgh.2022.02.007](https://doi.org/10.1016/j.jcmgh.2022.02.007) |
| `1b350d0a` | 4,487 | 1 | Multimodal profiling reveals tissue-directed signatures of human immun | [10.1038/s41590-025-02241-4](https://doi.org/10.1038/s41590-025-02241-4) |
| `4e74d843` | 4,329 | 2 | Single-cell transcriptome analysis reveals differential nutrient absor | [10.1084/jem.20191130](https://doi.org/10.1084/jem.20191130) |
| `80d7111d` | 3,797 | 2 | Single-cell transcriptome analysis reveals differential nutrient absor | [10.1084/jem.20191130](https://doi.org/10.1084/jem.20191130) |
| `1b9d8702` | 828 | 4 | Cross-tissue immune cell analysis reveals tissue-specific features in  | [10.1126/science.abl5197](https://doi.org/10.1126/science.abl5197) |

### liver — 33 个数据集 / 1,625,338 细胞

| 数据集 | 细胞数 | 供体 | 研究 | DOI |
|---|---:|---:|---|---|
| `34f5307e` | 250,739 | 19 | HTAN/HTAPP Broad - Spatio-molecular dissection of the breast cancer me | [10.1038/s41591-024-03215-z](https://doi.org/10.1038/s41591-024-03215-z) |
| `fd072bc3` | 219,296 | 14 | Mapping the developing human immune system across organs | [10.1126/science.abo0510](https://doi.org/10.1126/science.abo0510) |
| `e84f2780` | 167,598 | 19 | Spatial proteogenomics reveals distinct and evolutionarily conserved h | [10.1016/j.cell.2021.12.018](https://doi.org/10.1016/j.cell.2021.12.018) |
| `f7c1c579` | 113,138 | 9 | A human cell atlas of fetal gene expression | [10.1126/science.aba7721](https://doi.org/10.1126/science.aba7721) |
| `6c87755e` | 110,476 | 16 | HTAN/HTAPP Broad - Spatio-molecular dissection of the breast cancer me | [10.1038/s41591-024-03215-z](https://doi.org/10.1038/s41591-024-03215-z) |
| `1873a18a` | 79,265 | 9 | Single-cell and spatial transcriptomics characterisation of the immuno | [10.1016/j.jhep.2023.12.023](https://doi.org/10.1016/j.jhep.2023.12.023) |
| `e3ed2ba4` | 75,104 | 6 | Ex vivo modelling of human colorectal cancer liver metastasis by normo | [10.1186/s12943-025-02430-7](https://doi.org/10.1186/s12943-025-02430-7) |
| `6725ee8e` | 73,899 | 23 | Integrated human endoderm-derived organoids cell atlas (HEOCA) | [10.1038/s41588-025-02182-6](https://doi.org/10.1038/s41588-025-02182-6) |
| `ddb22b3d` | 73,295 | 4 | Single-Cell, Single-Nucleus, and Spatial RNA Sequencing of the Human L | [10.1002/hep4.1854](https://doi.org/10.1002/hep4.1854) |
| `84cfa5aa` | 72,710 | 15 | Single-Cell Atlas Of Human Pediatric Liver Reveals Age-Related Hepatic | [10.1097/HC9.0000000000000813](https://doi.org/10.1097/HC9.0000000000000813) |
| `ec6ceff8` | 67,110 | 3 | Identification of distinct tumor cell populations and key genetic mech | [10.1038/s42003-021-02562-8](https://doi.org/10.1038/s42003-021-02562-8) |
| `02792605` | 65,396 | 10 | Single-cell and spatial transcriptomics characterisation of the immuno | [10.1016/j.jhep.2023.12.023](https://doi.org/10.1016/j.jhep.2023.12.023) |
| `234c23d0` | 45,186 | 23 | Single-cell and spatial transcriptomics characterisation of the immuno | [10.1016/j.jhep.2023.12.023](https://doi.org/10.1016/j.jhep.2023.12.023) |
| `2adb1f8a` | 28,501 | 4 | Construction of a human cell landscape at single-cell level | [10.1038/s41586-020-2157-4](https://doi.org/10.1038/s41586-020-2157-4) |
| `afae4675` | 24,293 | 2 | HTAN CHOP - Longitudinal single-cell multiomic atlas of high-risk neur | [10.1038/s41588-025-02158-6](https://doi.org/10.1038/s41588-025-02158-6) |
| `53d208b0` | 22,214 | 4 | Tabula Sapiens | [10.1126/science.abl4896](https://doi.org/10.1126/science.abl4896) |
| `471647b3` | 16,099 | 3 | Blood and immune development in human fetal bone marrow and Down syndr | [10.1038/s41586-021-03929-x](https://doi.org/10.1038/s41586-021-03929-x) |
| `43560d7c` | 15,216 | 23 | Single-cell and spatial transcriptomics characterisation of the immuno | [10.1016/j.jhep.2023.12.023](https://doi.org/10.1016/j.jhep.2023.12.023) |
| `1b9d8702` | 13,368 | 7 | Cross-tissue immune cell analysis reveals tissue-specific features in  | [10.1126/science.abl5197](https://doi.org/10.1126/science.abl5197) |
| `76af615b` | 11,982 | 23 | Single-cell and spatial transcriptomics characterisation of the immuno | [10.1016/j.jhep.2023.12.023](https://doi.org/10.1016/j.jhep.2023.12.023) |
| `1d89d081` | 10,352 | 9 | Automatic cell-type harmonization and integration across Human Cell At | [10.1016/j.cell.2023.11.026](https://doi.org/10.1016/j.cell.2023.11.026) |
| `829a3cd1` | 10,332 | 25 | Progressive plasticity during colorectal cancer metastasis | [10.1038/s41586-024-08150-0](https://doi.org/10.1038/s41586-024-08150-0) |
| `576f193c` | 9,757 | 3 | HTAN MSK - Single cell profiling reveals novel tumor and myeloid subpo | [10.1016/j.ccell.2021.09.008](https://doi.org/10.1016/j.ccell.2021.09.008) |
| `b550fd94` | 9,515 | 23 | Single-cell and spatial transcriptomics characterisation of the immuno | [10.1016/j.jhep.2023.12.023](https://doi.org/10.1016/j.jhep.2023.12.023) |
| `b252b015` | 9,077 | 1 | MSK SPECTRUM – Ovarian cancer mutational processes drive site-specific | [10.1038/s41586-022-05496-1](https://doi.org/10.1038/s41586-022-05496-1) |
| `4ce8d8ec` | 8,219 | 23 | Single-cell and spatial transcriptomics characterisation of the immuno | [10.1016/j.jhep.2023.12.023](https://doi.org/10.1016/j.jhep.2023.12.023) |
| `1b350d0a` | 6,859 | 3 | Multimodal profiling reveals tissue-directed signatures of human immun | [10.1038/s41590-025-02241-4](https://doi.org/10.1038/s41590-025-02241-4) |
| `b617ee1b` | 6,661 | 19 | Dissecting novel myeloid-derived cell states through single-cell RNA-S | [10.1038/s41467-024-49916-4](https://doi.org/10.1038/s41467-024-49916-4) |
| `9968be68` | 4,075 | 1 | Charting human development using a multi-endodermal organ atlas and or | [10.1016/j.cell.2021.04.028](https://doi.org/10.1016/j.cell.2021.04.028) |
| `1e6a6ef9` | 2,366 | 4 | High-resolution single-cell atlas reveals diversity and plasticity of  | [10.1016/j.ccell.2022.10.008](https://doi.org/10.1016/j.ccell.2022.10.008) |
| `e347396c` | 1,241 | 22 | Single-cell and spatial transcriptomics characterisation of the immuno | [10.1016/j.jhep.2023.12.023](https://doi.org/10.1016/j.jhep.2023.12.023) |
| `4993d61c` | 1,167 | 23 | Single-cell and spatial transcriptomics characterisation of the immuno | [10.1016/j.jhep.2023.12.023](https://doi.org/10.1016/j.jhep.2023.12.023) |
| `0c9a8cfb` | 832 | 23 | Single-cell and spatial transcriptomics characterisation of the immuno | [10.1016/j.jhep.2023.12.023](https://doi.org/10.1016/j.jhep.2023.12.023) |

### kidney — 24 个数据集 / 1,463,942 细胞

| 数据集 | 细胞数 | 供体 | 研究 | DOI |
|---|---:|---:|---|---|
| `a12ccb9b` | 304,989 | 80 | An atlas of healthy and injured cell states and niches in the human ki | — |
| `dea717d4` | 225,177 | 77 | An atlas of healthy and injured cell states and niches in the human ki | — |
| `5af90777` | 192,956 | 12 | Mapping single-cell transcriptomes in the intra-tumoral and associated | [10.1016/j.ccell.2022.11.001](https://doi.org/10.1016/j.ccell.2022.11.001) |
| `f7c1c579` | 155,386 | 9 | A human cell atlas of fetal gene expression | [10.1126/science.aba7721](https://doi.org/10.1126/science.aba7721) |
| `bd65a70f` | 131,902 | 6 | Single-cell sequencing links multiregional immune landscapes and tissu | [10.1016/j.ccell.2021.03.007](https://doi.org/10.1016/j.ccell.2021.03.007) |
| `0b4a15a7` | 48,783 | 18 | HCA kidney seed network: University of Michigan | [10.1172/jci.insight.204385](https://doi.org/10.1172/jci.insight.204385) |
| `867757c1` | 46,957 | 12 | snRNA-seq of human kidney with ureteral obstruction | [10.1038/s41467-025-59997-4](https://doi.org/10.1038/s41467-025-59997-4) |
| `2adb1f8a` | 45,630 | 7 | Construction of a human cell landscape at single-cell level | [10.1038/s41586-020-2157-4](https://doi.org/10.1038/s41586-020-2157-4) |
| `9ea768a2` | 40,268 | 13 | Spatiotemporal immune zonation of the human kidney | [10.1126/science.aat5031](https://doi.org/10.1126/science.aat5031) |
| `09b518f9` | 37,717 | 10 | Multimodal benchmarking dataset for renal cortex characterization (mBD | [10.1101/2025.03.06.637075](https://doi.org/10.1101/2025.03.06.637075) |
| `7b20c613` | 29,993 | 7 | Integrated cancer cell-specific single-cell RNA-seq datasets of immune | [10.1038/s41597-025-04381-6](https://doi.org/10.1038/s41597-025-04381-6) |
| `0bae7ebf` | 27,675 | 19 | Healthy living donor kidney | [10.1038/s41467-022-35297-z](https://doi.org/10.1038/s41467-022-35297-z) |
| `0b75c598` | 27,491 | 7 | An atlas of healthy and injured cell states and niches in the human ki | [10.1038/s41586-023-05769-3](https://doi.org/10.1038/s41586-023-05769-3) |
| `d7dcfd8f` | 27,197 | 6 | Spatiotemporal immune zonation of the human kidney | [10.1126/science.aat5031](https://doi.org/10.1126/science.aat5031) |
| `be39785b` | 20,509 | 7 | Single-cell analyses of renal cell cancers reveal insights into tumor  | [10.1073/pnas.2103240118](https://doi.org/10.1073/pnas.2103240118) |
| `e067e5ca` | 20,215 | 6 | Multimodal single cell sequencing implicates chromatin accessibility a | [10.1038/s41467-022-32972-z](https://doi.org/10.1038/s41467-022-32972-z) |
| `9df60c57` | 19,985 | 5 | Single cell transcriptional and chromatin accessibility profiling rede | [10.1038/s41467-021-22368-w](https://doi.org/10.1038/s41467-021-22368-w) |
| `8dc29642` | 17,962 | 5 | Single-cell transcriptomics reveals gene expression dynamics of human  | [10.1371/journal.pbio.3000152](https://doi.org/10.1371/journal.pbio.3000152) |
| `53d208b0` | 11,376 | 1 | Tabula Sapiens | [10.1126/science.abl4896](https://doi.org/10.1126/science.abl4896) |
| `fd072bc3` | 10,705 | 3 | Mapping the developing human immune system across organs | [10.1126/science.abo0510](https://doi.org/10.1126/science.abo0510) |
| `60a29d0b` | 7,813 | 4 | Multimodal benchmarking dataset for renal cortex characterization (mBD | [10.1101/2025.03.06.637075](https://doi.org/10.1101/2025.03.06.637075) |
| `f801b7a9` | 6,044 | 5 | Single-cell analyses of renal cell cancers reveal insights into tumor  | [10.1073/pnas.2103240118](https://doi.org/10.1073/pnas.2103240118) |
| `ea01c125` | 4,636 | 1 | Single cell derived mRNA signals across human kidney tumors | [10.1038/s41467-021-23949-5](https://doi.org/10.1038/s41467-021-23949-5) |
| `4c6f9f26` | 2,576 | 1 | Single-cell analyses of renal cell cancers reveal insights into tumor  | [10.1073/pnas.2103240118](https://doi.org/10.1073/pnas.2103240118) |

### heart — 20 个数据集 / 3,104,821 细胞

| 数据集 | 细胞数 | 供体 | 研究 | DOI |
|---|---:|---:|---|---|
| `65badd7a` | 665,955 | 67 | Pathogenic variants damage cell composition and single cell transcript | [10.1126/science.abo1984](https://doi.org/10.1126/science.abo1984) |
| `d4e69e01` | 486,134 | 14 | Cells of the adult human heart | [10.1038/s41586-020-2797-4](https://doi.org/10.1038/s41586-020-2797-4) |
| `364bd0c7` | 437,851 | 32 | Automatic cell-type harmonization and integration across Human Cell At | [10.1016/j.cell.2023.11.026](https://doi.org/10.1016/j.cell.2023.11.026) |
| `f1606894` | 392,819 | 54 | Transcriptomic and spatial dissection of human ex vivo right atrial ti | [10.1016/j.xcrm.2024.101556](https://doi.org/10.1016/j.xcrm.2024.101556) |
| `8f4f8502` | 221,756 | 47 | Transcriptomic and spatial dissection of human ex vivo right atrial ti | [10.1016/j.xcrm.2024.101556](https://doi.org/10.1016/j.xcrm.2024.101556) |
| `d567b692` | 211,061 | 11 | Spatially resolved multiomics of human cardiac niches | [10.1038/s41586-023-06311-1](https://doi.org/10.1038/s41586-023-06311-1) |
| `1c739a3e` | 191,795 | 20 | Spatial multi-omic map of human myocardial infarction | [10.1038/s41586-022-05060-x](https://doi.org/10.1038/s41586-022-05060-x) |
| `f7c1c579` | 101,749 | 9 | A human cell atlas of fetal gene expression | [10.1126/science.aba7721](https://doi.org/10.1126/science.aba7721) |
| `fe7aae33` | 62,705 | 24 | The cellular and molecular cardiac tissue responses in human inflammat | [10.1038/s44161-025-00612-6](https://doi.org/10.1038/s44161-025-00612-6) |
| `99e317c0` | 51,176 | 9 | Sex-Specific Control of Human Heart Maturation by the Progesterone Rec | [10.1161/CIRCULATIONAHA.120.051921](https://doi.org/10.1161/CIRCULATIONAHA.120.051921) |
| `2e9d2f32` | 49,359 | 3 | Single-cell genomics improves the discovery of risk variants and genes | [10.1038/s41467-023-40505-5](https://doi.org/10.1038/s41467-023-40505-5) |
| `87a06a2a` | 41,279 | 12 | Recruited macrophages elicit atrial fibrillation | [10.1126/science.abq3061](https://doi.org/10.1126/science.abq3061) |
| `53d208b0` | 40,006 | 6 | Tabula Sapiens | [10.1126/science.abl4896](https://doi.org/10.1126/science.abl4896) |
| `4ed927e9` | 36,574 | 3 | Single-nucleus cross-tissue molecular reference maps to decipher disea | [10.1126/science.abl4290](https://doi.org/10.1126/science.abl4290) |
| `72955cdb` | 35,395 | 16 | Plaque Atlas | [10.1038/s41467-025-63202-x](https://doi.org/10.1038/s41467-025-63202-x) |
| `5500c673` | 30,889 | 7 | Integrated adult and foetal heart single-cell RNA sequencing | [10.1038/s44161-022-00183-w](https://doi.org/10.1038/s44161-022-00183-w) |
| `f15e263b` | 19,722 | 3 | Spatial multi-omic map of human myocardial infarction | [10.1038/s41586-022-05060-x](https://doi.org/10.1038/s41586-022-05060-x) |
| `43245158` | 15,842 | 6 | A Coding and Non-Coding Atlas of the Human Arterial Cell | [10.1016/j.xgen.2025.101034](https://doi.org/10.1016/j.xgen.2025.101034) |
| `2adb1f8a` | 10,783 | 4 | Construction of a human cell landscape at single-cell level | [10.1038/s41586-020-2157-4](https://doi.org/10.1038/s41586-020-2157-4) |
| `f3ee7613` | 1,971 | 9 | High-resolution single-cell transcriptomic survey of cardiomyocytes fr | [10.1111/cpr.13557](https://doi.org/10.1111/cpr.13557) |

### lung — 45 个数据集 / 6,167,731 细胞

| 数据集 | 细胞数 | 供体 | 研究 | DOI |
|---|---:|---:|---|---|
| `9f222629` | 1,525,043 | 360 | The integrated Human Lung Cell Atlas | [10.1038/s41591-023-02327-2](https://doi.org/10.1038/s41591-023-02327-2) |
| `1e6a6ef9` | 823,285 | 163 | High-resolution single-cell atlas reveals diversity and plasticity of  | [10.1016/j.ccell.2022.10.008](https://doi.org/10.1016/j.ccell.2022.10.008) |
| `350237e0` | 654,947 | 44 | Early human lung immune cell development and its role in epithelial ce | [10.1126/sciimmunol.adf9988](https://doi.org/10.1126/sciimmunol.adf9988) |
| `f14bc322` | 440,854 | 109 | Single-cell RNA-seq analysis of Interstitial Lung Disease (ILD) subtyp | [10.1038/s41588-024-01702-0](https://doi.org/10.1038/s41588-024-01702-0) |
| `d68a8b48` | 255,204 | 21 | Bronchopulmonary Dysplasia | — |
| `f7c1c579` | 217,738 | 11 | A human cell atlas of fetal gene expression | [10.1126/science.aba7721](https://doi.org/10.1126/science.aba7721) |
| `1b350d0a` | 213,154 | 19 | Multimodal profiling reveals tissue-directed signatures of human immun | [10.1038/s41590-025-02241-4](https://doi.org/10.1038/s41590-025-02241-4) |
| `6725ee8e` | 142,486 | 40 | Integrated human endoderm-derived organoids cell atlas (HEOCA) | [10.1038/s41588-025-02182-6](https://doi.org/10.1038/s41588-025-02182-6) |
| `7f7faf6b` | 134,076 | 1 | A molecular cell atlas of COVID-19 initiation in human lung | [10.1084/jem.20232192](https://doi.org/10.1084/jem.20232192) |
| `7b55fe5c` | 125,530 | 1 | A molecular cell atlas of COVID-19 initiation in human lung | [10.1084/jem.20232192](https://doi.org/10.1084/jem.20232192) |
| `eb499fd8` | 122,445 | 21 | Human CellCards Multi-Study CellRef 1.0 Atlas | [10.1038/s41467-023-40173-5](https://doi.org/10.1038/s41467-023-40173-5) |
| `093d3bfe` | 122,399 | 9 | A spatially resolved atlas of the human lung characterizes a gland-ass | [10.1038/s41588-022-01243-4](https://doi.org/10.1038/s41588-022-01243-4) |
| `d8da613f` | 116,313 | 27 | A molecular single-cell lung atlas of lethal COVID-19 | [10.1038/s41586-021-03569-1](https://doi.org/10.1038/s41586-021-03569-1) |
| `429d21fd` | 108,613 | 12 | Ischemia Reperfusion Responses in Human Lung Transplants at the Single | [10.1016/j.ajt.2024.08.019](https://doi.org/10.1016/j.ajt.2024.08.019) |
| `576f193c` | 90,384 | 26 | HTAN MSK - Single cell profiling reveals novel tumor and myeloid subpo | [10.1016/j.ccell.2021.09.008](https://doi.org/10.1016/j.ccell.2021.09.008) |
| `d41f45c1` | 82,991 | 23 | HTAN MSK - Transcriptional connectivity of regulatory T cells in the t | [10.1038/s41590-023-01504-2](https://doi.org/10.1038/s41590-023-01504-2) |
| `7c497c66` | 76,298 | 1 | A molecular cell atlas of COVID-19 initiation in human lung | [10.1084/jem.20232192](https://doi.org/10.1084/jem.20232192) |
| `3dc61ca1` | 71,752 | 12 | A human fetal lung cell atlas uncovers proximal-distal gradients of di | [10.1016/j.cell.2022.11.005](https://doi.org/10.1016/j.cell.2022.11.005) |
| `53d208b0` | 65,847 | 4 | Tabula Sapiens | [10.1126/science.abl4896](https://doi.org/10.1126/science.abl4896) |
| `b617ee1b` | 62,209 | 9 | Dissecting novel myeloid-derived cell states through single-cell RNA-S | [10.1038/s41467-024-49916-4](https://doi.org/10.1038/s41467-024-49916-4) |
| `4023a2bc` | 61,399 | 9 | A human fetal lung cell atlas uncovers proximal-distal gradients of di | [10.1016/j.cell.2022.11.005](https://doi.org/10.1016/j.cell.2022.11.005) |
| `8c42cfd0` | 60,993 | 3 | A molecular cell atlas of the human lung from single cell RNA sequenci | [10.1038/s41586-020-2922-4](https://doi.org/10.1038/s41586-020-2922-4) |
| `8fbed309` | 57,918 | 16 | Anomalous Epithelial Variations and Ectopic Inflammatory Response in C | [10.1165/rcmb.2021-0555OC](https://doi.org/10.1165/rcmb.2021-0555OC) |
| `2672b679` | 57,019 | 5 | scRNA-seq assessment of the human lung, spleen, and esophagus tissue s | [10.1186/s13059-019-1906-x](https://doi.org/10.1186/s13059-019-1906-x) |
| `2ac76f1b` | 48,226 | 1 | A molecular cell atlas of COVID-19 initiation in human lung | [10.1084/jem.20232192](https://doi.org/10.1084/jem.20232192) |
| `9968be68` | 47,909 | 3 | Charting human development using a multi-endodermal organ atlas and or | [10.1016/j.cell.2021.04.028](https://doi.org/10.1016/j.cell.2021.04.028) |
| `3de0ad6d` | 46,500 | 9 | LungMAP — Human data from a broad age healthy donor group | [10.7554/eLife.62522](https://doi.org/10.7554/eLife.62522) |
| `2f132ec9` | 39,778 | 12 | SARS-CoV-2 receptor ACE2 and TMPRSS2 are primarily expressed in bronch | [10.15252/embj.20105114](https://doi.org/10.15252/embj.20105114) |
| `1e4214ce` | 35,954 | 4 | Transcriptome Analysis of Archived Tumor Tissues by Visium, GeoMx DSP, | [10.1101/2024.11.01.621259](https://doi.org/10.1101/2024.11.01.621259) |
| `1e5bd3b8` | 35,699 | 6 | Emphysema Cell Atlas | [10.1016/j.immuni.2023.01.032](https://doi.org/10.1016/j.immuni.2023.01.032) |
| `1b9d8702` | 35,419 | 10 | Cross-tissue immune cell analysis reveals tissue-specific features in  | [10.1126/science.abl5197](https://doi.org/10.1126/science.abl5197) |
| `4ed927e9` | 35,284 | 3 | Single-nucleus cross-tissue molecular reference maps to decipher disea | [10.1126/science.abl4290](https://doi.org/10.1126/science.abl4290) |
| `2adb1f8a` | 33,698 | 5 | Construction of a human cell landscape at single-cell level | [10.1038/s41586-020-2157-4](https://doi.org/10.1038/s41586-020-2157-4) |
| `b61a921b` | 21,144 | 5 | A Unique Cellular Organization of Human Distal Airways and Its Disarra | [10.1164/rccm.202207-1384OC](https://doi.org/10.1164/rccm.202207-1384OC) |
| `43a1506b` | 20,475 | 15 | Organoid modeling of human fetal lung alveolar development reveals mec | [10.1016/j.stem.2022.11.013](https://doi.org/10.1016/j.stem.2022.11.013) |
| `4b6af54a` | 18,386 | 6 | Emphysema Cell Atlas | [10.1016/j.immuni.2023.01.032](https://doi.org/10.1016/j.immuni.2023.01.032) |
| `01209dce` | 11,059 | 2 | Single-cell transcriptomics of human T cells reveals tissue and activa | [10.1038/s41467-019-12464-3](https://doi.org/10.1038/s41467-019-12464-3) |
| `781cc817` | 9,619 | 5 | A novel human fetal lung-derived alveolar organoid model reveals mecha | [10.1038/s44318-024-00328-6](https://doi.org/10.1038/s44318-024-00328-6) |
| `e04daea4` | 8,657 | 3 | A molecular cell atlas of the human lung from single cell RNA sequenci | [10.1038/s41586-020-2922-4](https://doi.org/10.1038/s41586-020-2922-4) |
| `8dc858b0` | 8,556 | 1 | Early human lung immune cell development and its role in epithelial ce | [10.1126/sciimmunol.adf9988](https://doi.org/10.1126/sciimmunol.adf9988) |
| `9dbab10c` | 7,441 | 4 | COVID-19 immune features revealed by a large-scale single-cell transcr | [10.1016/j.cell.2021.01.053](https://doi.org/10.1016/j.cell.2021.01.053) |
| `ef1a4b05` | 5,594 | 2 | Transcriptomic Analysis of Air-Liquid Interface Culture in Human Lung  | [10.3390/cells13231991](https://doi.org/10.3390/cells13231991) |
| `6e00ccf7` | 5,499 | 2 | Single cell transcriptomic profiling identifies molecular phenotypes o | [10.3390/genes15030298](https://doi.org/10.3390/genes15030298) |
| `6c87755e` | 2,962 | 1 | HTAN/HTAPP Broad - Spatio-molecular dissection of the breast cancer me | [10.1038/s41591-024-03215-z](https://doi.org/10.1038/s41591-024-03215-z) |
| `829a3cd1` | 975 | 1 | Progressive plasticity during colorectal cancer metastasis | [10.1038/s41586-024-08150-0](https://doi.org/10.1038/s41586-024-08150-0) |
