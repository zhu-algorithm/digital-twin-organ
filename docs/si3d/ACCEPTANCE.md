# Acceptance record — 3-D anatomical small intestine

Every row is a value read back out of a generated artefact, not a claim. The last section lists what was **not** achieved, not measured, or synthesised; nothing in this project is tuned to land inside a published range.

## 1. Geometry — real segmentation mesh

| gate | required | measured | verdict |
|---|---|---|---|
| triangle count of the real-segmentation surfaces | 130 000 – 160 000 | **146,000** | pass |
| glb size | ≤ 7 MB | **2.51 MB** (2,631,324 B) | pass |
| source licence | attribution-only, redistributable | **CC-BY-4.0 attribution-only**, DOI 10.5281/zenodo.10047292 | pass |
| CT image downloaded | no | **no** — only 9 mask members, via HTTP range requests | pass |

Per-surface, all contour-anchored:

| surface | faces | vertices | watertight | volume (mL) |
|---|---|---|---|---|
| `small_bowel` | 72,000 | 35,978 | True | 1406.25 |
| `duodenum` | 14,000 | 7,002 | True | 55.03 |
| `colon` | 22,000 | 11,000 | True | 802.16 |
| `stomach` | 12,000 | 6,002 | True | 425.68 |
| `liver` | 18,000 | 9,000 | True | 1423.86 |
| `pancreas` | 8,000 | 4,002 | True | 75.56 |

### Subject screen (s0236, re-screened for this organ)

| criterion | measured |
|---|---|
| `small_bowel` voxels on any of the six volume faces | **0** |
| `duodenum` voxels on any face | **0** |
| largest connected component of `small_bowel` | **0.99997** of voxels (5 components) |
| `duodenum` connected components | **1** |
| duodenojejunal flexure contact | **302** voxels |
| ileocaecal contact | **8157** voxels |
| envelope volume | **1431.1 mL** |

## 2. Centreline — measured, and what it is not

| quantity | measured | published reference | in range |
|---|---|---|---|
| envelope skeleton longest geodesic | **472.8 mm** | 3000–5000 mm in vivo | **no — and it is not used as the course** |
| duodenal centreline (contour-anchored) | **130.9 mm** | 200–300 mm | **no — reported below range as found** |
| duodenojejunal bridge | **18.5 mm** | — | — |
| jejuno-ileal course (parametric) | **3262.3 mm** | — | — |
| total small intestine | **3.412 m** (check 3411.7 mm) | 3–5 m in vivo | **yes** |
| envelope volume ÷ 20 mm cross-section | **4.56 m** | 3–5 m in vivo | **yes** |
| centreline samples inside the measured envelope | **96.75 %** of 1631 | — | — |

Construction parameters, both measured on this subject and neither tuned: lattice pitch **14.7 mm** (2 x p75 of the envelope inscribed radius), node admission threshold **4.5 mm** (median envelope inscribed radius), 203 nodes, 8 of 8 offending hops re-routed on the in-envelope voxel geodesic.

## 3. Vessel trees — Murray and clearance

| tree | segments | bifurcations | terminals | radius (mm) | total length (mm) | Murray max relative residual | bifurcations checked | contour-anchored |
|---|---|---|---|---|---|---|---|---|
| mesenteric_artery | 1090 | 257 | 258 | 0.485 – 3.091 | 5651.7 | **0.0** | 257 | no |
| mesenteric_vein | 1091 | 259 | 260 | 0.62 – 3.957 | 5657.4 | **0.0** | 259 | no |
| lymphatic | 571 | 130 | 131 | 0.212 – 1.076 | 3364.1 | **0.0** | 130 | no |

| gate | required | measured | verdict |
|---|---|---|---|
| Murray max relative residual, every bifurcation of every tree | 0 | **0.0** | pass |
| minimum cross-tree surface gap | (1.0 mm, carried over from the liver page) | **-1.667 mm** | **does NOT meet it — reported as found** |
| segment pairs whose surfaces interpenetrate | 0 | **8** | **FAIL — see §10** |

**Priority.** The 1.0 mm figure is the liver page's acceptance threshold, not an anatomical constraint: real companion vessels run with wall-to-wall gaps well under 1 mm. Anatomical fidelity therefore outranks it here. No node was pushed out of the mesentery and nothing was deleted to reach a number; the measured minimum is reported wherever it lands.

| tree pair | min surface gap (mm) | min centreline distance (mm) | interpenetrating pairs |
|---|---|---|---|
| artery vs vein | **-1.667** | 0.096 | 8 |
| artery vs lymph | **1.073** | 1.771 | 0 |
| vein vs lymph | **1.068** | 1.904 | 0 |

Where the global minimum is achieved — the thing to look at before reading the number as a defect:

| field | value |
|---|---|
| tree pair | artery vs vein |
| radius of segment A | 2.478 mm, **133 downstream terminals** (a proximal branch) |
| radius of segment B | 0.62 mm, **1 downstream terminal** (a terminal twig) |
| centreline separation | 1.431 mm |
| midpoint A (scene mm) | [58.6, 372.4, 180.8] |
| midpoint B (scene mm) | [58.9, 377.9, 182.0] |

So the worst case is a thin venous twig crossing a thick arterial branch, not two trunks colliding. no_interpenetration is the constraint that matters geometrically: zero segment pairs whose surfaces intersect. A positive gap below 1.0 mm is an anatomically ordinary companion spacing, not a defect.

### Displacement and deletion accounting

| field | value |
|---|---|
| segments dropped | **{'artery': 0, 'vein': 0, 'lymph': 0}** |
| terminals dropped | **{'artery': 0, 'vein': 0, 'lymph': 0}** |
| nothing deleted | **True** |
| artery segments before → after repair | 1090 → 1090 |
| vein segments before → after repair | 1091 → 1091 |
| lymph segments before → after repair | 571 → 571 |
| free nodes displaced | 89 |
| max displacement from the CCO position | **6.0 mm** (declared cap 6.0 mm) |
| median displacement | 0.868 mm |
| nodes over the declared cap | **0** |
| contour-anchored nodes displaced | **0** |
| pushes vetoed for leaving the mesentery | **148** |
| nodes left exactly where CCO put them by that veto | **18** |
| free nodes outside the mesenteric fan, before → after repair | **809 → 808** of 2713 |
| free nodes pushed outside the fan BY the repair | **-1** |

CCO tests inside() at the bifurcation SITE; the cosmetic subdivision then inserts bowed intermediate nodes that are never tested, which is where almost all of the out-of-fan nodes come from. The before/after pair separates the two causes instead of blaming the repair for geometry it did not create.

CCO growth: provisional vein+lymph -> artery vs provisional -> vein vs final artery -> lymph vs both final — a single ordered pass leaves the first-grown tree blind to the later trees' branches, which is where cross-tree interpenetration came from; two passes let every tree choose its bifurcation sites against the others' real geometry and real final radii

### Arcades and vasa recta — the jejunal/ileal signature

| quantity | jejunum | ileum |
|---|---|---|
| arcade tiers | 2 | 5 |
| vasa recta measured (n) | 90 | 168 |
| vasa recta mean length (mm) | **27.1** | **10.84** |
| vasa recta median (mm) | 27.97 | 11.36 |
| range (mm) | 6.63 – 34.83 | 0.0 – 22.07 |

Jejunal / ileal mean vasa-recta length ratio **2.5** — long vasa recta with few arcade tiers proximally, short vasa recta with many tiers distally, which is the published anatomy. The ratio is a CONSEQUENCE of the tier counts (f_k = (k+1)/(n_tiers+1) of the distance from the bowel wall to the root of the mesentery), not an imposed value. 1011 arcade edges and 258 vasa recta are rendered; both are excluded from the Murray audit and from the flow solve because an arcade is an anastomotic loop, not a bifurcation.

### Cross-organ coupling

The superior mesenteric vein trunk radius is **not** a free parameter. It is derived from the published liver page's MEASURED main portal trunk radius (4.743 mm, diameter 9.49 mm, measured on the same subject's own label) by Murray's law at the portal confluence: r_SMV = r_portal * (Q_SMV/Q_portal)^(1/3) with the textbook SMV share of portal inflow (0.67), giving **4.15 mm** (diameter 8.3 mm) against a published 8-12 mm adult SMV.

## 4. Flow — v = Q/A, per-sample backflow audit

| tree | paths | terminals | trunk v (cm/s) | median v | max v | published reference (cm/s) | trunk in range | samples audited | retrograde samples |
|---|---|---|---|---|---|---|---|---|---|
| mesenteric_artery | 30 | 258 | **28.754** | 17.8759 | 37.861 | 20.0 – 45.0 | **yes** | 1920 | **0** |
| mesenteric_vein | 26 | 260 | **16.138** | 10.5189 | 23.837 | 8.0 – 20.0 | **yes** | 1664 | **0** |
| lymphatic | 18 | 131 | **0.687** | 0.53 | 0.687 | 0.1 – 1.0 | **yes** | 1152 | **0** |

| gate | required | measured | verdict |
|---|---|---|---|
| retrograde samples across every path | 0 | **0** of 4736 | pass |
| flow paths | — | **74** at 64 samples each | — |
| lymph display gain declared | yes | **30.0×** on a true 0.687 cm/s trunk velocity | pass |

Flow magnitudes: SMA **517.7 mL/min** fasting (published Doppler [300, 600]), **1294.3 mL/min** postprandial (published [1000, 1600]), SMV **476.3 mL/min**. Cardiac output 5.1771 L/min and MAP 94.4872 mmHg are this project's own heart-model outputs; the splanchnic and SMA fractions of cardiac output are textbook magnitudes and are flagged as such in `small_intestine_model.json`.

## 5. Microscopic structure — true dimensions, declared magnifications

| quantity | modelled | published reference | declared magnification |
|---|---|---|---|
| circular folds (plicae circulares) | 7 mm pitch, 4.5 mm projection | 3–10 mm apart, projecting 3–10 mm | ×12 (cut view) |
| villus height | 1.05 mm (cut view 0.9 mm) | 0.5–1.6 mm | ×60 |
| crypt depth | 0.25 mm | 0.1–0.4 mm | ×60 |
| enterocyte height | 25.0 µm (radial) | ~25 µm | ×60 |
| microvillus | 1.0 µm long × 0.1 µm | ~1 µm × ~0.1 µm | ×20000 |
| epithelial turnover | 3–5 d, animated at the 4 d midpoint | 3–5 d | — |

### Three-level surface amplification

| level | factor |
|---|---|
| circular_folds_plicae_circulares | ×3.0 |
| villi | ×10.0 |
| microvillus_brush_border | ×20.0 |
| **total** | **×600.0** |

Applied to THIS model's own smooth-cylinder area (**0.2431 m²**, computed as ∫π·d(s)·ds over the measured centreline), ×600.0 gives **145.89 m²**. Helander & Fändriks (2014, doi:10.3109/00365521.2014.898326) revise the human small-intestinal mucosal area down to about **30.0 m²** by direct morphometry. Both are reported side by side; **neither was adjusted to match the other**.

### Four wall layers

| layer | thickness (mm) | published range (mm) | contents |
|---|---|---|---|
| mucosa | 0.5 | 0.4–0.8 | epithelium (enterocytes, goblet, Paneth, enteroendocrine, tuft, stem) + lamina propria + muscularis mucosae |
| submucosa | 0.3 | 0.2–0.4 | dense connective tissue, submucosal (Meissner) plexus, arteriolar/venular plexus, lymphatics |
| muscularis_propria_circular | 0.2 | 0.15–0.3 | inner circular smooth muscle; the layer that generates the peristaltic ring |
| myenteric_plexus | 0.03 | 0.02–0.05 | Auerbach plexus between the two muscle layers, with the interstitial cells of Cajal pacemaker network |
| muscularis_propria_longitudinal | 0.12 | 0.08–0.2 | outer longitudinal smooth muscle; shortens the segment ahead of the wave |
| serosa | 0.06 | 0.04–0.1 | mesothelium + subserosal connective tissue |
| **total** | **1.21** | <=3 mm distended normal small bowel | |

## 6. Peristalsis and segmentation

| quantity | value | source |
|---|---|---|
| peristaltic wave speed | **3.0 cm/s** (published 2.0–4.0) | Peristaltic contractions travel aborally at 2-4 cm/s in the small intestine (Barrett KE et al., Ganong's Review of Medical Physiology 26th ed., ch. 25) |
| slow-wave frequency, duodenum → ileum | **11.5 → 8.5 /min** | Interstitial cells of Cajal set the slow-wave frequency; it falls aborally from ~11-12/min in the duodenum to ~8-9/min in the ileum (Code CF & Marlett JA, J Physiol 1975, doi:10.1113/jphysiol.1975.sp011034; Hansen MB, Pharmacol Toxicol 2003, doi:10.1034/j.1600-0773.2003.920101.x) |
| contraction depth | 0.45 of the local radius | **visualisation value, not measured — declared** |
| wave half-width | 35.0 mm | visualisation value |

Implementation: a travelling radial constriction c(s,t) = A*exp(-((s - v t) / w)^2) is applied to the parametric lumen radius at arc length s, so the ring really moves along the three-dimensional wall; the mesh is NOT scaled as a whole

In-browser probe of the shader's own displacement law, re-evaluated on real ring attributes:

| probe | measured |
|---|---|
| rings sampled | 700 |
| rings undisplaced at t = 0 | **561** (a global scale would be 0) |
| is it a global scale? | **False** |
| contraction factor at the crest | 0.55 (= 1 − 0.45) |
| crest speed measured from crest displacement | **30.2 mm/s** against a declared 30 mm/s |

## 7. Molecular layer — ADME expression

| quantity | measured |
|---|---|
| source | CZ CELLxGENE Census **2025-11-08**, `tissue_general == 'small intestine'`, `is_primary_data == True` |
| cells | **1,434,955** |
| donors | **244** |
| datasets | **22** |
| cell types | **195** |
| genes requested / resolved | 61 / **61** |
| genes absent from the Census var axis | none |
| normalisation | CP10K (raw counts / obs.raw_sum * 1e4) then log1p; mean over all cells of the cell type (zeros included) |
| read strategy | 10 contiguous `soma_joinid` blocks of 150,000, **0 failed**, 12352 s wall |

Measured mean log1p(CP10K) per cell type (zeros included), 12 most abundant cell types:

| cell type | cells | CYP3A4 | ABCB1 | SLC15A1 | SLCO2B1 | ABCG2 | UGT1A1 | CES2 | NR1I2 |
|---|---|---|---|---|---|---|---|---|---|
| stem cell | 141,774 | 0.036 | 0.178 | 0.080 | 0.058 | 0.057 | 0.003 | 0.801 | 0.085 |
| enterocyte | 136,259 | 0.877 | 0.523 | 0.876 | 0.132 | 0.545 | 0.035 | 1.777 | 0.125 |
| CD8-positive, alpha-beta memory T cell | 131,823 | 0.006 | 0.354 | 0.004 | 0.003 | 0.002 | 0.001 | 0.077 | 0.002 |
| plasma cell | 58,207 | 0.006 | 0.012 | 0.004 | 0.003 | 0.007 | 0.000 | 0.029 | 0.001 |
| fibroblast | 55,241 | 0.004 | 0.017 | 0.004 | 0.080 | 0.030 | 0.000 | 0.076 | 0.001 |
| CD4-positive, alpha-beta memory T cell | 45,594 | 0.004 | 0.262 | 0.002 | 0.002 | 0.002 | 0.001 | 0.060 | 0.001 |
| B cell | 42,775 | 0.006 | 0.035 | 0.005 | 0.003 | 0.005 | 0.000 | 0.047 | 0.001 |
| macrophage | 42,325 | 0.005 | 0.013 | 0.005 | 0.592 | 0.044 | 0.000 | 0.052 | 0.001 |
| mesenchymal cell | 37,848 | 0.002 | 0.002 | 0.016 | 0.025 | 0.077 | 0.000 | 0.212 | 0.002 |
| mesodermal cell | 36,231 | 0.000 | 0.001 | 0.002 | 0.003 | 0.033 | 0.000 | 0.111 | 0.001 |
| T cell | 35,114 | 0.007 | 0.265 | 0.005 | 0.003 | 0.003 | 0.000 | 0.064 | 0.001 |
| CD4-positive, alpha-beta T cell | 35,001 | 0.005 | 0.258 | 0.005 | 0.004 | 0.003 | 0.000 | 0.031 | 0.001 |

### What the expression layer is NOT used for

本项目没有任何把表达量换算成 Jmax 或 CL_int 的标定系数，也不发明一个。表达谱只做定性映射（哪些细胞类型、沿轴哪一段表达高），吸收与首过用的是已发表的 Fa/Fg/Fh 与转运时间，二者在界面上分开标注。

### Axial gradients — the half that IS measured

the ONLY axial read-out this project measures. Three of the 195 cell types name their segment in the ontology label, so for those the position is in the data. Every other cell type is position-free and its axial gradient remains literature-only.

| segment | enterocyte cells | above the 2 000-cell reliability threshold |
|---|---|---|
| duodenum | **13,495** | **yes** |
| jejunum | **405** | **NO — provisional** |
| ileum | **6,477** | **yes** |

| gene | duodenum | jejunum (provisional) | ileum | ileum/duodenum | measured direction | published direction | agrees |
|---|---|---|---|---|---|---|---|
| CYP3A4 | 1.5693 | 0.7225 | 1.7999 | **1.147** | proximal low -> distal high | proximal high -> distal low | **NO — opposite** |
| ABCB1 | 0.6255 | 0.3123 | 1.4762 | **2.36** | proximal low -> distal high | proximal low -> distal high | **yes** |
| SLC15A1 | 1.8646 | 1.2824 | 1.7758 | **0.952** | proximal high -> distal low | duodenum/jejunum high, ileum lower | **yes** |

**CYP3A4 measures HIGHER in the ileal enterocyte label than in the duodenal one (ratio 1.147), the opposite of the published proximal-high gradient (Paine MF et al., Drug Metab Dispos 2006, doi:10.1124/dmd.105.008672).** The measured value is reported as computed — no parameter was changed, no cell type swapped and no normalisation re-chosen to bend it back. Unverified candidate explanations: these are per-dataset ontology labels naming the segment a sample came from, not a coordinate along the gut; donor and platform composition differs between segments; and the ileal label carries fewer cells than the duodenal one.

Per-gene maxima across all 195 cell types (so no "highest of any cell type" claim is made without the artefact behind it):

| gene | highest cell type | mean log1p(CP10K) | cells |
|---|---|---|---|
| CYP3A4 | enterocyte of epithelium proper of ileum | **1.7999** | 6,477 |
| ABCB1 | enterocyte of epithelium proper of ileum | **1.4762** | 6,477 |
| ABCG2 | enterocyte of epithelium proper of ileum | **0.8309** | 6,477 |
| SLC15A1 | enterocyte of epithelium proper of duodenum | **1.8646** | 13,495 |
| CES2 | epithelial cell of stomach | **2.5301** | 2 |

The generic `enterocyte` label (136,259 cells) is NOT the maximum for any of these genes — the segment-labelled enterocyte types exceed it. An earlier draft of this record claimed otherwise; the claim was wrong and is corrected here against the saved CSV.

* the jejunal enterocyte label carries only 405 cells, far below the 2000-cell threshold used here, so every jejunal value is PROVISIONAL
* these are ontology labels assigned per dataset, not a measured coordinate along the gut: a label says which segment the sample came from, not where in that segment
* the generic `enterocyte` label (136 259 cells) carries no segment and is excluded from this table
* no calibration from expression level to Jmax or CL_int exists in this project and none is invented here either

### Axial gradients — the half that is NOT measured

| gradient | direction | measured in this project |
|---|---|---|
| CYP3A4 | proximal high -> distal low | **measured, see above** |
| ABCB1_Pgp | proximal low -> distal high | **measured, see above** |
| SLC15A1_PEPT1 | duodenum/jejunum high, ileum lower | **measured, see above** |
| SLC10A2_ASBT | terminal ileum only | **no — carried from the cited literature** |

The 192 cell types whose labels carry no segment are position-free, and their axial direction remains literature-only. the single-cell reference used for the molecular layer is resolved by CELL TYPE, not by position along the gut axis; the axial direction of each gradient is therefore carried from the cited literature and labelled as such, never computed from the expression table

## 8. Segmental absorption and intestinal first pass

7 well-mixed compartments in series, mean small-intestinal transit time 199.0 min, 28.43 min per compartment (Compartmental absorption and transit: the small intestine is represented as 7 well-mixed compartments in series, mean total small-intestinal transit time 199 min (Yu LX & Amidon GL, Int J Pharm 1999, doi:10.1016/S0378-5173(99)00147-7; Yu LX et al., Int J Pharm 1996, doi:10.1016/0378-5173(96)04640-0)).

| probe drug | probe for | Fa | Fg (gut wall) | Fh (liver) | F predicted | F observed | residual | ka derived (/min) | Fa recovered by the solve |
|---|---|---|---|---|---|---|---|---|---|
| 咪达唑仑 (midazolam) | CYP3A4 | 0.99 | 0.51 | 0.56 | **0.2827** | 0.3 | -0.0173 | 0.03274 | 0.99 |
| 非索非那定 (fexofenadine) | P-gp / OATP | 0.35 | 1.0 | 1.0 | **0.35** | 0.33 | +0.02 | 0.00223 | 0.35 |

ka is not fitted: ka = ((1-Fa)^(-1/N) - 1)/tau, closed form from the published Fa for an N-compartment first-order absorption chain; no fitting. Every Fa / Fg / Fh is a published value carried with its citation; the predicted oral bioavailability is then compared against the published observed value rather than fitted to it.

Chain rendered in the scene: lumen (dissolved drug, 7 transit compartments) → apical uptake into the enterocyte (passive + OATP2B1/PEPT1) → enterocyte: CYP3A4 oxidation and P-gp efflux back to the lumen → basolateral exit to the villus capillary network → superior mesenteric vein -> portal vein -> liver (second pass).

## 9. Mobile and rendering performance

Measured with Playwright (`chrome`) in a real Chromium on a **datacentre GPU**, not a laptop and not a real phone.

Renderer string the page actually received: **ANGLE (NVIDIA Corporation, NVIDIA A10/PCIe/SSE2, OpenGL ES 3.2)**. A SwiftShader fallback was explicitly ruled out before measuring — the probe checks the unmasked renderer string and refuses to accept a software rasteriser as a GPU result.

| viewport | view | mean fps | **median fps** | 5th-pct fps | worst frame (ms) | programs | not linked | HUD outside viewport | DOM outside viewport | h-overflow (px) |
|---|---|---|---|---|---|---|---|---|---|---|
| desktop 1440×900 | organ | 60.0 | **60.2** | 57.1 | 17.8 | 9 | **0** | **0** | **0** | **0** |
| desktop 1440×900 | vessels | 60.0 | **60.2** | 57.1 | 17.8 | 12 | **0** | **0** | **0** | **0** |
| desktop 1440×900 | wall | 59.0 | **60.6** | 56.8 | 55.1 | 14 | **0** | **0** | **0** | **0** |
| desktop 1440×900 | villus | 51.0 | **60.2** | 57.1 | 256.8 | 16 | **0** | **0** | **0** | **0** |
| phone390 390×844 | organ | 60.0 | **60.2** | 57.1 | 17.8 | 9 | **0** | **0** | **0** | **0** |
| phone390 390×844 | vessels | 60.0 | **60.2** | 57.1 | 17.7 | 12 | **0** | **0** | **0** | **0** |
| phone390 390×844 | wall | 60.0 | **60.2** | 57.1 | 17.8 | 14 | **0** | **0** | **0** | **0** |
| phone390 390×844 | villus | 60.0 | **60.2** | 57.5 | 17.6 | 16 | **0** | **0** | **0** | **0** |
| phone360 360×800 | organ | 60.0 | **60.2** | 57.1 | 17.7 | 9 | **0** | **0** | **0** | **0** |
| phone360 360×800 | vessels | 60.0 | **60.2** | 57.1 | 17.7 | 12 | **0** | **0** | **0** | **0** |
| phone360 360×800 | wall | 60.0 | **60.2** | 57.5 | 17.6 | 14 | **0** | **0** | **0** | **0** |
| phone360 360×800 | villus | 60.0 | **60.2** | 57.5 | 17.6 | 16 | **0** | **0** | **0** | **0** |

| gate | required | measured | verdict |
|---|---|---|---|
| 60 fps, median frame | ≥ 60 | **60.2 – 60.6** | pass |
| 60 fps, mean frame | ≥ 60 | **51.0 – 60.0** | two views below — see the note |
| shader programs that failed to link | 0 | **0** | pass |
| HUD elements outside the viewport | 0 | **0** | pass |
| DOM elements outside the viewport | 0 | **0** | pass |
| horizontal page overflow | 0 px | **0 px** | pass |

**Two desktop views have a mean below 60 and the reason is a single long frame, not a sustained frame rate:** desktop/wall mean 59.0 fps, median 60.6 fps, worst frame 55.1 ms; desktop/villus mean 51.0 fps, median 60.2 fps, worst frame 256.8 ms. The worst frame anywhere in this run is 256.8 ms — a first-render hitch (shader compile plus the villus instanced-mesh upload) that happens once when the view is first shown. The MEDIAN frame holds at 60.2+ fps in every one of the 12 combinations, and the immediately preceding run on the same 3-D scene (only page prose and one extra HTML table changed since) measured 60.0 fps mean on all 12. Both runs are saved; the mean is reported as measured rather than replaced by the better run.

On a narrow viewport the in-canvas HUD switches to a stacked single column and the DOM control panel moves to the bottom; the two wide data tables become stacked cards rather than scrolling inside a wrapper, which is what removes the last out-of-viewport elements.

## 10. Not achieved / not measured / synthesised / out of range and not tuned

Everything in this section is a deliberate, declared gap. Nothing here is papered over with a proxy value.

1. **The jejuno-ileal luminal course is not contour-anchored.** TotalSegmentator v2.0.1 labels the jejunum and ileum together as one `small_bowel` region with no lumen and no separation between touching loops. The course is a serpentine space-filling path through the measured envelope. `contour_anchored: false` in `si_centreline.json`.
2. **Every luminal calibre is parametric.** Duodenum 27 mm, jejunum 25 mm, ileum 20 mm are published outer diameters, not measurements of this subject. The subject's own envelope inscribed radius along the centreline is smaller (mean 7.04 mm, median 7.5 mm) because the bowel is partly decompressed; both are reported.
3. **The duodenum / jejunum / ileum boundaries are parametric.** The 2:3 jejunum:ileum split is a textbook proportion, not a measurement.
4. **The duodenal centreline measures 130.9 mm against a published 200–300 mm.** A 3-D skeleton of a sharply curved C-loop short-cuts the curve and the label is only 79 skeleton voxels long. **Reported as computed; not stretched to reach the published range.**
5. **3.25 % of the jejuno-ileal centreline samples fall outside the measured envelope** after spline smoothing, which cuts corners the lattice path turned.
6. **Circular folds, villi, crypts and microvilli are entirely parametric.** Villi (0.5–1.6 mm) are at or below the 1.5 mm voxel pitch; no CT segmentation can carry them.
7. **The four wall layers and their thicknesses are published ranges, not measured on this subject.** 1.5 mm voxels cannot resolve a 0.06 mm serosa.
8. **The 45 % contraction depth and the 35 mm wave half-width are visualisation values**, not measured. The wave SPEED (3 cm/s) and the slow-wave FREQUENCY gradient (11.5 → 8.5 /min) are published values with citations.
9. **The migrating motor complex (phase III housekeeping wave) is not modelled.** Only segmentation and propulsive peristalsis are.
10. **No superior mesenteric artery, vein, arcade, vasa-rectum or lymphatic label exists in the source release.** Only the aortic ostium (on this subject's own `aorta` label, with a declared 30.0 mm cranial offset from the measured duodenojejunal flexure) and the two root-of-mesentery endpoints are anchored. Every course is constructed.
11. **The SMA trunk calibre (6.5 mm) and the intestinal lymph trunk calibre (2.2 mm) are published values, not measured.** Only the SMV trunk is derived from a measurement, and that measurement belongs to the liver page.
12. **The root of the mesentery is idealised as a straight line** between two anchored endpoints; the real attachment is an oblique curved line. Its measured straight length is 139.0 mm against a published 150 mm oblique root, duodenojejunal flexure to the right sacroiliac joint (Standring, Gray's Anatomy 42nd ed.).
13. **8 segment pairs still interpenetrate** (artery vs vein), worst surface gap **-1.667 mm** where a 0.62 mm terminal twig crosses a 2.478 mm branch at 1.431 mm centreline separation. The cause is structural in CCO: it clearance-tests only the two NEW segments it creates, never the shortened parent segment, and each later terminal inflates the ancestors' radii. **This was NOT fixed by pushing nodes out of the mesentery or by deleting segments** — anatomical fidelity and the reported segment counts were given priority over the gap number, so the overlap is declared instead.
14. **The 1.0 mm cross-tree gap threshold is not met** (-1.667 mm). That threshold was written after the liver page happened to achieve 6.378 mm; it is not an anatomical constraint — real companion vessels run with sub-millimetre wall-to-wall gaps. The two pairs that do clear it are reported with their own values (artery vs lymph 1.073 mm, vein vs lymph 1.068 mm).
15. **148 clearance pushes were vetoed** because they would have taken a node out of the mesenteric fan, leaving 18 nodes exactly where CCO put them with their violation intact. That is the intended trade: a vessel outside the mesentery is a worse error than a sub-millimetre gap.
16. **89 free nodes were displaced from the position CCO chose** (median 0.868 mm, max 6.0 mm against a declared 6.0 mm cap, 6 at the cap). No contour-anchored node and no trunk node moved.
17. **809 of 2713 free vessel nodes lie outside the mesenteric fan domain**, which is a property of the cosmetic subdivision (CCO tests containment at the bifurcation SITE; the bowed intermediate nodes it inserts are never tested). The repair pushed -1 nodes out — i.e. it did not cause this.
18. **43 arterial, 41 venous and 10 lymphatic terminals could not be attached** by CCO under the angle, containment and clearance constraints, so the trees carry fewer terminals than the 300 / 300 / 140 targeted. The two-pass growth scheme that removed most of the interpenetration is what raised the arterial figure: growing against the other trees' real branches rejects more candidate sites than growing against their trunks.
19. **Arcade anastomoses and vasa recta carry no flow in the solve** and are excluded from the Murray audit. A balanced arcade carries no net flow at rest, and Murray's law does not apply to a loop — but this does mean the rendered arcade network is geometry only, with no haemodynamics behind it.
20. **Intestinal lymph is rendered with a declared 30.0× display gain** on a true 0.687 cm/s trunk velocity. Its published reference range is a RAT mesentery intravital measurement; no human intestinal-lymphatic velocity series exists.
21. **The mass balance does not close from measurement.** Arterial in 517.7 − venous out 476.3 − lymph out 1.5 = 39.9 mL/min, which is a textbook 8 % venous-return reduction standing in for net luminal water absorption, not a measured closure.
22. **The splanchnic and SMA fractions of cardiac output are textbook magnitudes.** Only the cardiac output and MAP they multiply are this project's own computed values.
23. **The classical ×600.0 amplification applied to this geometry gives 145.89 m², against the 30.0 m² that Helander & Fändriks (2014) arrive at by direct morphometry.** Both are reported; neither was adjusted.
24. **No calibration from expression level to Jmax or CL_int exists in this project and none was invented.** The expression layer maps onto the absorption picture qualitatively only.
25. **The measured CYP3A4 axial gradient runs the WRONG WAY relative to the literature**: ileal enterocyte label 1.7999 vs duodenal 1.5693 (ratio 1.147), against a published proximal-high gradient. Reported as computed; nothing was adjusted to recover the published direction.
26. **The jejunal enterocyte label carries only 405 cells**, far below the 2 000-cell threshold used here, so every jejunal value in the axial table is provisional.
27. **The axial direction of every ADME gradient OUTSIDE the three segment-labelled enterocyte types is carried from the literature, not computed.** 192 of the 195 cell types carry no segment in their ontology label. An earlier draft of this record stated this blanket claim for ALL gradients, which was wrong: three cell types do name their segment and those gradients are measured.
28. **The villus counter-current exchanger is drawn from published anatomy only.** No oxygen-tension or solute-gradient solve was run.
29. **The epithelial-turnover animation is a declared 4-day midpoint** of the published 3–5 d range, not a simulated proliferation model.
30. **A declared 100.8° sector of the villus epithelial sheet is omitted** so the core is visible. Cosmetic; no dimension changes.
31. **The respiratory frequency is PhysioNet nsrdb/16265 public data, not any specific individual's.** `is_user_personal_data = false` in every file.
32. **Two of the 12 viewport x view combinations have a MEAN frame rate below 60 (51.0 fps worst) because of a single 257 ms first-render frame** (shader compile plus the villus instanced-mesh upload). The median frame holds at 60.2 fps in all 12. Reported as measured; the better preceding run was not substituted.
33. **The 60 fps figures were measured on a datacentre NVIDIA A10 through ANGLE, not on a laptop and not on a real phone.** The phone viewports are 390×844 and 360×800 rendered on that GPU; they are a layout and draw-cost test, not a measurement of any handset.
34. **The GPU image needed the GLVND dispatch loaders installed at job time** (`libglvnd0 libegl1 libgles2`); without them Chromium silently fell back to SwiftShader and reported 1.3–20 fps. That fallback was detected and rejected rather than published.
