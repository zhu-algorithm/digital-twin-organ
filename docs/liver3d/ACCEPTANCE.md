# Acceptance — 3-D anatomical liver

Every number below was measured, not asserted. Items that fall short are marked
as such rather than quietly relaxed. The source JSON for each figure is named so
it can be checked without re-running anything.

> **v2 (this revision) adds sections 10–15**: the molecular layer rebuilt from the
> full Census, a per-cell layer in which every rendered cell is one real sequenced
> cell, the protein layer and its resolution ladder, the interactive anatomical
> section with a verified solid cut face, the dynamic liver-function channels, and
> the continuous zoom to single cells. Sections 1–9 are **unchanged** and their
> nine listed limitations still stand; v2 adds its own, collected at the end.
> Where a v2 measurement **disagrees** with a v1 number it is shown next to it
> rather than replacing it (see §10.3, CYP2C19).

## 1. Geometry — real anatomy, redistributable licence

`assets/geometry_report.json`

| check | result |
|---|---|
| real anatomical mesh (not hand-modelled, not the old 2-D outline) | **yes** — one openly licensed CT segmentation release, see `GEOMETRY_PROVENANCE.md` |
| source | TotalSegmentator dataset v2.0.1, [10.5281/zenodo.10047292](https://doi.org/10.5281/zenodo.10047292), **CC-BY-4.0**, subject `s0236` |
| redistributable in a public repository | **yes**, attribution only — a CC-BY-SA alternative (BodyParts3D) and an ambiguous multi-terms release (Gibson et al.) were both rejected for exactly this reason |
| parts obtained | liver parenchyma, portal vein, inferior vena cava, gallbladder (aorta used for frame verification only) |
| parts absent from the source | hepatic veins, hepatic artery, bile ducts, Couinaud labels, falciform ligament, separate caudate/quadrate lobe labels |
| capsule volume | **1455.3 mL** (liver label 1434.1 mL + 6 287 voxels of vessel channel bridged by the envelope closing) |
| closed-surface volume of the exported mesh | 1443.4 mL (0.8 % below the voxel count — marching-cubes surface vs voxel volume) |
| face count, `.glb` | **96 180** (parenchyma 78 000 + portal vein 7 180 + IVC 7 000 + gallbladder 4 000) |
| scene total measured in-browser | **137 256** triangles (the four vessel trees add 41 076) |
| `.glb` size / total page weight | 3.85 MB / **6.66 MB** over 21 files |
| right / left hemiliver volume fraction | **0.653 / 0.347** (reference 0.55–0.70 right) ✔ |
| Cantlie plane tilt from true sagittal | **48.8°** |
| retrohepatic IVC centreline curvature | median 2.63 mm, p90 3.83 mm from a straight line (this is the real caval groove's curvature, not an error term) |
| frame verified, not assumed | aorta principal axis cranio-caudal (singular values 69.8 / 6.7 / 5.4 mm), aorta left of and posterior to the IVC, `kidney_right` at larger x ⇒ **RAS+** |

**Cross-source registration residual is not reported** — unlike the heart, this
liver comes from a *single* source, so there is no cross-source registration to
measure. That removes a whole class of error and also means no second source
corroborates the geometry.

## 2. Vascular trees

`assets/vessel_report.json`

| | portal vein | hepatic veins | hepatic artery | bile ducts |
|---|---|---|---|---|
| contour-anchored | **yes** — trunks skeletonised from this subject's own `portal_vein_and_splenic_vein` label (59 anchored nodes) | **no** — ostia on the real IVC, trunk courses constructed from measured landmark planes (53 anchored nodes) | **no** — synthetic portal-triad companion | **no** — synthetic, retrograde |
| segments (after curvature subdivision) | **1 188** | **973** | 372 | 370 |
| bifurcations | **311** | **261** | 48 | 48 |
| terminals | **313** | **263** | 50 | 50 |
| radius range (mm) | **0.699 – 4.743** | 0.906 – 4.750 | 0.611 – 2.250 | 0.543 – 2.000 |
| diameter range (mm) | 1.40 – 9.49 | 1.81 – 9.50 | 1.22 – 4.50 | 1.09 – 4.00 |
| total centreline length (mm) | 5 547 | 4 487 | 1 517 | 1 502 |
| **Murray's law max relative residual** | **0.0** over 310 bifurcations | **0.0** over 261 | **0.0** over 47 | **0.0** over 47 |

- Murray is imposed structurally as `r = r_term · n_downstream^(1/3)`, which
  satisfies `r_parent³ = Σ r_child³` exactly; the residual is a check that the
  construction held, not a fit.
- **Portal trunk calibre is measured, not chosen**: the maximal inscribed sphere
  along the anchored trunk gives radius **4.74 mm** (diameter **9.49 mm**,
  reference 9–13 mm ✔), and `r_term` is back-solved from it, then rescaled after
  growth so the trunk radius matches exactly (4.743 vs 4.74 mm).
- **Hepatic-vein calibre is NOT measured**: no hepatic-vein label exists, so
  `r_term` comes from the published 8–11 mm right-hepatic-vein ostium diameter.
  Measured caval radius at the hepatocaval confluence: 13.5 mm.
- **MHV and LHV are modelled as a common caval trunk** — the majority variant
  (~60–80 % of people). An assumption, not a measurement of this subject.
- Curvature: segments longer than 7 mm are subdivided (portal +514 nodes,
  hepatic +409, arterial +6, biliary +4) and bowed laterally by 12 % of their
  length. **Cosmetic only** — no new terminal, no new bifurcation, no radius
  change; every Murray residual and downstream count is unaffected, and the
  clearance audit below runs on the subdivided geometry.

### Non-intersection — the basic fact of liver vascular anatomy

| check | result |
|---|---|
| minimum **surface** gap, portal vs hepatic venous | **6.378 mm** |
| minimum centreline distance | 8.265 mm |
| required gap (hard constraint during growth) | 1.0 mm |
| segment pairs overlapping | **0** |
| **trees non-intersecting** | **yes** |

## 3. Flow

`assets/flow_report.json`, plus the in-browser probe in `assets/benchmark_report.json`

Particle speed is a real velocity, `v(s) = Q(s)/A(s)`, evaluated per sample
along every path: `Q(s) =` per-terminal flow × downstream terminal count,
`A(s) = πr(s)²`.

| path set | proximal velocity | published reference | within? | retrograde samples |
|---|---|---|---|---|
| portal vein (26 paths) | **22.89 cm/s** | 15–25 cm/s, main portal vein | ✔ | **0** |
| hepatic veins (22) | **16.66 cm/s** | 16–40 cm/s, phasic peak | ✔ | **0** |
| hepatic artery (14) | **33.91 cm/s** | 30–70 cm/s, peak systolic | ✔ | **0** |
| bile ducts (12) | **0.07 cm/s** | 0.03–0.2 cm/s, computed from bile output (no Doppler reference exists) | ✔ | **0** |

| check | result |
|---|---|
| retrograde samples, all 74 paths | **0** |
| in-browser monotonicity probe, all 74 paths over 30 frames | **0 retrograde, 74 advancing** |
| flow direction by tree | portal 26 antegrade, hepatic venous 22 antegrade (parenchyma → caval ostium), arterial 14 antegrade, **bile 12 retrograde** |
| path points inside the liver bounding box | **4 736 / 4 736** |
| mass balance, inflow vs outflow | 1 294.3 vs 1 294.3 mL/min, residual **−0.0** |

**Bile is shown at a declared 120× display gain** (true ≈ 0.07 cm/s, about 300×
slower than portal blood). The gain is printed in the HUD and in the page's
provenance footer; without it the direction would be invisible.

**Provenance carried through unchanged from `liver_model.json`:** cardiac output
5.177 L/min is this project's own cardiac model; the **hepatic fraction of
cardiac output (25 %) and the portal share (75 %) are textbook-magnitude values,
not measurements of this project**, and so are the sinusoidal transit time, the
bile output and the liver mass.

## 4. Organ motion

| check | result |
|---|---|
| frequency | **0.18323 Hz = 10.99 breaths/min**, carried over unchanged from the project's cardiac model (respiratory sinus arrhythmia of PhysioNet nsrdb/16265) |
| `is_user_personal_data` | **false** — public placeholder data, not the user's own recording |
| configured craniocaudal excursion | 2.0 cm = midpoint of the **textbook 1–3 cm** range for quiet breathing; **not measured on this subject** |
| measured excursion of a real mesh vertex over one cycle | **19.13 mm** (below 20 mm because the 2.2° nodding rotation subtracts from the vertical travel of that vertex) |
| **parenchymal deformation** | **max inter-point distance drift 2.8 × 10⁻¹⁴ mm** over the cycle — rigid to floating-point precision. There is deliberately no per-vertex deformation anywhere in this scene |
| waveform | inspiration : expiration = 1 : 1.5, so not a plain sinusoid |
| rotation amplitude (2.2°) | **an assumption**, not from a source in this project |

## 5. Hepatic lobule

| check | result |
|---|---|
| geometry | hexagonal prism, **1.0 mm diameter × 1.6 mm high** (classical human lobule 0.7–2.0 mm), displayed at **60×**, factor shown in the HUD |
| portal tracts | **6**, each carrying the full triad: portal venule (35 µm radius), hepatic arteriole (16 µm), bile ductule (20 µm) |
| central vein | 45 µm radius, on the axis |
| sinusoids rendered / canaliculi rendered | **90 / 30**, over 5 vertical levels |
| blood direction | portal tract → sinusoid → **central vein** (inward) |
| **bile direction** | canaliculus → bile ductule at the portal tract (**outward, counter-current to the blood**) — the structure the 2-D page could only hint at |
| sinusoidal velocity | **0.333 mm/s** = lobule radius / sinusoidal transit time (transit time is textbook-magnitude, flagged in `liver_model.json`) |
| canalicular velocity | **0.033 mm/s**, from bile output per lobule over the canalicular cross-section |
| implied lobule count in the liver | **1 466 329** (textbook order ~10⁶ ✔) — a consistency check on the two inputs, not an independent measurement |
| cutaway | one 60° sector (a classical "portal unit") removed so the central vein, sinusoids and a whole portal tract are visible from outside |
| triangles / particles | 17 072 / 10 800 desktop, 2 160 narrow screen |

### Zonal colouring — and the thing that must not be over-claimed

Driven by this project's own measured single-cell expression
(`liver_model.json .zonal_expression`): three hepatocyte populations annotated by
lobular position, **103 402 cells over 24 donors**.

**The zonal gradient is not the same for every gene.** The colour ramp is rebound
per gene and its *direction* follows the measured ratio:

| gene | central / periportal | ramp |
|---|---|---|
| CYP1A2 | **2.21×** | pericentral-high |
| ABCB11 | 2.16× | pericentral-high |
| SLCO1B1 | 1.46× | pericentral-high |
| CYP3A4 | **1.37×** | pericentral-high |
| CYP2C9 | **1.04×** | essentially **flat** |
| UGT2B7 | 0.98× | flat |
| CYP2D6 | **0.74×** | **reversed** — periportal-high |
| CYP2C19 | 0.67× | reversed |

The ramp is **linear** in the zone coordinate on purpose: shaping it with an
exponent would make a shallow gradient look steep, which is the one claim this
scene must not make. Magnitudes are stated numerically in the HUD and the page
instead. Nothing here asserts that "the zonal gradient was reproduced" as a
single fact.

## 6. Couinaud segments — done, but an approximation

`assets/liver_fields_report.json`, `assets/couinaud.json`

Method: portal-territory Voronoi over the CCO tree's terminals (Selle et al.
2002), named by the standard geometric planes (Soyer 1995 / Fasel 1998).
Agreement between the territory-based and the purely geometric labelling:
**0.809**.

| segment | measured fraction | published CT-volumetry range | inside? |
|---|---|---|---|
| I (caudate) | 4.3 % (63.0 mL) | 1–4 % | **no, 0.3 pp high** |
| II | 7.3 % (106.1 mL) | 4–10 % | ✔ |
| III | 12.4 % (179.9 mL) | 4–10 % | **no, 2.4 pp high** |
| IV | 13.3 % (194.0 mL) | 12–22 % | ✔ |
| V | 15.2 % (221.4 mL) | 8–15 % | **no, 0.2 pp high** |
| VI | 18.8 % (273.9 mL) | 7–15 % | **no, 3.8 pp high** |
| VII | 12.0 % (174.1 mL) | 10–18 % | ✔ |
| VIII | 16.7 % (242.9 mL) | 12–20 % | ✔ |

**4 of 8 segments fall outside the published population ranges.** These values
are reported as they came out and were **not tuned** to land inside the
reference. Two causes, both structural: this subject's own anatomy differs from
the population mean (its right hemiliver is 65.3 %), and the right and left
intersectoral planes rest on *constructed* hepatic-vein trunks because the
source release has no hepatic-vein label.

## 7. Surface

Entirely procedural; no imported textures.

| check | result |
|---|---|
| Glisson's capsule | clearcoat 0.85 serous sheen over the parenchyma, roughness reduced on convex exposed surface |
| parenchymal grain | Worley F2−F1 cell pattern in the tangent plane at ~0.95 mm cell size, i.e. the real lobular length scale, not generic fractal noise |
| **falciform / round ligament** | **measured, not drawn**: the groove is where the umbilical-fissure plane (offset from the Cantlie plane through the umbilical portion of the real left portal branch) meets the capsule **AND** the segmented surface is genuinely concave. **256 of 39 002** parenchymal vertices carry a non-zero mask; measured concavity range −2.07 to +1.71 mm. The hand-drawn hard line of the 2-D page is gone |
| Couinaud and zonal overlays | per-vertex, from the baked segment id and the baked zone field |

## 8. Performance — measured, with the method stated

`assets/benchmark_report.json`. Hardware-accelerated WebGL2, renderer string
reported by the page: `ANGLE (NVIDIA Corporation, NVIDIA A10/PCIe/SSE2,
OpenGL ES 3.2)`.

| viewport | fps mean | fps median | median / p95 frame | frames | scene ready | triangles | particles drawn | HUD boxes outside viewport |
|---|---|---|---|---|---|---|---|---|
| 1440 × 900 | **60.00** | 59.88 | 16.7 / 16.7 ms | 601 | 2 476 ms | 137 256 | 23 680 | **0** |
| **390 × 844** | **60.00** | 59.88 | 16.7 / 16.7 ms | 600 | 839 ms | 137 256 | 4 662 | **0** |
| **360 × 800** | **60.00** | 59.88 | 16.7 / 16.7 ms | 600 | 697 ms | 137 256 | 4 662 | **0** |

A 16.7 ms median *and* p95 is the 60 Hz presentation interval: the scene is
vsync-limited, so **60 fps is a floor on this hardware, not a saturation
measurement**.

**Shader-link assertion.** 17 shader programs compiled, **0 failed to link**, and
the benchmark refuses to report a frame rate until that holds — the heart build
once returned a clean 60 fps while its particle fragment shader had silently
failed to compile, removing its work from the frame budget. The particle draw
count is asserted non-zero for the same reason.

**Honest limitation on where this was measured.** This is **not** an ordinary
laptop. Chrome, Playwright's Chromium and Playwright's WebKit all abort with
SIGABRT inside this agent's sandbox on the user's macOS machine, so no browser
could be run there at all. The measurement was taken in **headless Chromium in a
Linux container on an NVIDIA A10G**, which is a genuine hardware-accelerated
WebGL2 measurement but on datacentre hardware. A phone GPU is a fraction of an
A10G, so these numbers do **not** demonstrate 60 fps on a real phone — what the
phone viewports *do* demonstrate is the layout and the downgrade path.
`tools/benchmark_liver3d.py` runs unmodified on a normal macOS shell:

```bash
python tools/benchmark_liver3d.py --engine chrome --seconds 12
```

## 9. Mobile — screenshot-checked, not bounding-box-checked

The heart page had to hide its HUD on narrow screens: it was laid out for a
desktop canvas, its right-hand half fell outside the phone viewport, and because
it is painted *inside* a canvas, DOM bounding-box checks reported zero overflow.
Only a screenshot caught it.

| check | result |
|---|---|
| viewports actually screenshotted | **1440×900, 390×844, 360×800** — 5 views each (organ, vessels, Couinaud, zones, lobule), 15 images per run |
| HUD elements outside the viewport, after exercising every toggle | **0 of 18 (desktop), 0 of 11 (both phones)** — the HUD reports its own element rectangles in CSS pixels, so this is a real assertion rather than a DOM check that cannot see inside a canvas |
| narrow-screen HUD | switches to a stacked single column, drops the flow table from four rows to two, shrinks the type — it is **not hidden** |
| narrow-screen downgrade | particles 320 → 140 per path (23 680 → 10 360, 4 662 drawn), sprite scale 0.8×, device pixel ratio capped at 1.75 instead of 2 |
| narrow-screen controls | the DOM control panel starts collapsed behind a "panel" button so it cannot cover the organ |

## Not met / limitations, collected

1. **No contour-anchored hepatic veins.** The one source that would have provided
   them (MSD Task08_HepaticVessel) is unreachable from this environment. The
   hepatic venous trunks are constructed from measured landmark planes and are
   marked `contour_anchored: no` everywhere they appear.
2. **Hepatic artery and bile ducts are synthetic.** Anatomically coupled to the
   portal tree as portal-triad companions, calibrated to published trunk
   diameters, but not derived from any image.
3. **4 of 8 Couinaud segment volumes fall outside published population ranges**
   (§6), reported untuned.
4. **Respiratory excursion amplitude and the 2.2° rotation are not measured** on
   this subject; the frequency is public placeholder data.
5. **Bile needs a 120× display gain** to be visible at all.
6. **The zonal overlay on the whole-organ surface is weak.** Zonation is a
   lobule-scale phenomenon; on the organ capsule the zone field is a derived
   quantity and the tint is subtle. The lobule view is where it reads properly.
7. **Flow particles are invisible from outside** while the parenchyma is shown,
   because the parenchyma is opaque — the same limitation the heart has. The
   `parenchyma` toggle exposes them.
8. **No liver clearance model exists**, so the scene shows no drug concentration
   of any kind. Turning expression into clearance needs absolute-abundance
   calibration and IVIVE extrapolation, neither of which this project has built.
9. **Frame rate measured on datacentre hardware**, not a laptop or a phone (§8).

---

# v2 — section plane, zoom to single cells, per-cell gene and protein layers

## 10. Molecular layer — rebuilt from the whole Census

`assets/liver_expression.json`, `assets/liver_expression_normal.json`,
`liver_census_provenance.json`, `liver_expression.csv`

The v1 molecular layer was a bottleneck, not a data limit: Census could not be
read through this sandbox's single local HTTP proxy, so v1 fell back to HTTP
range reads of one source H5AD and reached 11 cell types. Run on a
direct-egress remote container with the proxy unset and TileDB at its own
concurrency, the same contiguous-`soma_joinid`-block route reads the whole liver.

| | v1 (source H5AD route) | **v2 (Census 2025-11-08)** |
|---|---|---|
| cell types | 11 | **171** (all primary) / **155** (normal only) |
| cells | 126 244 | **1 625 338** (all) / **956 026** (normal) |
| donors | 24 | **242** (all) / **137** (normal) |
| datasets | 1 | **33** (all) / **23** (normal) |
| genes | 60 | **61** (the kidney file's list verbatim, `UGT1A3` added) |
| blocks read / failed | — | **11 / 0** |
| wall time | 67 min for one block, locally | **430.3 s for all 11**, remote |

`n_cells_used`, `n_donors` are **per-cell-type arrays**, matching the liver/heart/
kidney file shape rather than the organ-level scalars the oral/oesophagus/stomach/
small-intestine/colon files use. The site renderer branches on the shape.

### 10.1 The Census liver is not all healthy tissue — and the file says so

Pulling `tissue_general == 'liver' and is_primary_data == True` returns tumour
and cholangiopathy datasets as well as normal liver. This is reported, not
hidden, and both aggregates ship:

| disease | cells |
|---|---|
| normal | 956 026 |
| breast cancer (liver metastasis) | 361 215 |
| primary sclerosing cholangitis | 104 667 |
| blastoma | 57 445 |
| primary biliary cholangitis | 39 994 |
| colorectal carcinoma / metastasis | 31 536 |
| neuroblastoma | 24 293 |
| intestinal failure-associated liver disease | 11 969 |
| colorectal cancer | 10 332 |
| small cell lung carcinoma | 9 757 |
| malignant ovarian serous tumour | 9 077 |
| liver cancer | 6 661 |
| lung adenocarcinoma | 2 366 |

`malignant cell` (196 802) and `neoplastic cell` (52 431) are the 1st and 10th
largest cell types in the all-primary aggregate. **A healthy digital twin must
not read those rows as normal liver**, which is why the normal-only aggregate
exists and why every per-cell record carries its own `disease` field.

### 10.2 Per-cell layer — each rendered cell IS one sequenced cell

`assets/liver_cells.json` + 6 lazy shards, `assets/liver_cells_index.json`

| check | result |
|---|---|
| sampling | **per-cell-type reservoir sampling (Algorithm R)** over the full contiguous-block stream of every liver cell in the Census — not a per-block convenience sample |
| quota per cell type | `clip(round(80 + 0.004 × n_cells_of_type), 40, 600)` |
| **random seed** | **20261030**, recorded in every shard |
| cells sampled in total | **16 389** |
| eagerly loaded (`liver_cells.json`) | **1 991** cells, 13 cell types, 1.27 MB — exactly the types the single-cell tier binds, **normal donors only** |
| lazily loadable | **6 shards, 14 398 cells**, all cell types, all disease states |
| fields per cell | `soma_joinid`, `cell_type`, `donor_id`, `dataset_id`, `assay`, `tissue`, `disease`, `sex`, `development_stage`, `raw_sum`, sparse raw counts, sparse CP10K-log1p |
| **browser recomputation check** | the panel recomputes CP10K-log1p in-page from each cell's own counts and raw_sum; **547 values checked against the file, 0 mismatches, max absolute difference 5.0 × 10⁻⁵** (4-decimal rounding in the file) |

**What this buys, stated precisely.** A clicked cell shows `soma_joinid`,
`donor_id`, `dataset_id`, `assay` and 61 raw integer counts that belong to *that
one cell*. It is not a cell-type mean wearing a cell's shape. Cells whose Census
type has no sampled member render **without any gene values at all** rather than
borrowing the type mean.

### 10.3 Zonation — now measured for all 61 genes, and it disagrees with v1 once

`assets/liver_zonation.json`

Three zone-annotated hepatocyte populations in the **normal-only** subset:
periportal **57 756** cells, midzonal **20 168**, centrilobular **29 214**,
**33 donors** each.

| | genes |
|---|---|
| ratio computable | **61 / 61** |
| pericentral-high (> 1.15×) | **45** |
| periportal-high (< 0.87×, i.e. **reversed**) | **6** |
| flat (0.87–1.15×) | **10** |

| gene | v2 central/periportal | v1 (8-gene block) | agreement |
|---|---|---|---|
| CYP1A2 | **3.11×** | 2.21× | same direction, steeper |
| ABCB11 | 2.04× | 2.16× | same |
| CYP3A4 | **1.53×** | 1.37× | same |
| SLCO1B1 | 1.49× | 1.46× | same |
| UGT2B7 | 1.20× | 0.98× | v1 flat → v2 pericentral-high |
| CYP2E1 | 1.20× | — | new |
| CYP2C9 | **1.09×** | 1.04× | both **flat** |
| CYP2D6 | **0.78×** | 0.74× | both **reversed** |
| CYP2C19 | **1.38×** | **0.67×** | **DISAGREE — direction flips** |

CYP2C19 is the one real conflict: one dataset / 24 donors said periportal-high,
the full Census normal subset / 33 donors says pericentral-high. **Neither number
is tuned and neither is withdrawn.** `liver_model.json`'s v1 block ships
unchanged alongside `liver_zonation.json`, and `api.diag.zonation.legacy_eight_gene_block`
exposes it so the two can be compared in-page.

**There is still no global "CYPs are pericentral" claim anywhere.** Each gene's
ramp direction is the sign of its own measured log-ratio and its amplitude scales
with the magnitude, so a flat gene renders flat and a reversed gene ramps the
other way.

## 11. Protein layer — four tiers, and a corrected premise

`liver_protein_evidence.csv` (269 rows), `assets/liver_protein_evidence.json`

**Two corrections to the brief this layer was built under.**

1. The brief said `protein_abundance_anchors.csv` is "named protein, contents are
   actually HPA mRNA" and that the project therefore holds no liver protein data.
   Reading the file: of its **133 liver rows**, **56** are
   `relative_mRNA_abundance` (nTPM, HPA) — so the complaint is right about those —
   but **34** are `relative_protein_intensity`, **7** are dataset references, and
   **36 are `absolute_protein_abundance` with `has_absolute_quantification = True`**,
   covering **23 of the 61 genes** by validated LC-MS/MS with signature peptides.
   Those 36 rows are the only absolute protein numbers in the project and are
   carried through as their own tier. The 56 mRNA rows are **not** used as protein.
2. The brief said liver single-cell proteomics "basically does not exist". A PRIDE
   search contradicts this. It exists and is published; see tier 4.

| tier | what it is | resolution | absolute? | coverage |
|---|---|---|---|---|
| **本细胞实测** gene counts | Census raw counts + CP10K-log1p | **this one cell** | n/a | 61/61 genes, 1 991 cells loaded |
| **本细胞类型** HPA IHC | antibody immunohistochemistry staining level, intensity, % of cells stained, subcellular location, per annotated liver cell type (Hepatocytes, Cholangiocytes) | cell **type** | **no** — ordinal categories (high/medium/low/not detected), no copies per cell | **50 / 61 genes** |
| **整块组织** shotgun iBAQ | MaxQuant iBAQ, label-free, 29 healthy tissues | **none** — liver homogenate, every cell type mixed | no — within-sample molar proxy | **52 / 61 genes** detected |
| **整块组织（绝对）** targeted LC-MS/MS | pmol/mg protein in human hepatocytes, liver microsomes, liver total membrane | none — pooled preparations | **YES** | **23 / 61 genes**, 36 rows |
| **他人单细胞** scDVP coverage | protein quantified per individual laser-microdissected hepatocyte | **single cell — but a different cohort** | no | **44 / 61 genes** |

### 11.1 Identifying the liver run in PXD010154 — an inference, reported as one

The deposit ships 30 anonymous experiment ids and **no run-to-tissue key**.

- A plasma-secretome panel (ALB, SERPINA1, HP, TF, fibrinogen…) was tried and
  **rejected**: those proteins sit in MaxQuant's contaminant list and are abundant
  in any blood-rich tissue. It ranked `P015159` first, which is not the liver.
- The call is made on a **hepatocyte-restricted, non-secreted** panel —
  SLC10A1, SLCO1B1, SLCO1B3, ABCB11, CYP2E1, ADH1B, CYP2C9, CYP3A4, FMO3, UGT2B7,
  CYP1A2, CES1 plus ASGR1 and CPS1.

| | result |
|---|---|
| experiment assigned | **P012502** |
| panel genes peaking in it | **12 / 14** |
| separation from the runner-up | **3.25×** on the aggregate panel fraction |
| ASGR1 | 375.3 ppm of total iBAQ here vs 0.6 ppm in the next tissue |
| CYP2E1 | 1 562.9 ppm vs 8.4 ppm |
| CYP1A2 | 525.7 ppm vs 14.4 ppm |

This is an **inference from marker dominance, not a label from the depositors**,
and the full per-gene evidence is in
`assets/liver_protein_evidence.json .liver_experiment_identification` so it can be
checked or overridden.

### 11.2 Single-cell liver proteomics — it exists, and here is what it reaches

**PXD062231**, "Single-cell spatial proteomics maps human liver zonation patterns
and their vulnerability to disruption in tissue architecture" — Weiss CAM et al.
**2026 Nat Metab 8(3):741-756**, **PMID 41721135**, DOI 10.1038/s42255-026-01459-2,
Orbitrap Astral, label-free DIA. Also **PXD070185** (Zenodo DOI
10.5281/zenodo.18824148, no PubMed id).

Read out of PXD062231's DIA-NN gene-group matrix (`report.gg_matrix.tsv`,
range-extracted from the 1.97 GB archive):

| check | result |
|---|---|
| single hepatocytes | **792** |
| gene groups quantified | 4 519 |
| proteins per single cell | median **2 230**, range 0–3 691 |
| **61-gene ADME panel quantified in ≥1 single cell** | **44 / 61** |
| best-covered | GAPDH 96% of cells, ADH1B 94%, ALDH2 88%, CES1 87%, CYP2C9 86%, CYP3A4 75%, CYP1A2 71% |

**What this tier does NOT say.** It does not give a protein value for the cell you
clicked. Those are different cells from a different cohort. **No cell anywhere in
this project has both its transcriptome and its proteome measured.**

**Why no per-zone protein value was derived.** The deposit ships 18 per-individual
zonation score files (`Individual*_scores.txt`, **792 score rows in total, matching
the 792 sample columns**) but **does not document the key joining a score row to an
mzML sample column**. Joining them by row order would be an invented mapping, so
no periportal/pericentral protein number was taken from this dataset.

### 11.3 The resolution gap, on every row

The per-cell panel tags each column with its tier (本细胞实测 / 本细胞类型 /
整块组织 / 他人单细胞) and states in plain language what each is and is not.
**There is no expression→abundance and no abundance→activity coefficient anywhere
in this project**, so the panel never converts a gene reading into a protein
amount or the reverse, and where a tier has no value for a gene it prints 无数据
rather than borrowing from a neighbouring tier. The absolute tier's denominator is
*a milligram of fraction protein*, not *a cell*: turning pmol/mg into molecules per
hepatocyte would need a protein-per-cell factor this project does not have.

## 12. The anatomical section — and how the cut face was proved solid

`section_report.json`, `assets/benchmark_report.json`

| degree of freedom | provided |
|---|---|
| presets | sagittal, coronal, axial |
| free plane | two angle sliders (azimuth 0–2π, elevation ±1.5 rad) |
| position along the normal | slider, 0.03–0.97 of the span |
| side kept | flip button |
| travels with respiration | **yes** — the plane is defined in the breathing body's frame and the cap geometry is parented to it |

### 12.1 Capping — three independent area computations

Two-pass stencil: the clipped solid's **back faces increment** the stencil and its
**front faces decrement** it; a quad on the plane is painted **only where the
stencil is non-zero**, i.e. exactly where the plane is inside the solid. Nothing
is painted outside the cross-section, so the cap can never be a disc that
overshoots the organ.

| plane | CPU mesh-slice (Python) | CPU mesh-slice (browser) | CPU grid raster (browser) | **GPU stencil cap (pixels)** | GPU / Python |
|---|---|---|---|---|---|
| sagittal | 7 774.98 mm² | 7 774.98 mm² | 7 772.97 mm² | **7 797.42 mm²** (15 942 px) | **1.0029** |
| coronal | 11 582.41 mm² | 11 582.41 mm² | 11 578.26 mm² | **11 579.72 mm²** (23 675 px) | **0.9998** |
| axial | 14 496.68 mm² | 14 496.68 mm² | 14 491.41 mm² | **14 402.39 mm²** (29 446 px) | **0.9935** |
| free oblique (az 0.9, el 0.35) | — | 7 612.36 mm² | 7 618.41 mm² | **7 579.28 mm²** (15 496 px) | 0.9957 |

The Python column intersects the 78 000-face parenchyma with the plane and
applies an **oriented** shoelace sum; it never touches the GPU. It is the same
algorithm as the browser mesh-slice column, independently reimplemented, and the
two agree to every digit printed here (ratio 1.0000 on all three presets) — that
column pair is a port check, not an independent measurement. The independent
measurements are the grid raster and the GPU stencil cap. The GPU column
counts cap pixels under an orthographic camera of known scale with every other
renderable hidden; it never touches the CPU geometry. **The GPU cap and the Python slice agree to within 0.02–0.65 %** across the three presets and the oblique plane.
A missing cap would give 0 px; an oversized cap would give an area far above the
slice area.

> Two bugs were found by this audit and are recorded because they are the kind of
> thing that silently produces a convincing-looking wrong picture:
> an **unoriented** shoelace sum partially cancels and under-reported the sagittal
> cut as 1 676 mm² instead of 7 775 mm²; and taking the face-frame origin as the
> plane's closest point to the *world* origin instead of the *organ* centre put
> the cut polygon several hundred millimetres outside the classification grid, so
> the grid area read 0 while the mesh slice read a plausible number.

### 12.2 Is anything actually being cut?

Two readbacks of the same framebuffer, section on then off, with nothing else
touched:

| check | result |
|---|---|
| framebuffer pixels | 907 200 |
| organ pixels with no section | 157 091 |
| **pixels changed by enabling the section** | **110 793** |
| fraction of the organ that changes | **70.5 %** |

### 12.3 Vessel lumens on the cut face — calibres straight out of the file

Every tree segment straddling the plane contributes one opening, placed at the
interpolated crossing with the radius interpolated between the two node radii of
`liver_vessels.json`. **`radius_scale_applied: 1.0` everywhere.** The browser's
counts and diameters were recomputed independently in Python:

| plane | portal vein | hepatic veins | hepatic artery | bile ducts |
|---|---|---|---|---|
| sagittal | **23**, 1.40–6.67 mm | **4**, 1.81–9.50 mm | **6**, 1.22–3.32 mm | **5**, 1.09–2.95 mm |
| coronal | **18**, 1.40–6.67 mm | **19**, 1.81–8.91 mm | **6**, 1.22–3.32 mm | **6**, 1.09–2.95 mm |
| axial | **30**, 1.40–5.48 mm | **19**, 1.81–2.77 mm | **11**, 1.22–2.48 mm | **9**, 1.09–2.17 mm |

**12 of 12 tree×plane combinations match exactly** between the
browser and the independent Python recomputation, and every maximum diameter agrees
to **0.0 mm**. Each opening is drawn as an opaque wall annulus plus a
translucent lumen disc so the flow particles inside the kept half of the vessel
remain visible through it, and carries a direction mark read from the sign of
(segment direction · plane normal) — into-page or out-of-page comes from the
tree's own geometry.

### 12.4 Couinaud on the cut face

Classified by the **same rule** that baked `_field0.y` onto the capsule: the
geometric Couinaud label of the nearest portal terminal (portal-territory
Voronoi, Selle et al. 2002), evaluated on a 256 × 256 grid over the face.
`tools/build_section_assets.py` replays `geometric_segment()` from
`build_liver_fields.py` verbatim over the **313** portal terminals and ships them
as `assets/couinaud_terminals.json`, so the section and the capsule cannot
disagree. Mid-sagittal cut: **S1 1 391, S3 2 268, S4 2 890, S5 151, S6 65,
S7 513, S8 497 mm²**.

### 12.5 The cut face travels with the organ — and it is still rigid

Two points on the cut face, **97.7 mm apart**, tracked over one breath:

| check | result |
|---|---|
| cut-face centroid peak-to-peak travel | **19.58 mm** |
| configured craniocaudal excursion | 20 mm (textbook midpoint, **not measured on this subject**) |
| **inter-point distance drift** | **4.26 × 10⁻¹⁴ mm** — rigid to floating-point precision |
| two-point displacement-vector difference | **3.75 mm** |

The 3.75 mm difference is the **signature of the real 2.2° rotation** about the
hepatocaval confluence: a pure translation would make it exactly zero, and any
resizing would show up in the distance drift instead. The 19.58 mm centroid travel
falls 0.42 mm short of the 20 mm configured craniocaudal excursion for the same
reason — on this side of the hepatocaval pivot the 2.2° nodding rotation subtracts
from the vertical travel rather than adding to it, which is the same effect that
put a single mesh vertex at 19.13 mm in §4. Neither figure is evidence about a real
liver's excursion: 20 mm is the configured input, and these are checks that the cut
face obeys it.

## 13. Dynamic liver function — five channels, each graded

`liver_function_report.json`

All five are plotted against the **same** zone coordinate,
`zone01 = d(nearest portal terminal) / (d_portal + d_hepatic-vein terminal)`,
computed from this subject's own two vascular trees (313 + 263 terminals). What
differs is where the numbers on that axis come from.

| channel | grade | numbers | citation / value grade |
|---|---|---|---|
| **氧梯度** pO2 | **已发表** | 60–65 → 30–35 mmHg, **linear** in the zone coordinate (profile shape not measured) | Jungermann & Kietzmann, Hepatology 2000;31:255-260, DOI 10.1002/hep.510310201 · value **B** (endpoints from the field's standard values, **not re-read from the source in this session**) · paper **A** |
| **代谢分区** zone 1 vs zone 3 | **已发表** | categorical pathway assignment, no value plotted | Jungermann & Kietzmann, Annu Rev Nutr 1996;16:179-203 · Kietzmann, Redox Biol 2017;11:622-630 · paper **A** |
| **CYP 逐基因** | **本项目实测** | 61 per-gene ratios, 45 pericentral-high / 6 reversed / 10 flat | this session's Census pull; `liver_zonation.json` |
| **底物衰减** | **示意** | `C(zone) = exp(-zone · (-ln(1-E)))`, **E is a slider the reader sets** | **none** |
| **胆汁逆流** | **本项目实测** | pericentral → periportal, ~0.07 cm/s, **declared 120× display gain** | derived from bile output per lobule over the canalicular cross-section |

**The hard statement, repeated in the UI next to the channel and in the report:**
this project has **no expression→CLint calibration** and has **not** done IVIVE.
The substrate-decay curve is **not** computed from any expression value on this
page. Reading it as a clearance predicted from the CYP expression shown elsewhere
in the scene is wrong, and the extraction ratio *E* is an input the reader chooses.

**A declared gap in the zonation channel**: the 61-gene panel is ADME only. It
contains no CPS1/OTC/ARG1 (urea cycle), no GLUL (glutamine synthesis) and no
PCK1/G6PC (gluconeogenesis), so the **glucose and urea arms of the metabolic-zonation
channel have no supporting data of this project's own** and are shown as published
qualitative assignments.

## 14. Continuous zoom — organ to one cell

Six rungs on one continuous slider; the scale bar is computed **from the camera
every frame** (millimetres per CSS pixel at the orbit target) while the
magnification printed beside it is the rung's **declared** scene-units-per-mm, so
the two are independent readouts and a wrong magnification would show up as a
disagreement.

| rung | field of view | declared magnification | geometry |
|---|---|---|---|
| 器官 | 220 mm | 1× | real CT segmentation, 137 256 triangles |
| 肝段 | 80 mm | 1× | same |
| **解剖切面** | 200 mm | 1× | same + stencil cap + lumen openings |
| 肝小叶 | 1.9 mm | **60×** | 1.0 mm × 1.6 mm hexagonal prism (v1) |
| 门管区／血窦 | 0.40 mm | **3 000×** | hepatic-plate slab |
| **单个细胞** | 0.075 mm | **3 000×** | same slab, camera at one hepatocyte |

### 14.1 Single-cell tier — real sizes, declared magnification

| structure | rendered size | count |
|---|---|---|
| hepatocyte | **22 µm**, polyhedral, plates one cell thick | **352** (308 shown, front plate cut away) |
| sinusoid lumen | **10 µm** | 14 |
| space of Disse | **0.5 µm** | 14 |
| bile canaliculus | **0.75 µm**, between adjacent hepatocytes | 16 |
| microvilli on the sinusoidal face | **0.5 µm** tall | 1 600 |
| Kupffer cell | **13 µm**, in the sinusoidal lumen | 88 |
| hepatic stellate cell | **10 µm** body in the space of Disse, with processes | 29 |
| sinusoidal endothelium | 9 × 5 µm flattened, 1.2 µm thick | 117 |
| cholangiocyte | **9 µm** cuboidal, around a 7 µm-radius ductule | 24 |
| field of view | **494 × 46 × 264 µm** (portal→central × vertical × across plates) | |
| **clickable cells bound to a real Census cell** | **610 / 610** | |

| honesty item | statement |
|---|---|
| cell-number proportions | published human-liver shares **by cell number** (hepatocyte 60 %, sinusoidal endothelium 20 %, Kupffer 15 %, stellate 5 %) |
| **microvillus density** | **SCHEMATIC** — the 0.5 µm height is real, the surface density is not taken from any measurement |
| **endothelial fenestrae** | **NOT RENDERED** — 50–150 nm is below this tier's resolution; the lining is drawn as a continuous sheet |
| plate cutaway | the front hepatic plate (**44 of 352** hepatocytes) is hidden so the sinusoid, its lining, the Kupffer and stellate cells and the canaliculus are visible from outside — the same device as the lobule's 60° sector |
| **what a rendered cell represents** | **610 cells are rendered; an adult human liver holds of the order of 2 × 10¹¹ hepatocytes. A rendered cell does NOT stand for a block of real cells — it IS the one Census cell that was sampled, and its gene values are that cell's own measurements.** |
| zone mapping | a hepatocyte in the portal third of the field is bound to a `periportal region hepatocyte`, the middle third to `midzonal`, the central third to `centrilobular region hepatocyte` — the render position matches the Census annotation of lobular position; `api.diag.cellTier.census_type_mapping` records every fallback taken |

### 14.2 Per-cell gene tint

Each hepatocyte can be coloured by **its own** measured CP10K-log1p for the bound
gene, normalised by the **maximum over the rendered hepatocytes** — the maximum is
printed, and a cool-coloured cell means that cell's own count for that gene was
zero. Measured for CYP1A2 in the audit run: maximum CP10K-log1p over the rendered
hepatocytes **1.9677**, and **46 of 352** rendered hepatocytes have a non-zero
CYP1A2 count of their own — which is itself the per-cell sparsity of 10x data made
visible, and exactly the thing a cell-type mean hides. This is the per-cell measurement in the geometry, not a cell-type mean
painted onto cells.

## 15. Performance — v2, measured, with the method stated

`assets/benchmark_report.json`. Unmasked renderer string the page actually
received: **`ANGLE (NVIDIA, Vulkan 1.4.341 (NVIDIA NVIDIA A10G (0x00002237)), NVIDIA)`**.

On this GPU **only `--use-angle=vulkan` reaches the hardware**; the GL and EGL
ANGLE backends fall back to SwiftShader and report software-rasteriser numbers as
if they were GPU results. The benchmark asserts the renderer string is not a
software rasteriser before reporting anything.

| viewport | per-frame render cost | frame-time p95 | 60-fps budget | headroom | scene ready | triangles (organ / cell tier) | particles drawn | **HUD boxes outside viewport** |
|---|---|---|---|---|---|---|---|---|
| 1440 × 900 | **16.666 ms** | 16.70 ms | 16.67 ms | 1.00× | 1 846 ms | 137 256 / 58 628 | 23 680 | **0 of 21** |
| **390 × 844** | **16.666 ms** | 16.70 ms | 16.67 ms | 1.00× | 1 626 ms | 137 256 / 58 628 | 4 662 | **0 of 14** |
| **360 × 800** | **16.666 ms** | 16.80 ms | 16.67 ms | 1.00× | 1 622 ms | 137 256 / 58 628 | 4 662 | **0 of 14** |

A per-frame cost pinned at 16.666 ms — with a p95 of 16.70 ms on 1440 × 900 and
390 × 844 and **16.80 ms** on 360 × 800, one 0.1 ms histogram bin higher rather
than a slower frame — is the **60 Hz presentation interval**. This headless
Chromium *is* vsynced, so **60 fps is a floor on this hardware, not a saturation
measurement**, and the figures do not say how much headroom remains. (The general
warning that headless Chromium can run rAF uncapped does not apply to this build —
it was checked, not assumed; each viewport carries its own `vsync_detected` flag
in `benchmark_report.json`, set from that viewport's own per-frame cost.)

**17 shader programs compiled, 0 failed to link**, and the benchmark refuses to
report a frame rate until that holds.

**Where this was measured.** Headless Chromium in a Linux container on an **NVIDIA
A10G** — a datacentre GPU, **not a laptop and not a real phone**. The phone
viewports demonstrate **layout and the downgrade path**, not phone GPU
performance. Playwright's browsers cannot be installed in this agent's sandbox on
the user's macOS machine (`EPERM` on the browser cache directory), so no local
browser run was possible at all.

**Page weight**, counted from the requests the page actually issued during the
benchmark run (`assets_served`): **29 files, 8.28 MB** eagerly fetched, of which
`liver.glb` is **3.85 MB** (≤ 7 MB). The **6 per-cell shards (6.93 MB)** are lazy and
appear in **none** of the three runs' request lists.

## v2 — not met / unmeasured / synthetic / out-of-range, collected

The nine v1 limitations above all still stand. In addition:

1. **No cell has both a transcriptome and a proteome in this project.** The
   per-cell gene values are real measurements of real cells; the protein values on
   the same row belong to the cell *type*, to a tissue homogenate, to a pooled
   subcellular fraction, or to 792 cells of a different cohort. This gap is
   labelled on every row and cannot be closed with the data to hand.
2. **No per-zone protein values were derived** from the published single-cell
   liver proteome (PXD062231) because the deposit does not document the key
   joining its 792 zonation scores to its 792 sample columns (§11.2).
3. **The liver run in PXD010154 is an inference** (`P012502`, 12/14 panel genes,
   3.25× separation), not a label supplied by the depositors (§11.1).
4. **The oxygen channel's endpoints are value-grade B** — 60–65 / 30–35 mmHg are
   reported from the field's standard values and were **not re-read from the
   source text in this session**. The profile *shape* between them is not measured
   at all.
5. **The substrate-decay channel is schematic.** No parameter in it is measured or
   published; the extraction ratio is a slider.
6. **The glucose and urea arms of the metabolic-zonation channel have no data of
   this project's own** — the 61-gene panel is ADME only (§13).
7. **Microvillus surface density is schematic** and **endothelial fenestrae are not
   rendered at all** (§14.1).
8. **CYP2C19's zonal direction disagrees between the v1 and v2 measurements**
   (0.67× vs 1.38×). Both are reported untuned; neither is withdrawn (§10.3).
9. **The all-primary molecular aggregate mixes disease states**, including 196 802
   `malignant cell` and 52 431 `neoplastic cell` rows. The normal-only aggregate is
   shipped alongside it and the per-cell layer is normal-only, but the site's
   default aggregate is still the all-primary one for cross-organ comparability
   (§10.1).
10. **The per-cell sample is 16 389 of 1 625 338 cells (1.0 %)**, 1 991 of them
    loaded eagerly. It is a seeded reservoir sample, not the population.
11. **The free oblique plane is not recomputed in Python** — only the three presets
    are independently checked (§12.1); the oblique plane's three browser-side
    numbers agree with each other to 0.51 %.
12. **Respiratory excursion amplitude and the 2.2° rotation remain unmeasured on
    this subject**, so the 19.58 mm cut-face travel inherits that assumption.
13. **Frame rate still measured on datacentre hardware**, and the measurement is
    vsync-limited, so it reports a floor rather than the available headroom (§15).
