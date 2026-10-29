# Geometry provenance — 3-D anatomical small intestine

Every surface in `assets/small_intestine.glb` derives from one openly licensed CT
segmentation release. Nothing is hand-modelled. Where a structure could **not**
come from that release — and for the small intestine that is most of what you
see — this file says so explicitly and names what replaced it.

## Source

| field | value |
|---|---|
| dataset | **Dataset with segmentations of 117 important anatomical structures in 1228 CT images** (TotalSegmentator dataset), Wasserthal J., Breit H.-C., Meyer M. T., Pradella M., Hinck D., Sauter A. W., Heye T., Boll D., Cyriac J., Yang S., Bach M., Segeroth M. |
| version | v2.0.1, published 2023-10-27 |
| DOI | [10.5281/zenodo.10047292](https://doi.org/10.5281/zenodo.10047292) |
| licence | **CC-BY-4.0** (`license.id: "cc-by-4.0"` in the Zenodo record metadata) |
| redistributable in a public repository | **yes** — attribution only, no non-commercial and no share-alike clause |
| subject used | **s0236**, `study_type = ct abdomen`, 1.5 mm isotropic (the release is resampled); volume [252, 198, 308] |
| structures used | `small_bowel`, `duodenum`, `colon`, `stomach`, `liver`, `pancreas`, `aorta`, `portal_vein_and_splenic_vein`, `inferior_vena_cava` |
| access method | HTTP range requests against the 23.5 GB Zenodo zip: the ZIP64 central directory is read from the tail, then only the nine mask members are fetched. See `tools/remote_zip_reader.py`. |
| CT image used | **no** — `ct.nii.gz` was never downloaded. This project contains no image-intensity information of any kind, for any organ. |

Attribution required by CC-BY-4.0 is carried in this file, in
`assets/geometry_report.json`, in the scene's own provenance footer and in the
page the scene is embedded in.

### Subject selection — re-screened for this organ

s0236 is the subject the published liver page uses, carried over so the five
organs share one anatomy. Reusing it is not by itself a justification, so it was
re-screened against small-bowel-specific criteria. All five are measured:

| criterion | measured on s0236 | verdict |
|---|---|---|
| `small_bowel` entirely inside the field of view | **0** voxels on any of the six volume faces | pass |
| `duodenum` entirely inside the field of view | **0** voxels on any face | pass |
| `small_bowel` effectively one connected component | 5 components, largest = **0.99997** of voxels (the other four are 1–6 stray voxels) | pass |
| `duodenum` present and connected | **1** component | pass |
| both bowel junctions present | touches `duodenum` in **302** voxels (duodenojejunal flexure), touches `colon` in **8157** voxels (ileocaecal junction) | pass |
| envelope volume consistent with a published in-vivo length | **1431.1 mL**; divided by the cross-section of a 20 mm-calibre tube this is **4.56 m**, inside the published in-vivo 3–5 m | pass |

## Frame

Identical to the liver page, so the two scenes are in register. The release's
NIfTI sform is diagonal and **RAS+**; scene space is
`(patient-left, superior, anterior)`:

```
x_scene = -x_RAS   y_scene = z_RAS   z_scene = y_RAS      det = +1
```

## What the release can and cannot give for this organ

This is the central honesty point of the page, so it is stated before any number.

**`small_bowel` is a single region.** TotalSegmentator v2.0.1 labels the whole
jejunum and ileum together as one region with **no lumen** and **no separation
between touching loops**. `duodenum` *is* a separate label — verified against the
release's own member list for this subject rather than assumed — and the villi
(0.5–1.6 mm) and circular folds sit at or below the 1.5 mm voxel pitch.

Consequences, carried into every output file as `contour_anchored`:

| structure | anchored? | what it actually is |
|---|---|---|
| small-bowel occupancy envelope | **yes** | the real loop territory of a real abdomen; used as a hard containment constraint |
| duodenal centreline | **yes** | skeletonised from the real `duodenum` label |
| duodenojejunal flexure, ileocaecal junction | **yes** | taken from where the labels actually touch |
| local lumen inscribed radius profile | **yes** | distance transform of the envelope |
| jejuno-ileal luminal course | **no** | serpentine space-filling path through the measured envelope (below) |
| every luminal calibre | **no** | published segmental outer diameters |
| duodenum / jejunum / ileum boundaries | **no** | textbook 2:3 jejunum:ileum proportion |
| circular folds, villi, crypts, microvilli | **no** | parametric, below CT resolution |
| four wall layers and their thicknesses | **no** | published ranges; 1.5 mm voxels cannot resolve them |
| SMA / SMV / lymphatic courses | **no** | no label exists; only the aortic ostium is anchored |

## The jejuno-ileal course

### Why the envelope skeleton is not the answer

Skeletonising the closed, hole-filled envelope gives 1801 skeleton voxels in
4 components (largest 1797). Its **longest geodesic path is
472.8 mm** — 0.47 m, an order of magnitude short of the published
in-vivo small-intestinal length. The reason is anatomical, not algorithmic:
adjacent loops are in contact on CT, so the skeleton merges them into a single
branch. Total skeleton branch length is 3743.6 mm, but that is a branched
network, not a course. **The envelope skeleton measures loop territory, not
luminal course**, and is reported here as a measurement rather than used as one.

### The construction rule, with no fitted length

| step | rule | measured value |
|---|---|---|
| lattice pitch | 2 × p75 of the envelope's inscribed radius | **14.7 mm** |
| node admission | local inscribed radius ≥ median inscribed radius | **4.5 mm** |
| nodes admitted | lattice ∩ envelope ∩ threshold | **203** |
| ordering | open TSP, first node fixed at the duodenojejunal flexure, last at the ileocaecal junction; nearest-neighbour seed then 2-opt + Or-opt | 4300.3 mm → **3753.2 mm** |
| containment repair | any straight hop leaving the envelope re-routed on the in-envelope voxel geodesic | **8** of 8 offending hops |
| smoothing | cubic spline, 2 mm resample | **3262.3 mm** |

Both admission numbers are measured on this subject. Neither is tuned, and no
step targets a length.

### Measured lengths against the literature

| quantity | measured | published reference | in range? |
|---|---|---|---|
| duodenal centreline | **130.9 mm** (raw skeleton geodesic 151.7 mm) | 200–300 mm (Gray's Anatomy 42nd ed.) | **no — reported below range as found** |
| jejuno-ileal centreline | **3262.3 mm** | — | — |
| duodenojejunal bridge | **18.5 mm** | — | — |
| total small intestine | **3.412 m** = 130.9 + 18.5 + 3262.3 mm | 3–5 m in vivo | **yes** |
| envelope volume ÷ 20 mm cross-section | **4.56 m** | 3–5 m in vivo | **yes** |
| samples inside the measured envelope | **96.75 %** of 1631 | — | — |

The 18.5 mm bridge is the straight segment between the duodenal skeleton tip and
the first jejunal lattice node. They are not the same point: the duodenojejunal
flexure anchor is the centroid of the `small_bowel`/`duodenum` contact patch,
while the duodenal centreline ends at the skeleton tip inside the duodenum
label. It is part of the assembled arc length and is reported rather than
absorbed, so the three addends sum to the stated total
(`geometry_report.json` → `centreline.total_length_check_mm` = 3411.7 mm).

The duodenal shortfall is explained, not corrected: a 3-D skeleton of a sharply
curved C-loop short-cuts the curve, and at 1.5 mm voxels the duodenal label is
only 79 skeleton voxels long. The measured value is what the
label gives and is reported as such; it was **not** stretched to reach the
published range.

The last two rows are an independent consistency check: a path built from the
envelope's own packing geometry (3.412 m) and a path-free volume argument
(4.56 m) both land inside 3–5 m without either being fitted to the other.

## Surfaces

| surface | faces | watertight | volume (mL) | anchored |
|---|---|---|---|---|
| `small_bowel` | 72,000 | True | 1406.25 | yes |
| `duodenum` | 14,000 | True | 55.03 | yes |
| `colon` | 22,000 | True | 802.16 | yes |
| `stomach` | 12,000 | True | 425.68 | yes |
| `liver` | 18,000 | True | 1423.86 | yes |
| `pancreas` | 8,000 | True | 75.56 | yes |
| **total** | **146,000** | | | |

`small_intestine.glb` is **2.51 MB** (2,631,324 bytes). The parametric lumen tube is
**not** in the glb — it is generated in the browser from `si_centreline.json` so
the peristaltic wave can deform it per-ring on the GPU.

## Structures the release does not contain

Carried here because each one propagates into the vascular and microscopic
layers:

* **no superior mesenteric artery label**, no jejunal or ileal artery labels, no
  arcades, no vasa recta → the arterial course is constructed; only the ostium
  on the real `aorta` label is anchored;
* **no superior mesenteric vein label** → the venous course is constructed; its
  trunk calibre is derived from the liver page's measured portal radius by
  Murray's law rather than invented;
* **no lymphatic label in any open CT segmentation release** → the lacteal /
  collecting-lymphatic / intestinal-trunk tree is entirely synthetic;
* **no jejunum/ileum boundary**, no circular folds, no villi, no crypts, no wall
  layers.

## Sources considered and rejected

| source | why rejected |
|---|---|
| Gibson et al. multi-organ abdominal CT, [10.5281/zenodo.1169361](https://doi.org/10.5281/zenodo.1169361) | the record's own terms stack TCIA Pancreas-CT and BTCV challenge terms on top; redistribution rights unclear |
| any mesh carrying a non-commercial (NC) clause | not redistributable under this repository's MIT licence |
| BodyParts3D | unreachable from this environment and not attribution-only |
| 3D-IRCADb-01 | unreachable from this environment and not attribution-only |
