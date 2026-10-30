# Geometry provenance — 3-D anatomical liver

Every surface in `assets/liver.glb` derives from one openly licensed CT
segmentation release. Nothing here is hand-modelled, and nothing is carried over
from the previous 2-D page's Bézier outline.

## Source

| field | value |
|---|---|
| dataset | **Dataset with segmentations of 117 important anatomical structures in 1228 CT images** (TotalSegmentator dataset), Wasserthal J., Breit H.-C., Meyer M. T., Pradella M., Hinck D., Sauter A. W., Heye T., Boll D., Cyriac J., Yang S., Bach M., Segeroth M. |
| version | v2.0.1, published 2023-10-27 |
| DOI | [10.5281/zenodo.10047292](https://doi.org/10.5281/zenodo.10047292) |
| licence | **CC-BY-4.0** (`license.id: "cc-by-4.0"` in the Zenodo record metadata) |
| redistributable in a public repository | **yes** — attribution only, no non-commercial and no share-alike clause |
| subject used | **s0236**, `study_type = ct abdomen`, 1.5 mm isotropic (the release is resampled) |
| structures used | `liver`, `portal_vein_and_splenic_vein`, `inferior_vena_cava`, `gallbladder`, `aorta` |
| access method | HTTP range requests against the 23.5 GB Zenodo zip: the ZIP64 central directory is read from the tail, then only the five mask members (≈ 1 MB total) are fetched. See `tools/remote_zip_reader.py`. |

Attribution required by CC-BY-4.0 is carried in this file, in
`assets/geometry_report.json`, in the scene's own provenance footer and in the
page the scene is embedded in.

### Subject selection

60 abdominal-CT subjects were screened. Selection criteria, all measured, none
cosmetic:

1. liver entirely inside the field of view (no truncation at the volume boundary);
2. capsule volume near the adult mean;
3. largest connected intrahepatic `portal_vein_and_splenic_vein` label;
4. right-hemiliver volume fraction inside the published CT-volumetry range.

s0236 measured 1455.3 mL capsule volume, 1923 intrahepatic portal-vein voxels,
and a right-hemiliver fraction of 0.653.

## Frame

The release's NIfTI sform is diagonal and **RAS+**. This was verified rather than
assumed: the aorta's principal axis is cranio-caudal (singular values
69.8 / 6.7 / 5.4 mm), the aorta sits left of and posterior to the IVC, and
`kidney_right` has a larger x than `kidney_left`. Scene space is
`(patient-left, superior, anterior)` = `diag(-1, ·, ·)` applied as

```
x_scene = -x_RAS   y_scene = z_RAS   z_scene = y_RAS      det = +1
```

so a camera on +Z sees the anterior surface with the right lobe on the viewer's
left, the convention of anatomical illustration.

## From mask to mesh

The release's label maps are **mutually exclusive**: intrahepatic vessels are
carved out of the liver label as open channels that connect to the hilum, so a
3-D hole fill recovers nothing (measured: 0 voxels). The parenchymal envelope is
instead a morphological closing of 7 voxels (10.5 mm — wider than the largest
intrahepatic vessel, narrower than the hilar notch) followed by a hole fill. That
bridges the vessel channels while leaving the porta hepatis open, and adds
6 287 voxels (21.2 mL) to the 1434.1 mL liver label.

Each structure is then: Gaussian-blurred (σ 0.8–1.1 voxel), marching cubes at
the 0.5 isolevel, largest connected component only, Taubin-smoothed
(λ 0.53 / µ −0.54, 10–14 iterations — volume-preserving, unlike Laplacian
smoothing), quadric-decimated, normals fixed.

## Landmarks — all measured on this subject's own labels

| landmark | how it is obtained |
|---|---|
| retrohepatic IVC axis | per-slice centre of mass of the `inferior_vena_cava` label over the liver's cranio-caudal span (99 points); the axis is its first principal direction |
| gallbladder fossa | the `gallbladder` voxels in contact with the parenchymal envelope — the bed, not the centroid of a body that may protrude below the inferior edge |
| **Cantlie plane** | the plane through the retrohepatic IVC axis and the gallbladder fossa centroid |
| porta hepatis | the extrahepatic portal skeleton point closest to the intrahepatic portal skeleton |
| portal bifurcation | the first skeleton node with two children, nearest the porta |
| umbilical fissure | parallel to the Cantlie plane, through the umbilical portion of the real left portal branch (the ligamentum teres lies in that fissure) |
| right / left intersectoral planes | each contains the caval axis and the **azimuthal volume bisector** of the corresponding parenchymal mass — the median of the azimuth distribution about the caval axis, not the mean of unit vectors, which is pulled by distant points and does not halve the volume |

## What the source does NOT contain

Stated plainly, because these gaps propagate into the vascular trees and the
Couinaud approximation:

- **no hepatic-vein label.** TotalSegmentator v2 has no hepatic veins. The
  hepatic venous tree's caval ostia are anchored on the real IVC; its three trunk
  courses are *constructed* from the measured landmark planes, and every report
  marks it `contour_anchored: no`.
- **no hepatic-artery label.** The proper hepatic artery is not a class in this
  release. The arterial tree is synthetic, generated as a portal-triad companion
  of the portal tree (anatomically the artery does run with the portal vein), and
  is marked synthetic everywhere.
- **no bile-duct label.** Same treatment, with retrograde flow.
- **no Couinaud segment labels.** The segmentation is an approximation — see
  `assets/liver_fields_report.json`, which reports per-segment volumes against
  published population ranges and names the four segments that fall outside them.
- **no falciform ligament label.** The groove is detected geometrically, from the
  umbilical-fissure plane combined with measured surface concavity, rather than
  drawn. 256 of 39 002 parenchymal vertices carry a non-zero ligament mask.

## Sources considered and rejected

| candidate | outcome |
|---|---|
| Medical Segmentation Decathlon Task08_HepaticVessel (CC-BY-SA 4.0) | **not obtained.** Its S3 host `msd-for-monai.s3-us-west-2.amazonaws.com` is outside this environment's network allowlist. The archive is also a ~10 GB sequential tar with no central directory, so the range-extraction trick that worked for the Zenodo zip is not available. This is the reason the hepatic venous trunks are constructed rather than contour-anchored. |
| Gibson et al., *Multi-organ Abdominal CT Reference Standard Segmentations*, [10.5281/zenodo.1169361](https://doi.org/10.5281/zenodo.1169361) | **downloaded and rejected.** Zenodo records it as CC-BY-4.0 and it carries liver, gallbladder, IVC and portal/splenic vein for 90 cases, but its own terms-of-use file layers the TCIA Pancreas-CT and BTCV challenge terms on top ("may not be used to develop submissions for the challenge", BTCV requires registration). That ambiguity is unacceptable for redistribution in a public repository, so it was not used. |
| BodyParts3D (CC-BY-SA 2.1 JP) | **not obtained.** `lifesciencedb.jp` is outside this environment's network allowlist. It is also share-alike, which would impose a licence on the repository. |
| 3D-IRCADb-01 (has hepatic veins *and* Couinaud labels) | **not attempted.** `ircad.fr` is outside the allowlist and the redistribution terms are not an attribution-only licence. |

## Reproduce

```bash
python tools/build_liver_geometry.py     # masks -> surfaces, landmarks, frame
python tools/build_liver_vessels.py      # portal / hepatic / arterial / biliary trees
python tools/build_liver_fields.py       # Couinaud, zone field, baked fields, liver.glb
python tools/build_liver_flow.py         # flow paths with v = Q/A, flow audit
```

Only the first step touches the network, and it fetches ≈ 1 MB.
`tools/ts_zip_index.pkl` is the parsed zip central directory, so the 14 MB index
read is not repeated.

---

# v2 additions — what is new, and where each new thing came from

The CT geometry is **unchanged**: still TotalSegmentator v2.0.1 subject `s0236`,
still 137 256 triangles, still CC-BY-4.0. Everything below is layered on top of
it, and each layer states whether it is anchored to that image, derived from it,
or not image-derived at all.

## New assets and their provenance

| asset | derived from | contour-anchored? |
|---|---|---|
| `couinaud_terminals.json` | the 313 portal and 263 hepatic venous terminals of `liver_vessels.json`, each portal terminal given the geometric Couinaud label by replaying `geometric_segment()` from `build_liver_fields.py` verbatim | **portal terminals: yes** (the portal trunks are skeletonised from this subject's own label). **Hepatic venous terminals: no** — the hepatic venous tree is constructed, as in v1 |
| the section plane itself | not a structure — a user-controlled plane in the segmentation's own frame | n/a |
| the cut face (cap) | the **real parenchyma mesh**: the cap is painted only where the two-pass stencil says the plane is inside the segmented solid, so its outline is the segmentation's own cross-section | **yes** — it is the segmentation, sliced |
| vessel lumen openings on the cut face | segment/plane crossings of the four trees in `liver_vessels.json`, radius interpolated between node radii, `radius_scale_applied: 1.0` | **portal: yes. Hepatic venous, arterial, biliary: no** — those trees were already `contour_anchored: false` in v1 and remain so |
| the zone coordinate on the cut face | `d(nearest portal terminal) / (d_portal + d_hepatic-vein terminal)`, the same field baked into `_field0.z` in v1 | derived from the two trees, so it inherits the hepatic venous tree's `false` |
| **single-cell tier geometry** (`js/cells.js`) | **NOT image-derived at all.** A synthetic hepatic-plate slab built at real published human dimensions and then magnified 3 000× for display | **no — `contour_anchored: false`**, and there is no CT at this scale to anchor it to. The 1.5 mm isotropic voxels of the source scan are 68 times a hepatocyte diameter |
| `liver_expression.json`, `liver_expression_normal.json`, `liver_cells.json` + 6 shards, `liver_zonation.json` | CELLxGENE Census 2025-11-08, `tissue_general == 'liver' and is_primary_data == True` | n/a — expression data, not geometry. **From other people's livers**, not from subject `s0236` |
| `liver_protein_evidence.csv` / `.json` | HPA v25 per-gene XML (IHC); PRIDE PXD010154 `proteinGroups.txt`; the 36 `absolute_protein_abundance` rows of `protein_abundance_anchors.csv`; PRIDE PXD062231 `report.gg_matrix.tsv` (coverage only) | n/a — **from other people's livers** |

## The thing that must not be misread

Subject `s0236` supplied the **shape**. Nobody sequenced subject `s0236`. The
expression and protein layers come from entirely different donors — 242 of them in
the all-primary aggregate, 137 in the normal-only one, 18 in the single-cell
proteome. **A cell shown inside this liver's geometry is a real cell from a
different liver**, and the page says so rather than implying one subject was
measured end to end.

## New source files fetched

| source | how it was fetched | size fetched |
|---|---|---|
| CELLxGENE Census 2025-11-08, human liver | TileDB-SOMA over **direct egress on a remote container**, contiguous `soma_joinid` blocks of 150 000 rows, 11 blocks, 0 failures, 430.3 s. The sandbox-local HTTP proxy cannot carry TileDB's S3 fan-out — one block took 67 min through it in v1 | ~1.6 M cells × 61 genes |
| Human Protein Atlas v25, 61 per-gene XML records | `https://www.proteinatlas.org/ENSG*.xml`, one request per gene | 2.5 MB gzipped |
| PRIDE **PXD010154**, `txt/proteinGroups.txt` | **HTTP range reads into a 28.17 GB ZIP64 archive**: read the ZIP64 end-of-central-directory, then the 1 491-byte central directory, then the 71.4 MB deflate stream of that one member and inflated it (280 MB uncompressed). The other 18 members, including a 41.8 GB `allPeptides.txt`, were never transferred | 71.4 MB of 28.17 GB |
| PRIDE **PXD062231**, `search_engine_results/report.gg_matrix.tsv` | the same range trick into a 1.97 GB archive (16 members; the two full `report.tsv` files are 4.9 GB and 3.2 GB uncompressed and were not touched) | 6.1 MB of 1.97 GB |
| PRIDE PXD062231, `Individual*_scores.txt` | direct HTTPS, 18 files | 19 kB |

HPA's per-dataset downloads (`normal_tissue.tsv.zip` and friends) **404 now** —
the download page offers only the whole-atlas `proteinatlas.tsv.zip` / `.xml.gz` /
`.json.gz`. The per-gene XML route was used instead because it carries the
per-cell-type IHC levels, intensities, stained fractions and subcellular locations
that the summary table does not.

## Reproduce — v2 steps

```bash
# geometry (unchanged from v1)
python tools/build_liver_geometry.py
python tools/build_liver_vessels.py
python tools/build_liver_fields.py
python tools/build_liver_flow.py

# v2
python tools/build_section_assets.py            # portal/hepatic terminals for the cut face
python tools/pull_liver_census.py               # RUN THIS ON A DIRECT-EGRESS HOST
python tools/build_liver_molecular.py           # aggregates, per-cell shards, zonation
LIVER_PROTEIN_ANCHORS=<path to protein_abundance_anchors.csv> \
  python tools/build_liver_protein_evidence.py  # four protein tiers
python tools/benchmark_liver3d.py --engine vulkan --seconds 12 --out bench
python tools/build_section_function_reports.py assets bench/benchmark_report.json .
```

`pull_liver_census.py` unsets every proxy variable itself and leaves TileDB at its
own concurrency. Run through a proxy it will appear to hang; that is the bottleneck
v1 hit, not a Census limit.
