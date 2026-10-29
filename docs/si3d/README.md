# si3d — 3-D anatomical digital-twin small intestine

Drop this directory into `docs/` next to `si3d.html` and the page works with no
build step and no backend. All paths are relative; `assetBase` is
`./si3d/assets/` and the import map points at `./si3d/vendor/three/three.module.js`,
matching the published liver page.

```
docs/
  si3d.html          ← the page
  si3d/
    js/              scene modules (main, peristalsis, villus, vessels, flow, hud)
    vendor/three/    pinned three.js build + the jsm modules the scene imports
    assets/          every generated artefact the scene reads
    ACCEPTANCE.md    per-item acceptance record, including what was NOT achieved
    GEOMETRY_PROVENANCE.md   source, licence, subject screen, anchoring table
```

## Assets and what produced them

| file | produced by | what it is |
|---|---|---|
| `small_intestine.glb` | `build_si_geometry.py` | the real CT-segmentation surfaces (146 000 triangles, 2.51 MB) |
| `si_geometry.npz` | `build_si_geometry.py` | the same surfaces plus the packed envelope mask and the centreline, for downstream tools |
| `si_centreline.json` | `build_si_geometry.py` | duodenal + jejuno-ileal centreline, arc length, station labels, measured envelope radius, anchoring flags |
| `geometry_report.json` | `build_si_geometry.py` | subject screen, face counts, volumes, measured lengths vs published ranges |
| `small_intestine_vessels.json` | `build_si_vessels.py` | three CCO trees + arcade anastomoses + vasa recta |
| `vessel_report.json` | `build_si_vessels.py` | Murray residuals, clearance audit, vasa-recta statistics, cross-organ coupling |
| `small_intestine_model.json` | `build_si_model.py` | haemodynamics, motility, wall layers, surface amplification, ACAT absorption |
| `si_flow.json` | `build_si_flow.py` | flow paths with real v = Q/A per sample |
| `flow_report.json` | `build_si_flow.py` | per-sample backflow audit, velocities vs published ranges |
| `small_intestine_expression.json` | `build_small_intestine_expression.py` | measured single-cell ADME profile. **If this file is absent the page says so explicitly and shows no substitute.** |
| `small_intestine_census_provenance.json` | same | cells, donors, datasets, assays, read strategy |

## Reproducing

```bash
# geometry  (HTTP range requests against the Zenodo zip; no CT image is fetched)
SI_OUT=web/si3d/assets python tools/build_si_geometry.py
# vessels, physiology, flow  (each reads the previous step's JSON)
SI_OUT=web/si3d/assets python tools/build_si_vessels.py
SI_OUT=web/si3d/assets HEART_TS=ref/heart_timeseries.json python tools/build_si_model.py
SI_OUT=web/si3d/assets python tools/build_si_flow.py
# molecular layer  (CZ CELLxGENE Census, contiguous soma_joinid blocks)
python tools/build_small_intestine_expression.py
# browser audit  (needs a real GPU; see the note in benchmark_si3d.py)
SI_WEB_ROOT=$PWD/web python tools/benchmark_si3d.py --out mobile_report.json --shots shots
# acceptance record, assembled from the reports above
SI_BENCH=mobile_report.json python tools/write_acceptance.py
```

## Scene API

```js
const api = await createIntestineScene({ container, assetBase: './si3d/assets/' });
api.setView('organ' | 'vessels' | 'wall' | 'villus');
api.setPlaybackRate(1.0);          // time scale
api.setParticleDensity(1.0);
api.setStation(-1 | 0 | 1 | 2);    // highlight duodenum / jejunum / ileum
api.setTree('mesenteric_artery' | 'mesenteric_vein' | 'lymphatic', true);
api.setArcades(true);
api.stats();                        // everything the acceptance record quotes
api.hudOverflow();                  // HUD boxes outside the viewport (should be [])
```

`window.__si3d` is the same object, which is what the headless audit drives.

## Honesty contract

The envelope, the duodenal centreline, the two bowel junctions and the molecular
layer are measured. **The whole lumen, every calibre, the station boundaries,
the circular folds, the villi, the crypts, the wall layers and all three vessel
courses are parametric** and carry `contour_anchored: false`. There is no
calibration from expression level to Jmax or CL_int anywhere in this project and
none was invented. Read `ACCEPTANCE.md` §10 before quoting any number from this
page.
