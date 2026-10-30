# liver3d — 3-D anatomical liver

A self-contained Three.js scene: real CT-segmentation geometry, four intrahepatic
tubular systems, flow at measured velocities, diaphragmatic motion, and a 3-D
hepatic lobule driven by the project's own single-cell zonal expression.
No backend, no network at run time, three.js r186 vendored under `vendor/`.

## Open it

```
# any static server rooted at web/
python -m http.server 8080 --directory web
# -> http://localhost:8080/liver3d/liver3d_demo.html     standalone demo
# -> http://localhost:8080/liver3d.html                  site page (Chinese)
```

`liver3d_demo.html` opens on its own, `liver3d.html` is the page for the public
site; both import the same `js/main.js`.

## API

```js
import { createLiverScene } from './js/main.js';

const api = await createLiverScene({
  container: document.getElementById('app'),
  assetBase: './assets/',        // default
  view: 'organ',                 // or 'lobule'
  gene: 'CYP1A2',                // any key of liver_model.json .zonal_expression
  particlesPerPath: 320,         // 140 is the narrow-screen default
  excursionCm: 2.0               // respiratory excursion; default = midpoint of 1-3 cm
});

api.setView('organ' | 'lobule');
api.setGene('CYP2D6');                  // rebinds the zonal ramp, direction included
api.setOverlay('segments' | 'zones' | 'ligament', bool);
api.setPartVisible('parenchyma' | 'portal_vein' | 'hepatic_vein' | ..., bool);
api.setPlaybackRate(1.0);
api.setParticleDensity(0.5);
api.setClip(normalVec3, offset, bool);  // cut plane through parenchyma + particles
api.hudOverflow();                      // HUD element boxes outside the viewport
api.frameStats();                       // measured frame times
api.dispose();
```

The signature deliberately mirrors `createHeartScene({container, assetBase, ...})`
from the heart module.

## Modules

| file | what it does |
|---|---|
| `js/main.js` | scene assembly, view switching, zonal gene binding, public API |
| `js/respiration.js` | rigid-body diaphragmatic motion (no parenchymal deformation) |
| `js/materials.js` | procedural Glisson capsule, parenchymal grain, measured falciform groove, Couinaud and zonal overlays |
| `js/vessels.js` | swept-tube meshes for the four trees at their true Murray calibres |
| `js/flow.js` | GPU particle advection, velocity baked per path sample |
| `js/lobule.js` | 3-D hexagonal lobule, portal triads, sinusoids, counter-current canaliculi |
| `js/hud.js` | viewport-responsive in-canvas HUD that reports its own element boxes |

## Assets

| file | contents |
|---|---|
| `assets/liver.glb` | parenchyma / portal vein / IVC / gallbladder with baked `_FIELD0`, `_FIELD1` |
| `assets/liver_vessels.json` | four trees: nodes, parents, Murray radii, labels, planes |
| `assets/liver_flow.json` | 74 flow paths, `v = Q/A` per sample |
| `assets/liver_model.json` | haemodynamics + respiration contract (upstream) |
| `assets/liver_expression.json` | measured single-cell expression, 60 genes × 11 cell types |
| `assets/couinaud.json` | segment planes and per-segment volumes |
| `assets/*_report.json` | the measured numbers behind `ACCEPTANCE.md` |

## Rebuild

```bash
python tools/build_liver_geometry.py   # CT masks -> surfaces + landmarks
python tools/build_liver_vessels.py    # portal / hepatic / arterial / biliary trees
python tools/build_liver_fields.py     # Couinaud, zone field, baked fields, liver.glb
python tools/build_liver_flow.py       # flow paths + flow audit
python tools/benchmark_liver3d.py      # browser measurement + screenshots
```

Only the first step touches the network (≈ 1 MB of range requests to Zenodo).

## Read these before trusting a number

- `GEOMETRY_PROVENANCE.md` — source, DOI, licence, frame verification, and a
  blunt list of what the source does *not* contain.
- `ACCEPTANCE.md` — every acceptance item measured, with the ones that fall
  short marked as such.
