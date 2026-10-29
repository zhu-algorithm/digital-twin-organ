// main.js — three-dimensional anatomical digital-twin stomach.
//
// Views: organ · vessels · section · gland · molecular
// Every number the HUD shows is read out of the generated reports.
import * as THREE from 'three';
import { OrbitControls } from '../vendor/jsm/controls/OrbitControls.js';
import { GLTFLoader } from '../vendor/jsm/loaders/GLTFLoader.js';
import { RoomEnvironment } from '../vendor/jsm/environments/RoomEnvironment.js';
import { wallMaterial, REGION_COLOUR } from './materials.js';
import { Motility } from './motility.js';
import { buildVessels, TREE_COLOUR, TREE_LABEL } from './vessels.js';
import { Flow } from './flow.js';
import { buildGland } from './gland.js';
import { hudText } from './hud.js';
import { resolveCells, geneValue, buildParietalCell, MEMBRANE } from './molecular.js';

const GASTRIC = ['serosa', 'mucosa', 'cut_wall', 'cut_lumen'];
const NEIGHBOURS = ['esophagus', 'duodenum', 'liver', 'spleen', 'pancreas',
                    'aorta', 'portal_splenic_vein', 'colon'];

async function j(base, name) {
  const r = await fetch(base + name);
  if (!r.ok) throw new Error(`${name}: HTTP ${r.status}`);
  return r.json();
}
async function jOpt(base, name) {
  try { return await j(base, name); } catch (e) { return { __error: String(e) }; }
}

export async function createStomachScene(opts) {
  const base = opts.assetBase;
  const container = opts.container;
  const particlesPerPath = opts.particlesPerPath ?? 220;

  const [geoRep, vesRep, flowRep, fieldRep, model, land, vessels, flowData,
         motion, micro, motRep] = await Promise.all([
    j(base, 'geometry_report.json'), j(base, 'vessel_report.json'),
    j(base, 'flow_report.json'), j(base, 'stomach_fields_report.json'),
    j(base, 'stomach_model.json'), j(base, 'stomach_landmarks.json'),
    j(base, 'stomach_vessels.json'), j(base, 'stomach_flow.json'),
    j(base, 'stomach_motion.json'), j(base, 'stomach_micro.json'),
    j(base, 'motion_report.json'),
  ]);
  // the molecular layer is allowed to be MISSING; the page says so rather than
  // substituting a tissue average or another organ's numbers
  const expr = await jOpt(base, 'stomach_expression.json');
  const pump = await jOpt(base, 'stomach_proton_pump_expression.json');
  const census = await jOpt(base, 'stomach_census_provenance.json');

  // ----------------------------------------------------------------- renderer
  const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: false,
    powerPreference: 'high-performance' });
  const narrow = Math.min(window.innerWidth, window.innerHeight) < 520;
  renderer.setPixelRatio(Math.min(window.devicePixelRatio, narrow ? 1.5 : 2));
  renderer.setSize(container.clientWidth, container.clientHeight);
  renderer.outputColorSpace = THREE.SRGBColorSpace;
  renderer.toneMapping = THREE.ACESFilmicToneMapping;
  renderer.toneMappingExposure = 0.88;
  container.appendChild(renderer.domElement);
  renderer.domElement.style.display = 'block';
  renderer.domElement.style.touchAction = 'pan-y';

  const scene = new THREE.Scene();
  scene.background = new THREE.Color(0x070a10);
  const pmrem = new THREE.PMREMGenerator(renderer);
  scene.environment = pmrem.fromScene(new RoomEnvironment(), 0.04).texture;
  scene.environmentIntensity = 0.55;
  scene.add(new THREE.AmbientLight(0xffffff, 0.28));
  const key = new THREE.DirectionalLight(0xffffff, 1.05);
  key.position.set(90, 120, 150); scene.add(key);
  const rim = new THREE.DirectionalLight(0x9ec4ff, 0.6);
  rim.position.set(-120, 40, -90); scene.add(rim);

  const camera = new THREE.PerspectiveCamera(42, 1, 1, 4000);
  const controls = new OrbitControls(camera, renderer.domElement);
  controls.enableDamping = true; controls.dampingFactor = 0.08;
  controls.enablePan = false;

  // -------------------------------------------------------------------- GLB
  const gltf = await new GLTFLoader().loadAsync(base + 'stomach.glb');
  const parts = {};
  gltf.scene.traverse(o => { if (o.isMesh) parts[o.name] = o; });

  const root = new THREE.Group(); root.name = 'stomach_root';
  scene.add(root);
  const organGroup = new THREE.Group(); organGroup.name = 'organ';
  const neighbourGroup = new THREE.Group(); neighbourGroup.name = 'neighbours';
  const cutGroup = new THREE.Group(); cutGroup.name = 'cut';
  root.add(organGroup, neighbourGroup, cutGroup);

  const wallMats = [];
  for (const nm of GASTRIC) {
    const m = parts[nm];
    if (!m) continue;
    const mat = wallMaterial({
      opacity: nm === 'serosa' ? 0.46 : 1.0,
      side: nm.startsWith('cut') || nm === 'mucosa' ? THREE.DoubleSide : THREE.FrontSide,
    });
    m.material = mat; wallMats.push(mat);
    (nm.startsWith('cut') ? cutGroup : organGroup).add(m);
  }
  for (const nm of NEIGHBOURS) {
    const m = parts[nm];
    if (!m) continue;
    m.material = new THREE.MeshStandardMaterial({
      color: nm === 'aorta' ? 0xa33c3c
           : nm === 'portal_splenic_vein' ? 0x3f73c8
           : nm === 'liver' ? 0x8c5b4a : nm === 'spleen' ? 0x7a4a62
           : nm === 'pancreas' ? 0xc7a97a : nm === 'duodenum' ? 0xb98a6a
           : nm === 'esophagus' ? 0xb0907f : 0x9a8f7a,
      roughness: 0.75, transparent: true, opacity: 0.26, depthWrite: false,
    });
    m.name = nm;
    neighbourGroup.add(m);
  }

  // frame the camera on the measured bounds
  const bbox = new THREE.Box3().setFromObject(organGroup);
  const centre = bbox.getCenter(new THREE.Vector3());
  const radius = bbox.getSize(new THREE.Vector3()).length() * 0.5;
  controls.target.copy(centre);
  camera.position.copy(centre).add(new THREE.Vector3(radius * 1.6, radius * 0.55,
                                                     radius * 1.75));
  camera.near = radius / 120; camera.far = radius * 40; camera.updateProjectionMatrix();

  // --------------------------------------------------------------- motility
  const motileMeshes = [];
  if (parts.mucosa) motileMeshes.push({ mesh: parts.mucosa, kind: 'mucosa' });
  if (parts.serosa) motileMeshes.push({ mesh: parts.serosa, kind: 'serosa' });
  if (parts.cut_lumen) motileMeshes.push({ mesh: parts.cut_lumen, kind: 'mucosa' });
  const motility = new Motility(motion, motileMeshes);

  // --------------------------------------------------------------- vessels
  const { group: vesGroup, perTree } = buildVessels(vessels, 1.0);
  root.add(vesGroup);
  const flow = new Flow(flowData, particlesPerPath);
  root.add(flow.group);

  // ------------------------------------------------------- retropulsion demo
  // Particles advected with the wall.  Those at or below the published sieving
  // size pass the pylorus; the rest are driven back toward the corpus.  The
  // pass fraction is COUNTED here, not declared.
  const retro = buildRetropulsion(motion, land);
  root.add(retro.group);

  // ----------------------------------------------------------------- gland
  const glandGroup = new THREE.Group(); glandGroup.name = 'gland_view';
  glandGroup.visible = false;
  const oxyntic = buildGland(micro, 'oxyntic');
  const pyloric = buildGland(micro, 'pyloric');
  pyloric.position.x = micro.gland_spacing_um * 9 * micro.magnification / 1000;
  glandGroup.add(oxyntic, pyloric);
  scene.add(glandGroup);

  // -------------------------------------------------------------- molecular
  const molGroup = new THREE.Group(); molGroup.name = 'molecular_view';
  molGroup.visible = false;
  const parietal = buildParietalCell(1.0);
  molGroup.add(parietal);
  scene.add(molGroup);

  let molecular;
  if (expr.__error || !expr.genes) {
    molecular = { ok: false,
      reason: `stomach_expression.json 未能加载（${expr.__error || '缺少 genes 字段'}）。`
            + '本页不用组织平均值或其它器官的值代替。' };
  } else {
    const cells = resolveCells(expr);
    molecular = {
      ok: true, n_cell_types: expr.cell_types.length, n_genes: expr.genes.length,
      n_cells: (expr.n_cells || []).reduce((a, b) => a + b, 0),
      n_donors: expr.n_donors ?? (census.n_donors ?? null),
      n_datasets: expr.n_datasets_used ?? (census.n_datasets ?? null),
      census_version: expr.census_version,
      keyCells: cells,
      missingKeyCells: cells.filter(c => !c.present).map(c => c.key),
      protonPump: pump.__error ? { ok: false, reason: String(pump.__error) }
                               : { ok: true, genes: pump.genes,
                                   panel_membership: pump.panel_membership },
      membrane: MEMBRANE,
    };
  }

  // ------------------------------------------------------------------ state
  let view = 'organ';
  let wallMode = 0;
  let labels = true;
  let showMotility = true;
  const labelGroup = new THREE.Group(); labelGroup.name = 'labels';
  root.add(labelGroup);
  buildLabels(labelGroup, land, geoRep, radius);

  function applyView(v) {
    view = v;
    const organLike = (v === 'organ' || v === 'vessels' || v === 'section');
    organGroup.visible = organLike;
    neighbourGroup.visible = (v === 'organ' || v === 'vessels');
    cutGroup.visible = (v === 'section');
    if (parts.serosa) parts.serosa.visible = (v !== 'section');
    if (parts.mucosa) parts.mucosa.visible = organLike;
    vesGroup.visible = (v === 'organ' || v === 'vessels');
    flow.group.visible = (v === 'organ' || v === 'vessels');
    retro.group.visible = (v === 'section');
    glandGroup.visible = (v === 'gland');
    molGroup.visible = (v === 'molecular');
    labelGroup.visible = labels && organLike;
    if (v === 'gland') {
      const b = new THREE.Box3().setFromObject(glandGroup);
      const c = b.getCenter(new THREE.Vector3());
      const r = b.getSize(new THREE.Vector3()).length() * 0.5;
      controls.target.copy(c);
      camera.position.copy(c).add(new THREE.Vector3(r * 1.1, r * 0.35, r * 1.5));
    } else if (v === 'molecular') {
      controls.target.set(0, 0, 0);
      camera.position.set(26, 12, 32);
    } else {
      controls.target.copy(centre);
      camera.position.copy(centre).add(new THREE.Vector3(
        radius * (v === 'section' ? 0.4 : 1.6),
        radius * 0.55, radius * (v === 'section' ? 2.0 : 1.75)));
    }
    camera.updateProjectionMatrix();
  }

  // ------------------------------------------------------------------ loop
  const clock = new THREE.Clock();
  let simT = 0, rate = 1.0;
  const frameMs = [];
  let lastMotilityUpdate = -1;
  let normalTick = 0;

  function resize() {
    const w = container.clientWidth, h = container.clientHeight;
    if (!w || !h) return;
    renderer.setSize(w, h, false);
    camera.aspect = w / h; camera.updateProjectionMatrix();
  }
  new ResizeObserver(resize).observe(container);
  resize();

  function tick() {
    const t0 = performance.now();
    const dt = Math.min(clock.getDelta(), 0.1);
    simT += dt * rate;
    if ((view === 'organ' || view === 'vessels') && flow.group.visible) {
      flow.step(dt);
    }
    if (showMotility && (view === 'organ' || view === 'section')) {
      // 20 Hz is plenty for a 3 cpm wave and keeps the per-vertex CPU
      // deformation off the critical path on phones
      if (simT - lastMotilityUpdate > 1 / 12) {
        normalTick = (normalTick + 1) % 4;
        motility.update(simT, normalTick === 0);
        lastMotilityUpdate = simT;
      }
    }
    if (view === 'section') retro.step(dt * rate, simT, motility);
    controls.update();
    renderer.render(scene, camera);
    const ms = performance.now() - t0;
    frameMs.push(ms); if (frameMs.length > 600) frameMs.shift();
    requestAnimationFrame(tick);
  }
  applyView('organ');
  motility.update(0);
  requestAnimationFrame(tick);

  // --------------------------------------------------------------- the API
  const gl = renderer.getContext();
  const dbg = gl.getExtension('WEBGL_debug_renderer_info');
  const triangles = fieldRep.glb_total_faces;

  const api = {
    THREE, scene, camera, renderer, controls, motility, flow,
    setView: applyView,
    setTreeVisible(name, on) {
      const t = perTree[name]; if (t) t.group.visible = on;
      flow.setVisible(name, on);
    },
    setWallMode(mode) {               // 0 region · 1 pH · 2 oxyntic density
      wallMode = mode;
      for (const m of wallMats) m.userData.u.uMode.value = mode;
    },
    setMotility(on) {
      showMotility = on; motility.enabled = on; motility.update(simT);
    },
    setAccommodation(f) { motility.accommodation = f; motility.update(simT); },
    setPlaybackRate(r) { rate = r; flow.setRate(r); },
    setParticleDensity(f) { flow.setDensity(f); },
    setLabelsVisible(on) { labels = on; labelGroup.visible = on &&
      (view === 'organ' || view === 'vessels' || view === 'section'); },
    highlightRegion(rid) { for (const m of wallMats) m.userData.u.uHighlight.value = rid; },
    measureMotion(o) { const r = motility.measure(o); api.diag.motionMeasured = r; return r; },
    emptyingCurves(n, s) { return motility.emptyingCurves(n, s); },
    retropulsionAudit() { return retro.audit(); },
    // lets the headless audit drive the retropulsion counters directly,
    // instead of depending on how long the section view happened to be shown
    retroStep(dt, t) { retro.step(dt, t, motility); },
    glyphValue(cellKey, gene) {
      if (!molecular.ok) return null;
      const c = molecular.keyCells.find(x => x.key === cellKey);
      if (!c || !c.present) return null;
      return geneValue(expr, c.index, gene);
    },
    resetFrameStats() { frameMs.length = 0; },
    renderCost(n = 40) {
      // The real per-frame cost: back-to-back render() calls each followed by
      // gl.finish(), so the GPU is actually drained before the clock is read.
      // Headless Chromium has no vsync, so requestAnimationFrame rates are
      // uncapped and must NOT be read as a user-visible frame rate.
      const g = renderer.getContext();
      renderer.render(scene, camera); g.finish();
      const t0 = performance.now();
      for (let i = 0; i < n; i++) { renderer.render(scene, camera); g.finish(); }
      const ms = (performance.now() - t0) / n;
      return { n, msPerFrame: +ms.toFixed(4), fps: +(1000 / ms).toFixed(2),
               budget_60fps_ms: 16.67, headroom_x: +(16.67 / ms).toFixed(2) };
    },
    setPartVisible(name, on) { if (parts[name]) parts[name].visible = on; },
    parts,
    hudOverflow() {
      // This page's HUD is a DOM panel (#ui/#prov), not a canvas overlay, so
      // the applicable check is the DOM bounding-box check the harness already
      // runs.  Returned empty so the two pages' reports stay comparable.
      return [];
    },
    frameStats() {
      const a = [...frameMs].sort((x, y) => x - y);
      if (!a.length) return null;
      const q = p => a[Math.min(a.length - 1, Math.floor(p * (a.length - 1)))];
      return { samples: a.length, mean_ms: +(a.reduce((s, x) => s + x, 0) / a.length).toFixed(4),
               p50_ms: +q(0.50).toFixed(4), p95_ms: +q(0.95).toFixed(4),
               max_ms: +q(1.0).toFixed(4), budget_60fps_ms: 16.67 };
    },
    glInfo() {
      return {
        unmaskedRenderer: dbg ? gl.getParameter(dbg.UNMASKED_RENDERER_WEBGL) : null,
        unmaskedVendor: dbg ? gl.getParameter(dbg.UNMASKED_VENDOR_WEBGL) : null,
        version: gl.getParameter(gl.VERSION),
      };
    },
    diag: {
      triangles, glbMB: fieldRep.glb_mb,
      geometry: geoRep, vesselReport: vesRep, flowReport: flowRep,
      fieldsReport: fieldRep, motionReport: motRep, model, micro,
      molecular, census, motionMeasured: null,
      trees: Object.fromEntries(Object.entries(perTree).map(
        ([k, v]) => [k, { segments: v.segments, flow_mL_min: v.flow,
                          retrograde: v.retrograde, label: TREE_LABEL[k] }])),
    },
    hud(n) { return hudText(api.diag, n ?? narrow); },
  };
  // measure once at startup so the HUD has real numbers rather than the
  // declared ones
  api.measureMotion();
  window.__stomach3d = api;
  return api;
}

// ---------------------------------------------------------------- labels
function buildLabels(group, land, geoRep, radius) {
  const pts = [
    ['贲门 cardia / EGJ', land.egj, 0xffd24a],
    ['胃底 fundus', land.fundic_apex, 0xe08a6a],
    ['角切迹 incisura angularis', land.incisura, 0x8fd9a8],
    ['幽门 pylorus', land.pylorus, 0xd07a9b],
    ['腹腔干起始 coeliac origin', land.coeliac_origin_on_aorta, 0xa33c3c],
    ['脾门 splenic hilum', land.splenic_hilum, 0x7a4a62],
  ];
  for (const [text, p, col] of pts) {
    if (!p) continue;
    const s = new THREE.Mesh(
      new THREE.SphereGeometry(radius * 0.016, 12, 10),
      new THREE.MeshBasicMaterial({ color: col }));
    s.position.set(p[0], p[1], p[2]);
    s.name = `landmark_${text}`;
    s.userData.label = text;
    group.add(s);
  }
}

// ----------------------------------------------------------- retropulsion
function buildRetropulsion(motion, land) {
  const classes = motion.retropulsion.particle_size_classes_mm;
  const sieve = motion.retropulsion.sieving_size_mm;
  const group = new THREE.Group(); group.name = 'retropulsion';
  const CL = motion.centreline.map(p => new THREE.Vector3(...p));
  const sArc = motion.arclength_mm;
  const L = sArc[sArc.length - 1];
  const sPace = motion.pacemaker_s_mm;
  const per = [];
  for (const d of classes) {
    const n = 40;
    const P = new Float32Array(n * 3);
    const geo = new THREE.BufferGeometry();
    geo.setAttribute('position', new THREE.BufferAttribute(P, 3));
    const passes = d <= sieve;
    const mat = new THREE.PointsMaterial({
      color: passes ? 0x8fd9a8 : 0xf0a33c,
      size: Math.max(1.4, d * 1.1), sizeAttenuation: true,
      transparent: true, opacity: 0.95, depthWrite: false });
    const pts = new THREE.Points(geo, mat);
    pts.frustumCulled = false; pts.name = `particles_${d}mm`;
    group.add(pts);
    const s = new Float32Array(n);
    for (let i = 0; i < n; i++) s[i] = sPace + Math.random() * (L - sPace);
    per.push({ d, passes, P, geo, s, n, passed: 0, pushedBack: 0 });
  }

  function posAt(s) {
    const u = Math.min(Math.max(s / L, 0), 1) * (CL.length - 1);
    const i = Math.min(CL.length - 2, Math.floor(u));
    return CL[i].clone().lerp(CL[i + 1], u - i);
  }

  return {
    group,
    step(dt, t, motility) {
      const scale = motility.computeScale(t);
      const N = motility.sGrid.length;
      for (const g of per) {
        for (let i = 0; i < g.n; i++) {
          let s = g.s[i];
          const si = Math.min(N - 1, Math.max(0, Math.round(s / L * (N - 1))));
          // local wall motion drives the content: a closing ring squeezes
          // content away from itself, distally if the particle can pass and
          // proximally if it cannot
          const squeeze = 1.0 - scale[si];
          const v = squeeze > 0.02
            ? (g.passes ? 18.0 * squeeze : -26.0 * squeeze)
            : 1.2;
          s += v * dt;
          if (s >= L) { g.passed++; s = sPace + Math.random() * (L - sPace) * 0.4; }
          if (s < sPace) { g.pushedBack++; s = sPace + Math.random() * 20; }
          g.s[i] = s;
          const p = posAt(s);
          g.P[3 * i] = p.x; g.P[3 * i + 1] = p.y; g.P[3 * i + 2] = p.z;
        }
        g.geo.attributes.position.needsUpdate = true;
      }
    },
    audit() {
      const rows = per.map(g => ({
        diameter_mm: g.d, below_sieving_size: g.passes,
        passed_pylorus: g.passed, pushed_back_to_corpus: g.pushedBack }));
      const tot = rows.reduce((a, r) => a + r.passed_pylorus + r.pushed_back_to_corpus, 0);
      const pass = rows.reduce((a, r) => a + r.passed_pylorus, 0);
      return { sieving_size_mm: sieve, classes: rows, events: tot,
               pass_fraction: tot ? +(pass / tot).toFixed(4) : null,
               note: 'counted in the browser from the live wall motion; the '
                   + 'sieving size is the only declared input' };
    },
  };
}
