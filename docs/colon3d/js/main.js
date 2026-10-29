import * as THREE from 'three';
import { OrbitControls } from '../vendor/jsm/controls/OrbitControls.js';
import { GLTFLoader } from '../vendor/jsm/loaders/GLTFLoader.js';
import { ColonTube, taeniaCurves, foldRings, appendicesMesh, vec } from './colon.js';
import { Motility, MODES } from './motility.js';
import { buildTree, FlowParticles, TREE_STYLE, arcadeLine } from './vessels.js';
import { buildCryptUnit, buildVillusGhost } from './crypt.js';
import { buildColonocyte, MICROBIAL_LINES } from './molecular.js';

const BASE = './colon3d/assets/';
const errBox = document.getElementById('err');
function fail(e) {
  errBox.style.display = 'block';
  errBox.textContent = (e && e.stack) ? e.stack : String(e);
  window.__colonError = String(e);
}
window.addEventListener('error', ev => fail(ev.error || ev.message));
window.addEventListener('unhandledrejection', ev => fail(ev.reason));

const J = async n => (await fetch(BASE + n)).json();

const state = {
  view: 'organ', rate: 1, dens: 1,
  show: { art: true, ven: true, arc: true, seg: false, lab: true,
          taenia: true, haust: true, epi: true, fold: true },
  mode: 'haustral'
};

(async function init() {
  const [sc, micro, geomRep, vesRep, flowRep, taeRep, vessels, flow] = await Promise.all([
    J('colon_scene.json'), J('colon_micro.json'), J('geometry_report.json'),
    J('vessel_report.json'), J('flow_report.json'), J('colon_taeniae.json'),
    J('colon_vessels.json'), J('colon_flow.json')]);
  window.__colonData = { sc, micro, geomRep, vesRep, flowRep, taeRep };

  const app = document.getElementById('app');
  const renderer = new THREE.WebGLRenderer({ antialias: true, alpha: false,
                                             powerPreference: 'high-performance' });
  renderer.setPixelRatio(Math.min(devicePixelRatio, 2));
  renderer.setSize(app.clientWidth, app.clientHeight);
  renderer.outputColorSpace = THREE.SRGBColorSpace;
  app.appendChild(renderer.domElement);
  const scene = new THREE.Scene();
  scene.background = new THREE.Color(0x070a10);
  const camera = new THREE.PerspectiveCamera(42, app.clientWidth / app.clientHeight, 1, 6000);
  const controls = new OrbitControls(camera, renderer.domElement);
  controls.enableDamping = true; controls.dampingFactor = 0.08;
  scene.add(new THREE.HemisphereLight(0xcfe0ff, 0x101620, 1.15));
  const key = new THREE.DirectionalLight(0xffffff, 1.5); key.position.set(1, 1.2, 1.1);
  scene.add(key);
  const fill = new THREE.DirectionalLight(0x9fc0ff, 0.55); fill.position.set(-1, -0.4, -0.8);
  scene.add(fill);

  // ---------------- groups -------------------------------------------------
  const gOrgan = new THREE.Group(), gHaustra = new THREE.Group(),
        gMotility = new THREE.Group(), gCrypt = new THREE.Group(), gMol = new THREE.Group();
  [gOrgan, gHaustra, gMotility, gCrypt, gMol].forEach(g => scene.add(g));

  const centre = vec(sc.centre);

  // ---- real segmentation surface (organ view) ------------------------------
  let glbMesh = null;
  await new Promise(res => {
    new GLTFLoader().load(BASE + 'colon.glb', gl => {
      gl.scene.traverse(o => {
        if (o.isMesh) {
          o.material = new THREE.MeshStandardMaterial({
            color: 0xb98a74, roughness: 0.82, metalness: 0.02,
            transparent: true, opacity: 0.42, side: THREE.DoubleSide,
            depthWrite: false });
          glbMesh = o;
        }
      });
      gOrgan.add(gl.scene);
      res();
    }, undefined, e => { fail(e); res(); });
  });

  // ---- segment colouring overlay ------------------------------------------
  const SEGCOL = { caecum: 0xff6b6b, ascending: 0xffa94d, hepatic_flexure: 0xffd43b,
                   transverse: 0x8ce99a, splenic_flexure: 0x63e6be, descending: 0x4dabf7,
                   sigmoid: 0xb197fc, rectum: 0xf783ac };
  const segTubes = new THREE.Group();
  for (const sg of sc.segments) {
    const pts = [];
    for (let i = sg.i0; i <= sg.i1; i++) pts.push(vec(sc.axis[i]));
    if (pts.length < 2) continue;
    const curve = new THREE.CatmullRomCurve3(pts);
    const t = new THREE.Mesh(
      new THREE.TubeGeometry(curve, Math.max(8, pts.length), 3.2, 8, false),
      new THREE.MeshStandardMaterial({ color: SEGCOL[sg.name] || 0xffffff, roughness: 0.5,
        emissive: SEGCOL[sg.name] || 0xffffff, emissiveIntensity: 0.25 }));
    segTubes.add(t);
  }
  segTubes.visible = false;
  gOrgan.add(segTubes);

  // ---- taeniae -------------------------------------------------------------
  const TAE = { mesocolica: 0x6ee7b7, omentalis: 0xfcd34d, libera: 0xfb7185 };
  function taeniaGroup(radiusScale) {
    const g = new THREE.Group();
    for (const w of ['mesocolica', 'omentalis', 'libera']) {
      const pts = taeniaCurves(sc, w);
      const curve = new THREE.CatmullRomCurve3(pts);
      const m = new THREE.Mesh(
        new THREE.TubeGeometry(curve, pts.length, (micro ? 4.0 : 4.0) * radiusScale, 7, false),
        new THREE.MeshStandardMaterial({ color: TAE[w], roughness: 0.45,
          emissive: TAE[w], emissiveIntensity: 0.3 }));
      m.name = 'taenia_' + w;
      g.add(m);
    }
    return g;
  }
  const taeOrgan = taeniaGroup(1.0); gOrgan.add(taeOrgan);
  const taeHaustra = taeniaGroup(1.0); gHaustra.add(taeHaustra);

  // ---- vessels -------------------------------------------------------------
  const treeGroups = {};
  for (const name of Object.keys(TREE_STYLE)) {
    if (!vessels[name]) continue;
    const g = buildTree(name, vessels[name], flow[name]);
    treeGroups[name] = g; gOrgan.add(g);
  }
  const arcade = arcadeLine(vessels._marginal_arcade, 0xffe066, 3);
  const riolan = arcadeLine(vessels._arc_of_riolan, 0xff9f1c, 3);
  gOrgan.add(arcade); gOrgan.add(riolan);
  const particles = new FlowParticles(
    Object.fromEntries(Object.keys(TREE_STYLE).filter(k => vessels[k]).map(k => [k, vessels[k]])),
    flow, n => TREE_STYLE[n].colour);
  gOrgan.add(particles.points);

  // ---- haustra view ---------------------------------------------------------
  const tube = new ColonTube(sc, { ncirc: 72 });
  const tubeMat = new THREE.MeshStandardMaterial({ color: 0xc98f76, roughness: 0.78,
    metalness: 0.02, side: THREE.DoubleSide, flatShading: false });
  const tubeMesh = new THREE.Mesh(tube.geom, tubeMat);
  gHaustra.add(tubeMesh);
  const foldGroup = new THREE.Group();
  for (const pts of foldRings(sc)) {
    const g = new THREE.BufferGeometry().setFromPoints(pts);
    foldGroup.add(new THREE.Line(g, new THREE.LineBasicMaterial({
      color: 0x5b3a2e, transparent: true, opacity: 0.9 })));
  }
  gHaustra.add(foldGroup);
  const epiMat = new THREE.MeshStandardMaterial({ color: 0xf5d98b, roughness: 0.85 });
  const epiMesh = appendicesMesh(sc, epiMat);
  gHaustra.add(epiMesh);
  // appendix base marker (constructed — the appendix is NOT in the label set)
  const appx = new THREE.Mesh(new THREE.SphereGeometry(5, 14, 10),
    new THREE.MeshStandardMaterial({ color: 0xff5d55, emissive: 0x6b1510 }));
  appx.position.copy(vec(sc.appendix_base));
  gHaustra.add(appx);

  // ---- motility view ---------------------------------------------------------
  const tube2 = new ColonTube(sc, { ncirc: 48 });
  const tube2Mesh = new THREE.Mesh(tube2.geom, new THREE.MeshStandardMaterial({
    color: 0xcf9880, roughness: 0.8, side: THREE.DoubleSide }));
  gMotility.add(tube2Mesh);
  const taeMot = taeniaGroup(0.8); gMotility.add(taeMot);
  const motility = new Motility(tube2, sc, micro);

  // luminal content markers, advanced at the DECLARED display gain
  const NC = 160;
  const contentGeom = new THREE.BufferGeometry();
  const cpos = new Float32Array(NC * 3);
  contentGeom.setAttribute('position', new THREE.BufferAttribute(cpos, 3));
  const content = new THREE.Points(contentGeom, new THREE.PointsMaterial({
    color: 0x9b7653, size: 5, sizeAttenuation: false, transparent: true, opacity: 0.9 }));
  gMotility.add(content);
  const contentU = new Float32Array(NC);
  for (let i = 0; i < NC; i++) contentU[i] = sc.arclen[sc.arclen.length - 1] * i / NC;

  // ---- crypt view -------------------------------------------------------------
  const MAG = micro.crypt.display_magnification;
  const crypt = buildCryptUnit(micro, MAG);
  gCrypt.add(crypt);
  gCrypt.add(buildVillusGhost(MAG));

  // ---- molecular view -----------------------------------------------------------
  let exprData = null;
  try {
    exprData = await J('colon_expression.json');
  } catch (e) { exprData = null; }
  const molGroup = new THREE.Group(); gMol.add(molGroup);
  function buildMol() {
    molGroup.clear();
    if (!exprData) {
      const warn = sprite('分子层数据缺失 · 未用任何代理值填补 · 见下表', '#ff9a8a');
      warn.scale.multiplyScalar(0.9);
      molGroup.add(warn);
      return null;
    }
    const want = ['colonocyte', 'enterocyte of epithelium of large intestine',
                  'intestinal epithelial cell', 'epithelial cell of large intestine'];
    let ct = exprData.cell_types.find(c => want.includes(c));
    if (!ct) ct = exprData.cell_types[exprData.n_cells.indexOf(Math.max(...exprData.n_cells))];
    const exprMap = {};
    exprData.cell_types.forEach((c, i) => exprMap[c] = exprData.expr[i]);
    const g = buildColonocyte(exprMap, exprData.genes, ct);
    molGroup.add(g);
    return { cell_type: ct, placed: g.userData.placed };
  }
  const molInfo = buildMol();
  window.__colonMolecular = molInfo;

  // ---------------- view switching ------------------------------------------
  const ROWS = { organ: 'r-organ', haustra: 'r-haustra', motility: 'r-motility' };
  function frameView(v) {
    state.view = v;
    gOrgan.visible = v === 'organ';
    gHaustra.visible = v === 'haustra';
    gMotility.visible = v === 'motility';
    gCrypt.visible = v === 'crypt';
    gMol.visible = v === 'mol';
    for (const [k, id] of Object.entries(ROWS)) {
      const el = document.getElementById(id);
      if (el) el.style.display = (k === v) ? 'flex' : 'none';
    }
    const box = new THREE.Box3();
    const g = { organ: gOrgan, haustra: gHaustra, motility: gMotility, crypt: gCrypt, mol: gMol }[v];
    box.setFromObject(g);
    const c = box.getCenter(new THREE.Vector3());
    const r = Math.max(box.getSize(new THREE.Vector3()).length() * 0.5, 10);
    controls.target.copy(c);
    camera.position.copy(c.clone().add(new THREE.Vector3(0.62, 0.42, 1).normalize()
      .multiplyScalar(r * 2.25)));
    camera.near = r / 120; camera.far = r * 40; camera.updateProjectionMatrix();
    controls.update();
  }

  const on = (id, fn) => { const e = document.getElementById(id); if (e) e.onclick = fn; };
  const vbtn = { organ: 'v-organ', haustra: 'v-haustra', motility: 'v-motility',
                 crypt: 'v-crypt', mol: 'v-mol' };
  for (const [v, id] of Object.entries(vbtn)) {
    on(id, () => {
      Object.values(vbtn).forEach(x => document.getElementById(x).classList.remove('on'));
      document.getElementById(id).classList.add('on');
      frameView(v);
    });
  }
  const tog = (id, k, apply) => on(id, () => {
    state.show[k] = !state.show[k];
    document.getElementById(id).classList.toggle('on', state.show[k]);
    apply();
  });
  tog('b-art', 'art', () => ['sma_colic', 'ima_colic', 'rectal_systemic']
    .forEach(n => { if (treeGroups[n]) treeGroups[n].visible = state.show.art; }));
  tog('b-ven', 'ven', () => ['smv_portal', 'imv_portal', 'rectal_vein_systemic']
    .forEach(n => { if (treeGroups[n]) treeGroups[n].visible = state.show.ven; }));
  tog('b-arc', 'arc', () => { arcade.visible = state.show.arc; riolan.visible = state.show.arc; });
  tog('b-seg', 'seg', () => { segTubes.visible = state.show.seg; });
  tog('b-lab', 'lab', () => { labelGroup.visible = state.show.lab; });
  tog('b-taenia', 'taenia', () => { taeHaustra.visible = state.show.taenia; });
  tog('b-haust', 'haust', () => { tube.build(0, state.show.haust); });
  tog('b-epi', 'epi', () => { epiMesh.visible = state.show.epi; });
  tog('b-fold', 'fold', () => { foldGroup.visible = state.show.fold; });
  const mbtn = { haustral: 'm-haustral', mass: 'm-mass', defecation: 'm-defec' };
  for (const [m, id] of Object.entries(mbtn)) {
    on(id, () => {
      Object.values(mbtn).forEach(x => document.getElementById(x).classList.remove('on'));
      document.getElementById(id).classList.add('on');
      state.mode = m; motility.setMode(m);
    });
  }
  const rate = document.getElementById('rate'), dens = document.getElementById('dens');
  rate.oninput = () => { state.rate = +rate.value;
    document.getElementById('rate-v').innerHTML = `<code>${state.rate.toFixed(2)}×</code>`; };
  dens.oninput = () => { state.dens = +dens.value;
    document.getElementById('dens-v').innerHTML = `<code>${Math.round(state.dens * 100)}%</code>`;
    particles.points.material.size = 2.6 * (0.5 + state.dens * 0.5); };
  const pb = document.getElementById('panelbtn');
  if (pb) pb.onclick = () => document.getElementById('ui').classList.toggle('collapsed');

  // ---------------- labels ----------------------------------------------------
  const labelGroup = new THREE.Group(); scene.add(labelGroup);
  function sprite(text, colour) {
    const c = document.createElement('canvas');
    const ctx = c.getContext('2d');
    ctx.font = '600 30px -apple-system,PingFang SC,sans-serif';
    const w = Math.ceil(ctx.measureText(text).width) + 20;
    c.width = w; c.height = 42;
    const g = c.getContext('2d');
    g.font = '600 30px -apple-system,PingFang SC,sans-serif';
    g.fillStyle = 'rgba(8,12,20,.78)'; g.fillRect(0, 0, w, 42);
    g.fillStyle = colour || '#e9eef6'; g.textBaseline = 'middle';
    g.fillText(text, 10, 22);
    const t = new THREE.CanvasTexture(c);
    t.colorSpace = THREE.SRGBColorSpace;
    const s = new THREE.Sprite(new THREE.SpriteMaterial({ map: t, transparent: true,
      depthTest: false, depthWrite: false }));
    s.scale.set(w * 0.38, 16, 1);
    return s;
  }
  const mid = sg => vec(sc.axis[Math.round((sg.i0 + sg.i1) / 2)]);
  // labels are sprites inside the canvas, so the DOM overflow checks cannot see them
  // colliding.  Stagger them vertically and push each one out along the local
  // circumferential direction so adjacent segment captions do not stack on each other.
  const lift = (i, base) => base + (i % 3) * 20;
  const outward = i => {
    const N = sc.normal[i], B = sc.binormal[i], th = sc.theta_free[i];
    return new THREE.Vector3(Math.cos(th) * N[0] + Math.sin(th) * B[0],
                             Math.cos(th) * N[1] + Math.sin(th) * B[1],
                             Math.cos(th) * N[2] + Math.sin(th) * B[2]);
  };
  const NAMES = { caecum: '盲肠', ascending: '升结肠', hepatic_flexure: '肝曲',
                  transverse: '横结肠', splenic_flexure: '脾曲', descending: '降结肠',
                  sigmoid: '乙状结肠', rectum: '直肠' };
  sc.segments.forEach((sg, k) => {
    const i = Math.round((sg.i0 + sg.i1) / 2);
    const s = sprite(`${NAMES[sg.name] || sg.name} ${sg.length_mm.toFixed(0)} mm`,
                     '#' + (SEGCOL[sg.name] || 0xffffff).toString(16).padStart(6, '0'));
    s.position.copy(mid(sg)).addScaledVector(outward(i), 26)
     .add(new THREE.Vector3(0, lift(k, 18), 0));
    labelGroup.add(s);
  });
  const iWs = vesRep.territories.watershed_axis_index;
  const ws = sprite('脾曲分水岭 (Griffiths) · SMA|IMA', '#ffe066');
  ws.position.copy(vec(sc.axis[iWs])).addScaledVector(outward(iWs), 34)
    .add(new THREE.Vector3(0, 86, 0));
  labelGroup.add(ws);
  const iRs = sc.segments.find(x => x.name === 'rectum').i0;
  const sud = sprite('直乙交界分水岭 (Sudeck) · 结肠带在此散开', '#ff9f1c');
  sud.position.copy(vec(sc.axis[iRs])).addScaledVector(outward(iRs), 34)
     .add(new THREE.Vector3(0, -58, 0));
  labelGroup.add(sud);
  const icv = sprite('回盲瓣 · 与小肠页相接 →', '#8ce99a');
  icv.position.copy(vec(sc.axis[0])).addScaledVector(outward(0), 40)
     .add(new THREE.Vector3(-30, 62, 0));
  labelGroup.add(icv);
  const apxL = sprite('阑尾根部（构造）· 三带汇聚', '#fb7185');
  apxL.position.copy(vec(sc.appendix_base)).add(new THREE.Vector3(-20, -40, 0));
  labelGroup.add(apxL);

  // ---------------- animation -------------------------------------------------
  const clock = new THREE.Clock();
  const pulseArr = flow.sma_colic.pulsatile_ratio;
  let pulseI = 0;
  const frames = [];
  let lastT = performance.now();
  const gainLumen = flowRep.luminal_transport.display_gain;
  const vBase = flowRep.luminal_transport.baseline_content_velocity_mm_s;

  function tick() {
    requestAnimationFrame(tick);
    const t0 = performance.now();
    const dt = Math.min(clock.getDelta(), 0.05);
    pulseI = (pulseI + 1) % pulseArr.length;
    const pulse = pulseArr[pulseI];
    if (state.view === 'organ') particles.step(dt, state.rate, pulse);
    if (state.view === 'motility') {
      motility.step(dt * state.rate);
      const L = sc.arclen[sc.arclen.length - 1];
      const v = (state.mode === 'mass') ? 0 : vBase * gainLumen;
      for (let i = 0; i < NC; i++) {
        contentU[i] = (contentU[i] + v * dt * state.rate) % L;
        let k = 1; while (k < sc.arclen.length - 1 && sc.arclen[k] < contentU[i]) k++;
        const a = sc.axis[k - 1], b = sc.axis[k];
        const f = (contentU[i] - sc.arclen[k - 1]) /
                  Math.max(1e-6, sc.arclen[k] - sc.arclen[k - 1]);
        const rr = sc.radius[k] * tube2.radial[k] * 0.45;
        const ph = i * 2.399963;
        const N = sc.normal[k], B = sc.binormal[k];
        const ex = Math.cos(ph) * N[0] + Math.sin(ph) * B[0];
        const ey = Math.cos(ph) * N[1] + Math.sin(ph) * B[1];
        const ez = Math.cos(ph) * N[2] + Math.sin(ph) * B[2];
        cpos[i * 3] = a[0] + (b[0] - a[0]) * f + ex * rr;
        cpos[i * 3 + 1] = a[1] + (b[1] - a[1]) * f + ey * rr;
        cpos[i * 3 + 2] = a[2] + (b[2] - a[2]) * f + ez * rr;
      }
      contentGeom.attributes.position.needsUpdate = true;
    }
    if (state.view === 'crypt') gCrypt.rotation.y += dt * 0.12 * state.rate;
    if (state.view === 'mol') gMol.rotation.y += dt * 0.18 * state.rate;
    controls.update();
    renderer.render(scene, camera);
    const t1 = performance.now();
    frames.push({ cpu: t1 - t0, gap: t0 - lastT });
    if (frames.length > 600) frames.shift();
    lastT = t0;
  }

  addEventListener('resize', () => {
    renderer.setSize(app.clientWidth, app.clientHeight);
    camera.aspect = app.clientWidth / app.clientHeight;
    camera.updateProjectionMatrix();
  });

  frameView('organ');
  tick();

  // ---------------- HUD + report tables ----------------------------------------
  const { renderTables } = await import('./hud.js');
  renderTables({ sc, micro, geomRep, vesRep, flowRep, taeRep, vessels,
                 expr: exprData, molInfo, MICROBIAL_LINES });

  // ---------------- audit hooks for the headless benchmark ---------------------
  window.__motionAudit = (mode) => {
    if (mode && mode !== motility.mode) motility.setMode(mode);
    return motility.audit();
  };
  window.__flowAudit = () => Object.assign({}, particles.audit,
    { paths_rendered: particles.paths.length,
      reflux_audit_python: flowRep.reflux_audit });
  window.__renderAudit = () => {
    const gl = renderer.getContext();
    const dbg = gl.getExtension('WEBGL_debug_renderer_info');
    const cpu = frames.map(f => f.cpu).sort((a, b) => a - b);
    const gaps = frames.map(f => f.gap).sort((a, b) => a - b);
    const q = (a, p) => a.length ? +a[Math.min(a.length - 1, Math.floor(a.length * p))].toFixed(3) : null;
    return {
      unmasked_renderer: dbg ? gl.getParameter(dbg.UNMASKED_RENDERER_WEBGL) : null,
      unmasked_vendor: dbg ? gl.getParameter(dbg.UNMASKED_VENDOR_WEBGL) : null,
      frames: frames.length,
      render_cost_ms_median: q(cpu, 0.5), render_cost_ms_p95: q(cpu, 0.95),
      frame_gap_ms_median: q(gaps, 0.5), frame_gap_ms_p95: q(gaps, 0.95),
      budget_60fps_ms: 16.67,
      programs: renderer.info.programs ? renderer.info.programs.length : null,
      draw_calls: renderer.info.render.calls, triangles: renderer.info.render.triangles,
      view: state.view, mode: state.mode
    };
  };
  window.__setView = v => { frameView(v); return state.view; };
  // Audit hook: drives the model clock faster than real time so a 3 cycles/min rhythm and
  // a 7.1 mm/s wave can both be measured inside a few seconds of wall time.  The motion
  // audit measures against MODEL time, so this multiplier cancels out of every frequency
  // and velocity it reports.  The user-facing slider is unaffected and still caps at 3x.
  window.__setRate = r => { state.rate = +r; return state.rate; };
  window.__colonReady = true;
})().catch(fail);
