/**
 * main.js — self-contained Three.js beating-heart scene.
 *
 * No backend. Assets consumed:
 *   assets/heart.glb                 anatomical surfaces with baked cardiac coords
 *   assets/deform_calibration.json   volume -> deformation calibration
 *   assets/coronary_tree.json        CCO coronary tree (Murray's law)
 *   assets/flow_paths.json           blood-flow centrelines + contract bindings
 *   assets/heart_timeseries_<state>.json   the cardiac-model-track contract
 *
 * Public API (also used by the automated frame-rate probe):
 *   const api = await createHeartScene({ container, assetBase, state })
 *   api.setState('rest'|'active'), api.setPlaybackRate(x), api.dispose()
 *   window.__heart3d  diagnostics for the headless benchmark
 */

import * as THREE from 'three';
import { GLTFLoader } from '../vendor/jsm/loaders/GLTFLoader.js';
import { OrbitControls } from '../vendor/jsm/controls/OrbitControls.js';
import { RoomEnvironment } from '../vendor/jsm/environments/RoomEnvironment.js';

import { loadTimeseries } from './contract.js';
import { DeformField } from './deform.js';
import {
  makeMyocardiumMaterial, makeBloodPoolMaterial, makeVesselWallMaterial
} from './materials.js';
import { buildCoronaryMesh, coronaryDiastolicFraction } from './coronary.js';
import { FlowSystem } from './flow.js';
import { HUD } from './hud.js';

const VESSEL_TINT = {
  aorta: 0xd9b7ae,
  pulmonary_artery: 0xb9c4d6,
  svc: 0xaab6c8,
  ivc: 0xaab6c8,
  pulmonary_veins: 0xc9b2b2
};

export async function createHeartScene (opts = {}) {
  const container = opts.container || document.body;
  const base = opts.assetBase || './assets/';
  const states = opts.states || ['rest', 'active'];
  let stateName = opts.state || 'rest';

  const stage = (s) => { window.__heart3dStage = s; };
  stage('renderer:begin');
  // ---------------------------------------------------------------- renderer
  const renderer = new THREE.WebGLRenderer({
    antialias: true,
    powerPreference: 'high-performance',
    alpha: false
  });
  renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
  renderer.setSize(container.clientWidth, container.clientHeight);
  renderer.toneMapping = THREE.ACESFilmicToneMapping;
  renderer.toneMappingExposure = 1.08;
  renderer.outputColorSpace = THREE.SRGBColorSpace;
  container.appendChild(renderer.domElement);

  stage('renderer:ok');
  const scene = new THREE.Scene();
  scene.background = new THREE.Color(0x070a10);

  stage('pmrem:begin');
  const pmrem = new THREE.PMREMGenerator(renderer);
  scene.environment = pmrem.fromScene(new RoomEnvironment(), 0.04).texture;
  scene.environmentIntensity = 0.55;

  stage('pmrem:ok');
  const camera = new THREE.PerspectiveCamera(
    38, container.clientWidth / container.clientHeight, 1, 4000);
  camera.position.set(190, 60, 185);

  const controls = new OrbitControls(camera, renderer.domElement);
  controls.enableDamping = true;
  controls.dampingFactor = 0.07;
  controls.minDistance = 70;
  controls.maxDistance = 900;

  // Three-point lighting tuned for a wet serous surface.
  const key = new THREE.DirectionalLight(0xfff2ea, 2.4);
  key.position.set(1, 1.1, 0.8);
  scene.add(key);
  const fill = new THREE.DirectionalLight(0x9fc4ff, 0.65);
  fill.position.set(-1.1, 0.2, -0.6);
  scene.add(fill);
  const rim = new THREE.DirectionalLight(0xffd0c0, 1.1);
  rim.position.set(-0.4, -0.8, -1.0);
  scene.add(rim);
  scene.add(new THREE.AmbientLight(0x30405a, 0.6));

  // ------------------------------------------------------------------ assets
  stage('fetch:begin');
  const [calib, coroTree, flowSpec] = await Promise.all([
    fetch(base + 'deform_calibration.json').then(r => r.json()),
    fetch(base + 'coronary_tree.json').then(r => r.json()),
    fetch(base + 'flow_paths.json').then(r => r.json())
  ]);

  stage('fetch:ok');
  const seriesByState = {};
  for (const s of states) {
    seriesByState[s] = await loadTimeseries(base + `heart_timeseries_${s}.json`);
  }
  let ts = seriesByState[stateName];

  stage('contract:ok');
  const field = new DeformField(calib);

  stage('glb:begin');
  const gltf = await new GLTFLoader().loadAsync(base + 'heart.glb');

  // Display orientation lives on the group; object space == heart space.
  const heartGroup = new THREE.Group();
  heartGroup.rotation.set(-Math.PI / 2, 0, 0);     // long axis up for display
  scene.add(heartGroup);

  const myoMat = makeMyocardiumMaterial(field);
  let triangles = 0;
  const parts = {};

  gltf.scene.traverse((o) => {
    if (!o.isMesh) return;
    const name = o.name;
    const g = o.geometry;
    if (!g.getAttribute('_field0')) {
      console.warn(`[heart3d] mesh "${name}" has no _field0 attribute; skipped`);
      return;
    }
    let mat;
    if (name.startsWith('myo_') || name.startsWith('chamber_')) {
      mat = myoMat;
    } else if (name.startsWith('vessel_')) {
      const tintKey = name.replace('vessel_', '');
      mat = makeVesselWallMaterial(field, VESSEL_TINT[tintKey] || 0xc0b0b0);
    } else if (name.startsWith('pool_')) {
      mat = makeBloodPoolMaterial(field,
        name.includes('_r') ? 0x3a62a8 : 0xa8243a, 0.14);
    } else {
      mat = myoMat;
    }
    const mesh = new THREE.Mesh(g, mat);
    mesh.name = name;
    mesh.frustumCulled = false;
    heartGroup.add(mesh);
    parts[name] = mesh;
    triangles += (g.index ? g.index.count : g.getAttribute('position').count) / 3;
  });

  stage('glb:ok');
  const coro = buildCoronaryMesh(coroTree, field);
  heartGroup.add(coro.mesh);
  triangles += coro.stats.triangles;

  stage('coronary:ok');
  const flow = new FlowSystem(flowSpec, field, opts.particlesPerPath || 900);
  flow.setPixelRatio(renderer.getPixelRatio());
  heartGroup.add(flow.points);

  // --------------------------------------------------------------------- HUD
  stage('flow:ok');
  const hudCanvas = document.createElement('canvas');
  hudCanvas.style.cssText =
    'position:absolute;inset:0;pointer-events:none;z-index:5';
  container.appendChild(hudCanvas);
  const hud = new HUD(hudCanvas);

  function resize () {
    const w = container.clientWidth; const h = container.clientHeight;
    renderer.setSize(w, h, false);
    camera.aspect = w / h;
    camera.updateProjectionMatrix();
    hud.resize(w, h, Math.min(window.devicePixelRatio, 2));
    flow.setPixelRatio(renderer.getPixelRatio());
  }
  window.addEventListener('resize', resize);
  resize();

  // Frame the heart from the actual mesh bounds.
  const box = new THREE.Box3().setFromObject(heartGroup);
  const centre = box.getCenter(new THREE.Vector3());
  const radius = box.getSize(new THREE.Vector3()).length() * 0.5;
  controls.target.copy(centre);
  camera.position.copy(centre).add(
    new THREE.Vector3(1.05, 0.45, 1.0).normalize().multiplyScalar(radius * 2.6));
  camera.near = radius * 0.02; camera.far = radius * 40;
  camera.updateProjectionMatrix();
  controls.update();

  // ------------------------------------------------------------------- clock
  stage('hud:ok');
  let simMs = 0;
  let playbackRate = 1.0;
  let running = true;
  let lastNow = performance.now();
  const frameTimes = [];
  let frameCount = 0;

  const diag = {
    triangles,
    particles: flow.nPath * flow.ppp,
    coronary: coro.stats,
    coronaryDiastolicFraction: coronaryDiastolicFraction(ts),
    partNames: Object.keys(parts),
    ready: false
  };

  function frame (now) {
    if (!running) return;
    const wallDt = Math.min(now - lastNow, 100);
    lastNow = now;
    frameTimes.push(wallDt);
    if (frameTimes.length > 600) frameTimes.shift();
    frameCount++;

    const dtSim = wallDt * playbackRate;
    simMs += dtSim;

    field.update(ts, simMs);
    flow.update(ts, simMs, dtSim);

    controls.update();
    renderer.render(scene, camera);

    hud.tickFps(now);
    hud.draw(ts, simMs, field, flow, diag);

    requestAnimationFrame(frame);
  }
  stage('loop:start');
  requestAnimationFrame(frame);
  diag.ready = true;
  stage('ready');

  const api = {
    renderer, scene, camera, controls, field, flow, hud, heartGroup, parts, diag,
    get timeseries () { return ts; },
    get simMs () { return simMs; },
    setState (s) {
      if (!seriesByState[s]) throw new Error(`no timeseries for state "${s}"`);
      stateName = s; ts = seriesByState[s];
      diag.coronaryDiastolicFraction = coronaryDiastolicFraction(ts);
    },
    stateName () { return stateName; },
    setPlaybackRate (x) { playbackRate = x; },
    setParticleDensity (f) { flow.setDensity(f); },
    frameStats () {
      const a = frameTimes.slice().sort((x, y) => x - y);
      const mean = a.reduce((s, v) => s + v, 0) / Math.max(a.length, 1);
      return {
        frames: frameCount,
        meanFrameMs: mean,
        medianFrameMs: a.length ? a[Math.floor(a.length / 2)] : NaN,
        p95FrameMs: a.length ? a[Math.floor(a.length * 0.95)] : NaN,
        fpsMean: 1000 / mean,
        fpsMedian: a.length ? 1000 / a[Math.floor(a.length / 2)] : NaN
      };
    },
    resetFrameStats () { frameTimes.length = 0; frameCount = 0; },
    dispose () {
      running = false;
      window.removeEventListener('resize', resize);
      renderer.dispose();
    }
  };

  window.__heart3d = api;
  return api;
}
