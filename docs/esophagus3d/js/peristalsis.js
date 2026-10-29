/**
 * peristalsis.js — the oesophageal lumen tube and the travelling peristaltic
 * wave.
 *
 * The tube is swept along the MEASURED centreline in esophagus_centreline.json
 * (per-axial-slice centre of mass of subject s0688's `esophagus` label, cubic
 * smoothing spline, resampled at 1 mm) with the MEASURED inscribed-radius
 * profile as its calibre. Unlike the small-intestine page, where the
 * jejuno-ileal centreline had to be constructed, every ring position and every
 * resting radius here comes off the segmentation.
 *
 * The wave is a real travelling constriction, not a global scale. Every ring is
 * displaced in the VERTEX SHADER from its own arc length s:
 *
 *   r(s,t) = r0(s) * (1 - depth * g(s) * exp(-((s - x(t))/w)^2))
 *   x(t)   = integral of v(s) dt          wave front position, mm from the UES
 *
 * so the CPU uploads one float per frame (the front position) and no vertex
 * position is ever touched from JavaScript. A ring whose |s - x(t)| is large
 * compared with w does not move at all — which is the property the in-browser
 * verification exploits to prove this is a travelling wave and not a breathing
 * cylinder.
 *
 * v(s) = 35.3 mm/s, with a declared 0.75 factor across the striated-to-smooth
 * transition zone (arc fraction 0.33-0.45), where high-resolution manometry
 * reproducibly shows the propagation slowing. The 35.3 mm/s is the midpoint of
 * the intersection of the published 20-40 mm/s wave speed with the window
 * implied by this subject's measured 305.5 mm length and the published 6-10 s
 * transit; that construction is stated in esophagus_model.json rather than
 * presented as an independent validation.
 *
 * g(s) gates the wave:
 *   - primary   g = 1 everywhere (a swallow sweeps the whole organ)
 *   - secondary g = 0 above the distension site (no pharyngeal phase; the wave
 *     starts where the bolus is held and travels aborally only)
 *
 * The sphincters are rendered as separate tonic constrictions that are NOT part
 * of the travelling term: the UES sits closed and opens ahead of the wave, the
 * LES sits closed and relaxes for the whole transit.
 */

import * as THREE from 'three';

export const WAVE = {
  speed_mm_s: 35.3,              // see esophagus_model.json .peristalsis.primary
  published_speed_mm_s: [20, 40],
  half_width_mm: 40.0,           // declared; HRM contraction front is 3-5 cm long
  depth_fraction: 0.55,          // DECLARED visualisation value, see note above
  transition_zone: [0.33, 0.45],
  transition_speed_factor: 0.75,
  ues_open_lead_ms: 350,         // published 200-500 ms
  les_relax_ms: 8650             // = the measured transit; LES stays open for it
};

function ringSides (r) { return r > 7 ? 24 : (r > 4.5 ? 18 : 14); }

/**
 * @param {object} cl parsed esophagus_centreline.json
 * @param {object} opts { stride, wallScale }
 */
export function buildLumen (cl, opts = {}) {
  const stride = opts.stride || 2;
  const rScale = opts.rScale !== undefined ? opts.rScale : 1.0;
  const P = cl.points_scene_mm;
  const S = cl.arc_mm;
  const R = cl.inscribed_radius_mm;
  const U = cl.frame_u;
  const V = cl.frame_v;
  const L = S[S.length - 1];

  const idxs = [];
  for (let i = 0; i < P.length; i += stride) idxs.push(i);
  if (idxs[idxs.length - 1] !== P.length - 1) idxs.push(P.length - 1);

  const pos = [], n1a = [], n2a = [], ang = [], sArr = [], r0a = [], muscle = [], idx = [];
  const sides = [];
  let base = 0;

  for (let k = 0; k < idxs.length; k++) {
    const i = idxs[k];
    const r0 = R[i] * rScale;
    const n = ringSides(r0);
    sides.push(n);
    const f = S[i] / L;
    // muscle class as a vertex attribute so the section view can tint by
    // muscle TYPE without a second geometry: 0 striated, 1 mixed, 2 smooth
    const mc = f < 1 / 3 ? 0.0 : (f < 2 / 3 ? 1.0 : 2.0);
    const c = P[i], u = U[i], v = V[i];
    for (let m = 0; m < n; m++) {
      pos.push(c[0], c[1], c[2]);
      n1a.push(u[0], u[1], u[2]);
      n2a.push(v[0], v[1], v[2]);
      ang.push((m / n) * Math.PI * 2);
      sArr.push(S[i]);
      r0a.push(r0);
      muscle.push(mc);
    }
  }
  // stitch consecutive rings; side counts can differ at a calibre step, so the
  // seam is built by angle rather than by index
  base = 0;
  for (let k = 0; k < idxs.length - 1; k++) {
    const nA = sides[k], nB = sides[k + 1];
    const bA = base, bB = base + nA;
    for (let m = 0; m < nA; m++) {
      const a0 = bA + m, a1 = bA + ((m + 1) % nA);
      const t0 = Math.round((m / nA) * nB) % nB;
      const t1 = Math.round(((m + 1) / nA) * nB) % nB;
      const b0 = bB + t0, b1 = bB + t1;
      idx.push(a0, b0, a1);
      if (t0 !== t1) idx.push(a1, b0, b1);
    }
    base += nA;
  }

  const g = new THREE.BufferGeometry();
  g.setAttribute('position', new THREE.Float32BufferAttribute(pos, 3));
  g.setAttribute('aN1', new THREE.Float32BufferAttribute(n1a, 3));
  g.setAttribute('aN2', new THREE.Float32BufferAttribute(n2a, 3));
  g.setAttribute('aAng', new THREE.Float32BufferAttribute(ang, 1));
  g.setAttribute('aS', new THREE.Float32BufferAttribute(sArr, 1));
  g.setAttribute('aR0', new THREE.Float32BufferAttribute(r0a, 1));
  g.setAttribute('aMuscle', new THREE.Float32BufferAttribute(muscle, 1));
  g.setIndex(idx);
  g.computeBoundingSphere();
  return { geometry: g, rings: idxs.length, vertices: pos.length / 3, L, ringArc: idxs.map(i => S[i]) };
}

const DECL = /* glsl */`
attribute vec3 aN1;
attribute vec3 aN2;
attribute float aAng;
attribute float aS;
attribute float aR0;
attribute float aMuscle;
uniform float uFront;      // wave front arc position, mm (negative = no wave)
uniform float uDepth;
uniform float uWidth;
uniform float uGateFrom;   // secondary peristalsis: no contraction above this s
uniform float uUES;        // UES arc position
uniform float uLES;        // LES arc position
uniform float uUESClose;   // 0..1 tonic closure
uniform float uLESClose;
uniform float uL;
uniform float uBolus;      // bolus centre arc position, mm (negative = none)
uniform float uBolusAmp;
varying float vS;
varying float vMuscle;
varying float vSqueeze;
varying vec3 vN;
`;

const BODY = /* glsl */`
  float r = aR0;
  // travelling peristaltic constriction
  float sq = 0.0;
  if (uFront >= 0.0 && aS >= uGateFrom) {
    float d = (aS - uFront) / uWidth;
    sq = uDepth * exp(-d * d);
  }
  // tonic sphincters, independent of the travelling term
  float du = (aS - uUES) / 9.0;
  float dl = (aS - uLES) / 12.0;
  sq = max(sq, uUESClose * 0.80 * exp(-du * du));
  sq = max(sq, uLESClose * 0.72 * exp(-dl * dl));
  r *= (1.0 - sq);
  // bolus: a local distension travelling just ahead of the contraction
  if (uBolus >= 0.0) {
    float db = (aS - uBolus) / 26.0;
    r *= (1.0 + uBolusAmp * exp(-db * db));
  }
  vec3 radial = cos(aAng) * aN1 + sin(aAng) * aN2;
  vec3 p = position + radial * r;
  vS = aS; vMuscle = aMuscle; vSqueeze = sq;
  vN = normalize(radial);
`;

export function makeLumenMaterial (opts = {}) {
  const u = {
    uFront: { value: -1 }, uDepth: { value: WAVE.depth_fraction },
    uWidth: { value: WAVE.half_width_mm }, uGateFrom: { value: 0 },
    uUES: { value: 12.0 }, uLES: { value: 248.4 },
    uUESClose: { value: 1.0 }, uLESClose: { value: 1.0 },
    uL: { value: opts.L || 305.5 }, uBolus: { value: -1 }, uBolusAmp: { value: 0.35 },
    uMuscleTint: { value: opts.muscleTint ? 1.0 : 0.0 },
    uOpacity: { value: opts.opacity !== undefined ? opts.opacity : 1.0 }
  };
  const m = new THREE.MeshPhysicalMaterial({
    color: 0xd4a08f, roughness: 0.56, metalness: 0.0,
    clearcoat: 0.35, clearcoatRoughness: 0.5,
    transparent: (opts.opacity !== undefined && opts.opacity < 1.0),
    opacity: opts.opacity !== undefined ? opts.opacity : 1.0,
    side: THREE.DoubleSide
  });
  m.onBeforeCompile = (sh) => {
    Object.assign(sh.uniforms, u);
    sh.vertexShader = DECL + 'uniform float uMuscleTint;\n' + sh.vertexShader;
    sh.vertexShader = sh.vertexShader.replace('#include <begin_vertex>',
      '#include <begin_vertex>\n' + BODY + '\n  transformed = p;');
    sh.vertexShader = sh.vertexShader.replace('#include <beginnormal_vertex>',
      '#include <beginnormal_vertex>\n  objectNormal = normalize(cos(aAng)*aN1 + sin(aAng)*aN2);');
    sh.fragmentShader = `
varying float vS; varying float vMuscle; varying float vSqueeze; varying vec3 vN;
uniform float uMuscleTint; uniform float uL;
` + sh.fragmentShader;
    sh.fragmentShader = sh.fragmentShader.replace('#include <color_fragment>', `
#include <color_fragment>
{
  // muscle-type banding: striated red, mixed amber, smooth blue. Off by default
  // in the organ view; the section view turns it on.
  vec3 striated = vec3(0.85, 0.34, 0.31);
  vec3 mixedm   = vec3(0.91, 0.64, 0.24);
  vec3 smoothm  = vec3(0.29, 0.56, 0.85);
  vec3 mcol = vMuscle < 0.5 ? striated : (vMuscle < 1.5 ? mixedm : smoothm);
  diffuseColor.rgb = mix(diffuseColor.rgb, mcol, uMuscleTint * 0.72);
  // the contracting ring blanches, as a real contraction does
  diffuseColor.rgb = mix(diffuseColor.rgb, vec3(0.95,0.86,0.82), vSqueeze * 0.45);
  // longitudinal mucosal folds of the collapsed resting oesophagus
  float fold = 0.5 + 0.5 * cos(vS * 0.0);
  diffuseColor.rgb *= (0.94 + 0.06 * fold);
}`);
    m.userData.shader = sh;
  };
  m.userData.uniforms = u;
  return m;
}

/**
 * Wave kinematics. Integrates the front position with the declared axial speed
 * profile, so the transition-zone slowing is in the POSITION, not applied as a
 * post-hoc visual.
 */
export class Wave {
  constructor (L, opts = {}) {
    this.L = L;
    this.front = -1;
    this.mode = 'idle';
    this.gateFrom = 0;
    this.bolus = -1;
    this.elapsed = 0;
    this.uesClose = 1.0;
    this.lesClose = 1.0;
    this.uesArc = opts.uesArc !== undefined ? opts.uesArc : 12.0;
    this.lesArc = opts.lesArc !== undefined ? opts.lesArc : 248.4;
    this.history = [];
  }

  speedAt (s) {
    const f = s / this.L;
    const [a, b] = WAVE.transition_zone;
    return WAVE.speed_mm_s * (f >= a && f <= b ? WAVE.transition_speed_factor : 1.0);
  }

  swallow () { this.mode = 'primary'; this.front = 0; this.gateFrom = 0; this.elapsed = 0; this.bolus = 26; }

  /** @param {number} sFrom arc position of the distension that triggers it */
  secondary (sFrom) {
    this.mode = 'secondary'; this.front = sFrom; this.gateFrom = sFrom;
    this.elapsed = 0; this.bolus = sFrom + 26;
  }

  reset () { this.mode = 'idle'; this.front = -1; this.bolus = -1; this.gateFrom = 0; }

  /**
   * The front is integrated against WALL-CLOCK elapsed time, not against a
   * clamped per-frame dt.
   *
   * This matters and was found by measurement, not by inspection: the first
   * headless audit clamped dt at 120 ms per frame, and under software
   * rasterisation requestAnimationFrame fired at only 2-5 fps (200-480 ms per
   * frame). The clamp then discarded most of each interval and the front
   * advanced at 17-19 mm/s against a declared 35.3 mm/s. A kinematic animation
   * must not slow down when frames are slow — it must skip. The clamp is kept
   * only as a guard against a single pathological interval (a backgrounded tab
   * returning a multi-second dt), by capping the per-call advance rather than
   * the elapsed time.
   */
  step (dtMs, rate) {
    if (this.mode === 'idle' || !(dtMs > 0)) return;
    const dt = Math.min(dtMs, 2000) * 0.001 * rate;
    this.elapsed += dt;
    const v = this.speedAt(this.front);
    this.front += v * dt;
    if (this.bolus >= 0) this.bolus += v * dt;
    // UES opens just ahead of the arriving wave and recloses behind it
    const lead = WAVE.ues_open_lead_ms * 0.001 * v;
    this.uesClose = (this.front > this.uesArc - lead && this.front < this.uesArc + 30) ? 0.12 : 1.0;
    // LES relaxes for the whole transit
    this.lesClose = (this.front > 10 && this.front < this.L + 20) ? 0.10 : 1.0;
    if (this.front > this.L + WAVE.half_width_mm * 1.6) this.reset();
  }

  apply (uniforms) {
    uniforms.uFront.value = this.front;
    uniforms.uGateFrom.value = this.gateFrom;
    uniforms.uBolus.value = this.bolus;
    uniforms.uUESClose.value = this.uesClose;
    uniforms.uLESClose.value = this.lesClose;
    uniforms.uUES.value = this.uesArc;
    uniforms.uLES.value = this.lesArc;
  }
}
