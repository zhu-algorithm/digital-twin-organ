/**
 * peristalsis.js — the parametric lumen tube and the travelling contraction.
 *
 * The tube is swept along the centreline in si_centreline.json. The duodenal
 * part of that centreline is skeletonised from the real `duodenum` label; the
 * jejuno-ileal part is a serpentine path through the measured `small_bowel`
 * envelope and is NOT contour-anchored — neither is any calibre. Both facts are
 * carried in the asset and shown in the page.
 *
 * The wave is a real travelling constriction, not a global scale. Every ring is
 * displaced in the vertex shader from its own arc length s:
 *
 *   r(s,t) = r0(s)
 *            * (1 - depth * exp(-((frac(s - v t, L)) / w)^2))     peristalsis
 *            * (1 - segAmp * (1 + cos(2 pi f(s) t)) / 2)          segmentation
 *            * (1 + foldAmp(s) * cos(2 pi s / foldPitch))         plicae circulares
 *
 * v = 30 mm/s (published 20-40 mm/s), f(s) falls linearly from 11.5/min in the
 * duodenum to 8.5/min in the ileum as the interstitial cells of Cajal dictate,
 * and foldAmp(s) tapers to zero in the distal ileum where plicae are absent.
 * Doing it per-ring in the shader is what keeps this at 60 fps on a phone: the
 * CPU uploads one float (time) per frame, not 14 000 positions.
 */

import * as THREE from 'three';

export const WAVE = {
  speed_mm_s: 30.0,              // 2-4 cm/s, Ganong's 26th ed. ch.25
  half_width_mm: 35.0,
  depth_fraction: 0.45,          // declared visualisation value, not measured
  segmentation_amp: 0.12,
  fold_pitch_mm: 7.0,            // plicae circulares 3-10 mm apart
  fold_amp_jejunum: 0.17,
  f_duodenum_per_min: 11.5,
  f_ileum_per_min: 8.5
};

function ringSides (r) { return r > 9 ? 20 : (r > 6 ? 16 : 12); }

/**
 * @param {object} cl parsed si_centreline.json
 * @param {number} stride take every n-th centreline sample (2 mm each)
 */
export function buildLumen (cl, opts = {}) {
  const stride = opts.stride || 2;
  const rScale = opts.rScale !== undefined ? opts.rScale : 1.0;
  const P = cl.points_scene_mm, S = cl.arc_mm, D = cl.outer_diameter_mm_parametric;
  const ST = cl.station;
  const L = S[S.length - 1];

  const idxs = [];
  for (let i = 0; i < P.length; i += stride) idxs.push(i);
  if (idxs[idxs.length - 1] !== P.length - 1) idxs.push(P.length - 1);

  const pos = [], n1a = [], n2a = [], ang = [], sArr = [], r0a = [], folda = [],
        freqa = [], stat = [], idx = [];
  const a = new THREE.Vector3(), b = new THREE.Vector3(), dir = new THREE.Vector3(),
        up = new THREE.Vector3(), t1 = new THREE.Vector3(), t2 = new THREE.Vector3(),
        prevT1 = new THREE.Vector3(0, 1, 0);
  let base = 0;
  const sides = [];

  for (let k = 0; k < idxs.length; k++) {
    const i = idxs[k];
    const j = idxs[Math.min(k + 1, idxs.length - 1)];
    const h = idxs[Math.max(k - 1, 0)];
    a.fromArray(P[h]); b.fromArray(P[j]);
    dir.subVectors(b, a);
    if (dir.lengthSq() < 1e-9) dir.set(0, 1, 0);
    dir.normalize();
    // parallel transport the frame so the fold crests do not corkscrew
    t1.copy(prevT1).addScaledVector(dir, -prevT1.dot(dir));
    if (t1.lengthSq() < 1e-6) {
      up.set(0, 1, 0);
      if (Math.abs(dir.dot(up)) > 0.94) up.set(1, 0, 0);
      t1.crossVectors(up, dir);
    }
    t1.normalize();
    t2.crossVectors(dir, t1).normalize();
    prevT1.copy(t1);

    const r0 = D[i] * 0.5 * rScale;
    const n = ringSides(r0);
    sides.push(n);
    // plicae fade out over the distal third of the ileum
    const fAxial = S[i] / L;
    const fold = WAVE.fold_amp_jejunum * Math.max(0, Math.min(1, (0.86 - fAxial) / 0.5));
    const fr = (WAVE.f_duodenum_per_min +
      (WAVE.f_ileum_per_min - WAVE.f_duodenum_per_min) * fAxial) / 60.0;   // Hz
    const c = P[i];
    for (let m = 0; m < n; m++) {
      pos.push(c[0], c[1], c[2]);
      n1a.push(t1.x, t1.y, t1.z);
      n2a.push(t2.x, t2.y, t2.z);
      ang.push((m / n) * Math.PI * 2);
      sArr.push(S[i]); r0a.push(r0); folda.push(fold); freqa.push(fr);
      stat.push(ST[i]);
    }
  }
  // stitch consecutive rings (they can differ in side count at a calibre step;
  // bridge by nearest angle)
  let off = 0;
  for (let k = 0; k < sides.length - 1; k++) {
    const nA = sides[k], nB = sides[k + 1], oA = off, oB = off + nA;
    for (let m = 0; m < nA; m++) {
      const m2 = (m + 1) % nA;
      const q = Math.round(m * nB / nA) % nB, q2 = Math.round(m2 * nB / nA) % nB;
      idx.push(oA + m, oB + q, oB + q2);
      idx.push(oA + m, oB + q2, oA + m2);
    }
    off += nA;
  }

  const g = new THREE.BufferGeometry();
  g.setAttribute('position', new THREE.Float32BufferAttribute(pos, 3));
  g.setAttribute('_n1', new THREE.Float32BufferAttribute(n1a, 3));
  g.setAttribute('_n2', new THREE.Float32BufferAttribute(n2a, 3));
  g.setAttribute('_ang', new THREE.Float32BufferAttribute(ang, 1));
  g.setAttribute('_s', new THREE.Float32BufferAttribute(sArr, 1));
  g.setAttribute('_r0', new THREE.Float32BufferAttribute(r0a, 1));
  g.setAttribute('_fold', new THREE.Float32BufferAttribute(folda, 1));
  g.setAttribute('_freq', new THREE.Float32BufferAttribute(freqa, 1));
  g.setAttribute('_stat', new THREE.Float32BufferAttribute(stat, 1));
  g.setAttribute('normal', new THREE.Float32BufferAttribute(new Float32Array(pos.length), 3));
  g.setIndex(idx);
  g.computeBoundingSphere();

  return { geometry: g, rings: sides.length, triangles: idx.length / 3,
           vertices: pos.length / 3, arcLength_mm: L };
}

/** Patch a MeshStandardMaterial so the ring is built and deformed on the GPU. */
export function applyPeristalsis (mat, L) {
  mat.userData.u = {
    uTime: { value: 0 }, uWaveSpeed: { value: WAVE.speed_mm_s },
    uHalfW: { value: WAVE.half_width_mm }, uDepth: { value: WAVE.depth_fraction },
    uSegAmp: { value: WAVE.segmentation_amp }, uFoldPitch: { value: WAVE.fold_pitch_mm },
    uArc: { value: L }, uNWave: { value: 4.0 }, uStationHi: { value: -1.0 }
  };
  mat.onBeforeCompile = (sh) => {
    Object.assign(sh.uniforms, mat.userData.u);
    sh.vertexShader = sh.vertexShader
      .replace('#include <common>', `#include <common>
        attribute vec3 _n1; attribute vec3 _n2;
        attribute float _ang; attribute float _s; attribute float _r0;
        attribute float _fold; attribute float _freq; attribute float _stat;
        uniform float uTime, uWaveSpeed, uHalfW, uDepth, uSegAmp, uFoldPitch,
                      uArc, uNWave, uStationHi;
        varying float vStat; varying float vSquash;`)
      .replace('#include <begin_vertex>', `
        float peri = 0.0;
        for (int w = 0; w < 4; w++) {
          float phase = uArc * (float(w) / uNWave);
          float d = _s - mod(uTime * uWaveSpeed + phase, uArc);
          d = d - uArc * floor(d / uArc + 0.5);          // shortest wrap
          peri = max(peri, exp(-(d * d) / (uHalfW * uHalfW)));
        }
        float segm = uSegAmp * 0.5 * (1.0 + cos(6.2831853 * _freq * uTime));
        float fold = _fold * cos(6.2831853 * _s / uFoldPitch);
        float r = _r0 * (1.0 - uDepth * peri) * (1.0 - segm) * (1.0 + fold);
        vec3 radial = _n1 * cos(_ang) + _n2 * sin(_ang);
        vec3 transformed = position + radial * r;
        vStat = _stat; vSquash = uDepth * peri;`)
      .replace('#include <beginnormal_vertex>', `
        vec3 objectNormal = normalize(_n1 * cos(_ang) + _n2 * sin(_ang));`);
    sh.fragmentShader = sh.fragmentShader
      .replace('#include <common>', `#include <common>
        varying float vStat; varying float vSquash;
        uniform float uStationHi;`)
      .replace('#include <color_fragment>', `#include <color_fragment>
        // station tint: duodenum / jejunum / ileum
        vec3 tintD = vec3(1.00, 0.76, 0.62);
        vec3 tintJ = vec3(0.95, 0.60, 0.52);
        vec3 tintI = vec3(0.86, 0.52, 0.50);
        vec3 tnt = vStat < 0.5 ? tintD : (vStat < 1.5 ? tintJ : tintI);
        diffuseColor.rgb *= tnt;
        if (uStationHi >= 0.0 && abs(vStat - uStationHi) < 0.5)
          diffuseColor.rgb = mix(diffuseColor.rgb, vec3(1.0, 0.92, 0.55), 0.35);
        // the contracted ring blanches, as a real constricted bowel does
        diffuseColor.rgb = mix(diffuseColor.rgb, vec3(0.98, 0.86, 0.80), vSquash * 0.45);`);
    mat.userData.shader = sh;
  };
  mat.needsUpdate = true;
  return mat;
}
