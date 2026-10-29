// colon.js — derived haustral tube, the three taeniae, semilunar folds, appendices
// epiploicae.  The sacculation is generated here from the SAME relation the Python side
// derived it from (rho = taenia/wall length ratio, lambda = haustral pitch), so the page
// is not showing a baked mesh with hand-made bumps.
import * as THREE from 'three';

export const NCIRC = 72;

export function vec(a) { return new THREE.Vector3(a[0], a[1], a[2]); }

// circular-arc sagitta profile: chord lambda, sagitta h, hence arc/chord = 1/rho exactly
export function bulgeProfile(sv, s0, s1, lambda, h) {
  if (sv < s0 || sv > s1) return 0;
  const Rc = (lambda * lambda / 4 + h * h) / (2 * h);
  let u = (sv - s0) % lambda;
  const x = u - lambda / 2;
  const y = Math.sqrt(Math.max(Rc * Rc - x * x, 0)) - (Rc - h);
  return Math.max(y, 0);
}

// circumferential weight: 0 on each of the three bands, 1 midway between them
export function bandWeight(psi, thFree) {
  const d = psi - thFree;
  const w = Math.sin(1.5 * d);
  return w * w;
}

export class ColonTube {
  /** sc: colon_scene.json */
  constructor(sc, opts = {}) {
    this.sc = sc;
    this.n = sc.axis.length;
    this.lambda = sc.derivation.haustral_pitch_mm;
    this.h = sc.derivation.derived_sagitta_mm;
    this.sFold0 = sc.s_folds[0];
    this.sFold1 = sc.s_folds[sc.s_folds.length - 1];
    this.ncirc = opts.ncirc || NCIRC;
    this.geom = new THREE.BufferGeometry();
    this.pos = new Float32Array(this.n * this.ncirc * 3);
    this.nor = new Float32Array(this.n * this.ncirc * 3);
    const idx = [];
    for (let i = 0; i < this.n - 1; i++) {
      for (let k = 0; k < this.ncirc; k++) {
        const k2 = (k + 1) % this.ncirc;
        const a = i * this.ncirc + k, b = i * this.ncirc + k2;
        const c = (i + 1) * this.ncirc + k, d = (i + 1) * this.ncirc + k2;
        idx.push(a, c, d, a, d, b);
      }
    }
    this.geom.setAttribute('position', new THREE.BufferAttribute(this.pos, 3));
    this.geom.setAttribute('normal', new THREE.BufferAttribute(this.nor, 3));
    this.geom.setIndex(idx);
    this.radial = new Float32Array(this.n);     // per-ring radial scale (motility writes here)
    this.shift = new Float32Array(this.n);      // per-ring axial shift (always 0 here)
    this.radial.fill(1);
    this.build(0);
  }

  /** phase: longitudinal phase shift of the haustral pattern (mm) */
  build(phase = 0, haustra = true) {
    const sc = this.sc, n = this.n, nc = this.ncirc;
    const P = sc.axis, N = sc.normal, B = sc.binormal, R = sc.radius;
    const s = sc.arclen, env = sc.envelope, thF = sc.theta_free;
    for (let i = 0; i < n; i++) {
      const bl = haustra
        ? bulgeProfile(s[i] + phase, this.sFold0, this.sFold1, this.lambda, this.h) * env[i]
        : 0;
      const px = P[i][0], py = P[i][1], pz = P[i][2];
      const nx = N[i][0], ny = N[i][1], nz = N[i][2];
      const bx = B[i][0], by = B[i][1], bz = B[i][2];
      const rr = R[i] * this.radial[i];
      for (let k = 0; k < nc; k++) {
        const psi = 2 * Math.PI * k / nc;
        const r = rr + bl * bandWeight(psi, thF[i]) * this.radial[i];
        const ca = Math.cos(psi), sa = Math.sin(psi);
        const ex = ca * nx + sa * bx, ey = ca * ny + sa * by, ez = ca * nz + sa * bz;
        const o = (i * nc + k) * 3;
        this.pos[o] = px + ex * r; this.pos[o + 1] = py + ey * r; this.pos[o + 2] = pz + ez * r;
        this.nor[o] = ex; this.nor[o + 1] = ey; this.nor[o + 2] = ez;
      }
    }
    this.geom.attributes.position.needsUpdate = true;
    this.geom.attributes.normal.needsUpdate = true;
    this.geom.computeBoundingSphere();
  }

  /** world position of a point on the derived surface — used by the motion audit */
  ringPoint(i, k) {
    const o = (i * this.ncirc + k) * 3;
    return new THREE.Vector3(this.pos[o], this.pos[o + 1], this.pos[o + 2]);
  }
  ringCentre(i) {
    const P = this.sc.axis[i];
    return new THREE.Vector3(P[0], P[1], P[2]);
  }
  ringRadius(i) {
    // mean radius of the ring as actually written into the buffer
    let sum = 0;
    const c = this.ringCentre(i);
    for (let k = 0; k < this.ncirc; k++) sum += this.ringPoint(i, k).distanceTo(c);
    return sum / this.ncirc;
  }
}

export function taeniaCurves(sc, which) {
  const key = { mesocolica: 'theta_mesocolic', libera: 'theta_free', omentalis: 'theta_omental' }[which];
  const th = sc[key];
  const pts = [];
  for (let i = 0; i < sc.axis.length; i++) {
    const P = sc.axis[i], N = sc.normal[i], B = sc.binormal[i], r = sc.radius[i] * 1.012;
    const c = Math.cos(th[i]), s = Math.sin(th[i]);
    pts.push(new THREE.Vector3(P[0] + (c * N[0] + s * B[0]) * r,
                               P[1] + (c * N[1] + s * B[1]) * r,
                               P[2] + (c * N[2] + s * B[2]) * r));
  }
  return pts;
}

export function foldRings(sc) {
  // one semilunar fold per declared fold position, drawn as a partial ring in each of the
  // three inter-taenial sectors (the fold does not cross a band)
  const out = [];
  for (const sf of sc.s_folds) {
    let i = 0;
    while (i < sc.arclen.length - 1 && sc.arclen[i] < sf) i++;
    if (sc.envelope[i] <= 0.02) continue;
    const P = sc.axis[i], N = sc.normal[i], B = sc.binormal[i], r = sc.radius[i] * 1.004;
    for (let sector = 0; sector < 3; sector++) {
      const a0 = sc.theta_free[i] + sector * 2 * Math.PI / 3 + 0.18;
      const a1 = a0 + 2 * Math.PI / 3 - 0.36;
      const pts = [];
      for (let t = 0; t <= 14; t++) {
        const a = a0 + (a1 - a0) * t / 14;
        const c = Math.cos(a), sn = Math.sin(a);
        pts.push(new THREE.Vector3(P[0] + (c * N[0] + sn * B[0]) * r,
                                   P[1] + (c * N[1] + sn * B[1]) * r,
                                   P[2] + (c * N[2] + sn * B[2]) * r));
      }
      out.push(pts);
    }
  }
  return out;
}

export function appendicesMesh(sc, mat) {
  const g = new THREE.SphereGeometry(1, 8, 6);
  const m = new THREE.InstancedMesh(g, mat, sc.appendices.length);
  const dummy = new THREE.Object3D();
  sc.appendices.forEach((a, i) => {
    const b = vec(a.base), d = vec(a.dir).normalize();
    dummy.position.copy(b.clone().addScaledVector(d, a.length_mm * 0.5));
    dummy.scale.set(a.width_mm * 0.5, a.length_mm * 0.5, a.width_mm * 0.45);
    dummy.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), d);
    dummy.updateMatrix();
    m.setMatrixAt(i, dummy.matrix);
  });
  m.instanceMatrix.needsUpdate = true;
  return m;
}
