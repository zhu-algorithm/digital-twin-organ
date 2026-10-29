// flow.js — particles advected along the baked flow paths at the velocities
// computed in build_stomach_flow.py.  Speed is never cosmetic: a particle's
// along-path speed is v(s) = Q(s)/A(s) in mm/s, interpolated from the path.
// There is NO display gain in the vascular layer.
import * as THREE from 'three';
import { particleMaterial } from './materials.js';
import { TREE_COLOUR } from './vessels.js';

export class Flow {
  constructor(flow, perPath = 220) {
    this.paths = flow.paths;
    this.group = new THREE.Group();
    this.group.name = 'flow';
    this.byTree = {};
    this.rate = 1.0;
    this.density = 1.0;

    const byTree = {};
    for (const p of this.paths) (byTree[p.tree] ||= []).push(p);
    for (const [tree, ps] of Object.entries(byTree)) {
      const total = ps.length * perPath;
      const P = new Float32Array(total * 3);
      const geo = new THREE.BufferGeometry();
      geo.setAttribute('position', new THREE.BufferAttribute(P, 3));
      const mat = particleMaterial(TREE_COLOUR[tree] ?? 0xffffff, 0.95);
      const pts = new THREE.Points(geo, mat);
      pts.frustumCulled = false;
      pts.name = tree;
      this.group.add(pts);
      // cumulative arc length per path, so u in [0,1] maps to real mm
      const meta = ps.map(p => {
        const n = p.points.length;
        const cum = new Float32Array(n);
        for (let i = 1; i < n; i++) {
          const a = p.points[i - 1], b = p.points[i];
          cum[i] = cum[i - 1] + Math.hypot(b[0] - a[0], b[1] - a[1], b[2] - a[2]);
        }
        return { p, cum, L: cum[n - 1] };
      });
      const seed = new Float32Array(total);
      for (let i = 0; i < total; i++) seed[i] = Math.random();
      this.byTree[tree] = { pts, geo, P, ps, meta, perPath, seed, s: new Float32Array(total) };
      for (let i = 0; i < total; i++) {
        this.byTree[tree].s[i] = seed[i] * meta[i % ps.length].L;
      }
    }
  }

  setVisible(tree, on) { const t = this.byTree[tree]; if (t) t.pts.visible = on; }
  setRate(r) { this.rate = r; }
  setDensity(f) {
    this.density = f;
    for (const t of Object.values(this.byTree)) {
      t.geo.setDrawRange(0, Math.max(1, Math.floor(t.P.length / 3 * f)));
    }
  }

  step(dt) {
    const d = dt * this.rate;
    for (const t of Object.values(this.byTree)) {
      const { P, meta, ps, perPath } = t;
      const n = P.length / 3;
      for (let i = 0; i < n; i++) {
        const mi = i % ps.length;
        const { p, cum, L } = meta[mi];
        let s = t.s[i] + d * sampleV(p, cum, t.s[i]);
        if (s >= L) s -= L;
        t.s[i] = s;
        // locate and interpolate
        let k = 1;
        while (k < cum.length - 1 && cum[k] < s) k++;
        const c0 = cum[k - 1], c1 = cum[k];
        const u = c1 > c0 ? (s - c0) / (c1 - c0) : 0;
        const a = p.points[k - 1], b = p.points[k];
        P[3 * i]     = a[0] + (b[0] - a[0]) * u;
        P[3 * i + 1] = a[1] + (b[1] - a[1]) * u;
        P[3 * i + 2] = a[2] + (b[2] - a[2]) * u;
      }
      t.geo.attributes.position.needsUpdate = true;
    }
  }
}

function sampleV(p, cum, s) {
  let k = 1;
  while (k < cum.length - 1 && cum[k] < s) k++;
  const c0 = cum[k - 1], c1 = cum[k];
  const u = c1 > c0 ? (s - c0) / (c1 - c0) : 0;
  return p.v_mm_s[k - 1] + (p.v_mm_s[k] - p.v_mm_s[k - 1]) * u;
}
