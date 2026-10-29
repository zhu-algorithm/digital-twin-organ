// vessels.js — draw the colonic trees and advance particles at the speed baked into
// colon_flow.json.  Particle speed is v(s) = Q(s)/A(s) from the Python side; nothing is
// re-derived here and no gain is applied to blood flow.
import * as THREE from 'three';

export const TREE_STYLE = {
  sma_colic:            { colour: 0xff5d55, kind: 'artery', label: '肠系膜上动脉系（近端结肠）' },
  ima_colic:            { colour: 0xff8a3d, kind: 'artery', label: '肠系膜下动脉系（远端结肠）' },
  rectal_systemic:      { colour: 0xffd166, kind: 'artery', label: '直肠中/下动脉（髂内，体循环回流）' },
  smv_portal:           { colour: 0x4f8cff, kind: 'vein',   label: '肠系膜上静脉 → 门静脉' },
  imv_portal:           { colour: 0x59b7ff, kind: 'vein',   label: '肠系膜下静脉 → 门静脉' },
  rectal_vein_systemic: { colour: 0x8bd3ff, kind: 'vein',   label: '直肠中/下静脉 → 髂内静脉 → 下腔静脉' }
};

export function buildTree(name, data, flow, opts = {}) {
  const group = new THREE.Group();
  group.name = name;
  const st = TREE_STYLE[name];
  const { parent, nodes, radius_mm } = data;
  const segs = [];
  for (let i = 0; i < parent.length; i++) if (parent[i] >= 0) segs.push(i);

  // one merged line set, width-coded by radius via three line groups
  const buckets = [[], [], []];
  const rmax = Math.max(...radius_mm);
  for (const i of segs) {
    const p = nodes[parent[i]], q = nodes[i];
    const b = radius_mm[i] > rmax * 0.33 ? 0 : radius_mm[i] > rmax * 0.1 ? 1 : 2;
    buckets[b].push(p[0], p[1], p[2], q[0], q[1], q[2]);
  }
  const widths = [2.6, 1.5, 0.8];
  buckets.forEach((arr, b) => {
    if (!arr.length) return;
    const g = new THREE.BufferGeometry();
    g.setAttribute('position', new THREE.Float32BufferAttribute(arr, 3));
    const m = new THREE.LineBasicMaterial({
      color: st.colour, transparent: true,
      opacity: b === 0 ? 0.96 : b === 1 ? 0.78 : 0.5, linewidth: widths[b]
    });
    group.add(new THREE.LineSegments(g, m));
  });

  // trunk drawn as a tube so the big vessels read as vessels, not as hairlines
  const trunkIdx = segs.filter(i => radius_mm[i] > rmax * 0.55);
  for (const i of trunkIdx) {
    const p = new THREE.Vector3(...nodes[parent[i]]), q = new THREE.Vector3(...nodes[i]);
    const len = p.distanceTo(q);
    if (len < 1e-3) continue;
    const cyl = new THREE.Mesh(
      new THREE.CylinderGeometry(radius_mm[i], radius_mm[i], len, 10, 1, true),
      new THREE.MeshStandardMaterial({ color: st.colour, roughness: 0.55, metalness: 0.05,
                                       transparent: true, opacity: 0.92 }));
    cyl.position.copy(p.clone().add(q).multiplyScalar(0.5));
    cyl.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0),
                                      q.clone().sub(p).normalize());
    group.add(cyl);
  }
  return group;
}

export class FlowParticles {
  /** paths: array of {nodes:[idx...], v_mm_s:[...]} ; nodePos: node coordinate array */
  constructor(trees, flow, colourOf, maxParticles = 2600) {
    this.items = [];
    const all = [];
    for (const name of Object.keys(trees)) {
      const f = flow[name]; if (!f) continue;
      const d = trees[name];
      const kids = {};
      d.parent.forEach((p, i) => { if (p >= 0) (kids[p] = kids[p] || []).push(i); });
      const terms = [];
      for (let i = 0; i < d.terminal.length; i++) if (d.terminal[i]) terms.push(i);
      for (const t of terms) {
        const chain = [t]; let c = t;
        while (d.parent[c] >= 0) { c = d.parent[c]; chain.push(c); }
        chain.reverse();
        const pts = chain.map(i => d.nodes[i]);
        const vs = chain.slice(1).map(i => f.v_mm_s[i]);
        if (!vs.length) continue;
        all.push({ name, pts, vs, colour: colourOf(name),
                   reverse: TREE_STYLE[name].kind === 'vein' });
      }
    }
    // even sample of paths so every tree is represented
    const stride = Math.max(1, Math.ceil(all.length / maxParticles));
    this.paths = all.filter((_, i) => i % stride === 0);
    const n = this.paths.length;
    const geom = new THREE.BufferGeometry();
    this.pos = new Float32Array(n * 3);
    this.col = new Float32Array(n * 3);
    geom.setAttribute('position', new THREE.BufferAttribute(this.pos, 3));
    geom.setAttribute('color', new THREE.BufferAttribute(this.col, 3));
    const c = new THREE.Color();
    this.paths.forEach((p, i) => {
      c.setHex(p.colour);
      this.col[i * 3] = c.r; this.col[i * 3 + 1] = c.g; this.col[i * 3 + 2] = c.b;
      p.cum = [0];
      for (let k = 1; k < p.pts.length; k++) {
        const a = p.pts[k - 1], b = p.pts[k];
        p.cum.push(p.cum[k - 1] + Math.hypot(b[0] - a[0], b[1] - a[1], b[2] - a[2]));
      }
      p.len = p.cum[p.cum.length - 1];
      p.u = (i * 0.6180339887) % 1 * p.len;
    });
    this.points = new THREE.Points(geom, new THREE.PointsMaterial({
      size: 2.6, vertexColors: true, transparent: true, opacity: 0.95,
      sizeAttenuation: false, depthWrite: false }));
    this.audit = { samples: 0, nonpositive: 0, nonfinite: 0 };
  }

  step(dt, rate, pulse) {
    const P = this.paths;
    for (let i = 0; i < P.length; i++) {
      const p = P[i];
      // local speed at the current arc position
      let k = 1;
      while (k < p.cum.length - 1 && p.cum[k] < p.u) k++;
      const v = p.vs[Math.min(k - 1, p.vs.length - 1)] * pulse;
      this.audit.samples++;
      if (!isFinite(v)) this.audit.nonfinite++;
      else if (v <= 0) this.audit.nonpositive++;
      p.u += v * dt * rate;
      if (p.u > p.len) p.u -= p.len;
      // interpolate
      let j = 1;
      while (j < p.cum.length - 1 && p.cum[j] < p.u) j++;
      const t = (p.u - p.cum[j - 1]) / Math.max(1e-6, p.cum[j] - p.cum[j - 1]);
      let a = p.pts[j - 1], b = p.pts[j];
      if (p.reverse) { const tmp = a; a = b; b = tmp; }
      this.pos[i * 3] = a[0] + (b[0] - a[0]) * t;
      this.pos[i * 3 + 1] = a[1] + (b[1] - a[1]) * t;
      this.pos[i * 3 + 2] = a[2] + (b[2] - a[2]) * t;
    }
    this.points.geometry.attributes.position.needsUpdate = true;
  }
}

export function arcadeLine(pts, colour, width) {
  const g = new THREE.BufferGeometry();
  g.setAttribute('position', new THREE.Float32BufferAttribute(pts.flat(), 3));
  return new THREE.Line(g, new THREE.LineBasicMaterial({
    color: colour, linewidth: width || 2, transparent: true, opacity: 0.95 }));
}
