// vessels.js — the five gastric arterial trees and the gastric venous tree,
// drawn as tapered tubes at their Murray radii.
import * as THREE from 'three';
import { vesselMaterial } from './materials.js';

export const TREE_COLOUR = {
  left_gastric:         0xd84a3a,
  right_gastric:        0xe0714a,
  left_gastroepiploic:  0xc2403f,
  right_gastroepiploic: 0xe8884e,
  short_gastric:        0xb8464f,
  gastric_veins:        0x3f73c8,
};
export const TREE_LABEL = {
  left_gastric: '胃左动脉 left gastric',
  right_gastric: '胃右动脉 right gastric',
  left_gastroepiploic: '胃网膜左动脉 left gastroepiploic',
  right_gastroepiploic: '胃网膜右动脉 right gastroepiploic',
  short_gastric: '胃短动脉 short gastric',
  gastric_veins: '胃静脉→门静脉 gastric veins → portal',
};

export function buildVessels(vessels, radiusGain = 1.0) {
  const group = new THREE.Group();
  group.name = 'vessels';
  const perTree = {};
  for (const tree of vessels.trees) {
    const g = new THREE.Group();
    g.name = tree.name;
    const pos = tree.nodes, par = tree.parent, rad = tree.radius_mm;
    const segs = [];
    for (let i = 0; i < pos.length; i++) {
      const p = par[i];
      if (p < 0) continue;
      const a = new THREE.Vector3(...pos[p]);
      const b = new THREE.Vector3(...pos[i]);
      const L = a.distanceTo(b);
      if (L < 1e-4) continue;
      segs.push({ a, b, L, ra: rad[p] * radiusGain, rb: rad[i] * radiusGain });
    }
    // one merged mesh per tree keeps the draw-call count down on phones
    const geoms = [];
    for (const s of segs) {
      const cyl = new THREE.CylinderGeometry(s.rb, s.ra, s.L, 6, 1, true);
      const mid = s.a.clone().add(s.b).multiplyScalar(0.5);
      const dir = s.b.clone().sub(s.a).normalize();
      const q = new THREE.Quaternion().setFromUnitVectors(new THREE.Vector3(0, 1, 0), dir);
      cyl.applyQuaternion(q); cyl.translate(mid.x, mid.y, mid.z);
      geoms.push(cyl);
    }
    const merged = mergeSimple(geoms);
    const mesh = new THREE.Mesh(merged, vesselMaterial(TREE_COLOUR[tree.name],
      tree.retrograde ? 0.10 : 0.20));
    mesh.name = tree.name;
    g.add(mesh);
    group.add(g);
    perTree[tree.name] = { group: g, segments: segs.length, flow: tree.flow_mL_min,
                           retrograde: tree.retrograde };
    for (const q of geoms) q.dispose();
  }
  return { group, perTree };
}

function mergeSimple(geoms) {
  let nv = 0, ni = 0;
  for (const g of geoms) { nv += g.attributes.position.count; ni += g.index.count; }
  const P = new Float32Array(nv * 3), N = new Float32Array(nv * 3);
  const I = nv < 65536 ? new Uint16Array(ni) : new Uint32Array(ni);
  let vo = 0, io = 0;
  for (const g of geoms) {
    const p = g.attributes.position.array, n = g.attributes.normal.array;
    P.set(p, vo * 3); N.set(n, vo * 3);
    const idx = g.index.array;
    for (let k = 0; k < idx.length; k++) I[io + k] = idx[k] + vo;
    vo += g.attributes.position.count; io += idx.length;
  }
  const out = new THREE.BufferGeometry();
  out.setAttribute('position', new THREE.BufferAttribute(P, 3));
  out.setAttribute('normal', new THREE.BufferAttribute(N, 3));
  out.setIndex(new THREE.BufferAttribute(I, 1));
  out.computeBoundingSphere();
  return out;
}
