/**
 * vessels.js — tube meshes for the seven oesophageal vascular trees plus the
 * palisade-zone anastomotic network.
 *
 * Each tree arrives as (nodes, parent, radius_mm) from esophagus_vessels.json.
 * Those radii satisfy Murray's law exactly (r_parent^3 = sum r_child^3) because
 * they are generated as r_terminal * n_downstream^(1/3) over strictly binary
 * bifurcations. The mesh does not alter a single radius, so the calibre you see
 * is the calibre in vessel_report.json — there is no display rScale on the
 * vessel geometry on this page.
 *
 * The palisade zone is NOT a tree. It is an explicit longitudinal network of 16
 * channels in the distal 25 mm of submucosa, where the PORTAL territory (left
 * gastric vein) meets the SYSTEMIC territory (azygos). It is drawn from its own
 * edge list at a declared calibre and carries no Murray residual, because
 * Murray's law does not apply to a non-dichotomous plexus.
 */

import * as THREE from 'three';

export const TREE_COLOR = {
  A1_inferior_thyroid_artery: 0xf7664c,
  A2_bronchial_artery: 0xee4d42,
  A3_aortic_oesophageal_branches: 0xdc3838,
  A4_left_gastric_ascending: 0xff8547,
  V1_inferior_thyroid_vein_systemic: 0x5c8cd9,
  V2_azygos_systemic: 0x426bc7,
  V3_left_gastric_portal: 0x76c29a,
  PZ_palisade_normal: 0x8cd2b3,
  PZ_palisade_varix: 0xb84ec7
};

export const TREE_LABEL = {
  A1_inferior_thyroid_artery: '甲状腺下动脉 inferior thyroid a.',
  A2_bronchial_artery: '支气管动脉 bronchial a.',
  A3_aortic_oesophageal_branches: '主动脉食管支 aortic oesophageal brs',
  A4_left_gastric_ascending: '胃左动脉升支 left gastric ascending brs',
  V1_inferior_thyroid_vein_systemic: '甲状腺下静脉 inf. thyroid v. (体循环)',
  V2_azygos_systemic: '奇静脉 azygos v. (体循环)',
  V3_left_gastric_portal: '胃左静脉 left gastric v. (门脉)',
  PZ_palisade_normal: '栅栏区静脉 palisade zone',
  PZ_palisade_varix: '食管静脉曲张 varices'
};

function ringSides (r) {
  return r > 1.4 ? 14 : (r > 0.7 ? 10 : (r > 0.35 ? 8 : 5));
}

function frameFor (dir, prev) {
  const t1 = prev.clone().addScaledVector(dir, -prev.dot(dir));
  if (t1.lengthSq() < 1e-8) {
    const up = Math.abs(dir.y) > 0.94 ? new THREE.Vector3(1, 0, 0) : new THREE.Vector3(0, 1, 0);
    t1.crossVectors(up, dir);
  }
  t1.normalize();
  const t2 = new THREE.Vector3().crossVectors(dir, t1).normalize();
  return [t1, t2];
}

/**
 * Swept tubes for one tree. Radius is interpolated along each parent->child
 * segment so a parent tapers into its children instead of stepping.
 * @returns {{geometry: THREE.BufferGeometry, segments: number, triangles: number}}
 */
export function buildVesselMesh (tree, opts = {}) {
  const N = tree.nodes, PA = tree.parent, R = tree.radius_mm;
  const pos = [], nrm = [], idx = [], rad = [];
  let prev = new THREE.Vector3(0, 1, 0);
  let segs = 0;
  const a = new THREE.Vector3(), b = new THREE.Vector3(), d = new THREE.Vector3();

  for (let i = 0; i < N.length; i++) {
    const p = PA[i];
    if (p < 0) continue;
    a.fromArray(N[p]); b.fromArray(N[i]);
    d.subVectors(b, a);
    const len = d.length();
    if (len < 1e-5) continue;
    d.divideScalar(len);
    const [t1, t2] = frameFor(d, prev);
    prev = t1;
    const rA = R[p], rB = R[i];
    const n = ringSides(Math.max(rA, rB));
    const base = pos.length / 3;
    for (let e = 0; e < 2; e++) {
      const c = e === 0 ? a : b;
      const r = e === 0 ? rA : rB;
      for (let m = 0; m < n; m++) {
        const th = (m / n) * Math.PI * 2;
        const nx = Math.cos(th) * t1.x + Math.sin(th) * t2.x;
        const ny = Math.cos(th) * t1.y + Math.sin(th) * t2.y;
        const nz = Math.cos(th) * t1.z + Math.sin(th) * t2.z;
        pos.push(c.x + nx * r, c.y + ny * r, c.z + nz * r);
        nrm.push(nx, ny, nz);
        rad.push(r);
      }
    }
    for (let m = 0; m < n; m++) {
      const m1 = (m + 1) % n;
      idx.push(base + m, base + n + m, base + m1);
      idx.push(base + m1, base + n + m, base + n + m1);
    }
    segs++;
  }
  const g = new THREE.BufferGeometry();
  g.setAttribute('position', new THREE.Float32BufferAttribute(pos, 3));
  g.setAttribute('normal', new THREE.Float32BufferAttribute(nrm, 3));
  g.setAttribute('aRad', new THREE.Float32BufferAttribute(rad, 1));
  g.setIndex(idx);
  g.computeBoundingSphere();
  return { geometry: g, segments: segs, triangles: idx.length / 3 };
}

/** The palisade-zone plexus: an explicit edge list, not a tree. */
export function buildPalisadeMesh (pz, diameter_mm) {
  const r = diameter_mm * 0.5;
  const pos = [], nrm = [], idx = [];
  let prev = new THREE.Vector3(0, 1, 0);
  const a = new THREE.Vector3(), b = new THREE.Vector3(), d = new THREE.Vector3();
  const n = 6;
  for (const [i, j] of pz.edges) {
    a.fromArray(pz.nodes[i]); b.fromArray(pz.nodes[j]);
    d.subVectors(b, a);
    const len = d.length();
    if (len < 1e-5) continue;
    d.divideScalar(len);
    const [t1, t2] = frameFor(d, prev);
    prev = t1;
    const base = pos.length / 3;
    for (let e = 0; e < 2; e++) {
      const c = e === 0 ? a : b;
      for (let m = 0; m < n; m++) {
        const th = (m / n) * Math.PI * 2;
        const nx = Math.cos(th) * t1.x + Math.sin(th) * t2.x;
        const ny = Math.cos(th) * t1.y + Math.sin(th) * t2.y;
        const nz = Math.cos(th) * t1.z + Math.sin(th) * t2.z;
        pos.push(c.x + nx * r, c.y + ny * r, c.z + nz * r);
        nrm.push(nx, ny, nz);
      }
    }
    for (let m = 0; m < n; m++) {
      const m1 = (m + 1) % n;
      idx.push(base + m, base + n + m, base + m1);
      idx.push(base + m1, base + n + m, base + n + m1);
    }
  }
  const g = new THREE.BufferGeometry();
  g.setAttribute('position', new THREE.Float32BufferAttribute(pos, 3));
  g.setAttribute('normal', new THREE.Float32BufferAttribute(nrm, 3));
  g.setIndex(idx);
  g.computeBoundingSphere();
  return { geometry: g, triangles: idx.length / 3 };
}

export function makeVesselMaterial (colour, opts = {}) {
  return new THREE.MeshPhysicalMaterial({
    color: colour, roughness: 0.38, metalness: 0.0,
    clearcoat: 0.55, clearcoatRoughness: 0.35,
    emissive: colour, emissiveIntensity: opts.emissive !== undefined ? opts.emissive : 0.10,
    transparent: opts.opacity !== undefined && opts.opacity < 1,
    opacity: opts.opacity !== undefined ? opts.opacity : 1.0
  });
}
