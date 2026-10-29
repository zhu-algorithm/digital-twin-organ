/**
 * coronary.js — renders the epicardial coronary tree produced by
 * tools/build_coronary_tree.py (constrained constructive optimisation anchored
 * on the real LM / LAD / LCx / RCA contours of the source atlas).
 *
 * The tree arrives as a node/segment list with per-segment radii obeying
 * Murray's law. Here it is turned into a single merged tube mesh whose vertices
 * carry the same `_field0`/`_field1` cardiac-coordinate attributes as the
 * myocardium, so the vessels are deformed by exactly the same GPU kinematics
 * and stay welded to the epicardial surface through the cardiac cycle.
 */

import * as THREE from 'three';
import { makeCoronaryMaterial } from './materials.js';

const RADIAL_SEG = 7;

/**
 * @param {object} tree parsed coronary_tree.json
 * @param {DeformField} field
 * @returns {{mesh: THREE.Mesh, stats: object}}
 */
export function buildCoronaryMesh (tree, field) {
  const nodes = tree.nodes;           // [{p:[x,y,z] heart space, zeta, theta, chamber_id}]
  const segs = tree.segments;         // [{a, b, radius_mm, territory, name}]

  const positions = [];
  const normals = [];
  const field0 = [];
  const field1 = [];
  const indices = [];

  const up = new THREE.Vector3();
  const tangent = new THREE.Vector3();
  const e1 = new THREE.Vector3();
  const e2 = new THREE.Vector3();
  const A = new THREE.Vector3();
  const B = new THREE.Vector3();

  let base = 0;
  for (const s of segs) {
    const na = nodes[s.a]; const nb = nodes[s.b];
    A.fromArray(na.p); B.fromArray(nb.p);
    tangent.subVectors(B, A);
    const len = tangent.length();
    if (len < 1e-6) continue;
    tangent.divideScalar(len);
    up.set(0, 0, 1);
    if (Math.abs(tangent.z) > 0.9) up.set(1, 0, 0);
    e1.crossVectors(tangent, up).normalize();
    e2.crossVectors(tangent, e1).normalize();

    const rA = s.radius_mm;
    const rB = s.radius_tip_mm !== undefined ? s.radius_tip_mm : s.radius_mm;

    for (let end = 0; end < 2; end++) {
      const node = end === 0 ? na : nb;
      const P = end === 0 ? A : B;
      const r = end === 0 ? rA : rB;
      for (let j = 0; j < RADIAL_SEG; j++) {
        const a = (j / RADIAL_SEG) * Math.PI * 2;
        const nx = e1.x * Math.cos(a) + e2.x * Math.sin(a);
        const ny = e1.y * Math.cos(a) + e2.y * Math.sin(a);
        const nz = e1.z * Math.cos(a) + e2.z * Math.sin(a);
        positions.push(P.x + nx * r, P.y + ny * r, P.z + nz * r);
        normals.push(nx, ny, nz);
        // epicardial: transmural = 1, so the vessel tracks the outer surface
        field0.push(1.0, node.zeta, 1.0, node.theta);
        field1.push(node.chamber_id, 0.0, s.territory_id, 0.0);
      }
    }
    for (let j = 0; j < RADIAL_SEG; j++) {
      const j2 = (j + 1) % RADIAL_SEG;
      const a0 = base + j; const a1 = base + j2;
      const b0 = base + RADIAL_SEG + j; const b1 = base + RADIAL_SEG + j2;
      indices.push(a0, b0, b1, a0, b1, a1);
    }
    base += RADIAL_SEG * 2;
  }

  const g = new THREE.BufferGeometry();
  g.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
  g.setAttribute('normal', new THREE.Float32BufferAttribute(normals, 3));
  g.setAttribute('_field0', new THREE.Float32BufferAttribute(field0, 4));
  g.setAttribute('_field1', new THREE.Float32BufferAttribute(field1, 4));
  g.setIndex(indices);
  g.computeBoundingSphere();

  // The tree is authored in heart space; the scene graph holds the heart frame.
  // Object space IS heart space (the calibration ships identity frame matrices
  // and the anatomical display orientation lives on the parent group), so the
  // tube mesh needs no extra transform here.
  const mesh = new THREE.Mesh(g, makeCoronaryMaterial(field));
  mesh.name = 'coronary_tree';
  mesh.renderOrder = 5;

  return {
    mesh,
    stats: {
      segments: segs.length,
      nodes: nodes.length,
      triangles: indices.length / 3
    }
  };
}

/**
 * Coronary flow is phase-shifted relative to ejection. The contract supplies
 * `Q_coronary_mL_s` directly; this helper exposes the diastolic fraction so the
 * HUD can state it rather than assert it.
 */
export function coronaryDiastolicFraction (ts) {
  const q = ts.sig.Q_coronary_mL_s;
  const ao = ts.valve.aortic;
  let sysQ = 0; let diaQ = 0;
  for (let i = 0; i < q.length; i++) {
    const v = Math.max(q[i], 0);
    if (ao[i]) sysQ += v; else diaQ += v;
  }
  return diaQ / Math.max(sysQ + diaQ, 1e-9);
}
