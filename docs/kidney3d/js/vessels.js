/**
 * vessels.js — tube meshes for the kidney's four tubular systems.
 *
 * Each tree arrives as (nodes, parent, radius_mm) from kidney_vessels.json. For the
 * renal artery, renal vein and collecting system those radii satisfy Murray's
 * law exactly (r_parent^3 = sum r_child^3) because they are generated as
 * r_term * n_downstream^(1/3); for the glomerular units they are published
 * per-level calibres. The mesh does not change a single radius, so the calibre
 * you see is the calibre in the audit — EXCEPT where an explicit, declared
 * `rScale` is applied, which is how the post-glomerular network is made visible
 * at organ scale. The scale in force is reported in the HUD.
 *
 * Segment geometry is a swept ring with radius interpolated along the segment,
 * so a parent tapers smoothly into its children rather than stepping. Ring
 * resolution falls with radius (14 sides on the trunks, 5 on terminal twigs),
 * which is what keeps the whole vascular bed inside a few tens of thousands of
 * triangles.
 */

import * as THREE from 'three';

const TREE_COLOR = {
  renal_artery: 0xcf3328,        // arterial scarlet
  renal_vein: 0x3b64ab,          // venous blue
  glomerular: 0xf2a62c,          // the glomerulus and the second bed: amber
  collecting_system: 0xd9d98a    // urine: pale straw
};

// Per-level tint inside a tree, so the cascade is readable without a legend
const LEVEL_TINT = {
  renal_artery: 0xcf3328, segmental: 0xd9452f, interlobar: 0xe05a35,
  arcuate: 0xf07045, cortical_radiate: 0xf78a5c,
  renal_vein: 0x3b64ab, segmental_vein: 0x3f6fb8, interlobar_vein: 0x4a7cc4,
  arcuate_vein: 0x5a8bd0, cortical_radiate_vein: 0x6c9bdb,
  afferent: 0xff6b4a, glomerulus: 0xffc24d, efferent: 0xff8d3f,
  peritubular: 0xffb36b, vasa_recta_descending: 0xff9a5c,
  vasa_recta_ascending: 0x7fb0e0,
  ureter: 0xd9d98a, major_calyx: 0xe3e39a, minor_calyx: 0xeaeaa8
};

function levelOf (label) {
  const s = String(label || '').replace('_left', '').replace('_right', '');
  const keys = ['cortical_radiate_vein', 'arcuate_vein', 'interlobar_vein',
    'segmental_vein', 'renal_vein', 'cortical_radiate', 'arcuate', 'interlobar',
    'renal_artery', 'vasa_recta_descending', 'vasa_recta_ascending', 'peritubular',
    'glomerulus', 'efferent', 'afferent', 'minor_calyx', 'major_calyx', 'ureter'];
  for (const k of keys) if (s.startsWith(k)) return k;
  if (s.startsWith('segmental')) return 'segmental';
  return s;
}

function ringSides (r) {
  if (r > 3.0) return 14;
  if (r > 1.6) return 10;
  if (r > 0.9) return 7;
  return 5;
}

export { levelOf, LEVEL_TINT };

/**
 * @param {object} tree one entry of liver_vessels.json .trees
 * @param {number} rScale cosmetic radius multiplier (1.0 = true calibre)
 */
export function buildVesselMesh (tree, opts = {}) {
  const rScale = opts.rScale !== undefined ? opts.rScale : 1.0;
  const minR = opts.minRadiusMm !== undefined ? opts.minRadiusMm : 0.0;
  const nodes = tree.nodes;
  const parent = tree.parent;
  const radius = tree.radius_mm;

  const pos = []; const nrm = []; const f0 = []; const f1 = []; const idx = [];
  const col = [];
  const base3 = new THREE.Color(TREE_COLOR[tree.name] || 0x888888);
  const byLevel = {};
  const up = new THREE.Vector3();
  const a = new THREE.Vector3(); const b = new THREE.Vector3();
  const dir = new THREE.Vector3(); const t1 = new THREE.Vector3(); const t2 = new THREE.Vector3();
  let segCount = 0;

  for (let i = 0; i < nodes.length; i++) {
    const p = parent[i];
    if (p < 0) continue;
    const rA = radius[p] * rScale; const rB = radius[i] * rScale;
    if (Math.max(rA, rB) < minR) continue;
    a.fromArray(nodes[p]); b.fromArray(nodes[i]);
    dir.subVectors(b, a);
    const len = dir.length();
    if (len < 1e-4) continue;
    dir.multiplyScalar(1 / len);
    up.set(0, 1, 0);
    if (Math.abs(dir.dot(up)) > 0.94) up.set(1, 0, 0);
    t1.crossVectors(up, dir).normalize();
    t2.crossVectors(dir, t1).normalize();

    const lv = levelOf(tree.label ? tree.label[i] : '');
    const tint = new THREE.Color(LEVEL_TINT[lv] !== undefined ? LEVEL_TINT[lv] : base3);
    byLevel[lv] = (byLevel[lv] || 0) + 1;
    const n = ringSides(Math.max(rA, rB));
    const base = pos.length / 3;
    for (let e = 0; e < 2; e++) {
      const c = e === 0 ? a : b;
      const r = e === 0 ? rA : rB;
      for (let k = 0; k < n; k++) {
        const th = (k / n) * Math.PI * 2;
        const ct = Math.cos(th); const st = Math.sin(th);
        const nx = t1.x * ct + t2.x * st;
        const ny = t1.y * ct + t2.y * st;
        const nz = t1.z * ct + t2.z * st;
        pos.push(c.x + nx * r, c.y + ny * r, c.z + nz * r);
        nrm.push(nx, ny, nz);
        col.push(tint.r, tint.g, tint.b);
        // _field0 = (0, 0, zone-ish placeholder, 0); _field1.x = 0 (no ligament)
        f0.push(0.0, 0.0, e === 0 ? 0.15 : 0.75, 0.0);
        f1.push(0.0, 0.0, 0.0, 0.0);
      }
    }
    for (let k = 0; k < n; k++) {
      const k2 = (k + 1) % n;
      idx.push(base + k, base + n + k, base + n + k2);
      idx.push(base + k, base + n + k2, base + k2);
    }
    segCount++;
  }

  const g = new THREE.BufferGeometry();
  g.setAttribute('position', new THREE.Float32BufferAttribute(pos, 3));
  g.setAttribute('normal', new THREE.Float32BufferAttribute(nrm, 3));
  g.setAttribute('_field0', new THREE.Float32BufferAttribute(f0, 4));
  g.setAttribute('_field1', new THREE.Float32BufferAttribute(f1, 4));
  g.setAttribute('color', new THREE.Float32BufferAttribute(col, 3));
  g.setIndex(idx);
  g.computeBoundingSphere();

  return {
    geometry: g,
    color: TREE_COLOR[tree.name] || 0x888888,
    stats: {
      name: tree.name,
      segments_rendered: segCount,
      segments_in_tree: parent.filter(x => x >= 0).length,
      triangles: idx.length / 3,
      vertices: pos.length / 3,
      radius_mm: [Math.min(...radius), Math.max(...radius)],
      radius_scale_applied: rScale,
      segments_by_level: byLevel
    }
  };
}

export { TREE_COLOR };
