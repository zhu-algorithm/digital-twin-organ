/**
 * wall.js — the four-layer oesophageal wall, built at TRUE PROPORTIONS and
 * rendered at a DECLARED magnification.
 *
 * Two nested scales, both stated in the HUD and in ACCEPTANCE.md:
 *
 *   WALL_MAG  = 60x   the four-layer stack. The wall is 4.80 mm thick at the
 *                     midpoints of its published layer thicknesses, so at 60x
 *                     it is 288 mm in scene units — comparable with the 305 mm
 *                     organ beside it.
 *   EPI_MAG   = 600x  the stratified squamous epithelium. A squamous cell is
 *                     10-20 um; at 600x a 15 um cell is 9 mm, which is the
 *                     smallest thing that stays legible on a 360 px viewport.
 *
 * LAYER THICKNESSES ARE PUBLISHED, NOT MEASURED. A 1.5 mm isotropic CT label
 * cannot resolve a 0.85 mm mucosa, so there is nothing in subject s0688's
 * segmentation to anchor them to. What IS measured and used here: the lumen
 * radius at the display station, taken from the inscribed-radius profile.
 *
 * The two facts this view exists to make unmissable:
 *   1. the epithelium is NON-KERATINISED STRATIFIED SQUAMOUS — the stomach,
 *      small intestine and colon pages all show SIMPLE COLUMNAR. The abrupt
 *      change at the Z-line is drawn at the distal end.
 *   2. the outer coat is ADVENTITIA, NOT SEROSA. There is no mesothelium over
 *      the cervical or thoracic oesophagus. That is why oesophageal
 *      anastomoses leak, why T staging has no "through serosa" step, and why
 *      tumour reaches periesophageal tissue early. Only the short
 *      intra-abdominal segment carries peritoneum on its anterior face.
 */

import * as THREE from 'three';

export const WALL_MAG = 60.0;
export const EPI_MAG = 600.0;

// published thicknesses (mm): [min, max] -> the midpoint is what is drawn
export const LAYERS = [
  { key: 'mucosa', zh: '黏膜', t: [0.5, 1.2], colour: 0xe8a7a0,
    sub: 'stratified squamous epithelium · lamina propria · muscularis mucosae' },
  { key: 'submucosa', zh: '黏膜下层', t: [0.6, 1.5], colour: 0xf0cf9a,
    sub: 'oesophageal glands proper · submucosal venous plexus · Meissner plexus' },
  { key: 'muscularis_propria', zh: '肌层', t: [1.5, 3.0], colour: 0xc9736c,
    sub: 'inner circular · outer longitudinal' },
  { key: 'adventitia', zh: '外膜 (非浆膜)', t: [0.3, 1.0], colour: 0xcfc3a6,
    sub: 'loose areolar connective tissue — NO mesothelium' }
];

// epithelial strata, published (um)
export const EPI_STRATA = [
  { key: 'basal', zh: '基底层', t_um: [15, 40], cells: 2, colour: 0x6f8fc7,
    note: 'proliferative compartment; the only mitotic layer. KRT5/KRT14, TP63.' },
  { key: 'spinous', zh: '棘层 (副基底+棘细胞)', t_um: [250, 450], cells: 12, colour: 0xa9b6d4,
    note: 'differentiating; desmosome-rich. KRT4/KRT13. Dilated intercellular spaces here are the earliest reflux lesion.' },
  { key: 'superficial', zh: '表层', t_um: [80, 200], cells: 4, colour: 0xd9dde8,
    note: 'flattened, non-keratinised, tight-junction bearing; the permeability barrier.' }
];

const mid = (ab) => (ab[0] + ab[1]) / 2;

export function wallThicknessMm () {
  return LAYERS.reduce((s, l) => s + mid(l.t), 0);
}

/**
 * A 240-degree annular wedge through the wall at a chosen station, so the lumen
 * and all four layers are visible at once.
 * @param {number} lumenRadiusMm measured inscribed radius at the station
 */
export function buildWallWedge (lumenRadiusMm, opts = {}) {
  const mag = opts.mag || WALL_MAG;
  const group = new THREE.Group();
  const arc = opts.arcDeg !== undefined ? opts.arcDeg : 250;
  const a0 = -arc * Math.PI / 360, a1 = arc * Math.PI / 360;
  const depth = (opts.depthMm || 6.0) * mag;
  let r = lumenRadiusMm;
  const meta = [];
  for (const L of LAYERS) {
    const t = mid(L.t);
    const g = new THREE.RingGeometry(r * mag, (r + t) * mag, 96, 1, a0, a1 - a0);
    const m = new THREE.MeshPhysicalMaterial({
      color: L.colour, roughness: 0.62, metalness: 0.0, side: THREE.DoubleSide,
      clearcoat: L.key === 'mucosa' ? 0.5 : 0.15
    });
    const mesh = new THREE.Mesh(g, m);
    mesh.name = 'wall_' + L.key;
    group.add(mesh);
    // a thin extruded band gives the wedge visible depth without a solid volume
    const band = new THREE.Mesh(
      new THREE.CylinderGeometry((r + t) * mag, (r + t) * mag, depth, 96, 1, true, a0 + Math.PI / 2, a1 - a0),
      new THREE.MeshPhysicalMaterial({ color: L.colour, roughness: 0.66, side: THREE.DoubleSide,
        transparent: true, opacity: 0.55 }));
    band.rotation.x = Math.PI / 2;
    band.name = 'wallband_' + L.key;
    group.add(band);
    meta.push({ key: L.key, zh: L.zh, inner_mm: +r.toFixed(3), outer_mm: +(r + t).toFixed(3),
      thickness_mm: +t.toFixed(3), published_mm: L.t, sub: L.sub });
    r += t;
  }
  group.userData = {
    magnification: mag, lumen_radius_mm: lumenRadiusMm,
    wall_thickness_mm: +wallThicknessMm().toFixed(3),
    outer_radius_mm: +r.toFixed(3), layers: meta,
    contour_anchored: false,
    declaration: `layer thicknesses are published midpoints, drawn at ${mag}x; only the ` +
                 `${lumenRadiusMm.toFixed(2)} mm lumen radius is measured from subject s0688`
  };
  return group;
}

/**
 * The stratified squamous epithelium at EPI_MAG, as discrete cell bodies so the
 * basal-to-superficial flattening is visible rather than asserted.
 */
export function buildEpithelium (opts = {}) {
  const mag = opts.mag || EPI_MAG;
  const group = new THREE.Group();
  const widthUm = opts.widthUm || 420;
  let yUm = 0;
  const meta = [];
  const rng = (() => { let s = 20251030; return () => ((s = (s * 1103515245 + 12345) & 0x7fffffff) / 0x7fffffff); })();
  for (const S of EPI_STRATA) {
    const tUm = (S.t_um[0] + S.t_um[1]) / 2;
    const rows = S.cells;
    const rowUm = tUm / rows;
    const mat = new THREE.MeshPhysicalMaterial({ color: S.colour, roughness: 0.55, metalness: 0.0,
      clearcoat: 0.3 });
    for (let rI = 0; rI < rows; rI++) {
      // basal cells are cuboidal/columnar, superficial cells flattened: the
      // aspect ratio is the published morphology, drawn to scale
      const f = (yUm + rI * rowUm) / (opts.totalUm || 400);
      const cellW = 9 + 9 * f;          // um, 9 -> 18 across the stack
      const cellH = Math.max(2.0, rowUm * (S.key === 'superficial' ? 0.55 : 0.9));
      const nAcross = Math.max(3, Math.round(widthUm / cellW));
      const geo = new THREE.BoxGeometry(cellW * 0.86 * mag / 1000, cellH * 0.82 * mag / 1000,
                                        cellW * 0.86 * mag / 1000);
      const inst = new THREE.InstancedMesh(geo, mat, nAcross);
      const o = new THREE.Object3D();
      for (let c = 0; c < nAcross; c++) {
        o.position.set((c - (nAcross - 1) / 2) * cellW * mag / 1000,
                       (yUm + rI * rowUm + cellH / 2) * mag / 1000,
                       (rng() - 0.5) * 2.0);
        o.rotation.set(0, (rng() - 0.5) * 0.3, (rng() - 0.5) * 0.12);
        o.updateMatrix();
        inst.setMatrixAt(c, o.matrix);
      }
      inst.name = `epi_${S.key}_row${rI}`;
      group.add(inst);
    }
    meta.push({ key: S.key, zh: S.zh, thickness_um: +tUm.toFixed(1), published_um: S.t_um,
      cell_rows: rows, note: S.note });
    yUm += tUm;
  }
  group.userData = {
    magnification: mag, total_thickness_um: +yUm.toFixed(1),
    published_total_um: [345, 690], strata: meta, contour_anchored: false,
    declaration: `non-keratinised stratified squamous epithelium at ${mag}x; cell sizes and ` +
                 'stratum thicknesses are published histology, not measured from this subject'
  };
  return group;
}
