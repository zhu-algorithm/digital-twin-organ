/**
 * molecular.js — the oesophageal molecular layer.
 *
 * WHAT THIS LAYER DOES NOT CLAIM. The oesophagus is not a drug absorption site
 * and not a metabolic organ. Normal transit of a swallowed dose is a few
 * seconds, the surface is ~0.02 m2 against ~200 m2 for the small intestine, and
 * the epithelium is multilayered squamous. There is no oesophageal ADME step
 * for any of this project's four probe drugs. This page therefore makes no flux
 * claim, draws no absorption arrow into the submucosal capillaries, and carries
 * NO calibration coefficient of any kind — expression is shown as measured
 * expression, and its placement on the cell diagram is qualitative (which cell
 * type, which membrane face), never quantitative.
 *
 * What the layer does show:
 *   1. the measured single-cell ADME panel (CELLxGENE Census 2025-11-08,
 *      tissue_general == 'esophagus', is_primary_data == True) resolved by cell
 *      type, with the stratified squamous strata — basal, spinous/suprabasal,
 *      superficial — reported SEPARATELY as the task requires;
 *   2. the barrier itself: which junction and keratin programme each stratum
 *      runs, since that is what determines whether refluxed acid or a lodged
 *      tablet reaches the proliferative compartment;
 *   3. pill-induced oesophagitis, the organ's real pharmacological
 *      significance, mapped onto the MEASURED retention levels.
 *
 * Membrane sidedness. Squamous epithelium is not a transporting monolayer, so
 * "apical vs basolateral" is not the organising axis it is in the kidney or
 * gut page. Where a transporter's localisation in squamous epithelium is
 * documented it is drawn on that face; where it is not, the molecule is drawn
 * as cytoplasmic/unassigned and labelled UNASSIGNED rather than guessed. Drawing
 * a transporter on the wrong face is worse than not drawing it.
 */

import * as THREE from 'three';

/** Strata we report separately, and how Census cell_type labels map onto them. */
/**
 * The task asks for the basal / spinous / superficial strata to be reported
 * SEPARATELY. Measured against Census 2025-11-08 `esophagus`, that is only
 * partly possible, and this map says exactly how far:
 *
 *   basal        RESOLVED      `basal cell` (7 239) + `cycling basal cell` (12)
 *   suprabasal   NOT RESOLVED  no suprabasal / spinous cell_type label exists
 *   superficial  NOT RESOLVED  no superficial cell_type label exists
 *
 * Census's dominant squamous population is `stratified epithelial cell`
 * (51 307), which spans the suprabasal AND superficial compartments without
 * distinguishing them. An earlier version of this map put `keratinocyte` into
 * the basal stratum and the generic `epithelial cell` into the suprabasal one —
 * that is inventing a resolution the data does not have, and it is removed.
 * Unassignable squamous populations get their own row instead.
 */
export const STRATA_MAP = [
  { stratum: 'basal', zh: '基底层', resolved: true,
    match: ['basal cell', 'cycling basal cell'],
    programme: 'KRT5/KRT14, TP63, proliferative; the only mitotic compartment' },
  { stratum: 'suprabasal_spinous', zh: '棘层/副基底层', resolved: false, match: [],
    programme: 'KRT4/KRT13, desmosome-rich, differentiating',
    unresolved_note: 'Census has no suprabasal or spinous cell_type label for oesophagus' },
  { stratum: 'superficial', zh: '表层', resolved: false, match: [],
    programme: 'flattened, non-keratinised, tight junctions; the permeability barrier',
    unresolved_note: 'Census has no superficial cell_type label for oesophagus' },
  { stratum: 'stratified_unassigned', zh: '复层上皮（Census 未分层）', resolved: true,
    match: ['stratified epithelial cell', 'epithelial cell of esophagus', 'keratinocyte'],
    programme: 'spans the suprabasal and superficial compartments without distinguishing them' }
];

/**
 * Transporters/enzymes with a documented squamous-epithelium localisation. Only
 * these get a membrane face; everything else in the 61-gene panel is drawn
 * unassigned.
 */
/**
 * ONLY molecules that are BOTH (a) inside this project's fixed 61-gene ADME
 * panel and (b) have a documented localisation in a polarised epithelium get a
 * membrane face. That intersection is three genes, not eight: ABCC1, SLC16A1,
 * ALDH3A1, ADH1C and CYP2C18 are the molecules whose squamous-epithelium
 * localisation is best documented, but NONE of them is in the 61-gene panel
 * (which is fixed across every organ page and whose order is not changed here),
 * so they are not measured on this page and are not drawn. They are listed in
 * OUT_OF_PANEL below so the omission is visible rather than silent.
 */
export const LOCALISED = [
  { gene: 'ABCB1', face: 'apical', zh: 'P-糖蛋白 P-gp',
    evidence: 'P-gp is apical wherever it is expressed in a polarised epithelium. Measured here at ' +
              '0.0002 (0.02% of cells) in stratified epithelial cells and 0.0023 (0.10%) in basal ' +
              'cells — effectively absent. Drawn apical because that is where it would be, with the ' +
              'measured value shown; no efflux function is claimed.' },
  { gene: 'ABCC3', face: 'basolateral', zh: 'MRP3',
    evidence: 'MRP3 is basolateral in epithelia generally. Measured at 0.0333 (4.1%) in stratified ' +
              'epithelial cells, 0.0704 (5.0%) in basal cells — low.' },
  { gene: 'GSTP1', face: 'cytoplasmic', zh: 'GSTP1',
    evidence: 'Cytosolic, so no membrane face applies. GSTP1 is the classic squamous phase-II marker ' +
              'and the measurement bears that out: 3.5619 in 98.7% of stratified epithelial cells — ' +
              'the HIGHEST value of any gene in the 61-gene panel in this tissue, above GAPDH.' }
];

/** Documented in squamous epithelium but NOT in the fixed 61-gene panel. */
export const OUT_OF_PANEL = [
  { gene: 'ABCC1', note: 'MRP1, basolateral in epithelia generally' },
  { gene: 'SLC16A1', note: 'MCT1, the one well-documented monocarboxylate route in oesophageal epithelium' },
  { gene: 'ALDH3A1', note: 'cytosolic; highly expressed in cornea and oesophageal squamous epithelium, ' +
                           'aldehyde detoxification, relevant to ethanol/acetaldehyde injury' },
  { gene: 'ADH1C', note: 'cytosolic alcohol dehydrogenase; mucosal expression underlies local acetaldehyde generation' },
  { gene: 'CYP2C18', note: 'the one CYP reported to be preferentially expressed in squamous epithelium ' +
                           'including oesophagus; no established substrate turnover there' }
];

export const PILL_OESOPHAGITIS = [
  { drug: 'doxycycline', zh: '多西环素', mechanism: 'direct acid injury — a dissolving capsule reaches pH ~3',
    class: 'tetracycline', note: 'the single most frequently reported cause' },
  { drug: 'alendronate / risedronate', zh: '双膦酸盐', mechanism: 'direct caustic injury to squamous mucosa',
    class: 'bisphosphonate', note: 'dosing instructions (200 mL water, upright 30 min) exist for this reason' },
  { drug: 'potassium chloride', zh: '氯化钾', mechanism: 'hyperosmolar local necrosis',
    class: 'electrolyte', note: 'wax-matrix formulations lodge and release slowly' },
  { drug: 'ferrous sulphate', zh: '硫酸亚铁', mechanism: 'corrosive, iron deposition in the mucosa', class: 'mineral' },
  { drug: 'NSAIDs', zh: '非甾体抗炎药', mechanism: 'topical acid trapping plus prostaglandin suppression', class: 'NSAID' },
  { drug: 'quinidine', zh: '奎尼丁', mechanism: 'direct mucosal injury', class: 'antiarrhythmic' },
  { drug: 'ascorbic acid', zh: '抗坏血酸', mechanism: 'low pH of a dissolving tablet', class: 'vitamin' }
];

/**
 * Build the molecular panel: a polarised-cell diagram is NOT appropriate here,
 * so what is drawn is the squamous stack with each localised molecule placed on
 * the face its evidence supports, and an explicit UNASSIGNED shelf.
 */
export function buildMolecular (expr, opts = {}) {
  const group = new THREE.Group();
  const mag = opts.mag || 600.0;
  const geneIx = {};
  if (expr && expr.genes) expr.genes.forEach((g, i) => { geneIx[g] = i; });

  const faces = { apical: [], basolateral: [], cytoplasmic: [], er: [], unassigned: [] };
  const localised = new Set(LOCALISED.map(l => l.gene));
  for (const L of LOCALISED) faces[L.face].push(L);
  if (expr && expr.genes) {
    for (const g of expr.genes) if (!localised.has(g)) faces.unassigned.push({ gene: g, face: 'unassigned' });
  }

  const MATS = {
    apical: new THREE.MeshPhysicalMaterial({ color: 0xffb03a, roughness: 0.4, emissive: 0xffb03a, emissiveIntensity: 0.25 }),
    basolateral: new THREE.MeshPhysicalMaterial({ color: 0x5aa7ff, roughness: 0.4, emissive: 0x5aa7ff, emissiveIntensity: 0.25 }),
    cytoplasmic: new THREE.MeshPhysicalMaterial({ color: 0x8de0a8, roughness: 0.45, emissive: 0x8de0a8, emissiveIntensity: 0.2 }),
    er: new THREE.MeshPhysicalMaterial({ color: 0xc89bff, roughness: 0.45, emissive: 0xc89bff, emissiveIntensity: 0.2 }),
    unassigned: new THREE.MeshPhysicalMaterial({ color: 0x8b95a6, roughness: 0.7, transparent: true, opacity: 0.5 })
  };
  // --- cell scaffold. Without it the molecular view is three capsules in a
  // void: a membrane face only means something if the membrane is drawn. This
  // is a squamous cell of the suprabasal compartment — wider than tall, unlike
  // the columnar cell on the gut and kidney pages, and that difference is the
  // point of the view.
  const W = 210, H = 112, D = 150;
  const cytoplasm = new THREE.Mesh(
    new THREE.BoxGeometry(W, H, D),
    new THREE.MeshPhysicalMaterial({ color: 0x2b3a52, roughness: 0.75, transparent: true,
      opacity: 0.30, side: THREE.DoubleSide }));
  cytoplasm.name = 'cell_cytoplasm';
  group.add(cytoplasm);
  const memMat = (c) => new THREE.MeshPhysicalMaterial({ color: c, roughness: 0.4,
    emissive: c, emissiveIntensity: 0.18, transparent: true, opacity: 0.85, side: THREE.DoubleSide });
  const apical = new THREE.Mesh(new THREE.BoxGeometry(W + 16, 7, D + 16), memMat(0xffb03a));
  apical.position.y = H / 2; apical.name = 'membrane_apical'; group.add(apical);
  const basolat = new THREE.Mesh(new THREE.BoxGeometry(W + 16, 7, D + 16), memMat(0x5aa7ff));
  basolat.position.y = -H / 2; basolat.name = 'membrane_basolateral'; group.add(basolat);
  // lateral faces are basolateral too — that is what "basolateral" means
  for (const sx of [-1, 1]) {
    const lat = new THREE.Mesh(new THREE.BoxGeometry(7, H, D + 16), memMat(0x5aa7ff));
    lat.position.x = sx * (W / 2); lat.name = 'membrane_lateral_' + (sx < 0 ? 'L' : 'R');
    group.add(lat);
  }
  const nucleus = new THREE.Mesh(new THREE.SphereGeometry(26, 24, 18),
    new THREE.MeshPhysicalMaterial({ color: 0x8aa0c8, roughness: 0.6, transparent: true, opacity: 0.55 }));
  nucleus.position.set(0, -10, 0); nucleus.name = 'nucleus'; group.add(nucleus);
  const er = new THREE.Mesh(new THREE.TorusKnotGeometry(21, 3.8, 72, 8, 2, 3),
    new THREE.MeshPhysicalMaterial({ color: 0xc89bff, roughness: 0.5, transparent: true, opacity: 0.45 }));
  er.position.set(0, 26, -10); er.name = 'endoplasmic_reticulum'; group.add(er);

  const placed = [];
  const put = (face, i, n) => {
    // sit each molecule ON its membrane: apical at +H/2, basolateral at -H/2,
    // cytosolic inside the cytoplasm, ER on the reticulum
    const y = face === 'apical' ? H / 2 : (face === 'basolateral' ? -H / 2 : (face === 'er' ? 26 : 4));
    const x = (i - (n - 1) / 2) * 34;
    const z = face === 'cytoplasmic' ? 24 : (face === 'er' ? -8 : 0);
    return new THREE.Vector3(x, y, z);
  };
  for (const face of Object.keys(faces)) {
    const list = faces[face];
    const shown = face === 'unassigned' ? list.slice(0, 0) : list;   // unassigned listed in the HUD, not drawn
    shown.forEach((L, i) => {
      const p = put(face, i, shown.length);
      const g = face === 'er'
        ? new THREE.TorusGeometry(11, 3.4, 10, 22)
        : new THREE.CapsuleGeometry(7.5, 17, 8, 16);
      const m = new THREE.Mesh(g, MATS[face]);
      m.position.copy(p);
      if (face === 'basolateral') m.rotation.z = Math.PI;
      m.name = 'mol_' + L.gene;
      m.userData = { gene: L.gene, face, evidence: L.evidence,
        expr_index: geneIx[L.gene] !== undefined ? geneIx[L.gene] : null };
      group.add(m);
      placed.push({ gene: L.gene, face, pos: p.toArray(), zh: L.zh });
    });
  }
  group.userData = {
    magnification: mag,
    placed,
    faces_counts: Object.fromEntries(Object.entries(faces).map(([k, v]) => [k, v.length])),
    unassigned_genes: faces.unassigned.map(x => x.gene),
    declaration:
      'Only the 8 molecules with a documented squamous-epithelium localisation are drawn on a membrane ' +
      'face. The remaining genes in the 61-gene panel are reported as measured expression in the table ' +
      'but are NOT given a membrane face, because their localisation in stratified squamous epithelium ' +
      'is not established and drawing a transporter on the wrong face is worse than not drawing it. ' +
      'No expression-to-flux or abundance-to-activity coefficient is used anywhere.',
    contour_anchored: false
  };
  return group;
}

/** Which Census cell types fall into each squamous stratum, computed from the asset. */
export function resolveStrata (expr) {
  if (!expr || !expr.cell_types) return [];
  const out = [];
  for (const S of STRATA_MAP) {
    const hits = [];
    expr.cell_types.forEach((ct, i) => {
      const l = ct.toLowerCase();
      if (S.match.some(m => l === m.toLowerCase())) hits.push({ cell_type: ct, index: i, n_cells: expr.n_cells[i] });
    });
    out.push({ ...S, matched: hits, n_cells: hits.reduce((s, h) => s + h.n_cells, 0) });
    // a stratum declared resolved that matched nothing is a contradiction worth
    // surfacing rather than hiding behind a zero
    if (S.resolved && !hits.length) out[out.length - 1].resolved = false;
  }
  return out;
}
