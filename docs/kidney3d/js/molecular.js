/**
 * molecular.js — renal drug handling on a POLARISED proximal tubule cell.
 *
 * The whole point of this view is the membrane. A renal tubular epithelial cell
 * has two different faces:
 *
 *   BASOLATERAL   outward, against the peritubular capillary — i.e. the BLOOD
 *   APICAL        inward, carrying the brush border, against the tubular
 *                 LUMEN — i.e. the future urine
 *
 * A secretory transporter drawn on the wrong face would move the drug the wrong
 * way, so each marker here is placed from the `membrane` field of
 * kidney_molecular.json and nothing is positioned by eye.
 *
 * The worked route is metformin, which is excreted entirely unchanged:
 *
 *   capillary blood --OCT2 (SLC22A2, basolateral)--> tubular cell
 *                   --MATE1 / MATE2-K (SLC47A1/2, apical)--> lumen --> urine
 *
 * Marker size and colour are an ORDINAL encoding of the measured mean
 * expression. There is no calibrated expression-to-flux conversion anywhere in
 * this project and none is invented here: the particles move at a single
 * declared display rate, not at a transport rate derived from expression.
 */

import * as THREE from 'three';

export const LUMEN_R = 14.0;          // scene units for this schematic view
export const CELL_H = 11.0;
export const SEG_LEN = 54.0;
const N_CELLS = 10;

const C = {
  cell: 0xb9754a,
  nucleus: 0x6a4a7a,
  brush: 0xe0b878,
  lumen: 0x2a2f1c,
  capillary: 0xc0392b,
  basolateral: 0x4f9bdd,
  apical: 0xf0a33c,
  metformin: 0x66e0b0
};

function label (pos, text, color, pri) {
  return { pos: pos.clone(), text, color, pri };
}

/**
 * @param {object} mol parsed kidney_molecular.json
 */
export function buildMolecular (mol) {
  const group = new THREE.Group();
  const labels = [];
  const paths = [];
  let triangles = 0;

  const pt = mol.transporters_by_cell_type['epithelial cell of proximal tubule'];
  const rows = pt ? pt.transporters : [];
  const maxE = Math.max(1e-6, ...rows.map(r => r.mean_expr_log1p_cp10k));

  const addMesh = (geo, mat) => {
    const m = new THREE.Mesh(geo, mat);
    group.add(m);
    triangles += (geo.index ? geo.index.count : geo.attributes.position.count) / 3;
    return m;
  };
  // Shared materials. Giving every microvillus and every transporter pin its
  // own material meant hundreds of shader compilations on first render, which
  // stalled the scene; colour now varies through a vertex-colour attribute or a
  // small cache instead.
  const brushMat = new THREE.MeshStandardMaterial({
    color: new THREE.Color(C.brush), roughness: 0.6 });
  const pinMat = new THREE.MeshStandardMaterial({
    vertexColors: true, roughness: 0.4 });
  const nucMat = new THREE.MeshStandardMaterial({
    color: new THREE.Color(C.nucleus), roughness: 0.5 });

  // ---------------------------------------------------------------- lumen
  addMesh(new THREE.CylinderGeometry(LUMEN_R, LUMEN_R, SEG_LEN, 48, 1, true),
    new THREE.MeshPhysicalMaterial({
      color: new THREE.Color(C.lumen), roughness: 0.9, metalness: 0,
      side: THREE.BackSide }))
    .rotation.z = Math.PI / 2;

  // brush border: the apical microvillous carpet. Real microvilli are about
  // 1 µm tall on a 60 µm tubule; here they are a declared schematic texture
  // ring, not a measurement.
  {
    const proto = new THREE.CylinderGeometry(0.22, 0.12, 1.5, 4);
    const inst = new THREE.InstancedMesh(proto, brushMat, 160);
    const mtx = new THREE.Matrix4();
    const dummy = new THREE.Object3D();
    for (let i = 0; i < 160; i++) {
      const a = (i / 160) * Math.PI * 2;
      dummy.position.set((i % 20) / 20 * SEG_LEN - SEG_LEN / 2,
        Math.cos(a) * (LUMEN_R - 0.6), Math.sin(a) * (LUMEN_R - 0.6));
      dummy.lookAt(new THREE.Vector3(dummy.position.x, 0, 0));
      dummy.rotateX(Math.PI / 2);
      dummy.updateMatrix();
      mtx.copy(dummy.matrix);
      inst.setMatrixAt(i, mtx);
    }
    inst.instanceMatrix.needsUpdate = true;
    group.add(inst);
    triangles += 160 * (proto.index ? proto.index.count : 0) / 3;
  }

  // -------------------------------------------------- epithelial cell ring
  const cellMat = new THREE.MeshPhysicalMaterial({
    color: new THREE.Color(C.cell), roughness: 0.52, metalness: 0,
    transparent: true, opacity: 0.93, side: THREE.DoubleSide });
  const cellMatHi = cellMat.clone();
  cellMatHi.color = new THREE.Color(0xd48f5c);
  cellMatHi.opacity = 0.97;

  const HERO = 0;                       // the cell the route is drawn on
  for (let i = 0; i < N_CELLS; i++) {
    const a0 = (i / N_CELLS) * Math.PI * 2;
    const g = new THREE.CylinderGeometry(
      LUMEN_R + CELL_H, LUMEN_R + CELL_H, SEG_LEN * 0.92, 36, 1, true,
      a0 + 0.035, (Math.PI * 2 / N_CELLS) - 0.07);
    const gi = new THREE.CylinderGeometry(
      LUMEN_R, LUMEN_R, SEG_LEN * 0.92, 36, 1, true,
      a0 + 0.035, (Math.PI * 2 / N_CELLS) - 0.07);
    for (const gg of [g, gi]) {
      const m = addMesh(gg, i === HERO ? cellMatHi : cellMat);
      m.rotation.z = Math.PI / 2;
    }
    const nuc = addMesh(new THREE.SphereGeometry(2.6, 12, 9), nucMat);
    const am = a0 + Math.PI / N_CELLS;
    // the nucleus of a proximal tubule cell sits basally
    nuc.position.set(0, Math.cos(am) * (LUMEN_R + CELL_H * 0.68),
      Math.sin(am) * (LUMEN_R + CELL_H * 0.68));
  }

  const heroA = Math.PI / N_CELLS;      // mid-angle of the hero cell
  const dirOut = new THREE.Vector3(0, Math.cos(heroA), Math.sin(heroA));
  const rBaso = LUMEN_R + CELL_H;
  const rApic = LUMEN_R;

  // ------------------------------------------------- peritubular capillary
  const cap = addMesh(new THREE.CylinderGeometry(4.2, 4.2, SEG_LEN * 1.1, 24, 1, true),
    new THREE.MeshPhysicalMaterial({ color: new THREE.Color(C.capillary),
      roughness: 0.45, metalness: 0, transparent: true, opacity: 0.5,
      side: THREE.FrontSide, depthWrite: false }));
  cap.rotation.z = Math.PI / 2;
  cap.position.copy(dirOut.clone().multiplyScalar(rBaso + 7.5));
  labels.push(label(cap.position.clone().add(new THREE.Vector3(0, 5, 0)),
    'peritubular capillary — BLOOD (2nd capillary bed)', '#ff8d7a', 12));
  labels.push(label(new THREE.Vector3(SEG_LEN * 0.36, 0, 0),
    'tubular lumen — future URINE', '#d9d98a', 12));
  labels.push(label(dirOut.clone().multiplyScalar(rBaso + 1.5)
    .add(new THREE.Vector3(-SEG_LEN * 0.34, 0, 0)),
    'BASOLATERAL membrane', '#6fb6ec', 11));
  labels.push(label(dirOut.clone().multiplyScalar(rApic - 1.5)
    .add(new THREE.Vector3(-SEG_LEN * 0.34, 0, 0)),
    'APICAL membrane + brush border', '#f5b85a', 11));

  // ------------------------------------------------------ transporter pins
  // position comes from the `membrane` field in the data, never by eye
  const placed = { basolateral: 0, apical: 0, other: 0 };
  const pins = [];
  const order = rows.filter(r => r.membrane === 'basolateral' || r.membrane === 'apical');
  order.sort((a, b) => b.mean_expr_log1p_cp10k - a.mean_expr_log1p_cp10k);
  order.forEach((r) => {
    const apical = r.membrane === 'apical';
    const k = apical ? placed.apical++ : placed.basolateral++;
    const rr = apical ? rApic + 0.9 : rBaso - 0.9;
    const spanA = 0.42;
    const a = heroA + (apical ? -1 : 1) * spanA * 0.34 +
      (k % 2 ? spanA * 0.30 : -spanA * 0.30);
    const x = -SEG_LEN * 0.30 + Math.floor(k / 2) * (SEG_LEN * 0.21);
    const p = new THREE.Vector3(x, Math.cos(a) * rr, Math.sin(a) * rr);
    const e = r.mean_expr_log1p_cp10k / maxE;
    const rad = 1.0 + 1.9 * Math.sqrt(Math.max(e, 0.02));
    const col = new THREE.Color(apical ? C.apical : C.basolateral)
      .lerp(new THREE.Color(0xffffff), 0.45 * e);
    const pg = new THREE.SphereGeometry(rad, 12, 9);
    const cols = new Float32Array(pg.attributes.position.count * 3);
    for (let v = 0; v < pg.attributes.position.count; v++) {
      cols[v * 3] = col.r; cols[v * 3 + 1] = col.g; cols[v * 3 + 2] = col.b;
    }
    pg.setAttribute('color', new THREE.BufferAttribute(cols, 3));
    const pin = addMesh(pg, pinMat);
    pin.position.copy(p);
    pins.push({ ...r, pos: p, radius: rad });
    const detTxt = r.detected
      ? `${r.mean_expr_log1p_cp10k.toFixed(3)} · ${r.pct_expressing.toFixed(1)}%`
      : 'NO DETECTION';
    labels.push(label(p, `${r.protein} (${r.gene}) ${apical ? 'apical' : 'basolat'} ${detTxt}`,
      apical ? '#f5b85a' : '#6fb6ec',
      (r.gene === 'SLC22A2' || r.gene === 'SLC47A1' || r.gene === 'SLC47A2') ? 14 : 5));
  });
  rows.filter(r => r.membrane !== 'basolateral' && r.membrane !== 'apical')
    .forEach(() => { placed.other++; });

  // ------------------------------------------------ the metformin route path
  const find = (g) => pins.find(p => p.gene === g);
  const oct2 = find('SLC22A2');
  const mate1 = find('SLC47A1');
  const mate2 = find('SLC47A2');
  const routePts = [];
  const push = (p) => routePts.push(p.clone());
  const bloodStart = cap.position.clone().add(new THREE.Vector3(-SEG_LEN * 0.34, 0, 0));
  push(bloodStart);
  if (oct2) {
    push(oct2.pos.clone().add(dirOut.clone().multiplyScalar(3.0)));
    push(oct2.pos);
    push(oct2.pos.clone().sub(dirOut.clone().multiplyScalar(2.6)));
  }
  const mid = dirOut.clone().multiplyScalar((rBaso + rApic) * 0.5);
  push(new THREE.Vector3(oct2 ? oct2.pos.x : 0, mid.y, mid.z));
  const exit = mate1 || mate2;
  if (exit) {
    push(new THREE.Vector3(exit.pos.x, mid.y, mid.z));
    push(exit.pos);
    push(exit.pos.clone().sub(dirOut.clone().multiplyScalar(3.5)));
  }
  push(new THREE.Vector3(SEG_LEN * 0.30, 0, 0));

  // resample so the particle speed is even along the route
  const rs = [];
  for (let i = 1; i < routePts.length; i++) {
    const n = Math.max(2, Math.round(routePts[i].distanceTo(routePts[i - 1]) / 1.4));
    for (let k = 0; k < n; k++) {
      rs.push(new THREE.Vector3().lerpVectors(routePts[i - 1], routePts[i], k / n));
    }
  }
  rs.push(routePts[routePts.length - 1].clone());
  paths.push({
    tree: 'metformin', retrograde: false, level: 'metformin secretion',
    points: rs.map(p => [p.x, p.y, p.z]),
    radius_mm: rs.map(() => 0.8),
    v_mm_s: rs.map(() => 26.0),
    length_mm: 1, display_gain: 1.0
  });

  labels.push(label(bloodStart.clone().add(new THREE.Vector3(0, -5, 0)),
    'metformin: excreted UNCHANGED, no CYP step', '#7ff0c0', 13));

  return {
    group,
    labels,
    flowSpec: { paths, urine_display_gain: 1.0, pulse: [1.0], rr_interval_ms: 1000 },
    info: {
      triangles,
      cell_type: 'epithelial cell of proximal tubule',
      n_cells: pt ? pt.n_cells : null,
      n_donors: pt ? pt.n_donors : null,
      transporters_placed: placed,
      transporters_unassigned_membrane: placed.other,
      encoding: 'marker radius and tint are an ORDINAL encoding of the measured mean '
                + 'expression (CP10K, log1p). They are NOT a transport rate: no '
                + 'calibrated expression-to-flux conversion exists in this project.',
      particle_speed: 'a single declared display rate, identical for every step of '
                      + 'the route; it does not encode any measured or inferred flux.',
      brush_border: 'the microvillous ring is a declared schematic. Real microvilli '
                    + 'are about 1 µm tall, far below the resolution of anything this '
                    + 'project measures.',
      membrane_rule: 'every marker is positioned from the membrane field of '
                     + 'kidney_molecular.json; basolateral markers sit on the '
                     + 'blood-facing surface and apical markers on the lumen-facing '
                     + 'surface, so a transporter cannot be drawn on the wrong face.'
    }
  };
}
