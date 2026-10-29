/**
 * molecular.js — the molecular layer of the oral cavity.
 *
 * Three things, in one offset frame:
 *
 * A. THE TWO ROUTES, drawn side by side. Left: a sublingual tablet dissolves,
 *    crosses the thin non-keratinised epithelium, and its venous drainage runs
 *    lingual/facial vein -> internal jugular vein -> brachiocephalic vein ->
 *    superior vena cava -> right atrium. It never meets the liver. Right: the
 *    same tablet swallowed goes stomach -> small intestine (enterocyte CYP3A4
 *    and P-gp) -> portal vein -> LIVER -> hepatic vein -> inferior vena cava.
 *    The liver block sits on the right chain only. That one difference is the
 *    whole pharmacological point of the oral cavity.
 *
 * B. A POLARISED ORAL MUCOSAL CELL with its transporters on the correct face:
 *    efflux pumps on the APICAL (saliva-facing) membrane push drug back into
 *    saliva; uptake and basolateral efflux face the lamina propria and its
 *    capillaries. Drawing a transporter on the wrong membrane would be worse
 *    than not drawing it, so each one carries its membrane in the label and the
 *    HUD states which assignments are extrapolated rather than measured in oral
 *    tissue.
 *
 * C. A POLARISED SALIVARY SECRETORY UNIT: acinar cell (basolateral NKCC1, apical
 *    AQP5) and striated duct cell (apical ENaC and CFTR, basolateral Na+/K+-
 *    ATPase, apical ABCB1), which is why saliva leaves the acinus isotonic and
 *    arrives in the mouth hypotonic.
 *
 * Expression values come from oral_cavity_expression.json (CZ CELLxGENE Census)
 * when it is present. If it is absent or empty the panel says so in plain words
 * and draws nothing - no tissue average and no other organ's values are ever
 * substituted.
 */

import * as THREE from 'three';

const APICAL = 'apical (saliva-facing)';
const BASO = 'basolateral (lamina propria / capillary-facing)';
const INTRA = 'intracellular (smooth ER / cytosol)';

// gene -> membrane face, with an explicit provenance grade for the assignment
export const MEMBRANE = {
  ABCB1: [APICAL, 'measured in human intestinal and renal epithelium; in '
          + 'stratified squamous oral epithelium the apical/superficial '
          + 'localisation is EXTRAPOLATED', 'extrapolated'],
  ABCG2: [APICAL, 'apical in intestine, placenta and mammary gland; oral '
          + 'assignment EXTRAPOLATED', 'extrapolated'],
  ABCC2: [APICAL, 'apical in hepatocyte and enterocyte; oral assignment '
          + 'EXTRAPOLATED', 'extrapolated'],
  ABCC1: [BASO, 'basolateral in most epithelia including squamous; '
          + 'well established', 'established'],
  ABCC3: [BASO, 'basolateral in hepatocyte and enterocyte', 'established'],
  ABCC4: [BASO, 'basolateral in several epithelia', 'established'],
  SLCO2B1: [BASO, 'basolateral uptake; reported in oral/oesophageal mucosa',
            'reported'],
  SLC22A3: [BASO, 'OCT3, basolateral uptake, broadly expressed', 'established'],
  SLC15A1: [APICAL, 'PEPT1, apical in enterocyte; oral assignment EXTRAPOLATED',
            'extrapolated'],
  CYP3A4: [INTRA, 'smooth endoplasmic reticulum', 'established'],
  CYP3A5: [INTRA, 'smooth endoplasmic reticulum', 'established'],
  CYP2D6: [INTRA, 'smooth endoplasmic reticulum', 'established'],
  CES1: [INTRA, 'ER lumen, retained by an HXEL motif', 'established'],
  CES2: [INTRA, 'ER lumen', 'established'],
  UGT1A1: [INTRA, 'ER membrane, lumen-facing active site', 'established'],
  UGT1A3: [INTRA, 'ER membrane', 'established'],
  UGT1A9: [INTRA, 'ER membrane', 'established'],
  UGT2B7: [INTRA, 'ER membrane', 'established'],
  SULT1A1: [INTRA, 'cytosolic', 'established'],
  GSTP1: [INTRA, 'cytosolic', 'established'],
  NAT1: [INTRA, 'cytosolic', 'established'],
  ALDH2: [INTRA, 'mitochondrial matrix', 'established'],
  ADH1B: [INTRA, 'cytosolic', 'established']
};

const ROUTE_SUBLINGUAL = [
  ['舌下黏膜 sublingual mucosa', 0xe06a5f],
  ['舌静脉 / 面静脉 lingual + facial vein', 0x4a7cc4],
  ['颈内静脉 internal jugular vein', 0x3b64ab],
  ['头臂静脉 brachiocephalic vein', 0x33579b],
  ['上腔静脉 superior vena cava', 0x2d4d8c],
  ['右心房 → 体循环 systemic', 0xcf3328]
];
const ROUTE_SWALLOWED = [
  ['吞咽 swallow', 0xb9b27a],
  ['胃 stomach', 0xc08a55],
  ['小肠 small intestine (CYP3A4 + P-gp)', 0xcf9a4a],
  ['门静脉 portal vein', 0x5c7f4a],
  ['肝脏 LIVER — 首过代谢', 0x8a4a2f],
  ['肝静脉 → 下腔静脉 IVC', 0x33579b],
  ['右心房 → 体循环 systemic', 0xcf3328]
];

function label (text, color) {
  const c = document.createElement('canvas');
  const pad = 10; const fs = 30;
  const ctx = c.getContext('2d');
  ctx.font = `600 ${fs}px -apple-system, "PingFang SC", sans-serif`;
  const w = Math.ceil(ctx.measureText(text).width) + pad * 2;
  c.width = w; c.height = fs + pad * 2;
  const g = c.getContext('2d');
  g.font = `600 ${fs}px -apple-system, "PingFang SC", sans-serif`;
  g.fillStyle = 'rgba(8,12,20,0.78)';
  g.fillRect(0, 0, c.width, c.height);
  g.fillStyle = '#' + new THREE.Color(color).getHexString();
  g.fillRect(0, 0, 5, c.height);
  g.fillStyle = '#e9eef6';
  g.textBaseline = 'middle';
  g.fillText(text, pad, c.height / 2);
  const tex = new THREE.CanvasTexture(c);
  tex.colorSpace = THREE.SRGBColorSpace;
  tex.minFilter = THREE.LinearFilter;
  const sp = new THREE.Sprite(new THREE.SpriteMaterial({
    map: tex, transparent: true, depthTest: false, depthWrite: false
  }));
  sp.scale.set(c.width * 0.19, c.height * 0.19, 1);
  sp.renderOrder = 920;
  return sp;
}

function chain (items, x0, dir) {
  const g = new THREE.Group();
  const pts = [];
  items.forEach((it, i) => {
    const y = 120 - i * 34;
    const p = new THREE.Vector3(x0 + dir * (i % 2 === 0 ? 0 : 6), y, 0);
    pts.push(p);
    const isLiver = /LIVER/.test(it[0]);
    const node = new THREE.Mesh(
      isLiver ? new THREE.BoxGeometry(34, 20, 16)
        : new THREE.SphereGeometry(7.5, 20, 14),
      new THREE.MeshPhysicalMaterial({
        color: it[1], roughness: isLiver ? 0.7 : 0.4, metalness: 0,
        clearcoat: 0.4, emissive: new THREE.Color(it[1]).multiplyScalar(0.16)
      }));
    node.position.copy(p);
    g.add(node);
    const lab = label(it[0], it[1]);
    lab.position.set(p.x + dir * 26, p.y, 2);
    g.add(lab);
  });
  const curve = new THREE.CatmullRomCurve3(pts);
  g.add(new THREE.Mesh(new THREE.TubeGeometry(curve, 120, 2.6, 10, false),
    new THREE.MeshPhysicalMaterial({ color: 0x5e6d85, roughness: 0.5,
      metalness: 0, transparent: true, opacity: 0.85 })));
  g.userData.curve = curve;
  return g;
}

/**
 * @param {object} model   parsed oral_model.json
 * @param {object|null} expr parsed oral_cavity_expression.json, or null
 */
export function buildMolecular (model, expr) {
  const root = new THREE.Group();
  const notes = [];

  // ----------------------------------------------------------- A. routes
  const routes = new THREE.Group();
  const left = chain(ROUTE_SUBLINGUAL, -150, -1);
  const right = chain(ROUTE_SWALLOWED, 90, 1);
  routes.add(left, right);
  const hA = label('舌下 / 颊黏膜：绕过肝首过 — bypasses hepatic first pass',
    0x5ad18a);
  hA.position.set(-150, 158, 0); routes.add(hA);
  const hB = label('吞咽：门静脉 → 肝脏首过 — portal vein, then the liver',
    0xf0a33c);
  hB.position.set(90, 158, 0); routes.add(hB);
  root.add(routes);

  // particles along both chains
  const mk = (curve, color, n) => {
    const pos = new Float32Array(n * 3);
    const geo = new THREE.BufferGeometry();
    geo.setAttribute('position', new THREE.BufferAttribute(pos, 3));
    const m = new THREE.Points(geo, new THREE.PointsMaterial({
      color, size: 4.2, sizeAttenuation: true, transparent: true,
      opacity: 0.95, depthWrite: false
    }));
    m.userData = { curve, n, phase: new Float32Array(n).map(() => Math.random()) };
    return m;
  };
  const pL = mk(left.userData.curve, 0x7fd0ff, 90);
  const pR = mk(right.userData.curve, 0xffc98a, 110);
  root.add(pL, pR);

  // --------------------------------------------- B. polarised mucosal cell
  const cell = new THREE.Group();
  cell.position.set(-130, -70, 0);
  const CW = 150; const CH = 78;
  cell.add(new THREE.Mesh(new THREE.BoxGeometry(CW, CH, 40),
    new THREE.MeshPhysicalMaterial({ color: 0xc9736e, roughness: 0.62,
      metalness: 0, transparent: true, opacity: 0.82 })));
  const nuc = new THREE.Mesh(new THREE.SphereGeometry(16, 24, 18),
    new THREE.MeshPhysicalMaterial({ color: 0x7e3f58, roughness: 0.5 }));
  nuc.position.set(-28, -6, 0); cell.add(nuc);
  const apic = label('顶侧 ' + APICAL, 0xf0a33c);
  apic.position.set(0, CH / 2 + 16, 0); cell.add(apic);
  const baso = label('基底侧 ' + BASO, 0x7fd0ff);
  baso.position.set(0, -CH / 2 - 16, 0); cell.add(baso);

  const genesShown = [];
  const placeT = (gene, i, nApi, nBas, nInt) => {
    const [face, why, grade] = MEMBRANE[gene] || [INTRA, 'unknown', 'unknown'];
    let p; let col;
    if (face === APICAL) {
      p = new THREE.Vector3(-CW / 2 + 18 + (i % nApi) * (CW - 36) / Math.max(nApi - 1, 1),
        CH / 2, 0); col = 0xf0a33c;
    } else if (face === BASO) {
      p = new THREE.Vector3(-CW / 2 + 18 + (i % nBas) * (CW - 36) / Math.max(nBas - 1, 1),
        -CH / 2, 0); col = 0x7fd0ff;
    } else {
      p = new THREE.Vector3(6 + (i % 3) * 26, -2 + ((i / 3) | 0) * 18 - 10, 6);
      col = 0xb0f08a;
    }
    const geo = face === INTRA ? new THREE.TorusGeometry(5.5, 2.0, 8, 18)
      : new THREE.CapsuleGeometry(3.4, 9.0, 4, 10);
    const m = new THREE.Mesh(geo, new THREE.MeshPhysicalMaterial({
      color: col, roughness: 0.35, metalness: 0.1,
      emissive: new THREE.Color(col).multiplyScalar(0.25)
    }));
    m.position.copy(p);
    if (face !== INTRA) m.rotation.z = 0;
    cell.add(m);
    const lab = label(gene, col);
    lab.position.set(p.x, p.y + (face === BASO ? -13 : 13), 10);
    lab.scale.multiplyScalar(0.62);
    cell.add(lab);
    genesShown.push({ gene, membrane: face, basis: why, grade });
  };
  const api = ['ABCB1', 'ABCG2', 'ABCC2'];
  const bas = ['ABCC1', 'ABCC3', 'SLCO2B1', 'SLC22A3'];
  const intr = ['CYP3A4', 'CYP3A5', 'CES1', 'CES2', 'UGT1A1', 'SULT1A1'];
  api.forEach((g, i) => placeT(g, i, api.length, 1, 1));
  bas.forEach((g, i) => placeT(g, i, 1, bas.length, 1));
  intr.forEach((g, i) => placeT(g, i, 1, 1, intr.length));
  root.add(cell);

  // ------------------------------------- C. polarised salivary secretory unit
  const sal = new THREE.Group();
  sal.position.set(150, -70, 0);
  const ac = new THREE.Mesh(new THREE.SphereGeometry(26, 26, 20),
    new THREE.MeshPhysicalMaterial({ color: 0xc9a87c, roughness: 0.55,
      transparent: true, opacity: 0.9 }));
  sal.add(ac);
  const lum = new THREE.Mesh(new THREE.SphereGeometry(9, 20, 14),
    new THREE.MeshBasicMaterial({ color: 0x7fd0ff, transparent: true,
      opacity: 0.55 }));
  sal.add(lum);
  const duct = new THREE.Mesh(new THREE.CylinderGeometry(7, 9, 70, 16, 1, true),
    new THREE.MeshPhysicalMaterial({ color: 0xd8dca0, roughness: 0.5,
      side: THREE.DoubleSide, transparent: true, opacity: 0.9 }));
  duct.rotation.z = Math.PI / 2; duct.position.set(56, 0, 0); sal.add(duct);
  const SALT = [
    ['NKCC1', BASO, -24, -16, 0x7fd0ff, '基底侧 Cl- 摄入，驱动等渗初级分泌'],
    ['AQP5', APICAL, -4, 12, 0xf0a33c, '顶侧水通道，水随 Cl- 进入腺泡腔'],
    ['ENaC', APICAL, 46, 10, 0xf0a33c, '纹状管顶侧回收 Na+'],
    ['CFTR', APICAL, 68, 10, 0xf0a33c, '纹状管顶侧 Cl- 通道'],
    ['ABCB1 (P-gp)', APICAL, 84, 10, 0xf0a33c, '导管顶侧外排 — 把底物泵回唾液'],
    ['ATP1A1 (Na/K-ATPase)', BASO, 60, -16, 0x7fd0ff, '基底侧钠泵']
  ];
  SALT.forEach(([g, face, x, y, col, why]) => {
    const m = new THREE.Mesh(new THREE.CapsuleGeometry(3.0, 8.0, 4, 10),
      new THREE.MeshPhysicalMaterial({ color: col, roughness: 0.35,
        emissive: new THREE.Color(col).multiplyScalar(0.25) }));
    m.position.set(x, y, 0); sal.add(m);
    const lab = label(g, col);
    lab.position.set(x, y + (face === BASO ? -12 : 12), 10);
    lab.scale.multiplyScalar(0.58);
    sal.add(lab);
    genesShown.push({ gene: g, membrane: face, basis: why, grade: 'established' });
  });
  const sh = label('唾液腺分泌单位：腺泡等渗 → 纹状管重吸收 Na+/Cl- → 低渗唾液',
    0xc9a87c);
  sh.position.set(140, 46, 0); sal.add(sh);
  root.add(sal);

  // ------------------------------------------------- expression bar panel
  const bars = new THREE.Group();
  bars.position.set(-250, -170, 0);
  let exprStatus;
  if (expr && expr.genes && expr.genes.length && expr.cell_types
      && expr.cell_types.length) {
    const want = ['CYP3A4', 'CYP3A5', 'CES1', 'CES2', 'UGT1A1', 'ABCB1',
      'ABCG2', 'SLCO2B1'];
    const gi = want.map((g) => expr.genes.indexOf(g));
    const nCT = Math.min(10, expr.cell_types.length);
    let vmax = 1e-6;
    for (let c = 0; c < nCT; c++) {
      gi.forEach((j) => { if (j >= 0) vmax = Math.max(vmax, expr.expr[c][j]); });
    }
    for (let c = 0; c < nCT; c++) {
      const lab = label(`${expr.cell_types[c]}  (n=${expr.n_cells[c]})`, 0x9fc4ff);
      lab.position.set(-96, -c * 30, 0); lab.scale.multiplyScalar(0.92);
      bars.add(lab);
      gi.forEach((j, k) => {
        const v = j >= 0 ? expr.expr[c][j] : 0;
        const h = Math.max(0.8, (v / vmax) * 24);
        const m = new THREE.Mesh(new THREE.BoxGeometry(9, h, 7),
          new THREE.MeshPhysicalMaterial({
            color: new THREE.Color().setHSL(0.58 - 0.55 * (v / vmax), 0.72, 0.52),
            roughness: 0.45, metalness: 0.05
          }));
        m.position.set(44 + k * 20, -c * 30 + h / 2 - 11, 0);
        bars.add(m);
      });
    }
    want.forEach((g, k) => {
      const lab = label(g, 0xb8c4d6);
      lab.position.set(44 + k * 20, 26, 0); lab.scale.multiplyScalar(0.75);
      bars.add(lab);
    });
    exprStatus = {
      ok: true, census_version: expr.census_version,
      cell_types: expr.cell_types.length, genes: expr.genes.length,
      cells: expr.n_cells.reduce((a, b) => a + b, 0),
      normalization: expr.normalization,
      genes_plotted: want, vmax: +vmax.toFixed(4)
    };
    // The measured profile is allowed to contradict the diagram, and when it
    // does the page says so rather than the diagram winning.
    const iABCB1 = expr.genes.indexOf('ABCB1');
    if (iABCB1 >= 0) {
      let epi = 0; let endo = 0; let epiName = ''; let endoName = '';
      expr.cell_types.forEach((ct, c) => {
        const v = expr.expr[c][iABCB1];
        if (/squamous|keratinocyte|epithelial/i.test(ct) && v > epi) {
          epi = v; epiName = ct;
        }
        if (/endothelial/i.test(ct) && v > endo) { endo = v; endoName = ct; }
      });
      if (endo > epi) {
        notes.push('实测与上面那张图不一致，照实说：ABCB1 (P-gp) 在本数据里'
          + `最高的是 ${endoName}（${endo.toFixed(3)}），`
          + `而不是 ${epiName || '鳞状上皮'}（${epi.toFixed(3)}）。`
          + '也就是说，口腔黏膜里的 P-gp 信号主要来自血管内皮与免疫细胞，'
          + '**本数据并不支持"顶侧上皮 P-gp 外排"这条叙述**。'
          + '图上该转运体仍标为 extrapolated。');
      }
    }
    notes.push(`单细胞：CELLxGENE Census ${expr.census_version}，`
      + `${expr.cell_types.length} 个细胞类型 × ${expr.genes.length} 个 ADME 基因，`
      + `${expr.n_cells.reduce((a, b) => a + b, 0).toLocaleString()} 个细胞。`
      + '颜色只能定性读：本项目没有任何"表达量→转运通量"标定。');
  } else {
    const msg = label('口腔单细胞表达数据缺失 — 未用任何替代值填补', 0xff9a8a);
    msg.position.set(-40, 0, 0); bars.add(msg);
    const msg2 = label('MISSING: no oral single-cell ADME profile in this build. '
      + 'No tissue average and no other organ was substituted.', 0xff9a8a);
    msg2.position.set(-40, -24, 0); msg2.scale.multiplyScalar(0.7); bars.add(msg2);
    exprStatus = { ok: false,
      reason: 'oral_cavity_expression.json absent or empty' };
    notes.push('口腔单细胞 ADME 谱缺失，页面不以任何替代数据填补。');
  }
  root.add(bars);

  root.userData = {
    particles: [pL, pR], genesShown, exprStatus, notes,
    first_pass: model.first_pass
  };
  return root;
}

export function stepMolecular (group, dt) {
  for (const p of group.userData.particles || []) {
    const { curve, n, phase } = p.userData;
    const arr = p.geometry.getAttribute('position').array;
    for (let i = 0; i < n; i++) {
      phase[i] = (phase[i] + dt * 0.16) % 1;
      const v = curve.getPointAt(phase[i]);
      arr[i * 3] = v.x + (Math.random() - 0.5) * 1.2;
      arr[i * 3 + 1] = v.y;
      arr[i * 3 + 2] = v.z + (Math.random() - 0.5) * 1.2;
    }
    p.geometry.getAttribute('position').needsUpdate = true;
  }
}
