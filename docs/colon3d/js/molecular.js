// molecular.js — the colonocyte membrane layer and the microbial drug-metabolism lines.
// Transporter side (apical vs basolateral) comes from the MEMBRANE field of the data, not
// from where it happens to look good.  There is NO expression-to-flux calibration anywhere
// in this project: colour and marker size are ORDINAL encodings of measured expression and
// are not transport rates.
import * as THREE from 'three';

// apical = luminal (brush border, facing the crypt lumen); basolateral = blood side
export const MEMBRANE = {
  ABCB1:  { apical: true,  protein: 'P-gp / MDR1',  role: 'efflux to lumen' },
  ABCG2:  { apical: true,  protein: 'BCRP',         role: 'efflux to lumen' },
  ABCC2:  { apical: true,  protein: 'MRP2',         role: 'efflux to lumen' },
  ABCC3:  { apical: false, protein: 'MRP3',         role: 'efflux to blood' },
  ABCC4:  { apical: false, protein: 'MRP4',         role: 'efflux to blood' },
  SLC15A1:{ apical: true,  protein: 'PEPT1',        role: 'uptake from lumen' },
  SLCO2B1:{ apical: true,  protein: 'OATP2B1',      role: 'uptake from lumen' },
  SLC22A1:{ apical: false, protein: 'OCT1',         role: 'uptake from blood' },
  SLC22A3:{ apical: false, protein: 'OCT3',         role: 'uptake from blood' },
  SLC47A1:{ apical: true,  protein: 'MATE1',        role: 'efflux to lumen' },
  CYP3A4: { apical: null,  protein: 'CYP3A4',       role: 'intracellular oxidation' },
  CYP3A5: { apical: null,  protein: 'CYP3A5',       role: 'intracellular oxidation' },
  CYP2C9: { apical: null,  protein: 'CYP2C9',       role: 'intracellular oxidation' },
  UGT1A1: { apical: null,  protein: 'UGT1A1',       role: 'intracellular glucuronidation' },
  UGT1A3: { apical: null,  protein: 'UGT1A3',       role: 'intracellular glucuronidation' },
  UGT1A9: { apical: null,  protein: 'UGT1A9',       role: 'intracellular glucuronidation' },
  UGT2B7: { apical: null,  protein: 'UGT2B7',       role: 'intracellular glucuronidation' },
  UGT2B15:{ apical: null,  protein: 'UGT2B15',      role: 'intracellular glucuronidation' },
  SULT1A1:{ apical: null,  protein: 'SULT1A1',      role: 'intracellular sulfation' },
  SULT1E1:{ apical: null,  protein: 'SULT1E1',      role: 'intracellular sulfation' },
  SULT2A1:{ apical: null,  protein: 'SULT2A1',      role: 'intracellular sulfation' },
  CES2:   { apical: null,  protein: 'CES2',         role: 'intracellular hydrolysis' },
  NAT1:   { apical: null,  protein: 'NAT1',         role: 'intracellular acetylation' },
  NAT2:   { apical: null,  protein: 'NAT2',         role: 'intracellular acetylation' }
};

export function buildColonocyte(expr, genes, cellType) {
  const g = new THREE.Group();
  const W = 34, H = 46, D = 22;
  const cellMat = new THREE.MeshStandardMaterial({ color: 0xe0b4a2, roughness: 0.72,
    transparent: true, opacity: 0.30, side: THREE.DoubleSide });
  const cell = new THREE.Mesh(new THREE.BoxGeometry(W, H, D), cellMat);
  g.add(cell);
  // apical brush border (microvilli — the colon's ONLY surface amplification, ~6.5x)
  const bb = new THREE.Group();
  for (let i = -5; i <= 5; i++) for (let k = -3; k <= 3; k++) {
    const mv = new THREE.Mesh(new THREE.CylinderGeometry(0.5, 0.5, 5.5, 5),
      new THREE.MeshStandardMaterial({ color: 0xffd9c6, roughness: 0.6 }));
    mv.position.set(i * 2.8, H / 2 + 2.8, k * 2.8);
    bb.add(mv);
  }
  g.add(bb);
  const nuc = new THREE.Mesh(new THREE.SphereGeometry(7.5, 20, 14),
    new THREE.MeshStandardMaterial({ color: 0x8aa8d8, roughness: 0.6,
      transparent: true, opacity: 0.6 }));
  nuc.position.y = -H / 2 + 11;
  g.add(nuc);

  const idx = {}; genes.forEach((gn, i) => idx[gn] = i);
  const row = expr[cellType] || [];
  const placed = [];
  const apicalKeys = Object.keys(MEMBRANE).filter(k => MEMBRANE[k].apical === true);
  const basoKeys = Object.keys(MEMBRANE).filter(k => MEMBRANE[k].apical === false);
  const cytoKeys = Object.keys(MEMBRANE).filter(k => MEMBRANE[k].apical === null);
  const maxV = Math.max(1e-6, ...Object.keys(MEMBRANE).map(k => row[idx[k]] || 0));

  const put = (keys, y, zsign, where) => keys.forEach((k, i) => {
    const v = row[idx[k]];
    if (v === undefined) return;
    const t = v / maxV;
    const col = new THREE.Color().setHSL(0.58 - 0.58 * Math.min(1, t), 0.75, 0.35 + 0.28 * Math.min(1, t));
    const size = 1.3 + 2.6 * Math.min(1, t);
    const m = new THREE.Mesh(new THREE.SphereGeometry(size, 12, 8),
      new THREE.MeshStandardMaterial({ color: col, roughness: 0.4,
        emissive: col.clone().multiplyScalar(0.3) }));
    const n = keys.length;
    if (where === 'cyto') m.position.set((i % 6 - 2.5) * 5, y + Math.floor(i / 6) * 7, 0);
    else m.position.set((i - (n - 1) / 2) * (W / Math.max(1, n)) * 0.92, y, zsign * D / 2);
    g.add(m);
    placed.push({ gene: k, protein: MEMBRANE[k].protein, membrane: where,
                  role: MEMBRANE[k].role, value: +(+v).toFixed(4) });
  });
  put(apicalKeys, H / 2 - 1.5, 0, 'apical');
  put(basoKeys, -H / 2 + 1.5, 0, 'basolateral');
  put(cytoKeys, -4, 0, 'cyto');
  g.userData.placed = placed;
  return g;
}

export const MICROBIAL_LINES = [
  { id: 'azo', title: '偶氮还原 · 结肠靶向给药的原理',
    chain: ['柳氮磺吡啶 (口服，小肠不吸收)', '结肠细菌偶氮还原酶裂解 —N=N—',
            '5-氨基水杨酸（局部活性） + 磺胺吡啶（被吸收）'],
    note: '本项目没有做任何宏基因组或酶动力学测定；这条线是文献描述的机制，不含任何速率常数。' },
  { id: 'glyc', title: '糖苷水解',
    chain: ['糖苷结合物到达结肠', '细菌糖苷酶脱糖基', '苷元被结肠吸收'],
    note: '同上：机制描述，无速率常数。' },
  { id: 'gus', title: 'β-葡萄糖醛酸酶 · 肠肝循环再激活（直接连到本项目肝脏页）',
    chain: ['肝脏 UGT 把原形药接上葡萄糖醛酸（→ 见肝脏页）', '结合物随胆汁排入肠道',
            '结肠细菌 β-葡萄糖醛酸酶水解回原形', '原形药在结肠被重吸收 → 暴露时间延长'],
    note: '临床例子：伊立替康的 SN-38G 在肠道被解离回 SN-38，是其迟发性腹泻的机制。' +
          '本项目不计算任何再吸收分数。' }
];
