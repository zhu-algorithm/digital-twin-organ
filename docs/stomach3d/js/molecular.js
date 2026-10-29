// molecular.js — the gastric cell/molecule layer.
//
// Transporters and enzymes are placed on the membrane face the DATA says, not
// on whichever face looks tidy: drawing the proton pump basolaterally would be
// worse than not drawing it at all.
import * as THREE from 'three';

// Membrane assignment is a literature fact, carried explicitly so the scene
// cannot silently put a protein on the wrong face.
export const MEMBRANE = {
  ATP4A: 'apical', ATP4B: 'apical',
  ABCB1: 'apical', ABCG2: 'apical', ABCC2: 'apical',
  ABCC3: 'basolateral', ABCC4: 'basolateral',
  SLC15A1: 'apical', SLCO2B1: 'apical',
  SLC22A1: 'basolateral', SLC22A3: 'basolateral',
  CYP3A4: 'intracellular (ER)', CYP3A5: 'intracellular (ER)',
  CYP2C19: 'intracellular (ER)', CYP2D6: 'intracellular (ER)',
  CES2: 'intracellular', SULT1A1: 'cytosolic', GSTP1: 'cytosolic',
  UGT1A1: 'intracellular (ER)',
};

// The gastric cell types the page highlights, with the Census labels they are
// matched against.  A type that is ABSENT from the Census stomach cohort is
// reported as absent, never substituted.
// The gastric cell types the page highlights, with the Census labels they are
// matched against.  Aliases were checked against the actual stomach cohort of
// Census 2025-11-08, not guessed:
//   * the chief cell appears as "peptic cell", NOT "chief cell"
//   * the surface mucous cell appears as "foveolar cell of stomach"
//   * there is NO "enterochromaffin-like cell" and NO D cell label in this
//     cohort.  "type EC enteroendocrine cell" is the ENTEROCHROMAFFIN cell
//     (serotonin), which is a DIFFERENT cell from the enterochromaffin-LIKE
//     (ECL, histamine) cell of the oxyntic gland.  Conflating them to fill the
//     row would be a fabrication, so ECL and D are reported ABSENT.
export const KEY_CELLS = [
  ['parietal cell',   ['parietal cell'],                        '壁细胞（泌酸 + 内因子）'],
  ['chief cell',      ['peptic cell', 'chief cell',
                       'gastric chief cell'],                   '主细胞（胃蛋白酶原）'],
  ['G cell',          ['type G enteroendocrine cell', 'G cell'], 'G 细胞（胃窦，胃泌素）'],
  ['surface mucous',  ['foveolar cell of stomach',
                       'surface mucous cell of stomach',
                       'mucous cell of stomach'],               '表面黏液细胞（胃小凹）'],
  ['mucous neck',     ['mucous neck cell of gastric gland'],     '颈黏液细胞（峡部/祖细胞区）'],
  ['ECL cell',        ['enterochromaffin-like cell'],            'ECL 细胞（组胺）'],
  ['D cell',          ['type D enteroendocrine cell', 'D cell'], 'D 细胞（生长抑素）'],
];

// Labels that exist in the cohort and are NEARBY but must not be used as
// stand-ins for the two missing types above.
export const NOT_SUBSTITUTES = {
  'ECL cell': ['type EC enteroendocrine cell is the ENTEROCHROMAFFIN cell '
             + '(serotonin), a different cell from the enterochromaffin-LIKE '
             + '(histamine) cell; it is NOT used as a stand-in',
               'type X enteroendocrine cell and P/D1 enteroendocrine cell are '
             + 'ghrelin cells, also not stand-ins'],
  'D cell': ['the generic "enteroendocrine cell" label is present but pools '
           + 'several lineages and is NOT used as a stand-in'],
};

export function resolveCells(expr) {
  const out = [];
  const lower = expr.cell_types.map(c => c.toLowerCase());
  for (const [key, aliases, zh] of KEY_CELLS) {
    let idx = -1;
    for (const a of aliases) {
      const j = lower.indexOf(a.toLowerCase());
      if (j >= 0) { idx = j; break; }
    }
    out.push({ key, zh, index: idx, present: idx >= 0,
               censusLabel: idx >= 0 ? expr.cell_types[idx] : null,
               n_cells: idx >= 0 ? expr.n_cells[idx] : 0 });
  }
  return out;
}

export function geneValue(expr, cellIndex, gene) {
  const gi = expr.genes.indexOf(gene);
  if (gi < 0 || cellIndex < 0) return null;
  return { expr: expr.expr[cellIndex][gi], pct: expr.pct[cellIndex][gi] };
}

/** A schematic parietal cell with the pump on the apical face. */
export function buildParietalCell(scale = 1.0) {
  const g = new THREE.Group();
  g.name = 'parietal_cell';
  const body = new THREE.Mesh(
    new THREE.SphereGeometry(10 * scale, 24, 18),
    new THREE.MeshStandardMaterial({ color: 0xf0a33c, roughness: 0.55,
      transparent: true, opacity: 0.30, depthWrite: false }));
  body.name = 'cell_body';
  g.add(body);
  // the secretory canaliculus: the invagination of the APICAL membrane that
  // the pump sits on, and where a PPI is protonated and trapped
  const can = new THREE.Mesh(
    new THREE.TorusKnotGeometry(4.4 * scale, 1.1 * scale, 96, 10, 2, 3),
    new THREE.MeshStandardMaterial({ color: 0xffd24a, roughness: 0.4,
      emissive: 0x5a4300 }));
  can.name = 'secretory_canaliculus';
  can.userData.membrane = 'apical';
  g.add(can);
  const nuc = new THREE.Mesh(new THREE.SphereGeometry(3.0 * scale, 16, 12),
    new THREE.MeshStandardMaterial({ color: 0x8aa0c0, roughness: 0.6 }));
  nuc.position.set(-4.5 * scale, -3.5 * scale, 0);
  nuc.name = 'nucleus';
  g.add(nuc);
  return g;
}
