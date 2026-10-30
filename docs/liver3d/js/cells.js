/**
 * cells.js — the single-cell tier: a plate of real hepatocytes you can click.
 *
 * Two things make this tier different from the lobule tier above it.
 *
 * (1) EVERY RENDERED CELL IS ONE REAL CELL FROM THE CENSUS. The instances are
 *     not "a hepatocyte" drawn 352 times; each instance is bound to one row of
 *     liver_cells.json, which is one `soma_joinid` that was actually sequenced,
 *     with its own raw counts over the 61 ADME genes, its own raw_sum, donor,
 *     dataset, assay and disease status. Clicking an instance shows THAT cell's
 *     numbers. A cell-type mean is never shown in a per-cell field — where only
 *     a cell-type-level or tissue-level number exists (all protein evidence) the
 *     panel says so on that row.
 *
 * (2) THE GEOMETRY IS AT TRUE SIZE AND THE MAGNIFICATION IS DECLARED. The
 *     dimensions below are real human values; the whole tier is then scaled by
 *     SCALE_PER_UM for display and that factor is printed in the HUD.
 *
 *       hepatocyte               22 um across, polyhedral, plates one cell thick
 *       sinusoid lumen           10 um
 *       space of Disse            0.5 um
 *       bile canaliculus          0.75 um, between adjacent hepatocytes
 *       microvilli                0.5 um tall on the sinusoidal face
 *       Kupffer cell             13 um, inside the sinusoidal lumen
 *       hepatic stellate cell    10 um body in the space of Disse, long processes
 *       cholangiocyte             9 um, cuboidal, around the bile ductule
 *
 *     Sinusoidal endothelial fenestrae (50-150 nm) are BELOW the rendered
 *     resolution of this tier and are not drawn; the lining is a continuous
 *     sheet here and the acceptance record says so.
 *
 * Cell-number proportions follow the published composition of human liver by
 * cell NUMBER (hepatocytes ~60%, sinusoidal endothelium ~20%, Kupffer ~15%,
 * stellate ~5%); the counts actually instanced are reported, together with the
 * statement that a rendered cell stands for itself and not for a block of real
 * cells.
 */

import * as THREE from 'three';

export const SCALE_PER_UM = 3.0;         // scene units per micrometre (3000x on mm)

export const CELL_DIM_UM = {
  hepatocyte: 22.0,
  sinusoid_lumen: 10.0,
  disse: 0.5,
  canaliculus: 0.75,
  microvillus: 0.5,
  kupffer: 13.0,
  stellate_body: 10.0,
  cholangiocyte: 9.0,
  lsec_thickness: 0.3
};

const PLATES = 8;                  // hepatic plates across the field
const LAYERS = 2;                  // plates stacked in the vertical direction
const CELLS_ALONG = 22;            // hepatocytes from portal to central end

const COMPOSITION = {              // published human-liver cell-number shares
  hepatocyte: 0.60, sinusoidal_endothelium: 0.20, kupffer: 0.15, stellate: 0.05
};

const COL = {
  hepatocyte: 0xc8704c,
  lsec: 0xa36a86,
  kupffer: 0x4e9a6a,
  stellate: 0xd2b14a,
  cholangiocyte: 0x6f8fc0,
  canaliculus: 0x9ab83a,
  disse: 0x7a6a90
};

/** Census cell_type candidates for each rendered class, in preference order. */
const CT_CANDIDATES = {
  hepatocyte_zone1: ['periportal region hepatocyte', 'hepatocyte'],
  hepatocyte_zone2: ['midzonal region hepatocyte', 'hepatocyte'],
  hepatocyte_zone3: ['centrilobular region hepatocyte', 'hepatocyte'],
  sinusoidal_endothelium: ['endothelial cell of hepatic sinusoid',
    'endothelial cell of pericentral hepatic sinusoid',
    'endothelial cell of periportal hepatic sinusoid',
    'liver sinusoidal endothelial cell', 'endothelial cell'],
  kupffer: ['Kupffer cell', 'macrophage'],
  stellate: ['hepatic stellate cell'],
  cholangiocyte: ['cholangiocyte', 'intrahepatic cholangiocyte', 'epithelial cell of bile duct']
};

function polyhedralCell (segments = 2) {
  // A box subdivided then pushed part-way towards its circumsphere: hepatocytes
  // are polyhedral with flat contact faces, not spheres.
  const g = new THREE.BoxGeometry(1, 1, 1, segments, segments, segments);
  const p = g.getAttribute('position');
  const v = new THREE.Vector3();
  for (let i = 0; i < p.count; i++) {
    v.fromBufferAttribute(p, i);
    const l = v.length();
    v.lerp(v.clone().multiplyScalar(0.866 / Math.max(l, 1e-6)), 0.45);
    p.setXYZ(i, v.x, v.y, v.z);
  }
  g.computeVertexNormals();
  return g;
}

function starCell () {
  // stellate cell: a small body with long perisinusoidal processes
  const parts = [];
  const body = new THREE.IcosahedronGeometry(0.42, 1);
  parts.push(body);
  for (let k = 0; k < 6; k++) {
    const th = (k / 6) * Math.PI * 2;
    const arm = new THREE.CylinderGeometry(0.055, 0.02, 1.05, 5);
    arm.translate(0, 0.52, 0);
    arm.rotateZ(Math.PI / 2);
    arm.rotateY(th);
    parts.push(arm);
  }
  return mergeGeoms(parts);
}

function ruffledCell () {
  // Kupffer cell: irregular body with pseudopodia reaching across the lumen
  const parts = [new THREE.IcosahedronGeometry(0.46, 1)];
  for (let k = 0; k < 5; k++) {
    const th = (k / 5) * Math.PI * 2 + 0.4;
    const arm = new THREE.ConeGeometry(0.1, 0.62, 5);
    arm.translate(0, 0.42, 0);
    arm.rotateZ(Math.PI / 2 * (k % 2 ? 1 : -1) * 0.7);
    arm.rotateY(th);
    parts.push(arm);
  }
  return mergeGeoms(parts);
}

function mergeGeoms (list) {
  let nv = 0; let ni = 0;
  for (const g of list) {
    nv += g.getAttribute('position').count;
    ni += g.index ? g.index.count : g.getAttribute('position').count;
  }
  const pos = new Float32Array(nv * 3);
  const nrm = new Float32Array(nv * 3);
  const idx = new Uint32Array(ni);
  let vo = 0; let io = 0;
  for (const g of list) {
    const p = g.getAttribute('position'); const n = g.getAttribute('normal');
    pos.set(p.array, vo * 3);
    if (n) nrm.set(n.array, vo * 3);
    const gi = g.index;
    if (gi) { for (let i = 0; i < gi.count; i++) idx[io + i] = gi.getX(i) + vo; io += gi.count; }
    else { for (let i = 0; i < p.count; i++) idx[io + i] = i + vo; io += p.count; }
    vo += p.count;
  }
  const out = new THREE.BufferGeometry();
  out.setAttribute('position', new THREE.BufferAttribute(pos, 3));
  out.setAttribute('normal', new THREE.BufferAttribute(nrm, 3));
  out.setIndex(new THREE.BufferAttribute(idx, 1));
  return out;
}

export class CellTier {
  /**
   * @param {object} o
   * @param {object} o.cellData   parsed liver_cells.json (merged shards)
   * @param {object} o.model      parsed liver_model.json
   */
  constructor (o) {
    this.cellData = o.cellData;
    this.model = o.model;
    this.group = new THREE.Group();
    this.group.name = 'cell_tier';
    this.S = SCALE_PER_UM;
    this.pick = [];                 // instance -> {cls, cell}
    this.meshes = {};
    this.mapping = {};
    this._index();
    this._build();
  }

  /** Bucket the sampled Census cells by cell_type, preserving file order. */
  _index () {
    const byType = {};
    for (const cell of this.cellData.cells) {
      (byType[cell.cell_type] = byType[cell.cell_type] || []).push(cell);
    }
    this.byType = byType;
    this.cursor = {};
    this.resolved = {};
    for (const [cls, cands] of Object.entries(CT_CANDIDATES)) {
      const hit = cands.find(c => byType[c] && byType[c].length);
      this.resolved[cls] = hit || null;
      this.mapping[cls] = {
        census_cell_type: hit || null,
        candidates_tried: cands,
        fallback_used: hit ? cands.indexOf(hit) > 0 : false,
        available: hit ? byType[hit].length : 0
      };
      this.cursor[hit] = 0;
    }
  }

  /** Next real cell of the Census type bound to this rendered class. */
  _take (cls) {
    const t = this.resolved[cls];
    if (!t) return null;
    const arr = this.byType[t];
    const i = (this.cursor[t] || 0) % arr.length;
    this.cursor[t] = i + 1;
    return arr[i];
  }

  _build () {
    const S = this.S;
    const D = CELL_DIM_UM;
    const hep = D.hepatocyte;
    const platePitchZ = hep + D.sinusoid_lumen + 2 * D.disse;   // plate + sinusoid
    const layerPitchY = hep * 1.04;
    const along = D.hepatocyte * 1.02;
    this.fieldUm = {
      x: CELLS_ALONG * along, y: LAYERS * layerPitchY, z: PLATES * platePitchZ
    };

    const mat = (col, opt = {}) => new THREE.MeshPhysicalMaterial(Object.assign({
      color: new THREE.Color(col), roughness: 0.52, metalness: 0.0,
      clearcoat: 0.25, clearcoatRoughness: 0.4
    }, opt));

    // ------------------------------------------------------- hepatocytes
    const hepGeo = polyhedralCell(2);
    const nHep = PLATES * LAYERS * CELLS_ALONG;
    const hepMesh = new THREE.InstancedMesh(hepGeo, mat(COL.hepatocyte), nHep);
    hepMesh.name = 'hepatocytes';
    // InstancedMesh multiplies instanceColor into the material colour, so a
    // zero-filled buffer renders every hepatocyte black. It starts at white and
    // is only ever rewritten by setGeneTint() below.
    hepMesh.instanceColor = new THREE.InstancedBufferAttribute(
      new Float32Array(nHep * 3).fill(1.0), 3);
    const m = new THREE.Matrix4();
    const q = new THREE.Quaternion();
    const sc = new THREE.Vector3();
    const pos = new THREE.Vector3();
    const rnd = mulberry32(20261030);
    let k = 0;
    this.hepatocyteSlots = [];
    for (let p = 0; p < PLATES; p++) {
      for (let l = 0; l < LAYERS; l++) {
        for (let i = 0; i < CELLS_ALONG; i++) {
          const zone01 = (i + 0.5) / CELLS_ALONG;        // 0 portal -> 1 central
          const cls = zone01 < 1 / 3 ? 'hepatocyte_zone1'
            : (zone01 < 2 / 3 ? 'hepatocyte_zone2' : 'hepatocyte_zone3');
          const cell = this._take(cls);
          pos.set((i - CELLS_ALONG / 2 + 0.5) * along * S,
                  (l - (LAYERS - 1) / 2) * layerPitchY * S,
                  (p - (PLATES - 1) / 2) * platePitchZ * S);
          // real hepatocytes are not a perfect lattice: jitter is a fraction of
          // a cell diameter, it never changes the declared cell size
          pos.x += (rnd() - 0.5) * 0.10 * hep * S;
          pos.y += (rnd() - 0.5) * 0.08 * hep * S;
          q.setFromEuler(new THREE.Euler((rnd() - 0.5) * 0.5, (rnd() - 0.5) * 0.5,
                                         (rnd() - 0.5) * 0.5));
          const sz = hep * S * (0.92 + 0.16 * rnd());
          sc.set(sz, sz * 0.98, sz);
          m.compose(pos, q, sc);
          hepMesh.setMatrixAt(k, m);
          this.hepatocyteSlots.push({ index: k, zone01 });
          this.pick.push({ mesh: 'hepatocytes', instance: k, cls, cell, zone01 });
          k++;
        }
      }
    }
    hepMesh.instanceMatrix.needsUpdate = true;
    this.group.add(hepMesh);
    this.meshes.hepatocytes = hepMesh;

    // ---------------------------------------------- sinusoids + Disse space
    const sinLenUm = CELLS_ALONG * along;
    const sinGeo = new THREE.CylinderGeometry(0.5, 0.5, 1, 14, 1, true);
    sinGeo.rotateZ(Math.PI / 2);
    const sinMat = new THREE.MeshPhysicalMaterial({
      color: 0x8a4a52, roughness: 0.42, transparent: true, opacity: 0.30,
      side: THREE.DoubleSide, depthWrite: false
    });
    const disseMat = new THREE.MeshPhysicalMaterial({
      color: COL.disse, roughness: 0.6, transparent: true, opacity: 0.22,
      side: THREE.DoubleSide, depthWrite: false
    });
    this.sinusoidPaths = [];
    const nSin = PLATES - 1;
    const sinMesh = new THREE.InstancedMesh(sinGeo, sinMat, nSin * LAYERS);
    const disseMesh = new THREE.InstancedMesh(sinGeo, disseMat, nSin * LAYERS);
    let si = 0;
    for (let p = 0; p < nSin; p++) {
      for (let l = 0; l < LAYERS; l++) {
        const z = (p - (PLATES - 1) / 2 + 0.5) * platePitchZ * S;
        const y = (l - (LAYERS - 1) / 2) * layerPitchY * S;
        pos.set(0, y, z);
        q.identity();
        sc.set(sinLenUm * S, D.sinusoid_lumen * S, D.sinusoid_lumen * S * 0.72);
        m.compose(pos, q, sc);
        sinMesh.setMatrixAt(si, m);
        sc.set(sinLenUm * S, (D.sinusoid_lumen + 2 * D.disse) * S,
               (D.sinusoid_lumen + 2 * D.disse) * S * 0.72);
        m.compose(pos, q, sc);
        disseMesh.setMatrixAt(si, m);
        const pts = [];
        for (let i = 0; i < 24; i++) {
          const t = i / 23;
          pts.push([(-0.5 + t) * sinLenUm * S, y, z]);
        }
        this.sinusoidPaths.push(pts);
        si++;
      }
    }
    sinMesh.instanceMatrix.needsUpdate = true;
    disseMesh.instanceMatrix.needsUpdate = true;
    sinMesh.renderOrder = 2; disseMesh.renderOrder = 1;
    this.group.add(disseMesh, sinMesh);
    this.meshes.sinusoid_lumen = sinMesh;

    // ------------------------------------------------- bile canaliculi
    const canGeo = new THREE.CylinderGeometry(0.5, 0.5, 1, 8);
    canGeo.rotateZ(Math.PI / 2);
    const canMesh = new THREE.InstancedMesh(canGeo, mat(COL.canaliculus, {
      roughness: 0.4, emissive: new THREE.Color(0x27300c), emissiveIntensity: 0.5
    }), PLATES * LAYERS);
    this.canaliculusPaths = [];
    let ci = 0;
    for (let p = 0; p < PLATES; p++) {
      for (let l = 0; l < LAYERS; l++) {
        const z = (p - (PLATES - 1) / 2) * platePitchZ * S;
        const y = (l - (LAYERS - 1) / 2) * layerPitchY * S;
        pos.set(0, y, z);
        q.identity();
        sc.set(sinLenUm * S, D.canaliculus * S, D.canaliculus * S);
        m.compose(pos, q, sc);
        canMesh.setMatrixAt(ci, m);
        const pts = [];
        for (let i = 0; i < 24; i++) {
          const t = i / 23;
          pts.push([(0.5 - t) * sinLenUm * S, y, z]);   // central -> portal
        }
        this.canaliculusPaths.push(pts);
        ci++;
      }
    }
    canMesh.instanceMatrix.needsUpdate = true;
    this.group.add(canMesh);
    this.meshes.canaliculi = canMesh;

    // --------------------------------------------------------- microvilli
    // Rendered at the real height (0.5 um) on the sinusoidal faces. The surface
    // DENSITY is schematic: no measurement of microvillus density was used.
    const mvGeo = new THREE.ConeGeometry(0.5, 1, 4);
    mvGeo.translate(0, 0.5, 0);
    const N_MV = 1600;
    const mvMesh = new THREE.InstancedMesh(mvGeo, mat(0xc06a4e, { roughness: 0.75 }), N_MV);
    for (let i = 0; i < N_MV; i++) {
      const p = Math.floor(rnd() * nSin);
      const l = Math.floor(rnd() * LAYERS);
      const side = rnd() < 0.5 ? -1 : 1;
      const z = (p - (PLATES - 1) / 2 + 0.5) * platePitchZ * S
              + side * (D.sinusoid_lumen * 0.5 + D.disse) * S;
      const y = (l - (LAYERS - 1) / 2) * layerPitchY * S + (rnd() - 0.5) * hep * S * 0.8;
      pos.set((rnd() - 0.5) * sinLenUm * S, y, z);
      q.setFromUnitVectors(new THREE.Vector3(0, 1, 0),
                           new THREE.Vector3(0, 0, -side));
      const h = D.microvillus * S;
      sc.set(h * 0.5, h, h * 0.5);
      m.compose(pos, q, sc);
      mvMesh.setMatrixAt(i, m);
    }
    mvMesh.instanceMatrix.needsUpdate = true;
    this.group.add(mvMesh);
    this.meshes.microvilli = mvMesh;

    // ------------------------------------------------ non-parenchymal cells
    const nLSEC = Math.round(nHep / COMPOSITION.hepatocyte * COMPOSITION.sinusoidal_endothelium);
    const nKup = Math.round(nHep / COMPOSITION.hepatocyte * COMPOSITION.kupffer);
    const nStel = Math.round(nHep / COMPOSITION.hepatocyte * COMPOSITION.stellate);
    const nChol = 24;

    // sinusoidal endothelium: thin flattened cells lining the lumen
    const lsecGeo = new THREE.SphereGeometry(0.5, 10, 6);
    const lsecMesh = new THREE.InstancedMesh(lsecGeo,
      mat(COL.lsec, { roughness: 0.35, clearcoat: 0.5 }), nLSEC);
    for (let i = 0; i < nLSEC; i++) {
      const p = Math.floor(rnd() * nSin);
      const l = Math.floor(rnd() * LAYERS);
      const side = rnd() < 0.5 ? -1 : 1;
      const z = (p - (PLATES - 1) / 2 + 0.5) * platePitchZ * S
              + side * (D.sinusoid_lumen * 0.5) * S;
      pos.set((rnd() - 0.5) * sinLenUm * S * 0.96,
              (l - (LAYERS - 1) / 2) * layerPitchY * S + (rnd() - 0.5) * hep * S * 0.5, z);
      q.identity();
      sc.set(9 * S, 5 * S, CELL_DIM_UM.lsec_thickness * 4 * S);
      m.compose(pos, q, sc);
      lsecMesh.setMatrixAt(i, m);
      this.pick.push({ mesh: 'lsec', instance: i, cls: 'sinusoidal_endothelium',
                       cell: this._take('sinusoidal_endothelium') });
    }
    lsecMesh.instanceMatrix.needsUpdate = true;
    lsecMesh.name = 'lsec';
    this.group.add(lsecMesh);
    this.meshes.lsec = lsecMesh;

    // Kupffer cells sit INSIDE the sinusoidal lumen, adherent to the endothelium
    const kupMesh = new THREE.InstancedMesh(ruffledCell(),
      mat(COL.kupffer, { roughness: 0.46 }), nKup);
    for (let i = 0; i < nKup; i++) {
      const p = Math.floor(rnd() * nSin);
      const l = Math.floor(rnd() * LAYERS);
      const z = (p - (PLATES - 1) / 2 + 0.5) * platePitchZ * S;
      pos.set((rnd() - 0.5) * sinLenUm * S * 0.94,
              (l - (LAYERS - 1) / 2) * layerPitchY * S + (rnd() - 0.5) * 3 * S, z);
      q.setFromEuler(new THREE.Euler(rnd() * 6, rnd() * 6, rnd() * 6));
      sc.setScalar(D.kupffer * S);
      m.compose(pos, q, sc);
      kupMesh.setMatrixAt(i, m);
      this.pick.push({ mesh: 'kupffer', instance: i, cls: 'kupffer',
                       cell: this._take('kupffer') });
    }
    kupMesh.instanceMatrix.needsUpdate = true;
    kupMesh.name = 'kupffer';
    this.group.add(kupMesh);
    this.meshes.kupffer = kupMesh;

    // stellate cells sit in the space of Disse, wrapping the sinusoid
    const stelMesh = new THREE.InstancedMesh(starCell(),
      mat(COL.stellate, { roughness: 0.55 }), nStel);
    for (let i = 0; i < nStel; i++) {
      const p = Math.floor(rnd() * nSin);
      const l = Math.floor(rnd() * LAYERS);
      const side = rnd() < 0.5 ? -1 : 1;
      const z = (p - (PLATES - 1) / 2 + 0.5) * platePitchZ * S
              + side * (D.sinusoid_lumen * 0.5 + D.disse * 1.4) * S;
      pos.set((rnd() - 0.5) * sinLenUm * S * 0.94,
              (l - (LAYERS - 1) / 2) * layerPitchY * S + (rnd() - 0.5) * 6 * S, z);
      q.setFromEuler(new THREE.Euler(0, rnd() * 6, rnd() * 6));
      sc.setScalar(D.stellate_body * S);
      m.compose(pos, q, sc);
      stelMesh.setMatrixAt(i, m);
      this.pick.push({ mesh: 'stellate', instance: i, cls: 'stellate',
                       cell: this._take('stellate') });
    }
    stelMesh.instanceMatrix.needsUpdate = true;
    stelMesh.name = 'stellate';
    this.group.add(stelMesh);
    this.meshes.stellate = stelMesh;

    // cholangiocytes: a cuboidal ring around the bile ductule at the portal end
    const cholMesh = new THREE.InstancedMesh(polyhedralCell(1),
      mat(COL.cholangiocyte, { roughness: 0.4 }), nChol);
    const ductR = 7.0;              // interlobular bile ductule radius, um
    for (let i = 0; i < nChol; i++) {
      const ring = Math.floor(i / 8);
      const th = ((i % 8) / 8) * Math.PI * 2;
      pos.set((-CELLS_ALONG / 2 - 2.4) * along * S,
              (ring - 1) * D.cholangiocyte * 1.05 * S + Math.cos(th) * ductR * S,
              Math.sin(th) * ductR * S + (PLATES / 2 + 0.6) * platePitchZ * S);
      q.setFromEuler(new THREE.Euler(th, 0, 0));
      sc.setScalar(D.cholangiocyte * S);
      m.compose(pos, q, sc);
      cholMesh.setMatrixAt(i, m);
      this.pick.push({ mesh: 'cholangiocyte', instance: i, cls: 'cholangiocyte',
                       cell: this._take('cholangiocyte') });
    }
    cholMesh.instanceMatrix.needsUpdate = true;
    cholMesh.name = 'cholangiocyte';
    this.group.add(cholMesh);
    this.meshes.cholangiocyte = cholMesh;

    // --------------------------------------------------------- highlight
    const hl = new THREE.Mesh(
      new THREE.SphereGeometry(1, 20, 14),
      new THREE.MeshBasicMaterial({ color: 0xffe28a, wireframe: true,
                                    transparent: true, opacity: 0.85 }));
    hl.visible = false;
    hl.renderOrder = 50;
    this.group.add(hl);
    this.highlight = hl;

    // ------------------------------------------------------------- flow
    const hmo = this.model.haemodynamics;
    const transit = (hmo.sinusoid_transit_s[0] + hmo.sinusoid_transit_s[1]) / 2;
    const vSinMmS = 0.5 / transit;                       // lobule radius / transit
    const vSinUmS = vSinMmS * 1000.0;
    // bile: same derivation as lobule.js, re-expressed per micrometre
    const liverMassG = (hmo.liver_mass_g[0] + hmo.liver_mass_g[1]) / 2;
    const lobuleVol = 2.598 * 0.5 * 0.5 * 1.6;
    const nLobules = (liverMassG / 1.05) * 1000.0 / lobuleVol;
    const bileMm3S = ((hmo.bile_flow_mL_min[0] + hmo.bile_flow_mL_min[1]) / 2) * 1000 / 60;
    const aCan = Math.PI * Math.pow(D.canaliculus / 2 / 1000, 2) * 6;   // mm^2, 6 canaliculi
    const vBileMmS = (bileMm3S / nLobules) / aCan;
    const vBileUmS = vBileMmS * 1000.0;
    const BILE_GAIN = 120.0;        // declared display gain, as on the lobule tier
    this.flowSpec = {
      bile_display_gain: BILE_GAIN,
      paths: [
        ...this.sinusoidPaths.map(p => ({
          tree: 'portal_vein', retrograde: false, points: p,
          v_mm_s: new Array(p.length).fill(vSinUmS * S), display_gain: 1.0
        })),
        ...this.canaliculusPaths.map(p => ({
          tree: 'bile_duct', retrograde: true, points: p,
          v_mm_s: new Array(p.length).fill(vBileUmS * S), display_gain: BILE_GAIN
        }))
      ]
    };

    // Declared cutaway: the front hepatic plate is hidden so the sinusoid behind
    // it, its endothelial lining, the Kupffer and stellate cells in and beside the
    // lumen and the canaliculus inside the plate are visible from outside. Same
    // device as the lobule tier's 60-degree sector, and counted here.
    this.plateCutawayCells = LAYERS * CELLS_ALONG;
    this.hepatocyteTotal = nHep;
    this.setPlateCutaway(true);

    const withCell = this.pick.filter(x => x.cell).length;
    this.info = {
      magnification: `${(SCALE_PER_UM * 1000).toFixed(0)}x (scene mm per tissue mm)`,
      scale_scene_units_per_um: SCALE_PER_UM,
      field_of_view_um: {
        portal_to_central: Number(this.fieldUm.x.toFixed(1)),
        vertical: Number(this.fieldUm.y.toFixed(1)),
        across_plates: Number(this.fieldUm.z.toFixed(1))
      },
      dimensions_um: CELL_DIM_UM,
      plates: PLATES, layers: LAYERS, cells_along_axis: CELLS_ALONG,
      rendered: {
        hepatocytes: nHep, sinusoidal_endothelium: nLSEC, kupffer: nKup,
        stellate: nStel, cholangiocytes: nChol,
        sinusoid_lumens: nSin * LAYERS, canaliculi: PLATES * LAYERS,
        microvilli: N_MV
      },
      clickable_cells: this.pick.length,
      plate_cutaway: {
        enabled: true,
        hepatocytes_hidden: LAYERS * CELLS_ALONG,
        hepatocytes_shown: nHep - LAYERS * CELLS_ALONG,
        reason: 'the front hepatic plate is hidden so the sinusoid, its endothelial '
              + 'lining, the Kupffer and stellate cells and the canaliculus are visible '
              + 'from outside; the same device as the lobule tier 60-degree sector'
      },
      clickable_bound_to_a_real_census_cell: withCell,
      composition_source: 'published human-liver cell-NUMBER shares '
        + '(hepatocyte 60%, sinusoidal endothelium 20%, Kupffer 15%, stellate 5%)',
      microvillus_density: 'SCHEMATIC — height is real (0.5 um), surface density is not '
        + 'taken from any measurement',
      fenestrae: 'NOT RENDERED — sinusoidal endothelial fenestrae are 50-150 nm, below '
        + 'this tier\'s rendered resolution; the lining is drawn as a continuous sheet',
      scaling_claim: `${this.pick.length} cells are rendered. An adult human liver holds `
        + 'of the order of 2e11 hepatocytes. A rendered cell does NOT stand for a block '
        + 'of real cells: it IS the one Census cell that was sampled, and its gene values '
        + 'are that cell\'s own measurements.',
      sinusoid_velocity_um_s: Number(vSinUmS.toFixed(3)),
      canalicular_velocity_um_s: Number(vBileUmS.toExponential(3)),
      bile_display_gain: BILE_GAIN,
      census_type_mapping: this.mapping
    };
  }

  /**
   * Tint every rendered hepatocyte by ITS OWN measured value for one gene.
   *
   * This is the per-cell measurement made visible in the geometry rather than
   * only in the click panel: instance i is coloured from the CP10K-log1p value
   * of the single Census cell bound to instance i, computed from that cell's own
   * raw count and its own raw_sum. It is NOT a cell-type mean painted onto a
   * cell. The ramp is normalised by the MAXIMUM over the rendered hepatocytes,
   * and that maximum is returned so the colour scale is a stated quantity.
   *
   * @param {number|null} geneIndex index into cellData.genes, or null to reset
   */
  setGeneTint (geneIndex) {
    const mesh = this.meshes.hepatocytes;
    const col = mesh.instanceColor;
    if (geneIndex === null || geneIndex === undefined) {
      col.array.fill(1.0);
      col.needsUpdate = true;
      this.geneTint = null;
      return null;
    }
    const vals = new Float32Array(this.hepatocyteTotal);
    let vmax = 0; let nWith = 0;
    for (const rec of this.pick) {
      if (rec.mesh !== 'hepatocytes' || !rec.cell) continue;
      let v = 0;
      for (const [i, x] of rec.cell.counts) {
        if (i === geneIndex) { v = Math.log1p(x / rec.cell.raw_sum * 1e4); break; }
      }
      vals[rec.instance] = v;
      if (v > vmax) vmax = v;
      if (v > 0) nWith++;
    }
    for (let i = 0; i < this.hepatocyteTotal; i++) {
      const t = vmax > 0 ? Math.min(vals[i] / vmax, 1) : 0;
      // cool slate where the gene is absent in that cell, warm amber where it is
      // highest among the rendered cells
      col.array[i * 3] = 0.30 + 1.30 * t;
      col.array[i * 3 + 1] = 0.36 + 0.52 * t;
      col.array[i * 3 + 2] = 0.62 - 0.42 * t;
    }
    col.needsUpdate = true;
    this.geneTint = {
      gene_index: geneIndex,
      max_cp10k_log1p_over_rendered_hepatocytes: Number(vmax.toFixed(4)),
      hepatocytes_with_nonzero_count: nWith,
      hepatocytes_rendered: this.hepatocyteTotal,
      scale: 'per-instance colour = that one Census cell\'s own CP10K-log1p for this '
           + 'gene, divided by the maximum over the rendered hepatocytes'
    };
    return this.geneTint;
  }

  /** Hide the front hepatic plate so the sinusoid behind it is visible. */
  setPlateCutaway (on) {
    this.cutaway = !!on;
    this.meshes.hepatocytes.count = on
      ? this.hepatocyteTotal - this.plateCutawayCells
      : this.hepatocyteTotal;
  }

  /** Raycast against the clickable meshes; returns the bound Census cell. */
  intersect (raycaster) {
    const targets = ['hepatocytes', 'lsec', 'kupffer', 'stellate', 'cholangiocyte']
      .map(n => this.meshes[n]).filter(Boolean);
    const hits = raycaster.intersectObjects(targets, false);
    if (!hits.length) return null;
    const h = hits[0];
    const meshName = h.object.name;
    const rec = this.pick.find(x => x.mesh === meshName && x.instance === h.instanceId);
    if (rec) {
      const m = new THREE.Matrix4();
      h.object.getMatrixAt(h.instanceId, m);
      const p = new THREE.Vector3().setFromMatrixPosition(m);
      const s = new THREE.Vector3().setFromMatrixScale(m);
      this.highlight.position.copy(p);
      this.highlight.scale.set(s.x * 0.72, s.y * 0.72, s.z * 0.72);
      this.highlight.visible = true;
    }
    return rec || null;
  }

  clearHighlight () { this.highlight.visible = false; }
}

function mulberry32 (a) {
  return function () {
    a |= 0; a = (a + 0x6D2B79F5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}
