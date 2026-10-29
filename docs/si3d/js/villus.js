/**
 * villus.js — the two magnified views: the cut wall segment and the
 * villus-crypt unit.
 *
 * Everything here is modelled at TRUE SIZE in millimetres from the published
 * dimensions in small_intestine_model.json and then scaled UNIFORMLY for the
 * screen, so proportions survive and only one number has to be declared per
 * view. The declared magnifications are
 *
 *   cut wall segment   x12     lumen radius 10 mm, wall 1.21 mm,
 *                              plicae 7 mm pitch, villi 0.9 mm
 *   villus-crypt unit  x60     villus 0.9 mm, crypt 0.25 mm,
 *                              enterocyte 25 um
 *   brush border       x20000  microvillus 1 um long, 0.1 um across;
 *                              the patch spans 13 x 13 microvilli, i.e. ~2.1 um
 *
 * None of it is contour-anchored: villi are an order of magnitude below the
 * 1.5 mm CT voxel, so no segmentation can carry them. That is stated in the
 * page and in ACCEPTANCE.md.
 */

import * as THREE from 'three';

export const MAG = { wall: 12, villus: 60, brush: 20000 };

export const LAYER_COLOR = {
  mucosa: 0xd9736a, lamina_propria: 0xc96a76, muscularis_mucosae: 0xa85a62,
  submucosa: 0xe0c48a, muscularis_circular: 0xb8585c,
  myenteric_plexus: 0xf0e07a, muscularis_longitudinal: 0xa04a50,
  serosa: 0xcfd6e2
};

function annulusSector (rIn, rOut, len, thetaStart, thetaLen, seg = 72) {
  // a tube wall sector with real end caps, so the cut face shows the layer
  const pos = [], nrm = [], idx = [];
  const push = (x, y, z, nx, ny, nz) => { pos.push(x, y, z); nrm.push(nx, ny, nz); };
  const z0 = -len / 2, z1 = len / 2;
  for (let s = 0; s <= seg; s++) {
    const th = thetaStart + thetaLen * (s / seg);
    const c = Math.cos(th), si = Math.sin(th);
    push(rOut * c, rOut * si, z0, c, si, 0); push(rOut * c, rOut * si, z1, c, si, 0);
    push(rIn * c, rIn * si, z0, -c, -si, 0); push(rIn * c, rIn * si, z1, -c, -si, 0);
  }
  for (let s = 0; s < seg; s++) {
    const a = s * 4, b = (s + 1) * 4;
    idx.push(a, b, b + 1, a, b + 1, a + 1);                 // outer
    idx.push(a + 2, b + 3, b + 2, a + 2, a + 3, b + 3);     // inner
  }
  // end caps at both ends of the sector plus the two cut faces
  const cap = (zi, sgn) => {
    for (let s = 0; s < seg; s++) {
      const a = s * 4 + zi, b = (s + 1) * 4 + zi;
      if (sgn > 0) idx.push(a, a + 2, b + 2, a, b + 2, b);
      else idx.push(a, b + 2, a + 2, a, b, b + 2);
    }
  };
  cap(0, -1); cap(1, 1);
  const face = (s) => {
    const o = s * 4;
    idx.push(o, o + 2, o + 3, o, o + 3, o + 1);
  };
  face(0); face(seg);
  const g = new THREE.BufferGeometry();
  g.setAttribute('position', new THREE.Float32BufferAttribute(pos, 3));
  g.setAttribute('normal', new THREE.Float32BufferAttribute(nrm, 3));
  g.setIndex(idx); g.computeVertexNormals(); g.computeBoundingSphere();
  return g;
}

function std (color, rough = 0.62, opts = {}) {
  return new THREE.MeshStandardMaterial({
    color, roughness: rough, metalness: 0.02, side: THREE.DoubleSide, ...opts
  });
}

/**
 * Cut segment of bowel wall: four layers (six sub-layers) in true relative
 * thickness, circular folds, a villus carpet and the mesenteric vessels
 * entering through the serosa.
 * @param {object} model parsed small_intestine_model.json
 */
export function buildWallSegment (model, opts = {}) {
  const mag = opts.mag || MAG.wall;
  const g = new THREE.Group();
  const lumenR = (opts.lumenRadiusMm || 10.0);
  const len = opts.lengthMm || 46.0;
  const theta0 = -Math.PI * 0.62, thetaLen = Math.PI * 1.24;    // open on the viewer side
  const layers = model.wall.layers;
  const stats = { layers: [], magnification: mag };

  // build outward from the lumen
  let r = lumenR;
  const subs = [
    ['mucosa', layers[0].thickness_mm * 0.62, LAYER_COLOR.mucosa, '黏膜上皮'],
    ['lamina_propria', layers[0].thickness_mm * 0.30, LAYER_COLOR.lamina_propria, '固有层'],
    ['muscularis_mucosae', layers[0].thickness_mm * 0.08, LAYER_COLOR.muscularis_mucosae, '黏膜肌层'],
    ['submucosa', layers[1].thickness_mm, LAYER_COLOR.submucosa, '黏膜下层（Meissner 丛）'],
    ['muscularis_circular', layers[2].thickness_mm, LAYER_COLOR.muscularis_circular, '内环行肌'],
    ['myenteric_plexus', layers[3].thickness_mm, LAYER_COLOR.myenteric_plexus, 'Auerbach 丛 + Cajal 细胞'],
    ['muscularis_longitudinal', layers[4].thickness_mm, LAYER_COLOR.muscularis_longitudinal, '外纵行肌'],
    ['serosa', layers[5].thickness_mm, LAYER_COLOR.serosa, '浆膜']
  ];
  for (const [name, th, col, cn] of subs) {
    const m = new THREE.Mesh(annulusSector(r * mag, (r + th) * mag, len * mag, theta0, thetaLen),
                             std(col, name.includes('muscular') ? 0.52 : 0.66));
    m.name = 'layer_' + name;
    m.userData = { label: cn, thickness_mm: th };
    g.add(m);
    stats.layers.push({ name, cn, thickness_mm: +th.toFixed(3),
                        inner_radius_mm: +r.toFixed(3) });
    r += th;
  }
  stats.wall_thickness_mm = +(r - lumenR).toFixed(3);
  stats.outer_radius_mm = +r.toFixed(3);

  // plicae circulares: inward-projecting rings on the mucosal surface
  const pitch = model.motility ? 7.0 : 7.0;
  const nFold = Math.floor(len / pitch);
  const foldGroup = new THREE.Group(); foldGroup.name = 'plicae';
  for (let i = 0; i <= nFold; i++) {
    const z = (-len / 2 + i * pitch) * mag;
    const h = 0.45 * lumenR;        // plicae project 3-10 mm into a 20 mm lumen
    const t = new THREE.Mesh(
      new THREE.TorusGeometry((lumenR - h * 0.5) * mag, (h * 0.38) * mag, 10, 56,
                              thetaLen).rotateZ(theta0),
      std(LAYER_COLOR.mucosa, 0.70));
    t.position.z = z;
    foldGroup.add(t);
  }
  g.add(foldGroup);
  stats.plicae = { count: nFold + 1, pitch_mm: pitch,
                   projection_mm: +(0.45 * lumenR).toFixed(2),
                   published_reference: 'plicae circulares 3-10 mm apart, ' +
                     'projecting 3-10 mm into the lumen, tallest in the distal ' +
                     'duodenum and proximal jejunum, absent in the distal ileum' };

  // villus carpet on the mucosal surface (instanced cones, true height)
  const vH = 0.9, vR = 0.11;
  const vg = new THREE.ConeGeometry(vR * mag, vH * mag, 7, 1);
  vg.translate(0, -vH * mag / 2, 0);
  const vm = std(0xe88a80, 0.74);
  const nAround = 46, nAlong = Math.floor(len / (vR * 3.2));
  const inst = new THREE.InstancedMesh(vg, vm, nAround * nAlong);
  const dummy = new THREE.Object3D(); let k = 0;
  for (let i = 0; i < nAlong; i++) {
    for (let j = 0; j < nAround; j++) {
      const th = theta0 + thetaLen * ((j + 0.5) / nAround);
      const z = (-len / 2 + (i + 0.5) * len / nAlong) * mag;
      dummy.position.set(Math.cos(th) * lumenR * mag, Math.sin(th) * lumenR * mag, z);
      dummy.lookAt(0, 0, z);
      dummy.rotateX(Math.PI / 2);
      dummy.scale.setScalar(0.78 + 0.44 * ((i * 7 + j * 13) % 11) / 11);
      dummy.updateMatrix();
      inst.setMatrixAt(k++, dummy.matrix);
    }
  }
  inst.count = k; inst.name = 'villi';
  g.add(inst);
  stats.villi = { instances: k, height_mm: vH, radius_mm: vR,
                  published_height_range_mm: [0.5, 1.6] };

  // mesenteric vessels crossing the serosa into the submucosal plexus
  const vess = new THREE.Group(); vess.name = 'wall_vessels';
  for (let i = 0; i < 7; i++) {
    const th = theta0 + thetaLen * (0.08 + 0.84 * i / 6);
    const z = (-len / 2 + len * ((i % 3) + 0.5) / 3) * mag;
    const a = new THREE.Vector3(Math.cos(th) * r * mag * 1.13, Math.sin(th) * r * mag * 1.13, z);
    const b = new THREE.Vector3(Math.cos(th) * (lumenR + subs[0][1] + subs[1][1]) * mag,
                                Math.sin(th) * (lumenR + subs[0][1] + subs[1][1]) * mag, z);
    const curve = new THREE.CatmullRomCurve3([a, a.clone().lerp(b, 0.45)
      .add(new THREE.Vector3(0, 0, 6 * mag * 0.08)), b]);
    vess.add(new THREE.Mesh(new THREE.TubeGeometry(curve, 10, 0.1 * mag, 7, false),
                            std(i % 2 ? 0x2f6fc0 : 0xc8322e, 0.4)));
  }
  g.add(vess);

  g.scale.setScalar(1);
  return { group: g, stats };
}

/**
 * One villus with its flanking crypts, at x60: central lacteal, capillary
 * network, the counter-current arteriole/venule pair, the enterocyte sheet and
 * the stem-cell compartment at the crypt base.
 */
export function buildVillusUnit (model, opts = {}) {
  const mag = opts.mag || MAG.villus;
  const g = new THREE.Group();
  const m = model.micro;
  const vH = 1.05;                    // mm, inside the published 0.5-1.6
  const vR = 0.135;
  const cD = 0.25;                    // crypt depth, published 0.1-0.4
  const encH = m.enterocyte_height_um / 1000.0;
  const stats = { magnification: mag, villus_height_mm: vH, villus_radius_mm: vR,
                  crypt_depth_mm: cD, enterocyte_height_mm: encH,
                  published_villus_height_mm: m.villus_height_mm,
                  published_crypt_depth_mm: m.crypt_depth_mm };

  // lamina propria core of the villus (a capsule) + the epithelial sheet
  const coreCurve = new THREE.CatmullRomCurve3([
    new THREE.Vector3(0, -cD * mag, 0), new THREE.Vector3(0, 0, 0),
    new THREE.Vector3(0, vH * 0.55 * mag, 0), new THREE.Vector3(0, vH * mag, 0)]);
  const core = new THREE.Mesh(
    new THREE.TubeGeometry(coreCurve, 26, (vR - encH) * mag, 20, false),
    std(0xc96a76, 0.66));
  core.name = 'lamina_propria'; g.add(core);
  core.material.transparent = true; core.material.opacity = 0.55;
  core.material.depthWrite = false;
  const tip = new THREE.Mesh(new THREE.SphereGeometry((vR - encH) * mag, 20, 14),
                             std(0xc96a76, 0.66, { transparent: true, opacity: 0.55,
                                                   depthWrite: false }));
  tip.position.y = vH * mag; g.add(tip);

  // epithelium: a translucent shell so the core stays visible, plus an
  // instanced columnar-cell band to make the 25 um enterocyte height readable
  const epi = new THREE.Mesh(new THREE.TubeGeometry(coreCurve, 26, vR * mag, 24, false),
    std(0xf0a89c, 0.52, { transparent: true, opacity: 0.20, depthWrite: false }));
  epi.name = 'epithelium'; g.add(epi);

  // Enterocytes are COLUMNAR: the published 25 um is their height measured
  // along the radius (apical brush border inward, basal pole on the basement
  // membrane), not a stack dimension. Each cell is therefore a box 25 um deep
  // radially and ~8 um across, tiled over the villus surface.
  const encW = 0.008;                       // ~8 um cell width
  const cellG = new THREE.BoxGeometry(encW * 0.86 * mag, encW * 0.86 * mag, encH * mag);
  cellG.translate(0, 0, encH * mag * 0.5);  // grow inward from the surface
  const cellM = std(0xffd0c0, 0.5);
  const circ = 2 * Math.PI * (vR - encH * 0.5);
  const nRow = Math.max(12, Math.round(circ / encW));
  const nCol = Math.max(4, Math.round((vH + cD) / encW));
  // A declared 100-degree sector of the epithelial sheet is omitted, the way an
  // anatomical plate cuts a villus open, so the lacteal, the capillary net and
  // the counter-current pair are visible. The omission is cosmetic and declared;
  // no dimension changes.
  const CUT0 = -Math.PI * 0.28, CUT1 = Math.PI * 0.28;
  const cells = new THREE.InstancedMesh(cellG, cellM, nRow * nCol);
  const dm = new THREE.Object3D(); let kc = 0;
  for (let i = 0; i < nCol; i++) {
    const y = (-cD + (i + 0.5) * (vH + cD) / nCol) * mag;
    for (let j = 0; j < nRow; j++) {
      const th = ((j + (i % 2) * 0.5) / nRow) * Math.PI * 2;
      const thw = th > Math.PI ? th - 2 * Math.PI : th;
      if (thw > CUT0 && thw < CUT1) continue;
      dm.position.set(Math.cos(th) * (vR - encH) * mag, y,
                      Math.sin(th) * (vR - encH) * mag);
      dm.lookAt(0, y, 0); dm.updateMatrix(); cells.setMatrixAt(kc++, dm.matrix);
    }
  }
  cells.count = kc; cells.name = 'enterocytes'; g.add(cells);
  stats.enterocyte_instances = kc;
  stats.enterocyte_width_mm = encW;
  stats.enterocyte_orientation = 'columnar: the 25 um dimension is radial '
    + '(apical brush border inward, basal pole on the basement membrane)';
  stats.epithelial_cutaway_deg = 100.8;
  stats.epithelial_cutaway_note = 'a declared sector of the epithelial sheet is '
    + 'omitted so the villus core is visible; cosmetic, no dimension changes';

  // central lacteal: blind-ended, wider than a capillary, ends below the tip
  const lac = new THREE.Mesh(
    new THREE.TubeGeometry(new THREE.CatmullRomCurve3([
      new THREE.Vector3(0, -cD * 0.4 * mag, 0),
      new THREE.Vector3(0, vH * 0.88 * mag, 0)]), 16, 0.030 * mag, 12, false),
    std(0xd9c66a, 0.35));
  lac.name = 'lacteal'; g.add(lac);
  const lacTip = new THREE.Mesh(new THREE.SphereGeometry(0.030 * mag, 14, 10),
                                std(0xd9c66a, 0.35));
  lacTip.position.y = vH * 0.88 * mag; lacTip.name = 'lacteal_tip'; g.add(lacTip);

  // counter-current pair: one arteriole up the core axis-adjacent, venules down
  const art = new THREE.Mesh(new THREE.TubeGeometry(new THREE.CatmullRomCurve3([
    new THREE.Vector3(vR * 0.42 * mag, -cD * mag, 0),
    new THREE.Vector3(vR * 0.46 * mag, vH * 0.5 * mag, 0),
    new THREE.Vector3(vR * 0.30 * mag, vH * 0.85 * mag, 0)]), 20, 0.017 * mag, 9, false),
    std(0xc8322e, 0.35));
  art.name = 'arteriole'; g.add(art);
  const ven = new THREE.Mesh(new THREE.TubeGeometry(new THREE.CatmullRomCurve3([
    new THREE.Vector3(-vR * 0.42 * mag, -cD * mag, 0),
    new THREE.Vector3(-vR * 0.48 * mag, vH * 0.5 * mag, 0),
    new THREE.Vector3(-vR * 0.32 * mag, vH * 0.85 * mag, 0)]), 20, 0.021 * mag, 9, false),
    std(0x2f6fc0, 0.35));
  ven.name = 'venule'; g.add(ven);

  // subepithelial capillary network: rings linked to the axial pair
  const caps = new THREE.Group(); caps.name = 'capillary_network';
  for (let i = 1; i <= 9; i++) {
    const y = (-cD * 0.2 + (vH + cD * 0.2) * i / 10) * mag;
    const rr = (vR - encH * 1.25) * mag * (1 - 0.22 * Math.pow(i / 10, 3));
    const t = new THREE.Mesh(new THREE.TorusGeometry(rr, 0.011 * mag, 7, 28),
                             std(0x9c4d7a, 0.4));
    t.rotation.x = Math.PI / 2; t.position.y = y; caps.add(t);
  }
  g.add(caps);

  // crypts of Lieberkuhn either side, with the stem-cell compartment at the base
  const crypts = new THREE.Group(); crypts.name = 'crypts';
  for (const sx of [-1, 1]) {
    for (const sz of [-1, 1]) {
      const x = sx * vR * 2.05 * mag, z = sz * vR * 1.55 * mag;
      const c = new THREE.Mesh(new THREE.CylinderGeometry(
        0.055 * mag, 0.075 * mag, cD * mag, 14, 1, true), std(0xb9606c, 0.6));
      c.position.set(x, -cD * mag / 2, z); crypts.add(c);
      const base = new THREE.Mesh(new THREE.SphereGeometry(0.075 * mag, 14, 10),
                                  std(0x7fd0a8, 0.45, { emissive: 0x1d5e40 }));
      base.position.set(x, -cD * mag, z); base.name = 'stem_niche'; crypts.add(base);
      const pan = new THREE.Mesh(new THREE.SphereGeometry(0.030 * mag, 10, 8),
                                 std(0xf2e07a, 0.5));
      pan.position.set(x + 0.05 * mag, -cD * 0.93 * mag, z); pan.name = 'paneth';
      crypts.add(pan);
    }
  }
  g.add(crypts);
  stats.crypts = 4;

  // migrating-cell markers: crypt base -> villus tip over the published 3-5 d
  const migG = new THREE.SphereGeometry(encH * 0.55 * mag, 10, 8);
  const migM = std(0xaef0cf, 0.4, { emissive: 0x2a7f5b });
  const migr = new THREE.InstancedMesh(migG, migM, 18);
  migr.name = 'migrating_cells'; g.add(migr);
  stats.turnover_days = m.epithelial_turnover_days;

  const dummy2 = new THREE.Object3D();
  const animate = (tSec, daysPerSec) => {
    for (let i = 0; i < 18; i++) {
      const u = ((tSec * daysPerSec / 4.0) + i / 18) % 1.0;      // 4 d transit
      const th = (i / 18) * Math.PI * 2;
      const y = (-cD + u * (vH + cD)) * mag;
      const rr = (vR - encH * 0.5) * mag;
      dummy2.position.set(Math.cos(th) * rr, y, Math.sin(th) * rr);
      dummy2.scale.setScalar(u > 0.96 ? Math.max(0.05, (1 - u) / 0.04) : 1);
      dummy2.updateMatrix(); migr.setMatrixAt(i, dummy2.matrix);
    }
    migr.instanceMatrix.needsUpdate = true;
  };

  return { group: g, stats, animate };
}

/** Brush-border inset: microvilli at true 1 um length, x2000. */
export function buildBrushBorder (model, opts = {}) {
  const mag = opts.mag || MAG.brush;
  const g = new THREE.Group();
  const Lmv = model.micro.microvillus_length_um / 1000.0;   // mm
  const Dmv = model.micro.microvillus_diameter_um / 1000.0;
  // The apical patch spans exactly the microvilli drawn on it, so the 1 um
  // microvillus is a readable fraction of the plate instead of a sliver.
  const n = 13;
  const pitch = Dmv * 1.6;                 // ~0.16 um centre-to-centre
  const span = n * pitch;                  // ~2.08 um of brush border
  const plate = new THREE.Mesh(
    new THREE.BoxGeometry(span * mag, Lmv * 0.22 * mag, span * mag), std(0xffd0c0, 0.5));
  plate.position.y = -Lmv * 0.11 * mag; g.add(plate);
  const mv = new THREE.CylinderGeometry(Dmv * 0.5 * mag, Dmv * 0.5 * mag, Lmv * mag, 7);
  mv.translate(0, Lmv * mag / 2, 0);
  const inst = new THREE.InstancedMesh(mv, std(0xffe2d2, 0.42), n * n);
  const d = new THREE.Object3D(); let k = 0;
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      d.position.set((i - (n - 1) / 2) * pitch * mag, 0,
                     (j - (n - 1) / 2) * pitch * mag);
      d.updateMatrix(); inst.setMatrixAt(k++, d.matrix);
    }
  }
  inst.count = k; g.add(inst);
  const gly = new THREE.Mesh(new THREE.BoxGeometry(span * mag, Lmv * 0.35 * mag, span * mag),
    std(0x8fd0e8, 0.3, { transparent: true, opacity: 0.22, depthWrite: false }));
  gly.position.y = Lmv * 0.92 * mag; g.add(gly);
  return { group: g, stats: { magnification: mag, microvilli_rendered: k,
    patch_span_um: +(span * 1000).toFixed(3), microvillus_pitch_um: +(pitch * 1000).toFixed(3),
    microvillus_length_um: model.micro.microvillus_length_um,
    microvillus_diameter_um: model.micro.microvillus_diameter_um,
    microvilli_per_enterocyte_published: model.micro.microvilli_per_enterocyte } };
}
