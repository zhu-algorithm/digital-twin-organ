// crypt.js — the colonic crypt unit, modelled at its real micrometre size and then shown
// at a declared magnification.  The point of this view is the CONTRAST with the small
// intestine: the colonic mucosa is flat.  There are no villi anywhere in the colon.
import * as THREE from 'three';

export function buildCryptUnit(micro, MAG) {
  const g = new THREE.Group();
  const c = micro.crypt;
  const depth = c.depth_um * 1e-3 * MAG;          // mm of scene
  const outer = c.outer_diameter_um * 1e-3 * MAG;
  const lumen = c.lumen_diameter_um * 1e-3 * MAG;
  const epi = c.epithelial_height_um * 1e-3 * MAG;
  const pitch = outer * 1.18;
  const NX = 5, NZ = 5;

  const mucosa = new THREE.Mesh(
    new THREE.BoxGeometry(pitch * NX, epi * 0.45, pitch * NZ),
    new THREE.MeshStandardMaterial({ color: 0xd9a898, roughness: 0.85 }));
  mucosa.position.y = epi * 0.22;
  g.add(mucosa);                                  // FLAT surface — the whole point

  const wallMat = new THREE.MeshStandardMaterial({
    color: 0xe8b8a6, roughness: 0.75, side: THREE.DoubleSide, transparent: true, opacity: 0.95 });
  const gobletMat = new THREE.MeshStandardMaterial({ color: 0x7fd8ff, roughness: 0.5,
    emissive: 0x1b4d66, emissiveIntensity: 0.4 });
  const stemMat = new THREE.MeshStandardMaterial({ color: 0xffd166, roughness: 0.5,
    emissive: 0x6b5112, emissiveIntensity: 0.5 });
  const eeMat = new THREE.MeshStandardMaterial({ color: 0xc792ea, roughness: 0.5 });
  const dcsMat = new THREE.MeshStandardMaterial({ color: 0x8fd39a, roughness: 0.5 });

  const cells = { colonocyte: 0, goblet: 0, stem: 0, enteroendocrine: 0, deep_crypt_secretory: 0 };
  const along = 18;                                // rendered cells along the crypt axis
  const around = 10;                               // rendered cells per circumference
  for (let ix = 0; ix < NX; ix++) {
    for (let iz = 0; iz < NZ; iz++) {
      const cx = (ix - (NX - 1) / 2) * pitch, cz = (iz - (NZ - 1) / 2) * pitch;
      const tube = new THREE.Mesh(
        new THREE.CylinderGeometry(outer / 2, outer / 2 * 0.86, depth, 18, 1, true), wallMat);
      tube.position.set(cx, -depth / 2, cz);
      g.add(tube);
      const base = new THREE.Mesh(new THREE.SphereGeometry(outer / 2 * 0.86, 16, 10,
        0, Math.PI * 2, Math.PI / 2, Math.PI / 2), wallMat);
      base.position.set(cx, -depth, cz);
      g.add(base);
      if (ix === Math.floor(NX / 2) && iz === Math.floor(NZ / 2)) continue;  // cut-away crypt
      for (let a = 0; a < along; a++) {
        const f = (a + 0.5) / along;                     // 0 at the mouth, 1 at the base
        const r = (outer / 2) * (1 - 0.14 * f) - epi * 0.18;
        for (let k = 0; k < around; k++) {
          const th = 2 * Math.PI * (k + (a % 2) * 0.5) / around;
          let mat = wallMat, kind = 'colonocyte', sc = 1.0;
          const hash = (a * 31 + k * 17 + ix * 7 + iz * 13) % 100;
          if (f > 0.86 && hash < 45) { mat = stemMat; kind = 'stem'; }
          else if (f > 0.78 && hash >= 45 && hash < 62) { mat = dcsMat; kind = 'deep_crypt_secretory'; }
          else if (hash < 25) { mat = gobletMat; kind = 'goblet'; sc = 1.25; }
          else if (hash < 28) { mat = eeMat; kind = 'enteroendocrine'; }
          else continue;                                   // plain colonocytes stay as wall
          cells[kind]++;
          const m = new THREE.Mesh(new THREE.SphereGeometry(epi * 0.19 * sc, 7, 5), mat);
          m.position.set(cx + Math.cos(th) * r, -depth * f, cz + Math.sin(th) * r);
          g.add(m);
        }
      }
    }
  }
  g.userData.cells = cells;
  g.userData.scale_note = `每个隐窝按实测 ${c.depth_um} µm 深、${c.outer_diameter_um} µm 外径建模，再统一放大 ${MAG}×`;
  return g;
}

export function buildVillusGhost(MAG) {
  // a ghosted small-intestinal villus at the SAME magnification, for scale comparison
  const g = new THREE.Group();
  const h = 0.8 * MAG;                        // 0.8 mm villus, mid published range
  const mat = new THREE.MeshStandardMaterial({ color: 0x5f8dd3, transparent: true,
    opacity: 0.22, roughness: 0.9, side: THREE.DoubleSide });
  for (let k = 0; k < 3; k++) {
    const v = new THREE.Mesh(new THREE.CapsuleGeometry(0.12 * MAG, h * 0.7, 4, 10), mat);
    v.position.set(2.6 * MAG + k * 0.42 * MAG, h / 2, -1.1 * MAG);
    g.add(v);
  }
  return g;
}
