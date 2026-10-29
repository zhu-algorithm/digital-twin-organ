// gland.js — the gastric gland column and the mucus-bicarbonate barrier, built
// at TRUE published dimensions in micrometres and then scaled by the DECLARED
// magnification.  Nothing here is at organ scale and the HUD says so.
import * as THREE from 'three';
import { phColour } from './materials.js';

export function buildGland(micro, variant = 'oxyntic') {
  const g = new THREE.Group();
  g.name = 'gland';
  const M = micro.magnification;            // mm of scene per mm of tissue
  const um = M / 1000.0;                    // mm of scene per micrometre
  const spec = variant === 'pyloric' ? micro.pyloric_gland : micro.oxyntic_gland;
  const depth = spec.true_total_depth_um;
  const fov = spec.true_foveolar_depth_um;
  const bore = micro.gland_bore_um;
  const spacing = micro.gland_spacing_um;
  const mucus = micro.mucus_barrier.true_thickness_um;

  // --- the mucus-bicarbonate barrier, with the pH field painted on it -------
  // Drawn as a stack of thin slabs so the gradient is a real spatial field
  // rather than a texture gradient; the slab count is the published profile's
  // sample count.
  const zs = micro.mucus_barrier.z_um;
  const phs = micro.mucus_barrier.pH_at_luminal_1_5;
  const barrier = new THREE.Group();
  barrier.name = 'mucus_barrier';
  const W = spacing * 7;
  for (let i = 0; i < zs.length - 1; i++) {
    const h = (zs[i + 1] - zs[i]) * um;
    const geo = new THREE.BoxGeometry(W * um, h, W * um);
    const mat = new THREE.MeshStandardMaterial({
      color: phColour(phs[i]), transparent: true, opacity: 0.5,
      roughness: 0.85, metalness: 0.0, depthWrite: false,
    });
    const m = new THREE.Mesh(geo, mat);
    // z = 0 is the luminal face; the epithelium is at z = thickness
    m.position.y = (mucus - 0.5 * (zs[i] + zs[i + 1])) * um;
    m.name = `mucus_slab_${i}`;
    barrier.add(m);
  }
  g.add(barrier);

  // --- the gland columns ----------------------------------------------------
  const glands = new THREE.Group();
  glands.name = 'glands';
  const nx = 5, nz = 5;
  for (let ix = 0; ix < nx; ix++) {
    for (let iz = 0; iz < nz; iz++) {
      const x = (ix - (nx - 1) / 2) * spacing * um;
      const z = (iz - (nz - 1) / 2) * spacing * um;
      // the pit (foveola) and the gland proper are different compartments and
      // different cell populations, so they are separate meshes
      const pit = new THREE.Mesh(
        new THREE.CylinderGeometry(bore * 0.62 * um, bore * 0.5 * um, fov * um, 12, 1, true),
        new THREE.MeshStandardMaterial({ color: 0xd9a89a, roughness: 0.7,
          side: THREE.DoubleSide }));
      pit.position.set(x, -fov * um * 0.5, z);
      pit.name = 'foveola';
      const body = new THREE.Mesh(
        new THREE.CylinderGeometry(bore * 0.5 * um, bore * 0.42 * um,
          (depth - fov) * um, 12, 1, true),
        new THREE.MeshStandardMaterial({ color: 0xc98f86, roughness: 0.72,
          side: THREE.DoubleSide }));
      body.position.set(x, -(fov + (depth - fov) * 0.5) * um, z);
      body.name = 'gland_body';
      glands.add(pit, body);

      // cells, placed on the published axial order.  Sizes are indicative and
      // the page says so; the AXIAL POSITION is the point being made.
      const cells = variant === 'pyloric'
        ? [[0.04, 0xe8cdbe, 'surface mucous'], [0.30, 0xd9b39c, 'mucous neck'],
           [0.72, 0x7fc1a8, 'G cell'], [0.86, 0x9b8fd0, 'D cell']]
        : [[0.04, 0xe8cdbe, 'surface mucous'], [0.26, 0xd9b39c, 'mucous neck'],
           [0.42, 0xf0a33c, 'parietal'], [0.58, 0xf0a33c, 'parietal'],
           [0.82, 0x6aa0e0, 'chief'], [0.93, 0xd05f9b, 'ECL']];
      for (const [f, col, nm] of cells) {
        const r = bore * 0.30 * um;
        const s = new THREE.Mesh(new THREE.SphereGeometry(r, 10, 8),
          new THREE.MeshStandardMaterial({ color: col, roughness: 0.55,
            emissive: new THREE.Color(col).multiplyScalar(0.18) }));
        s.position.set(x + bore * 0.42 * um, -depth * um * f, z);
        s.name = `cell_${nm.replace(/\s+/g, '_')}`;
        s.userData.cellType = nm;
        glands.add(s);
      }
    }
  }
  g.add(glands);

  // --- the proton pump, on the APICAL membrane of the parietal cell ---------
  // Drawing it on the wrong face would be worse than not drawing it, so its
  // position is taken from the data: APICAL = facing the gland lumen.
  const pumps = new THREE.Group();
  pumps.name = 'proton_pumps';
  const pumpMat = new THREE.MeshStandardMaterial({ color: 0xffd24a,
    emissive: 0x6a4e00, roughness: 0.35 });
  for (let k = 0; k < 24; k++) {
    const ang = (k / 24) * Math.PI * 2;
    const f = 0.42 + 0.16 * ((k % 3) / 3);
    const r = bore * 0.5 * um;
    const p = new THREE.Mesh(new THREE.SphereGeometry(bore * 0.09 * um, 8, 6), pumpMat);
    p.position.set(Math.cos(ang) * r, -depth * um * f, Math.sin(ang) * r);
    p.name = 'ATP4A_ATP4B';
    p.userData.membrane = 'apical';
    pumps.add(p);
  }
  g.add(pumps);

  g.userData.scale = {
    magnification: M,
    true_total_depth_um: depth,
    rendered_depth_mm: depth * um,
    mucus_true_um: mucus,
    mucus_rendered_mm: mucus * um,
    variant,
  };
  return g;
}
