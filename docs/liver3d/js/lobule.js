/**
 * lobule.js — the classical hepatic lobule in three dimensions.
 *
 * This is the structure that makes the liver the liver, and it is built here at
 * its real proportions and then magnified for display (the magnification factor
 * is reported in the HUD and in the acceptance record — nothing about the
 * geometry is invented to fit the screen):
 *
 *   lobule diameter      1.0 mm  (classical human lobule 0.7-2.0 mm)
 *   lobule height        1.6 mm
 *   central vein radius  0.045 mm
 *   portal venule        0.035 mm, hepatic arteriole 0.016 mm,
 *   bile ductule         0.020 mm   (interlobular calibres)
 *   sinusoid radius      0.0055 mm
 *
 * Six portal tracts sit at the hexagon's corners, each carrying the portal
 * TRIAD: a portal venule, a hepatic arteriole and a bile ductule. The central
 * vein runs up the axis. Sinusoids run from the corners inward.
 *
 * The two directions that matter, and which the 2-D page could only hint at:
 *   blood  portal tract -> sinusoid -> CENTRAL vein   (inward)
 *   bile   hepatocyte canaliculus -> bile ductule     (OUTWARD, counter-current)
 * Bile therefore flows against the blood inside the same plate of hepatocytes.
 *
 * Sinusoidal velocity is not chosen: it is lobule radius / sinusoidal transit
 * time, with the transit time taken from liver_model.json (flagged there as a
 * textbook-magnitude value, not a measurement of this project). Bile velocity
 * follows from bile output per lobule over the canalicular cross-section.
 *
 * Zonal colouring comes from the project's own measured single-cell expression
 * (liver_model.json .zonal_expression): three hepatocyte populations annotated
 * by lobular position, 103 402 cells over 24 donors. The gradient is NOT the
 * same for every gene — CYP1A2 is 2.21x higher pericentrally, CYP2C9 is flat at
 * 1.04x, CYP2D6 runs the OTHER WAY at 0.74x — so the colour ramp is rebound per
 * gene and the ramp direction follows the measured ratio.
 */

import * as THREE from 'three';

export const MAGNIFICATION = 60.0;      // display scale: 1 mm of tissue -> 60 mm
const R_LOBULE_MM = 0.50;
const H_LOBULE_MM = 1.60;
const R_CENTRAL_MM = 0.045;
const R_PORTAL_VENULE_MM = 0.035;
const R_ARTERIOLE_MM = 0.016;
const R_DUCTULE_MM = 0.020;
const R_SINUSOID_MM = 0.0055;
const SINUSOID_RINGS = 5;              // vertical levels of sinusoids
const SINUSOIDS_PER_CORNER = 3;

const COL = {
  portal: 0x4a6fa5,
  artery: 0xb23a3a,
  bile: 0x7d9a3c,
  central: 0x39527d,
  plate: 0x8d4036
};

const S = MAGNIFICATION;

function tube (a, b, r, sides, bend) {
  /** Straight or gently bowed tube from a to b (scene units). */
  const pts = [];
  const n = bend ? 10 : 2;
  for (let i = 0; i < n; i++) {
    const t = i / (n - 1);
    const p = a.clone().lerp(b, t);
    if (bend) {
      const mid = Math.sin(Math.PI * t);
      p.add(bend.clone().multiplyScalar(mid));
    }
    pts.push(p);
  }
  const curve = new THREE.CatmullRomCurve3(pts);
  return new THREE.TubeGeometry(curve, bend ? 14 : 1, r, sides, false);
}

/**
 * @param {object} model parsed liver_model.json
 * @returns {{group: THREE.Group, flowSpec: object, info: object}}
 */
export function buildLobule (model) {
  const group = new THREE.Group();
  const R = R_LOBULE_MM * S;
  const H = H_LOBULE_MM * S;
  const mats = [];

  // ---------------------------------------------------------- parenchyma
  // A 60-degree sector is cut away (one hexagonal sixth, i.e. one classical
  // "portal unit") so the central vein, the sinusoids and one whole portal tract
  // are visible from outside instead of being hidden inside an opaque block.
  const hexShape = new THREE.Shape();
  hexShape.moveTo(0, 0);
  for (let k = 1; k <= 6; k++) {
    const th = (k / 6) * Math.PI * 2 + Math.PI / 6;
    hexShape.lineTo(R * Math.cos(th), R * Math.sin(th));
  }
  hexShape.closePath();
  const prism = new THREE.ExtrudeGeometry(hexShape, { depth: H, bevelEnabled: false });
  prism.rotateX(-Math.PI / 2);
  prism.translate(0, 0, 0);
  const plateMat = makePlateMaterial(R, H);
  mats.push(plateMat);
  const plate = new THREE.Mesh(prism, plateMat);
  plate.name = 'lobule_parenchyma';
  group.add(plate);

  // ------------------------------------------------------- central vein
  const cvMat = new THREE.MeshPhysicalMaterial({
    color: COL.central, roughness: 0.3, clearcoat: 0.6, side: THREE.DoubleSide
  });
  const cv = new THREE.Mesh(
    new THREE.CylinderGeometry(R_CENTRAL_MM * S, R_CENTRAL_MM * S * 1.15, H * 1.04, 20),
    cvMat);
  cv.position.set(0, H / 2, 0);
  cv.name = 'central_vein';
  group.add(cv);

  // -------------------------------------------------- six portal triads
  const corners = [];
  const triadMats = {
    portal: new THREE.MeshPhysicalMaterial({ color: COL.portal, roughness: 0.32, clearcoat: 0.5 }),
    artery: new THREE.MeshPhysicalMaterial({ color: COL.artery, roughness: 0.3, clearcoat: 0.6 }),
    bile: new THREE.MeshPhysicalMaterial({ color: COL.bile, roughness: 0.38, clearcoat: 0.4 })
  };
  for (let k = 0; k < 6; k++) {
    const th = (k / 6) * Math.PI * 2 + Math.PI / 6;
    const c = new THREE.Vector3(R * Math.cos(th), 0, R * Math.sin(th));
    corners.push(c);
    const radial = c.clone().normalize();
    const tangent = new THREE.Vector3(-radial.z, 0, radial.x);
    // the three triad members sit side by side in the portal tract
    const members = [
      ['portal', R_PORTAL_VENULE_MM, radial.clone().multiplyScalar(-0.055 * S)],
      ['artery', R_ARTERIOLE_MM, tangent.clone().multiplyScalar(0.05 * S)],
      ['bile', R_DUCTULE_MM, tangent.clone().multiplyScalar(-0.05 * S)]
    ];
    for (const [kind, rmm, off] of members) {
      const m = new THREE.Mesh(
        new THREE.CylinderGeometry(rmm * S, rmm * S, H * 1.04, 14), triadMats[kind]);
      m.position.copy(c).add(off);
      m.position.y = H / 2;
      m.name = `triad_${kind}_${k}`;
      group.add(m);
    }
  }

  // ------------------------------------------------- sinusoids and canaliculi
  const sinMat = new THREE.MeshPhysicalMaterial({
    color: 0x9a5050, roughness: 0.5, transparent: true, opacity: 0.55
  });
  const canMat = new THREE.MeshPhysicalMaterial({
    color: COL.bile, roughness: 0.45, transparent: true, opacity: 0.8
  });
  const sinusoidPaths = [];
  const canaliculusPaths = [];
  for (let ring = 0; ring < SINUSOID_RINGS; ring++) {
    const y = H * (ring + 0.5) / SINUSOID_RINGS;
    for (let k = 0; k < 6; k++) {
      for (let j = 0; j < SINUSOIDS_PER_CORNER; j++) {
        const spread = (j - (SINUSOIDS_PER_CORNER - 1) / 2) * 0.17;
        const th = (k / 6) * Math.PI * 2 + Math.PI / 6 + spread;
        const a = new THREE.Vector3(R * 0.94 * Math.cos(th), y, R * 0.94 * Math.sin(th));
        const b = new THREE.Vector3(0, y, 0);
        const bend = new THREE.Vector3(-Math.sin(th), 0, Math.cos(th))
          .multiplyScalar(R * 0.1 * ((j % 2) ? 1 : -1));
        const g = tube(a, b.clone().setLength(R_CENTRAL_MM * S * 1.1).setY(y), R_SINUSOID_MM * S, 5, bend);
        const mesh = new THREE.Mesh(g, sinMat);
        mesh.name = `sinusoid_${ring}_${k}_${j}`;
        group.add(mesh);
        // sampled centreline for the particle paths (inward: tract -> centre)
        const curvePts = [];
        for (let i = 0; i < 24; i++) {
          const t = i / 23;
          const p = a.clone().lerp(b, t).add(bend.clone().multiplyScalar(Math.sin(Math.PI * t)));
          curvePts.push([p.x, p.y, p.z]);
        }
        sinusoidPaths.push(curvePts);
        if (j === 1) {
          // one canaliculus per corner per ring, running OUTWARD
          canaliculusPaths.push(curvePts.slice().reverse());
          const cg = tube(b.clone().setY(y), a, R_SINUSOID_MM * S * 0.55, 4,
                          bend.clone().multiplyScalar(-0.5));
          const cm = new THREE.Mesh(cg, canMat);
          cm.name = `canaliculus_${ring}_${k}`;
          group.add(cm);
        }
      }
    }
  }

  // ------------------------------------------------------------- velocities
  const hm = model.haemodynamics;
  const transit = (hm.sinusoid_transit_s[0] + hm.sinusoid_transit_s[1]) / 2;
  // real sinusoidal velocity: the path length the blood covers / transit time
  const v_sin_mm_s = R_LOBULE_MM / transit;
  // bile: output per lobule over the canalicular cross-section. Liver mass /
  // lobule volume gives the lobule count; both inputs are flagged as
  // textbook-magnitude in liver_model.json.
  const liverMassG = (model.haemodynamics.liver_mass_g[0] + model.haemodynamics.liver_mass_g[1]) / 2;
  const lobuleVol_mm3 = 2.598 * R_LOBULE_MM * R_LOBULE_MM * H_LOBULE_MM;  // hex prism
  // 1.05 g/cm^3 -> cm^3, times 1000 for mm^3. Dropping the cm^3->mm^3 factor
  // undercounts the lobules by 1000x and makes the canalicular velocity absurd.
  const nLobules = (liverMassG / 1.05) * 1000.0 / lobuleVol_mm3;
  const bile_mm3_s = ((hm.bile_flow_mL_min[0] + hm.bile_flow_mL_min[1]) / 2) * 1000 / 60;
  const bilePerLobule = bile_mm3_s / nLobules;
  const a_can = Math.PI * Math.pow(R_SINUSOID_MM * 0.55, 2) * 6;   // 6 canaliculi
  const v_bile_mm_s = bilePerLobule / a_can;

  const SIN_GAIN = 1.0;   // sinusoid particles run at true velocity
  const BILE_GAIN = v_sin_mm_s / Math.max(v_bile_mm_s, 1e-9) * 0.55;

  const flowSpec = {
    bile_display_gain: BILE_GAIN,
    paths: [
      ...sinusoidPaths.map(p => ({
        tree: 'portal_vein', retrograde: false, points: p,
        v_mm_s: new Array(p.length).fill(v_sin_mm_s * S), display_gain: SIN_GAIN,
        length_mm: R_LOBULE_MM * S
      })),
      ...canaliculusPaths.map(p => ({
        tree: 'bile_duct', retrograde: true, points: p,
        v_mm_s: new Array(p.length).fill(v_bile_mm_s * S), display_gain: BILE_GAIN,
        length_mm: R_LOBULE_MM * S
      }))
    ]
  };

  const info = {
    magnification: MAGNIFICATION,
    lobule_diameter_mm: R_LOBULE_MM * 2,
    lobule_height_mm: H_LOBULE_MM,
    portal_tracts: 6,
    triad_members_per_tract: ['portal venule', 'hepatic arteriole', 'bile ductule'],
    sinusoids_rendered: sinusoidPaths.length,
    canaliculi_rendered: canaliculusPaths.length,
    sinusoid_transit_s: hm.sinusoid_transit_s,
    sinusoid_transit_source: hm.sinusoid_transit_source,
    sinusoid_velocity_mm_s: Number(v_sin_mm_s.toFixed(5)),
    bile_velocity_mm_s: Number(v_bile_mm_s.toExponential(3)),
    bile_display_gain: Number(BILE_GAIN.toFixed(1)),
    estimated_lobules_in_liver: Math.round(nLobules),
    blood_direction: 'portal tract -> sinusoid -> central vein (inward)',
    bile_direction: 'canaliculus -> bile ductule at the portal tract (outward, '
                    + 'counter-current to the blood)'
  };

  group.userData = { info, plateMat, mats };
  return { group, flowSpec, info, plateMat, corners };
}

/**
 * Parenchyma of the lobule: hepatocyte plates radiating from the central vein,
 * tinted by the measured zonal expression of whichever gene is bound.
 */
function makePlateMaterial (R, H) {
  const mat = new THREE.MeshPhysicalMaterial({
    color: new THREE.Color(COL.plate),
    roughness: 0.58, metalness: 0.0, clearcoat: 0.30,
    transparent: true, opacity: 0.74, side: THREE.DoubleSide,
    depthWrite: false
  });
  mat.userData.uniforms = {
    uR: { value: R },
    uH: { value: H },
    uZoneExpr: { value: new THREE.Vector3(1, 1, 1) },   // periportal, mid, central
    uZoneRatio: { value: 1.0 },
    uZoneMix: { value: 1.0 },
    uShowCells: { value: 1.0 }
  };
  mat.onBeforeCompile = (shader) => {
    Object.assign(shader.uniforms, mat.userData.uniforms);
    shader.vertexShader = shader.vertexShader
      .replace('#include <common>', `#include <common>
        varying vec3 vLocal;`)
      .replace('#include <begin_vertex>', `#include <begin_vertex>
        vLocal = position;`);
    shader.fragmentShader = shader.fragmentShader
      .replace('#include <common>', `#include <common>
        varying vec3 vLocal;
        uniform float uR; uniform float uH;
        uniform vec3 uZoneExpr; uniform float uZoneRatio;
        uniform float uZoneMix; uniform float uShowCells;
        float h21 (vec2 p) {
          p = fract(p * vec2(123.34, 456.21)); p += dot(p, p + 45.32);
          return fract(p.x * p.y);
        }
        vec2 h22 (vec2 p) {
          return fract(sin(vec2(dot(p, vec2(127.1, 311.7)), dot(p, vec2(269.5, 183.3)))) * 43758.5453);
        }
        float cells (vec2 p) {
          vec2 n = floor(p), f = fract(p); float f1 = 8.0, f2 = 8.0;
          for (int j = -1; j <= 1; j++) for (int i = -1; i <= 1; i++) {
            vec2 g = vec2(float(i), float(j));
            float d = length(g + h22(n + g) - f);
            if (d < f1) { f2 = f1; f1 = d; } else if (d < f2) { f2 = d; }
          }
          return f2 - f1;
        }`)
      .replace('#include <map_fragment>', `#include <map_fragment>
        float rr = length(vLocal.xz) / uR;                 // 0 centre -> 1 corner
        // zone coordinate: 0 periportal (corner) -> 1 pericentral (axis)
        float z = clamp(1.0 - rr, 0.0, 1.0);
        // hepatocyte plates: radial cords, one cell wide, ~20 um cells
        float ang = atan(vLocal.z, vLocal.x);
        float cord = smoothstep(0.35, 0.75, abs(sin(ang * 46.0)));
        float cellpat = smoothstep(0.03, 0.30, cells(vec2(ang * 22.0, vLocal.y * 0.55)));
        vec3 base = diffuseColor.rgb * (0.80 + 0.30 * cellpat) * (0.88 + 0.18 * cord);
        // measured zonal expression: sample the three annotated zones
        float e = z < 0.5 ? mix(uZoneExpr.x, uZoneExpr.y, z * 2.0)
                          : mix(uZoneExpr.y, uZoneExpr.z, (z - 0.5) * 2.0);
        float emax = max(max(uZoneExpr.x, uZoneExpr.y), uZoneExpr.z);
        float en = emax > 1e-6 ? e / emax : 0.0;
        // Periportal -> pericentral ramp. The endpoints are a cool blue and a
        // warm amber with the MEASURED expression deciding where on the ramp each
        // zone sits, so a gene that is flat (CYP2C9, 1.04x) renders as a nearly
        // uniform band and a reversed gene (CYP2D6, 0.74x) ramps the other way.
        vec3 cold = vec3(0.16, 0.26, 0.50);
        vec3 hot  = vec3(0.95, 0.68, 0.16);
        vec3 zc = mix(cold, hot, clamp(en, 0.0, 1.0));
        base = mix(base, mix(base * 0.55, zc, 0.85), uZoneMix);
        diffuseColor.rgb = base;`);
  };
  mat.customProgramCacheKey = () => 'lobule-plate-v2';
  return mat;
}
