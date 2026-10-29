/**
 * deform.js — physiological myocardial contraction kinematics on the GPU.
 *
 * The ventricles are NOT uniformly scaled. The deformation is the composition
 * of four components, each with an independent, calibrated amplitude:
 *
 *   1. Longitudinal shortening with a fixed apex and a descending base.
 *      Axial strain is non-uniform: basal and mid segments shorten more than
 *      apical ones (local strain weight f(zeta) = a + b*zeta, normalised so the
 *      mean over the long axis equals 1). The apex stays at zeta = 0; the base
 *      descends by Lbase * eps_long * mean(f) towards it.
 *   2. Wall thickening: the endocardium moves inward further than the
 *      epicardium, so the transmural gap grows. Per-vertex interpolation uses
 *      the baked transmural coordinate (0 = endocardium, 1 = epicardium).
 *   3. Radial (circumferential) shortening about each chamber's own long axis.
 *   4. Torsion: apex rotates counter-clockwise viewed from the apex, base
 *      clockwise, with a linear gradient along the long axis.
 *
 * Amplitudes are driven, every frame, by the chamber volume time series of the
 * `heart_timeseries.json` contract. The mapping volume -> deformation is a
 * one-off calibration performed against the real mesh (see
 * tools/calibrate_deformation.py): for each chamber the endocardial radial
 * scale is solved numerically so that the enclosed cavity volume of the actual
 * surface matches the contract volume, under the imposed longitudinal
 * shortening. The calibration is shipped as cubic polynomials in the
 * contraction coordinate c, so any heart rate and any stroke volume drive the
 * geometry correctly without re-calibration.
 *
 * Chamber ids (baked into vertex attribute _field1.x):
 *   0 LV   1 RV   2 LA   3 RA   4 arterial great vessels   5 venous vessels
 */

import * as THREE from 'three';

export const NCHAMBER = 6;
export const CH = { LV: 0, RV: 1, LA: 2, RA: 3, ART: 4, VEN: 5 };

/** GLSL injected into every deformable shader (myocardium, coronaries, particles). */
export const DEFORM_GLSL = /* glsl */`
uniform mat4 uWorldToHeart;
uniform mat4 uHeartToWorld;
uniform float uLbase;
uniform float uZtop;
uniform float uStrainA;
uniform float uStrainB;
uniform float uEpsLong[${NCHAMBER}];
uniform float uScaleEndo[${NCHAMBER}];
uniform float uScaleEpi[${NCHAMBER}];
uniform float uTwist[${NCHAMBER}];
uniform vec2  uAxisXY[${NCHAMBER}];

// Integrated axial map. The local axial strain weight is f(u) = a + b*u,
// normalised so a + b/2 == 1, i.e. mean strain over the long axis is exactly
// eps. Integrating 1 - eps*f gives
//     z'/Lbase = zeta - eps*(a*zeta + b*zeta^2/2)
// which at zeta = 1 is 1 - eps: total long-axis shortening equals eps, as the
// calibration requires. (Integrating a*u + b*u^2 instead would deliver only
// 0.558*eps and silently under-contract the ventricle.)
float axialMap (float zeta, float eps) {
  return zeta - eps * (uStrainA * zeta + uStrainB * zeta * zeta * 0.5);
}

// Local axial stretch dz'/dz = 1 - eps*f(zeta), needed to transform normals.
float axialStretch (float zeta, float eps) {
  return 1.0 - eps * (uStrainA + uStrainB * zeta);
}

/**
 * Deform a heart-space point.
 *  pH          heart-space position
 *  transmural  0 = endocardium, 1 = epicardium
 *  ci          chamber id
 * Returns the deformed heart-space position; writes the local diagonal
 * Jacobian (radial, radial, axial) and the twist angle into outJ / outPhi so
 * callers can transform normals and tangents.
 */
vec3 deformHeart (vec3 pH, float transmural, int ci, out vec3 outJ, out float outPhi) {
  float eps  = uEpsLong[ci];
  float zetaRaw = pH.z / uLbase;
  float zeta = clamp(zetaRaw, 0.0, 1.0);

  // Below the valve plane the ventricle shortens with a fixed apex. Above it the
  // atria and great-vessel roots are tethered superiorly, so the descending
  // valve plane stretches them: this is the atrioventricular-plane-displacement
  // mechanism by which atrial filling is driven during ventricular systole.
  float zNew;
  float axialJ;
  if (zetaRaw <= 1.0) {
    zNew = uLbase * axialMap(zeta, eps);
    axialJ = axialStretch(zeta, eps);
  } else {
    float zBase = uLbase * axialMap(1.0, eps);
    float k = (uZtop - zBase) / max(uZtop - uLbase, 1e-3);
    zNew = zBase + (pH.z - uLbase) * k;
    axialJ = k;
  }

  float s = mix(uScaleEndo[ci], uScaleEpi[ci], clamp(transmural, 0.0, 1.0));

  float phi = uTwist[ci] * (1.0 - 2.0 * zeta);
  float cs = cos(phi), sn = sin(phi);

  vec2 axis = uAxisXY[ci];
  vec2 rel  = pH.xy - axis;
  vec2 relR = vec2(rel.x * cs - rel.y * sn, rel.x * sn + rel.y * cs);

  outJ   = vec3(s, s, axialJ);
  outPhi = phi;
  return vec3(axis + relR * s, zNew);
}

/** Convenience wrapper: world-space position in, world-space position out. */
vec3 deformWorld (vec3 pW, float transmural, int ci) {
  vec3 pH = (uWorldToHeart * vec4(pW, 1.0)).xyz;
  vec3 J; float phi;
  vec3 dH = deformHeart(pH, transmural, ci, J, phi);
  return (uHeartToWorld * vec4(dH, 1.0)).xyz;
}

/** Transform a heart-space direction by the inverse transpose of the local Jacobian. */
vec3 deformNormalHeart (vec3 nH, vec3 J, float phi) {
  vec3 n = vec3(nH.x / max(J.x, 1e-5), nH.y / max(J.y, 1e-5), nH.z / max(J.z, 1e-5));
  float cs = cos(phi), sn = sin(phi);
  return normalize(vec3(n.x * cs - n.y * sn, n.x * sn + n.y * cs, n.z));
}
`;

/**
 * Holds the shared uniform objects. Every deformable material references the
 * SAME uniform objects, so one update per frame moves myocardium, coronary tree
 * and blood particles coherently.
 */
export class DeformField {
  /** @param {object} calib parsed deform_calibration.json */
  constructor (calib) {
    this.calib = calib;
    const m4 = (arr) => new THREE.Matrix4().fromArray(arr);
    this.uniforms = {
      uWorldToHeart: { value: m4(calib.world_to_heart) },
      uHeartToWorld: { value: m4(calib.heart_to_world) },
      uLbase: { value: calib.L_base_mm },
      uZtop: { value: calib.Z_top_mm },
      uStrainA: { value: calib.strain_shape.a },
      uStrainB: { value: calib.strain_shape.b },
      uEpsLong: { value: new Array(NCHAMBER).fill(0) },
      uScaleEndo: { value: new Array(NCHAMBER).fill(1) },
      uScaleEpi: { value: new Array(NCHAMBER).fill(1) },
      uTwist: { value: new Array(NCHAMBER).fill(0) },
      uAxisXY: { value: calib.axis_xy.map(a => new THREE.Vector2(a[0], a[1])) }
    };
    this._poly = calib.chambers; // per-chamber polynomial coefficients
    this.state = { cLV: 0, cRV: 0, cLA: 0, cRA: 0 };
  }

  static _poly3 (co, c) {
    return co[0] + c * (co[1] + c * (co[2] + c * co[3]));
  }

  /**
   * Update all uniforms from the contract at wall-clock time tMs.
   * @param {HeartTimeseries} ts
   */
  update (ts, tMs) {
    const volKey = ['V_LV_mL', 'V_RV_mL', 'V_LA_mL', 'V_RA_mL'];
    for (let ci = 0; ci < 4; ci++) {
      const c = Math.min(Math.max(ts.contraction(volKey[ci], tMs), 0), 1);
      const p = this._poly[ci];
      this.uniforms.uScaleEndo.value[ci] = DeformField._poly3(p.scale_endo, c);
      this.uniforms.uScaleEpi.value[ci] = DeformField._poly3(p.scale_epi, c);
      this.uniforms.uEpsLong.value[ci] = p.eps_long_max * c;
      this.uniforms.uTwist.value[ci] = p.twist_max_rad * c;
    }
    this.state.cLV = Math.min(Math.max(ts.contraction('V_LV_mL', tMs), 0), 1);
    this.state.cRV = Math.min(Math.max(ts.contraction('V_RV_mL', tMs), 0), 1);
    this.state.cLA = Math.min(Math.max(ts.contraction('V_LA_mL', tMs), 0), 1);
    this.state.cRA = Math.min(Math.max(ts.contraction('V_RA_mL', tMs), 0), 1);

    // Great vessels: attached to the valve plane, so they follow basal descent
    // rigidly; their calibre pulses with transmural pressure (aortic root
    // distension on the arterial side, right-atrial pressure on the venous).
    const epsBase = this.uniforms.uEpsLong.value[CH.LV];
    const pAo = ts.sample('P_Ao_mmHg', tMs);
    const dAo = this.calib.arterial_distension;
    const sArt = 1 + dAo.k * (pAo - dAo.p_ref_mmHg);
    this.uniforms.uEpsLong.value[CH.ART] = epsBase;
    this.uniforms.uScaleEndo.value[CH.ART] = sArt;
    this.uniforms.uScaleEpi.value[CH.ART] = sArt;
    this.uniforms.uEpsLong.value[CH.VEN] = epsBase * 0.5;
    this.uniforms.uScaleEndo.value[CH.VEN] = 1.0;
    this.uniforms.uScaleEpi.value[CH.VEN] = 1.0;
  }

  /** Reproduce the GPU deformation on the CPU (used by the self-check probe). */
  deformCPU (pWorld, transmural, ci) {
    const u = this.uniforms;
    const pH = pWorld.clone().applyMatrix4(u.uWorldToHeart.value);
    const eps = u.uEpsLong.value[ci];
    const L = u.uLbase.value;
    const a = u.uStrainA.value; const b = u.uStrainB.value;
    const zetaRaw = pH.z / L;
    const zeta = Math.min(Math.max(zetaRaw, 0), 1);
    const amap = (zt) => zt - eps * (a * zt + b * zt * zt * 0.5);
    let zNew;
    if (zetaRaw <= 1) {
      zNew = L * amap(zeta);
    } else {
      const zBase = L * amap(1);
      const zTop = u.uZtop.value;
      zNew = zBase + (pH.z - L) * ((zTop - zBase) / Math.max(zTop - L, 1e-3));
    }
    const t = Math.min(Math.max(transmural, 0), 1);
    const s = u.uScaleEndo.value[ci] * (1 - t) + u.uScaleEpi.value[ci] * t;
    const phi = u.uTwist.value[ci] * (1 - 2 * zeta);
    const ax = u.uAxisXY.value[ci];
    const rx = pH.x - ax.x; const ry = pH.y - ax.y;
    const cs = Math.cos(phi); const sn = Math.sin(phi);
    const out = new THREE.Vector3(
      ax.x + (rx * cs - ry * sn) * s,
      ax.y + (rx * sn + ry * cs) * s,
      zNew);
    return out.applyMatrix4(u.uHeartToWorld.value);
  }
}
