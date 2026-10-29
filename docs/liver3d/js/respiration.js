/**
 * respiration.js — diaphragmatic liver motion.
 *
 * The liver is NOT a contractile organ. Unlike the myocardium it does not
 * deform of its own accord: during quiet breathing it is displaced bodily by the
 * diaphragm, mainly head-to-foot, with a small anterior-posterior rotation as
 * the right hemiliver rides over the dome. The parenchyma itself stays
 * essentially undeformed. This module therefore produces a RIGID-BODY transform
 * (rotation + translation) and nothing else — there is deliberately no
 * per-vertex deformation anywhere in this scene, so the shape you see is the
 * segmented shape at every instant of the cycle.
 *
 * Frequency: 0.18323 Hz (10.99 breaths/min), carried over unchanged from the
 * project's cardiac model, where it was fitted as the respiratory sinus
 * arrhythmia frequency of PhysioNet nsrdb record 16265. That is PUBLIC
 * placeholder data, not the user's own recording (is_user_personal_data = false).
 *
 * Amplitude: craniocaudal excursion 1-3 cm is a textbook range for quiet
 * breathing, NOT a measurement of this subject. The default is the midpoint of
 * that range and the scene reports which value it is using.
 *
 * Waveform: inspiration is shorter than expiration (I:E about 1:1.5 at rest), so
 * a plain sinusoid would be wrong in shape even with the right period. The
 * excursion follows a piecewise-smooth profile with that ratio; the profile is
 * normalised so peak-to-peak displacement equals the configured excursion
 * exactly.
 */

import * as THREE from 'three';

export const IE_RATIO = 1.5;            // expiratory : inspiratory duration

export class RespirationField {
  /**
   * @param {object} spec  the `respiration` block of liver_model.json
   * @param {THREE.Vector3} cranial  unit cranial direction in scene space
   * @param {THREE.Vector3} pivot    rotation centre (the hepatocaval confluence)
   */
  constructor (spec, cranial, pivot, opts = {}) {
    this.fHz = spec.f_hz;
    this.breathsPerMin = spec.breaths_per_min;
    this.source = spec.source;
    this.isPersonal = !!spec.is_user_personal_data;
    this.excursionRange = spec.craniocaudal_excursion_cm || [1, 3];
    this.excursionCm = opts.excursionCm !== undefined
      ? opts.excursionCm
      : (this.excursionRange[0] + this.excursionRange[1]) * 0.5;
    this.rotationDeg = opts.rotationDeg !== undefined ? opts.rotationDeg : 2.2;
    this.cranial = cranial.clone().normalize();
    this.pivot = pivot.clone();
    // rotation axis: transverse, i.e. perpendicular to the cranial axis and to
    // the anterior direction, so the motion is a nodding of the right hemiliver
    this.axis = new THREE.Vector3(1, 0, 0)
      .sub(this.cranial.clone().multiplyScalar(this.cranial.x)).normalize();
    this.phase01 = 0;
    this.displacementMm = 0;
    this.matrix = new THREE.Matrix4();
    this._q = new THREE.Quaternion();
    this._v = new THREE.Vector3();
  }

  /** Normalised excursion in [0,1]: 0 = end expiration, 1 = end inspiration. */
  static profile (u) {
    const tIn = 1 / (1 + IE_RATIO);
    if (u < tIn) {
      const x = u / tIn;                       // inspiration: liver pushed caudally
      return 0.5 - 0.5 * Math.cos(Math.PI * x);
    }
    const x = (u - tIn) / (1 - tIn);           // expiration: slower recoil
    return 0.5 + 0.5 * Math.cos(Math.PI * x);
  }

  /** @param {number} tMs simulation time */
  update (tMs) {
    const u = ((tMs * 1e-3 * this.fHz) % 1 + 1) % 1;
    this.phase01 = u;
    const a = RespirationField.profile(u);
    // inspiration drives the liver CAUDALLY (away from the cranial direction)
    this.displacementMm = -a * this.excursionCm * 10.0;
    this.rotationRad = -a * THREE.MathUtils.degToRad(this.rotationDeg);
    this._q.setFromAxisAngle(this.axis, this.rotationRad);
    this._v.copy(this.cranial).multiplyScalar(this.displacementMm);
    this.matrix.makeTranslation(this.pivot.x, this.pivot.y, this.pivot.z);
    this.matrix.multiply(new THREE.Matrix4().makeRotationFromQuaternion(this._q));
    this.matrix.multiply(new THREE.Matrix4().makeTranslation(
      -this.pivot.x, -this.pivot.y, -this.pivot.z));
    this.matrix.premultiply(new THREE.Matrix4().makeTranslation(
      this._v.x, this._v.y, this._v.z));
    return this.matrix;
  }

  /** CPU twin used by the automated probe: where does a point go at time t? */
  transformPoint (p, tMs) {
    const m = this.update(tMs);
    return p.clone().applyMatrix4(m);
  }

  describe () {
    return {
      f_hz: this.fHz,
      breaths_per_min: this.breathsPerMin,
      excursion_cm: this.excursionCm,
      excursion_reference_range_cm: this.excursionRange,
      rotation_deg: this.rotationDeg,
      rigid_body: true,
      parenchymal_deformation: 'none — the liver is displaced, not deformed',
      source: this.source,
      is_user_personal_data: this.isPersonal
    };
  }
}
