/**
 * motion.js — renal respiratory displacement.
 *
 * The kidney is not a contractile organ. During quiet breathing it is carried
 * up and down by the diaphragm, mainly head-to-foot, with a slight rotation as
 * it rides on the psoas; the parenchyma itself does not deform appreciably.
 * This module therefore produces a RIGID-BODY transform and nothing else —
 * there is deliberately no per-vertex deformation anywhere in this scene, so
 * the shape you see at any instant of the cycle is the segmented shape.
 *
 * Frequency: 0.18323 Hz (10.99 breaths/min), carried over unchanged from this
 * project's cardiac model, where it was fitted as the respiratory sinus
 * arrhythmia frequency of PhysioNet nsrdb record 16265. That is PUBLIC
 * placeholder data, not the user's own recording (is_user_personal_data=false).
 *
 * Amplitude: craniocaudal excursion 1-2 cm is a published range for the kidney
 * in quiet breathing, NOT a measurement of this subject. The default is the
 * midpoint and the scene reports which value it is using.
 *
 * Waveform: inspiration is shorter than expiration (I:E about 1:1.5 at rest),
 * so a plain sinusoid would be the wrong shape even at the right period.
 */

import * as THREE from 'three';

export const IE_RATIO = 1.5;

export class RespirationField {
  /**
   * @param {object} spec  the `respiration` block of kidney_model.json
   * @param {THREE.Vector3} cranial  unit cranial direction in scene space
   * @param {THREE.Vector3} pivot    rotation centre
   */
  constructor (spec, cranial, pivot, opts = {}) {
    this.fHz = spec.f_hz;
    this.breathsPerMin = spec.breaths_per_min;
    this.source = spec.source;
    this.isPersonal = !!spec.is_user_personal_data;
    this.excursionRange = spec.craniocaudal_excursion_cm || [1, 2];
    this.excursionCm = opts.excursionCm !== undefined
      ? opts.excursionCm
      : (this.excursionRange[0] + this.excursionRange[1]) * 0.5;
    this.rotationDeg = opts.rotationDeg !== undefined ? opts.rotationDeg : 1.6;
    this.cranial = cranial.clone().normalize();
    this.pivot = pivot.clone();
    this.axis = new THREE.Vector3(1, 0, 0)
      .sub(this.cranial.clone().multiplyScalar(this.cranial.x)).normalize();
    this.phase01 = 0;
    this.displacementMm = 0;
    this.matrix = new THREE.Matrix4();
    this._q = new THREE.Quaternion();
    this._v = new THREE.Vector3();
    this._m1 = new THREE.Matrix4();
    this._m2 = new THREE.Matrix4();
  }

  /** Normalised excursion in [0,1]: 0 = end expiration, 1 = end inspiration. */
  static profile (u) {
    const tIn = 1 / (1 + IE_RATIO);
    if (u < tIn) {
      const x = u / tIn;
      return 0.5 - 0.5 * Math.cos(Math.PI * x);
    }
    const x = (u - tIn) / (1 - tIn);
    return 0.5 + 0.5 * Math.cos(Math.PI * x);
  }

  /** @param {number} tMs simulation time */
  update (tMs) {
    const u = ((tMs * 1e-3 * this.fHz) % 1 + 1) % 1;
    this.phase01 = u;
    const a = RespirationField.profile(u);
    // inspiration drives the kidney CAUDALLY
    this.displacementMm = -a * this.excursionCm * 10.0;
    this.rotationRad = -a * THREE.MathUtils.degToRad(this.rotationDeg);
    this._q.setFromAxisAngle(this.axis, this.rotationRad);
    this._v.copy(this.cranial).multiplyScalar(this.displacementMm);
    this.matrix.makeTranslation(this.pivot.x, this.pivot.y, this.pivot.z);
    this.matrix.multiply(this._m1.makeRotationFromQuaternion(this._q));
    this.matrix.multiply(this._m2.makeTranslation(
      -this.pivot.x, -this.pivot.y, -this.pivot.z));
    this.matrix.premultiply(this._m1.makeTranslation(
      this._v.x, this._v.y, this._v.z));
    return this.matrix;
  }

  describe () {
    return {
      f_hz: this.fHz,
      breaths_per_min: this.breathsPerMin,
      excursion_cm: this.excursionCm,
      excursion_reference_range_cm: this.excursionRange,
      rotation_deg: this.rotationDeg,
      rigid_body: true,
      parenchymal_deformation: 'none — the kidney is displaced, not deformed',
      source: this.source,
      is_user_personal_data: this.isPersonal
    };
  }
}
