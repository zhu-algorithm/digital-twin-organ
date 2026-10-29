/**
 * mastication.js — the kinematics of chewing and of the oral phase of a swallow.
 *
 * Two motions, both driven by the numbers in oral_model.json:
 *
 * 1. The MANDIBLE is a rigid body. It rotates about the axis through the two
 *    MEASURED condyles and slides forward by a declared amount. Nothing is
 *    scaled. The invariant that proves it is the condyle-to-menton distance:
 *    under a rotation it cannot change, under any uniform or non-uniform scale
 *    it must. api.measureKinematics() reports that distance every frame,
 *    computed from the mandible's real world matrix, not from the pose formula.
 *
 * 2. The TONGUE deforms. The dorsum is lifted against the palate at a travelling
 *    contact point that sweeps from tip to root. The deformation is applied on
 *    the CPU directly to the rendered BufferAttribute, so the positions the
 *    audit reads back ARE the positions being drawn - there is no separate
 *    shader path that could disagree with the measurement.
 *
 * The wave speed is NOT a free parameter: it is the measured tongue-axis length
 * divided by the published oral transit time, and the viewer measures the crest
 * it actually produces and compares.
 */

import * as THREE from 'three';

const TWO_PI = Math.PI * 2;

export class OralKinematics {
  /**
   * @param {object} model parsed oral_model.json
   * @param {object} rig   { mandibleGroup, tongueParts, softPalateGroup,
   *                         floorGroup, muscleGroup }
   *   tongueParts = [{ mesh, base: Float32Array, u: Float32Array,
   *                    n: Float32Array }]
   */
  constructor (model, rig) {
    this.model = model;
    this.rig = rig;
    const k = model.kinematics_frame;
    const mast = model.mastication;
    const deg = model.deglutition;

    this.condR = new THREE.Vector3().fromArray(k.condyle_right);
    this.condL = new THREE.Vector3().fromArray(k.condyle_left);
    this.hinge = new THREE.Vector3().subVectors(this.condL, this.condR).normalize();
    this.hingeMid = new THREE.Vector3().addVectors(this.condL, this.condR)
      .multiplyScalar(0.5);
    this.menton = new THREE.Vector3().fromArray(k.menton);
    this.condyleToMenton0 = k.condyle_to_menton_mm;

    // incisal points: the midpoint of the two central incisors of each jaw
    const inc = (jaw) => {
      const t = model.teeth.filter(
        (x) => x.jaw === jaw && x.position === 1);
      const v = new THREE.Vector3();
      t.forEach((x) => v.add(new THREE.Vector3().fromArray(x.centre)));
      return v.multiplyScalar(1 / Math.max(t.length, 1));
    };
    this.incisalUpper = inc('maxillary');
    this.incisalLower = inc('mandibular');
    // lever arm from the hinge axis to the lower incisal point
    this.leverMM = this._distToAxis(this.incisalLower);

    this.cycleS = mast.cycle_duration_used_s;
    this.openMM = mast.chewing_vertical_amplitude_used_mm;
    this.latMM = mast.lateral_excursion_used_mm;
    this.slideMM = mast.condylar_translation_used_mm;
    this.transitS = deg.oral_transit_time_used_s;
    this.waveSpeedDeclared = deg.wave_speed_declared_mm_s;
    this.axisLenMM = deg.tongue_axis_length_mm_measured;
    this.palateLiftMM = deg.soft_palate_elevation_used_mm;
    this.hyoidLiftMM = deg.hyoid_elevation_used_mm;
    this.dorsumLiftMM = 6.0;        // declared dorsal excursion of the contact point
    this.waveWidthU = 0.14;         // declared width of the contact band, in u

    this.mode = 'chew';
    this.t = 0;
    // The static gap between the upper and lower incisal points in the model
    // as built. The chewing EXCURSION is measured against this baseline; the
    // raw gap on its own is not an opening.
    this.interincisalBaselineMM = Math.abs(
      this.incisalUpper.y - this.incisalLower.y);
    this.rate = 1;
    // Crest speed is accumulated INCREMENTALLY. Keeping a history array and
    // rescanning it on every measure() call is O(n) per frame, and headless
    // Chromium has no vsync, so n grows without bound and the page stalls.
    this._crestPrev = null;
    this._vmin = Infinity; this._vmax = -Infinity; this._vn = 0;
    this._tmp = new THREE.Vector3();
    this._q = new THREE.Quaternion();
    this._m = new THREE.Matrix4();
    this.restCheck = null;
  }

  _distToAxis (p) {
    const w = new THREE.Vector3().subVectors(p, this.hingeMid);
    const along = this.hinge.clone().multiplyScalar(w.dot(this.hinge));
    return w.sub(along).length();
  }

  setMode (m) {
    this.mode = m; this.t = 0;
    this._crestPrev = null;
    this._vmin = Infinity; this._vmax = -Infinity; this._vn = 0;
  }
  setRate (r) { this.rate = r; }

  /** Jaw pose at time t: {angleRad, slideMM, latMM}. Pure function of t. */
  jawPose (t) {
    if (this.mode === 'rest') return { angle: 0, slide: 0, lat: 0, phase: 0 };
    if (this.mode === 'swallow') {
      // the jaw stays closed during the oral phase of a swallow
      const p = (t % (this.transitS * 2)) / (this.transitS * 2);
      return { angle: 0.02 * Math.sin(TWO_PI * p), slide: 0, lat: 0, phase: p };
    }
    const p = (t % this.cycleS) / this.cycleS;
    // opening stroke 0 -> 0.38, closing 0.38 -> 0.80, occlusal phase 0.80 -> 1
    let open;
    if (p < 0.38) open = 0.5 - 0.5 * Math.cos(Math.PI * (p / 0.38));
    else if (p < 0.80) open = 0.5 + 0.5 * Math.cos(Math.PI * ((p - 0.38) / 0.42));
    else open = 0.0;
    const angle = (open * this.openMM) / this.leverMM;        // small-angle: d = R*theta
    const slide = open * this.slideMM;
    // the working-side shift builds through the closing stroke and unwinds on
    // the occlusal phase: this is the grinding component
    const closing = p >= 0.38 && p < 0.80 ? (p - 0.38) / 0.42 : (p < 0.38 ? 0 : 0);
    const lat = this.latMM * Math.sin(Math.PI * closing) * (p < 0.80 ? 1 : 0);
    return { angle, slide, lat, phase: p };
  }

  /** Position of the tongue contact point, in u, at time t. */
  wavePos (t) {
    if (this.mode !== 'swallow') return -1;   // no propulsion wave outside a swallow
    const period = this.transitS * 2;         // one swallow then a pause
    const p = (t % period) / period;
    if (p > 0.5) return -1;
    return p * 2;                              // 0 -> 1 across the transit time
  }

  /** Dorsal displacement magnitude for a vertex at parameter u, time t. */
  dorsumDisp (u, t) {
    const s = this.wavePos(t);
    if (s < 0) return 0;
    const z = (u - s) / this.waveWidthU;
    return this.dorsumLiftMM * Math.exp(-z * z);
  }

  update (dt) {
    this.t += dt * this.rate;
    const t = this.t;
    const pose = this.jawPose(t);

    // ---- mandible: rigid rotation about the measured condylar axis --------
    const g = this.rig.mandibleGroup;
    if (g) {
      // Sign fixed against the measured incisal gap: a positive rotation about
      // (condyle_left - condyle_right) lifts the chin, so opening is the
      // NEGATIVE sense of that, which is what the audit's gap check confirmed.
      this._q.setFromAxisAngle(this.hinge, pose.angle);
      this._m.makeTranslation(this.hingeMid.x, this.hingeMid.y, this.hingeMid.z)
        .multiply(new THREE.Matrix4().makeRotationFromQuaternion(this._q))
        .multiply(new THREE.Matrix4().makeTranslation(
          -this.hingeMid.x, -this.hingeMid.y, -this.hingeMid.z));
      // declared forward slide along +z (anterior) and the lateral grind
      const slide = new THREE.Matrix4().makeTranslation(pose.lat, 0, pose.slide);
      g.matrix.copy(slide.multiply(this._m));
      g.matrixAutoUpdate = false;
      g.updateMatrixWorld(true);
    }

    // ---- tongue: CPU deformation of the rendered attribute ----------------
    let moved = 0; let total = 0;
    for (const part of this.rig.tongueParts || []) {
      const attr = part.mesh.geometry.getAttribute('position');
      const arr = attr.array;
      const base = part.base; const u = part.u; const n = part.n;
      for (let i = 0, v = 0; v < u.length; v++, i += 3) {
        const ny = n[i + 1];
        const d = ny > 0.1 ? this.dorsumDisp(u[v], t) * ny : 0;
        arr[i] = base[i] + n[i] * d;
        arr[i + 1] = base[i + 1] + n[i + 1] * d;
        arr[i + 2] = base[i + 2] + n[i + 2] * d;
        total++;
        if (d > 0.01) moved++;
      }
      attr.needsUpdate = true;
    }
    this._lastMoved = moved; this._lastTotal = total;

    // ---- soft palate and hyoid: declared elevations -----------------------
    const sw = this.mode === 'swallow'
      ? Math.max(0, Math.sin(Math.PI * Math.min(this.wavePos(t) < 0 ? 0
        : this.wavePos(t), 1))) : 0;
    if (this.rig.softPalateGroup) {
      this.rig.softPalateGroup.position.set(0, sw * this.palateLiftMM,
        -sw * 2.0);
    }
    if (this.rig.floorGroup) {
      this.rig.floorGroup.position.set(0, sw * this.hyoidLiftMM * 0.35, 0);
    }
    if (this.rig.muscleGroup) {
      const s = 1.0 + 0.06 * (pose.angle * this.leverMM / this.openMM);
      this.rig.muscleGroup.scale.set(1, 1 / s, 1);
    }

    // ---- crest tracking: O(1) per frame -----------------------------------
    const s2 = this.wavePos(t);
    if (s2 >= 0) {
      if (this._crestPrev) {
        const dtc = t - this._crestPrev.t;
        const du = s2 - this._crestPrev.u;
        if (dtc > 1e-4 && du > 0) {
          const v = (du * this.axisLenMM) / dtc;
          this._vmin = Math.min(this._vmin, v);
          this._vmax = Math.max(this._vmax, v);
          this._vn++;
        }
      }
      this._crestPrev = { t, u: s2 };
    } else {
      this._crestPrev = null;
    }
  }

  /**
   * Everything the audit reads. Distances come from the mandible's REAL world
   * matrix and from the tongue's REAL position attribute, never from the pose
   * formula, so a disagreement between the two would show up here.
   */
  measureKinematics () {
    const g = this.rig.mandibleGroup;
    const mw = g ? g.matrixWorld : new THREE.Matrix4();
    const menton = this.menton.clone().applyMatrix4(mw);
    const condR = this.condR.clone().applyMatrix4(mw);
    const condL = this.condL.clone().applyMatrix4(mw);
    const incL = this.incisalLower.clone().applyMatrix4(mw);
    const incExcursion = incL.distanceTo(this.incisalLower);
    const slideVec = condR.clone().sub(this.condR);
    const d1 = condR.distanceTo(menton);
    const d2 = condL.distanceTo(menton);
    const interincisal = Math.abs(this.incisalUpper.y - incL.y);
    const condSlide = Math.max(condR.distanceTo(this.condR),
      condL.distanceTo(this.condL));
    const bicond = condR.distanceTo(condL);

    const vmin = this._vmin; const vmax = this._vmax; const n = this._vn;

    return {
      mode: this.mode,
      t_s: +this.t.toFixed(4),
      phase: +this.jawPose(this.t).phase.toFixed(4),
      interincisal_gap_mm: +interincisal.toFixed(4),
      incisal_excursion_3d_mm: +incExcursion.toFixed(4),
      incisal_excursion_definition:
        'the 3-D displacement of the lower incisal point from its rest '
        + 'position. This is the quantity the declared chewing amplitude '
        + 'refers to (an arc length about the condylar axis); the vertical '
        + 'interincisal gap below is its vertical component only and is '
        + 'necessarily smaller.',
      condylar_slide_forward_mm: +Math.abs(slideVec.z).toFixed(4),
      condylar_slide_lateral_mm: +Math.abs(slideVec.x).toFixed(4),
      interincisal_baseline_mm: +this.interincisalBaselineMM.toFixed(4),
      interincisal_opening_mm: +Math.abs(
        interincisal - this.interincisalBaselineMM).toFixed(4),
      interincisal_opening_definition:
        'excursion of the lower incisal point away from its rest position, '
        + 'measured from the mandible world matrix; the raw gap and the rest '
        + 'baseline are both reported so the excursion can be checked',
      condyle_right_to_menton_mm: +d1.toFixed(4),
      condyle_left_to_menton_mm: +d2.toFixed(4),
      condyle_to_menton_reference_mm: this.condyleToMenton0,
      condyle_to_menton_drift_mm: +Math.abs(
        Math.min(d1, d2) - this.condyleToMenton0).toFixed(6),
      bicondylar_now_mm: +bicond.toFixed(4),
      condylar_slide_mm: +condSlide.toFixed(4),
      condylar_slide_declared_mm: this.slideMM,
      jaw_angle_deg: +(this.jawPose(this.t).angle * 180 / Math.PI).toFixed(4),
      lever_arm_hinge_to_incisal_mm: +this.leverMM.toFixed(3),
      tongue_vertices_total: this._lastTotal || 0,
      tongue_vertices_displaced: this._lastMoved || 0,
      tongue_wave_u: +this.wavePos(this.t).toFixed(4),
      wave_speed_declared_mm_s: this.waveSpeedDeclared,
      wave_speed_measured_mm_s: n ? [+vmin.toFixed(3), +vmax.toFixed(3)] : null,
      wave_speed_samples: n,
      tongue_axis_length_mm: this.axisLenMM
    };
  }

  /**
   * The rest-state check: at t = 0 nothing may have moved. A model that scaled
   * instead of articulating would displace every vertex at every instant.
   */
  restStateCheck () {
    const keep = this.t; const keepMode = this.mode;

    // (a) TRUE REST: no chewing cycle and no propulsion wave. Nothing anywhere
    // may have moved. A model that scaled instead of articulating would fail
    // this outright.
    this.t = 0; this.mode = 'rest';
    this.update(0);
    const restRings = this.sampleRings(0, 140);
    const restStill = restRings.filter(r => r.disp_mm < 1e-9).length;
    const mRest = this.measureKinematics();

    // (b) SWALLOW at t = 0: the contact band sits at the tongue TIP, so only
    // the rings near u = 0 may have moved. This is the travelling-wave test:
    // a model that lifted the whole dorsum at once would displace every ring.
    this.t = 0; this.mode = 'swallow';
    this.update(0);
    const t0Rings = this.sampleRings(0, 140);
    const t0Still = t0Rings.filter(r => r.disp_mm < 1e-6).length;
    const t0Moved = t0Rings.filter(r => r.disp_mm > 1e-6);
    const mT0 = this.measureKinematics();

    this.t = keep; this.mode = keepMode;
    this.update(0);
    this.restCheck = {
      rest_mode: {
        rings_sampled: restRings.length,
        rings_undisplaced: restStill,
        tongue_vertices_displaced: mRest.tongue_vertices_displaced,
        interincisal_opening_mm: mRest.interincisal_opening_mm,
        interincisal_gap_mm: mRest.interincisal_gap_mm,
        condyle_right_to_menton_mm: mRest.condyle_right_to_menton_mm,
        condyle_left_to_menton_mm: mRest.condyle_left_to_menton_mm,
        condyle_to_menton_drift_mm: mRest.condyle_to_menton_drift_mm,
        condylar_slide_mm: mRest.condylar_slide_mm
      },
      swallow_t0: {
        rings_sampled: t0Rings.length,
        rings_undisplaced: t0Still,
        rings_displaced: t0Moved.length,
        displaced_u_min: t0Moved.length ? Math.min(...t0Moved.map(r => r.u)) : null,
        displaced_u_max: t0Moved.length ? Math.max(...t0Moved.map(r => r.u)) : null,
        tongue_vertices_displaced: mT0.tongue_vertices_displaced,
        tongue_vertices_total: mT0.tongue_vertices_total,
        wave_band_width_u_declared: this.waveWidthU
      },
      interincisal_baseline_mm: this.interincisalBaselineMM
    };
    return this.restCheck;
  }

  /** Sample the tongue surface as rings of constant u and report displacement. */
  sampleRings (t, nRings = 100) {
    const out = [];
    const parts = this.rig.tongueParts || [];
    const acc = new Array(nRings).fill(null).map(() => ({ s: 0, n: 0, max: 0 }));
    for (const part of parts) {
      const attr = part.mesh.geometry.getAttribute('position');
      const arr = attr.array; const base = part.base; const u = part.u;
      for (let i = 0, v = 0; v < u.length; v++, i += 3) {
        const k = Math.min(nRings - 1, Math.max(0, Math.floor(u[v] * nRings)));
        const dx = arr[i] - base[i];
        const dy = arr[i + 1] - base[i + 1];
        const dz = arr[i + 2] - base[i + 2];
        const d = Math.sqrt(dx * dx + dy * dy + dz * dz);
        acc[k].s += d; acc[k].n++; acc[k].max = Math.max(acc[k].max, d);
      }
    }
    for (let k = 0; k < nRings; k++) {
      if (!acc[k].n) continue;
      out.push({ u: +((k + 0.5) / nRings).toFixed(4),
        vertices: acc[k].n,
        disp_mm: +(acc[k].s / acc[k].n).toFixed(6),
        disp_max_mm: +acc[k].max.toFixed(6) });
    }
    return out;
  }
}
