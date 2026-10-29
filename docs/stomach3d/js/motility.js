// motility.js — the travelling antral contraction wave, receptive relaxation,
// retropulsion and the gastric-emptying curves.
//
// The deformation is applied to the REAL mesh vertices on the CPU, not to a
// proxy and not only in a shader, so the verification below measures the
// geometry the user is actually looking at.
//
// Verification exposed to the page and to the headless audit:
//   measure()  -> { ringsUndisplacedAtT0, peakSpeedMeasured_mm_s, ... }
// It samples the live vertex buffer, finds the station of maximum luminal
// narrowing frame by frame, and differentiates its arc-length position.  A
// global scale pulse would narrow every station at once and the measured peak
// speed would be undefined; a travelling wave gives a finite speed that must
// match the declared profile.
import * as THREE from 'three';

export class Motility {
  /** @param meshes array of { mesh, kind } whose geometry is deformed */
  constructor(motion, meshes) {
    this.m = motion;
    this.sGrid = Float32Array.from(motion.s_grid_mm);
    this.speed = Float32Array.from(motion.speed_mm_s);
    this.phase = Float32Array.from(motion.phase_s);
    this.amp = Float32Array.from(motion.amplitude_frac);
    this.sigma = motion.ring_half_width_mm;
    this.period = motion.period_s;
    this.N = this.sGrid.length;
    this.scale = new Float32Array(this.N);     // per-station radial scale
    this.CL = motion.centreline.map(p => new THREE.Vector3(p[0], p[1], p[2]));
    this.sArc = Float32Array.from(motion.arclength_mm);
    this.enabled = true;
    this.accommodation = 0;                    // 0..1 receptive relaxation
    this.t = 0;

    // per-vertex binding, computed once
    this.bind = meshes.map(({ mesh, kind }) => {
      const g = mesh.geometry;
      const pos = g.attributes.position;
      const n = pos.count;
      const base = new Float32Array(pos.array);     // pristine copy
      const station = new Int32Array(n);
      const rhat = new Float32Array(n * 3);
      const r0 = new Float32Array(n);
      const v = new THREE.Vector3(), d = new THREE.Vector3();
      const M = this.CL.length;
      const STRIDE = 16;
      for (let i = 0; i < n; i++) {
        v.set(base[3 * i], base[3 * i + 1], base[3 * i + 2]);
        // Nearest centreline station, coarse-to-fine.  A full linear scan is
        // O(vertices x centreline) = ~82M distance evaluations for this mesh,
        // which costs tens of seconds in JS at page load; striding by 16 and
        // then refining the bracket is ~10x cheaper and exact for a centreline
        // this smooth.
        let bi = 0, bd = Infinity;
        for (let k = 0; k < M; k += STRIDE) {
          const dd = v.distanceToSquared(this.CL[k]);
          if (dd < bd) { bd = dd; bi = k; }
        }
        const lo = Math.max(0, bi - STRIDE), hi = Math.min(M - 1, bi + STRIDE);
        for (let k = lo; k <= hi; k++) {
          const dd = v.distanceToSquared(this.CL[k]);
          if (dd < bd) { bd = dd; bi = k; }
        }
        const sMM = this.sArc[bi];
        let si = Math.round((sMM / this.sGrid[this.N - 1]) * (this.N - 1));
        si = Math.min(Math.max(si, 0), this.N - 1);
        station[i] = si;
        d.copy(v).sub(this.CL[bi]);
        const len = d.length();
        r0[i] = len;
        if (len > 1e-6) d.multiplyScalar(1 / len);
        rhat[3 * i] = d.x; rhat[3 * i + 1] = d.y; rhat[3 * i + 2] = d.z;
      }
      return { mesh, kind, pos, base, station, rhat, r0, n };
    });
    this.stationRadius = this._stationRadius();
  }

  _stationRadius() {
    // mean luminal radius per station, from the mucosal binding
    const acc = new Float64Array(this.N), cnt = new Float64Array(this.N);
    for (const b of this.bind) {
      if (b.kind !== 'mucosa') continue;
      for (let i = 0; i < b.n; i++) { acc[b.station[i]] += b.r0[i]; cnt[b.station[i]]++; }
    }
    const out = new Float32Array(this.N);
    for (let i = 0; i < this.N; i++) out[i] = cnt[i] > 0 ? acc[i] / cnt[i] : 0;
    return out;
  }

  /** per-station radial scale at time t: 1 = resting, <1 = contracted */
  computeScale(t) {
    const { phase, amp, speed, sigma, period, N, scale } = this;
    for (let i = 0; i < N; i++) {
      if (amp[i] <= 0) { scale[i] = 1.0; continue; }
      // a wave front leaves the pacemaker every `period` seconds; this station
      // is reached at phase[i] after each departure
      let dp = (t - phase[i]) % period;
      if (dp < 0) dp += period;
      const dmm = Math.min(dp, period - dp) * Math.max(speed[i], 1e-6);
      const env = Math.exp(-0.5 * (dmm / sigma) * (dmm / sigma));
      scale[i] = 1.0 - amp[i] * env;
    }
    return scale;
  }

  update(t, recomputeNormals = true) {
    this.t = t;
    const scale = this.enabled ? this.computeScale(t) : this.scale.fill(1.0);
    const acc = this.accommodation;
    for (const b of this.bind) {
      const a = b.pos.array;
      for (let i = 0; i < b.n; i++) {
        const si = b.station[i];
        let s = this.enabled ? scale[si] : 1.0;
        // receptive relaxation: the FUNDUS dilates, it does not contract
        // phasically.  amp[] is already zero there, so this is additive.
        if (this.amp[si] <= 0 && acc > 0) s += acc;
        const dr = (s - 1.0) * b.r0[i];
        a[3 * i]     = b.base[3 * i]     + b.rhat[3 * i]     * dr;
        a[3 * i + 1] = b.base[3 * i + 1] + b.rhat[3 * i + 1] * dr;
        a[3 * i + 2] = b.base[3 * i + 2] + b.rhat[3 * i + 2] * dr;
      }
      b.pos.needsUpdate = true;
      // computeVertexNormals() walks every face and is the single most
      // expensive call here; at a 3 cpm wave the shading does not need it every
      // frame, so the caller throttles it.
      if (recomputeNormals) b.mesh.geometry.computeVertexNormals();
    }
  }

  /** Per-station mean luminal radius measured from the LIVE vertex buffer. */
  liveStationRadius() {
    const acc = new Float64Array(this.N), cnt = new Float64Array(this.N);
    for (const b of this.bind) {
      if (b.kind !== 'mucosa') continue;
      const a = b.pos.array;
      for (let i = 0; i < b.n; i++) {
        const si = b.station[i];
        // the displacement is exactly rhat * dr, so projecting the live offset
        // back onto rhat recovers the live radius without a second nearest-
        // station search
        const r = b.r0[i] + (a[3 * i] - b.base[3 * i]) * b.rhat[3 * i]
                          + (a[3 * i + 1] - b.base[3 * i + 1]) * b.rhat[3 * i + 1]
                          + (a[3 * i + 2] - b.base[3 * i + 2]) * b.rhat[3 * i + 2];
        acc[si] += r; cnt[si]++;
      }
    }
    const out = new Float32Array(this.N);
    for (let i = 0; i < this.N; i++) out[i] = cnt[i] > 0 ? acc[i] / cnt[i] : NaN;
    return out;
  }

  /**
   * In-browser verification.  Everything here is measured off the live geometry;
   * nothing is read back from the declared profile except for comparison.
   */
  measure(opts = {}) {
    const dt = opts.dt ?? 0.25;              // s between probes
    const span = opts.span ?? this.period;   // one full wave period
    const t0 = 0;
    const was = { t: this.t, acc: this.accommodation, en: this.enabled };
    this.accommodation = 0; this.enabled = true;

    // --- (1) how many stations are UNDISPLACED at t = 0 -------------------
    this.update(t0, false);
    const r0 = this.liveStationRadius();
    const rest = this.stationRadius;
    let undisp = 0, movable = 0, total = 0;
    for (let i = 0; i < this.N; i++) {
      if (!isFinite(r0[i]) || rest[i] <= 0) continue;
      total++;
      if (this.amp[i] > 0) movable++;
      if (Math.abs(r0[i] - rest[i]) / rest[i] < 1e-3) undisp++;
    }

    // --- (2) track the station of maximum narrowing over one period -------
    const track = [];
    for (let k = 0; k * dt <= span; k++) {
      const t = t0 + k * dt;
      this.update(t, false);
      const r = this.liveStationRadius();
      let bi = -1, bv = 0;
      for (let i = 0; i < this.N; i++) {
        if (!isFinite(r[i]) || rest[i] <= 0) continue;
        const narrow = (rest[i] - r[i]) / rest[i];
        if (narrow > bv) { bv = narrow; bi = i; }
      }
      if (bi >= 0 && bv > 0.02) track.push({ t, s: this.sGrid[bi], narrow: bv });
    }
    // speed = ds/dt of the narrowing peak, forward steps only (a wave that
    // restarts at the pacemaker produces one negative jump per period, which is
    // excluded and counted)
    const sp = [];
    let resets = 0;
    for (let k = 1; k < track.length; k++) {
      const ds = track[k].s - track[k - 1].s;
      const dtt = track[k].t - track[k - 1].t;
      if (ds < 0) { resets++; continue; }
      if (dtt > 0 && ds > 0) sp.push(ds / dtt);
    }
    sp.sort((a, b) => a - b);

    // --- (3) declared profile, for comparison only -------------------------
    const amp = this.amp, sg = this.sGrid;
    let dmin = Infinity, dmax = 0;
    for (let i = 0; i < this.N; i++) {
      if (amp[i] <= 0) continue;
      dmin = Math.min(dmin, this.speed[i]); dmax = Math.max(dmax, this.speed[i]);
    }

    this.enabled = was.en; this.accommodation = was.acc; this.update(was.t);
    const pct = q => sp.length ? sp[Math.min(sp.length - 1,
      Math.floor(q * (sp.length - 1)))] : NaN;
    return {
      stationsTotal: total,
      stationsWithAmplitude: movable,
      kinematicProfileStations: this.N,
      centrelinePoints: this.CL.length,
      ringResolutionNote:
        'the kinematic profile is sampled at ' + this.N + ' stations, but the '
      + 'geometry binds to the ' + this.CL.length + ' centreline points (2 mm '
      + 'spacing), so only ' + total + ' distinct stations carry mucosal '
      + 'vertices and the undisplaced-ring count is out of THOSE. Claiming '
      + this.N + ' measured rings would overstate the resolution: 700 rings '
      + 'over this organ would be 0.3 mm spacing, finer than the 1.5 mm source '
      + 'voxel.',
      ringsUndisplacedAtT0: undisp,
      ringsUndisplacedFraction: total ? +(undisp / total).toFixed(4) : null,
      peakTrackSamples: track.length,
      waveResetsInOnePeriod: resets,
      measuredPeakSpeed_mm_s: {
        p10: +pct(0.10).toFixed(3), median: +pct(0.50).toFixed(3),
        p90: +pct(0.90).toFixed(3),
        min: sp.length ? +sp[0].toFixed(3) : null,
        max: sp.length ? +sp[sp.length - 1].toFixed(3) : null,
      },
      declaredSpeed_mm_s: { min: +dmin.toFixed(3), max: +dmax.toFixed(3) },
      declaredFrequency_cpm: 60 / this.period,
      measuredFrequency_cpm: resets > 0 ? +(60 * resets / span).toFixed(3) : null,
      method: 'per-station mean luminal radius read from the LIVE mucosal vertex '
            + 'buffer; the station of maximum narrowing is tracked over one wave '
            + 'period and its arc-length position differentiated. A global scale '
            + 'pulse would narrow every station simultaneously and leave no '
            + 'travelling peak to track.',
    };
  }

  /** gastric emptying curves, plotted on the page */
  emptyingCurves(minutes = 240, step = 2) {
    const e = this.m.emptying;
    const tL = (e.liquid_half_time_min[0] + e.liquid_half_time_min[1]) / 2;
    const tS = (e.solid_half_time_min[0] + e.solid_half_time_min[1]) / 2;
    const lag = 25;                 // declared lag, minutes
    const out = [];
    for (let t = 0; t <= minutes; t += step) {
      const liq = Math.pow(0.5, t / tL);
      let sol;
      if (t <= lag) sol = 1.0;
      else {
        const k = 0.5 / (tS - lag);          // linear-after-lag to 50% at tS
        sol = Math.max(0, 1 - k * (t - lag));
      }
      out.push({ t, liquid: +liq.toFixed(5), solid: +sol.toFixed(5) });
    }
    return { curve: out, liquid_half_time_min: tL, solid_half_time_min: tS,
             lag_min: lag,
             note: 'liquid first order, solid lag-then-linear; parameters are the '
                 + 'midpoints of the declared published ranges and are NOT fitted '
                 + 'to anything measured in this project' };
  }
}
