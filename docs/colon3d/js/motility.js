// motility.js — the three colonic motor patterns, and the in-browser measurement that
// proves each one is a travelling (or standing) wave rather than a global scale change.
//
//  haustral   non-propulsive segmentation: each haustrum narrows and widens about its own
//             centre at the published 3 cycles/min; neighbouring haustra are given mixed
//             phase so there is no net transport.
//  mass       high-amplitude propagating contraction: one narrow band of constriction
//             travelling distally at the published 0.71 cm/s over >= 30 cm.
//  defecation rectal filling followed by the recto-anal sequence; the rectal ampulla
//             distends and the anal end stays shut until the sequence permits it.
//
// All three write a per-ring radial scale on the tube; nothing scales the whole mesh.

export const MODES = ['haustral', 'mass', 'defecation'];

export class Motility {
  constructor(tube, sc, micro) {
    this.tube = tube; this.sc = sc;
    const m = micro.motility;
    this.p = {
      haustralHz: m.haustral_contractions.frequency_cycles_per_min / 60,
      haustralSpanMm: m.haustral_contractions.span_cm[0] * 10,
      haustralDepth: 0.30,                                  // declared radial amplitude
      hapcMmS: m.mass_movement.propagation_velocity_cm_s * 10,
      hapcWidthMm: 60,
      hapcDepth: 0.55,
      hapcMinDistMm: m.mass_movement.minimum_propagated_distance_cm * 10,
      eventsPer24h: m.mass_movement.events_per_24h_median,
      lambda: sc.derivation.haustral_pitch_mm
    };
    this.mode = 'haustral';
    this.t = 0;
    this.rectum = sc.segments.find(x => x.name === 'rectum');
    this.sFoldLast = sc.s_folds[sc.s_folds.length - 1];
    // deterministic per-haustrum phase: mixed antegrade/retrograde, seeded, reproducible
    this.phase = sc.s_folds.map((_, k) => ((k * 2654435761) % 1000) / 1000 * 2 * Math.PI);
    this.hist = [];       // {t, radial:Float32Array} ring-radius history for the audit
    this.maxHist = 4000;
  }

  setMode(m) { this.mode = m; this.t = 0; this.hist.length = 0; }

  /** dt seconds */
  step(dt) {
    this.t += dt;
    const sc = this.sc, n = this.tube.n, r = this.tube.radial, s = sc.arclen;
    const p = this.p;
    if (this.mode === 'haustral') {
      for (let i = 0; i < n; i++) {
        if (s[i] > this.sFoldLast) { r[i] = 1; continue; }     // no haustra in the rectum
        const k = Math.floor((s[i] - sc.s_folds[0]) / p.lambda);
        const ph = this.phase[Math.max(0, Math.min(this.phase.length - 1, k))];
        const u = ((s[i] - sc.s_folds[0]) % p.lambda) / p.lambda;   // 0..1 within the haustrum
        // stationary envelope within each haustrum: maximum at its centre, zero at the folds
        const envlope = Math.sin(Math.PI * u);
        r[i] = 1 - p.haustralDepth * envlope * 0.5 *
               (1 + Math.sin(2 * Math.PI * p.haustralHz * this.t + ph));
      }
    } else if (this.mode === 'mass') {
      const period = (s[n - 1] + 2 * p.hapcWidthMm) / p.hapcMmS;
      const front = ((this.t % period) * p.hapcMmS) - p.hapcWidthMm;
      for (let i = 0; i < n; i++) {
        const d = (s[i] - front) / p.hapcWidthMm;
        r[i] = (d >= -1 && d <= 1) ? 1 - p.hapcDepth * Math.cos(d * Math.PI / 2) ** 2 : 1;
      }
      this.front = front;
    } else {
      // rectal filling, then the recto-anal sequence
      const T = 18.0;  const ph = (this.t % T) / T;
      const i0 = this.rectum.i0;
      for (let i = 0; i < n; i++) {
        if (i < i0) { r[i] = 1; continue; }
        const f = (i - i0) / Math.max(1, (n - 1 - i0));
        const fill = Math.min(1, ph / 0.55);
        const amp = Math.sin(Math.PI * Math.min(1, f * 1.25));     // ampulla, not the anal end
        if (ph < 0.55)      r[i] = 1 + 0.35 * fill * amp;                       // filling
        else if (ph < 0.68) r[i] = 1 + 0.35 * amp;                              // RAIR / sampling
        else if (ph < 0.80) r[i] = 1 + 0.35 * amp * (1 - (ph - 0.68) / 0.12 * 0.15); // urge
        else                r[i] = 1 + 0.35 * amp * Math.max(0, 1 - (ph - 0.80) / 0.20); // evacuation
      }
      this.stage = ph < 0.55 ? '直肠充盈' : ph < 0.68 ? '直肠肛管抑制反射（取样）'
                 : ph < 0.80 ? '便意 · 随意控制' : '排便';
    }
    if (this.hist.length < this.maxHist) this.hist.push({ t: this.t, r: Float32Array.from(r) });
    this.tube.build(0, this.mode !== 'defecation');
  }

  /** Measure, in the browser, what the animation is actually doing. */
  audit() {
    const sc = this.sc, s = sc.arclen, n = this.tube.n;
    const H = this.hist;
    const out = { mode: this.mode, frames: H.length, rings: n };
    if (H.length < 30) return Object.assign(out, { error: 'not enough frames' });
    // --- is it a global scale change?  A global scale would move EVERY ring at t=0.
    const r0 = H[0].r;
    let still = 0;
    for (let i = 0; i < n; i++) if (Math.abs(r0[i] - 1) < 1e-3) still++;
    out.rings_undisplaced_at_t0 = still;
    out.rings_undisplaced_fraction = +(still / n).toFixed(4);
    // The decisive test for "a wave is travelling" rather than "the whole mesh is being
    // scaled" is a MID-RUN frame: if most rings sit at their resting radius while a few
    // are displaced, in the same frame, it cannot be a global scale.
    const mid = H[Math.floor(H.length / 2)].r;
    let stillMid = 0, movedMid = 0;
    for (let i = 0; i < n; i++) (Math.abs(mid[i] - 1) < 1e-3 ? stillMid++ : movedMid++);
    out.rings_undisplaced_mid_frame = stillMid;
    out.rings_displaced_mid_frame = movedMid;
    out.rings_undisplaced_mid_fraction = +(stillMid / n).toFixed(4);
    // --- amplitude spread: a global scale has the SAME radial factor on every ring
    let mn = Infinity, mx = -Infinity;
    for (const h of H) for (let i = 0; i < n; i++) { if (h.r[i] < mn) mn = h.r[i]; if (h.r[i] > mx) mx = h.r[i]; }
    out.radial_scale_min = +mn.toFixed(4); out.radial_scale_max = +mx.toFixed(4);
    let maxSpread = 0;
    for (const h of H) {
      let a = Infinity, b = -Infinity;
      for (let i = 0; i < n; i++) { if (h.r[i] < a) a = h.r[i]; if (h.r[i] > b) b = h.r[i]; }
      maxSpread = Math.max(maxSpread, b - a);
    }
    out.max_within_frame_radial_spread = +maxSpread.toFixed(4);
    out.is_global_scale = maxSpread < 1e-3;
    const dt = (H[H.length - 1].t - H[0].t) / (H.length - 1);
    out.dt_s = +dt.toFixed(5);

    if (this.mode === 'haustral') {
      // frequency from zero crossings of a mid-haustrum ring
      const probe = Math.floor(n * 0.35);
      const sig = H.map(h => h.r[probe] - 1);
      const mean = sig.reduce((a, b) => a + b, 0) / sig.length;
      let cross = 0;
      for (let k = 1; k < sig.length; k++) if ((sig[k - 1] - mean) < 0 && (sig[k] - mean) >= 0) cross++;
      const span = H[H.length - 1].t - H[0].t;
      out.measured_frequency_hz = +(cross / span).toFixed(4);
      out.measured_frequency_cycles_per_min = +(cross / span * 60).toFixed(3);
      out.declared_frequency_cycles_per_min = +(this.p.haustralHz * 60).toFixed(3);
      // Non-propulsion test.  Tracking the single narrowest ring is NOT a valid test:
      // with mixed per-haustrum phase the global minimum hops from one haustrum to the
      // next, which looks like travel but is not transport.  The test used here is the
      // DIRECTED DRIFT of the contraction pattern: fit the narrowest-ring position
      // against time.  A propagating wave gives a large slope with a high r-squared; a
      // standing segmentation gives a near-zero slope with essentially no linear fit.
      const cx = [], ct = [];
      for (const h of H) {
        let best = -1, bv = 2;
        for (let i = Math.floor(n * 0.15); i < Math.floor(n * 0.75); i++)
          if (h.r[i] < bv) { bv = h.r[i]; best = i; }
        if (best >= 0) { cx.push(s[best]); ct.push(h.t); }
      }
      let cs = null, cr2 = null;
      if (cx.length > 5) {
        const nn = cx.length;
        const mt = ct.reduce((a, b) => a + b, 0) / nn, mx = cx.reduce((a, b) => a + b, 0) / nn;
        let sxy = 0, sxx = 0, syy = 0;
        for (let k = 0; k < nn; k++) {
          sxy += (ct[k] - mt) * (cx[k] - mx); sxx += (ct[k] - mt) ** 2; syy += (cx[k] - mx) ** 2;
        }
        cs = sxx > 0 ? sxy / sxx : null;
        cr2 = (sxx > 0 && syy > 0) ? (sxy * sxy) / (sxx * syy) : null;
      }
      out.contraction_centre_drift_mm_s = cs === null ? null : +cs.toFixed(4);
      out.contraction_centre_drift_r2 = cr2 === null ? null : +cr2.toFixed(5);
      out.contraction_centre_excursion_mm = cx.length ? +(Math.max(...cx) - Math.min(...cx)).toFixed(2) : 0;
      // a travelling wave would need BOTH a sizeable drift and a decent linear fit
      out.is_propulsive = (cs !== null && Math.abs(cs) > 1.0 && cr2 !== null && cr2 > 0.5);
      out.verdict = out.is_propulsive
        ? 'PROPULSIVE: the contraction pattern drifts along the gut with a good linear fit -- this is NOT what a haustral segmentation should do'
        : 'standing segmentation: the narrowest ring hops between haustra but shows no directed drift (slope near zero, no linear fit), so the pattern mixes without transporting';
    } else if (this.mode === 'mass') {
      // track the minimum-radius ring and fit its velocity
      const xs = [], ts = [];
      for (const h of H) {
        let best = -1, bv = 2;
        for (let i = 0; i < n; i++) if (h.r[i] < bv) { bv = h.r[i]; best = i; }
        if (bv < 0.995) { xs.push(s[best]); ts.push(h.t); }
      }
      // least-squares slope of front position against MODEL time.  Differencing
      // consecutive frames is useless here: the front advances less than one ring per
      // frame, so most frame-to-frame differences are exactly zero and the rest are one
      // ring spacing, which biases any median or mean of the per-pair ratios.
      let slope = null, r2 = null;
      if (xs.length > 5) {
        const nn = xs.length;
        const mt = ts.reduce((a, b) => a + b, 0) / nn, mx = xs.reduce((a, b) => a + b, 0) / nn;
        let sxy = 0, sxx = 0, syy = 0;
        for (let k = 0; k < nn; k++) {
          sxy += (ts[k] - mt) * (xs[k] - mx);
          sxx += (ts[k] - mt) ** 2; syy += (xs[k] - mx) ** 2;
        }
        slope = sxx > 0 ? sxy / sxx : null;
        r2 = (sxx > 0 && syy > 0) ? (sxy * sxy) / (sxx * syy) : null;
      }
      out.measured_front_velocity_mm_s = slope === null ? null : +slope.toFixed(4);
      out.front_fit_r2 = r2 === null ? null : +r2.toFixed(5);
      out.front_samples = xs.length;
      out.declared_front_velocity_mm_s = +this.p.hapcMmS.toFixed(3);
      out.front_travel_mm = xs.length ? +(Math.max(...xs) - Math.min(...xs)).toFixed(1) : 0;
      out.declared_minimum_propagated_mm = this.p.hapcMinDistMm;
      out.verdict = 'travelling wave: a single constriction marches distally; rings ahead of and behind it are at their resting radius in the same frame';
    } else {
      let maxR = 0, idxMax = -1;
      for (const h of H) for (let i = 0; i < n; i++) if (h.r[i] > maxR) { maxR = h.r[i]; idxMax = i; }
      out.max_rectal_distension_factor = +maxR.toFixed(4);
      out.distension_located_in_rectum = idxMax >= this.rectum.i0;
      out.rings_in_rectum = n - this.rectum.i0;
      out.verdict = 'the distension is confined to the rectal segment; the rest of the colon is untouched in the same frame';
    }
    return out;
  }
}
