/**
 * contract.js — consumer of the `heart_timeseries.json` interface contract.
 *
 * The contract is produced by the cardiac-model track. This module does NOT
 * invent any field: it reads exactly the keys fixed in the contract and
 * exposes a time-indexed sampler plus the derived calibration constants that
 * drive the 3-D deformation.
 */

export const SIGNAL_KEYS = [
  't_ms',
  'V_LV_mL', 'V_RV_mL', 'V_LA_mL', 'V_RA_mL',
  'P_LV_mmHg', 'P_RV_mmHg', 'P_Ao_mmHg', 'P_PA_mmHg',
  'Q_mitral_mL_s', 'Q_aortic_mL_s', 'Q_tricuspid_mL_s',
  'Q_pulmonary_mL_s', 'Q_coronary_mL_s',
  'Vm_mV'
];

export const VALVE_KEYS = ['mitral', 'aortic', 'tricuspid', 'pulmonary'];

export class HeartTimeseries {
  constructor (json) {
    this.raw = json;
    this.meta = json.meta;
    this.derived = json.derived;
    this.perBeat = json.per_beat;

    this.dt = this.meta.dt_ms;
    this.duration = this.meta.duration_ms;
    this.n = Math.round(this.duration / this.dt);

    this.sig = {};
    for (const k of SIGNAL_KEYS) {
      const a = json.signals[k];
      if (!a) throw new Error(`heart_timeseries.json: missing signals.${k}`);
      if (a.length !== this.n) {
        throw new Error(
          `heart_timeseries.json: signals.${k} length ${a.length} != duration_ms/dt_ms (${this.n})`);
      }
      this.sig[k] = Float32Array.from(a);
    }
    this.valve = {};
    for (const k of VALVE_KEYS) {
      const a = json.signals.valve_open[k];
      if (!a) throw new Error(`heart_timeseries.json: missing signals.valve_open.${k}`);
      this.valve[k] = Uint8Array.from(a);
    }

    // Per-chamber volume extrema read from the contract signals themselves, so
    // beat-to-beat variation (short RR -> less filling -> smaller excursion)
    // emerges from the model rather than from cosmetic jitter.
    this.vRange = {};
    for (const k of ['V_LV_mL', 'V_RV_mL', 'V_LA_mL', 'V_RA_mL']) {
      const a = this.sig[k];
      let lo = Infinity; let hi = -Infinity;
      for (let i = 0; i < a.length; i++) { if (a[i] < lo) lo = a[i]; if (a[i] > hi) hi = a[i]; }
      this.vRange[k] = { min: lo, max: hi, span: Math.max(hi - lo, 1e-6) };
    }

    this.qPeak = {};
    for (const k of ['Q_mitral_mL_s', 'Q_aortic_mL_s', 'Q_tricuspid_mL_s',
      'Q_pulmonary_mL_s', 'Q_coronary_mL_s']) {
      let m = 0; const a = this.sig[k];
      for (let i = 0; i < a.length; i++) m = Math.max(m, Math.abs(a[i]));
      this.qPeak[k] = Math.max(m, 1e-6);
    }
  }

  index (tMs) {
    let t = tMs % this.duration;
    if (t < 0) t += this.duration;
    return t / this.dt;
  }

  sample (key, tMs) {
    const a = this.sig[key];
    const f = this.index(tMs);
    const i0 = Math.floor(f);
    const i1 = (i0 + 1) % a.length;
    const w = f - i0;
    return a[i0] * (1 - w) + a[i1] * w;
  }

  valveOpen (key, tMs) {
    const v = this.valve[key];
    return v[Math.floor(this.index(tMs)) % v.length];
  }

  /**
   * Normalised contraction coordinate c in [0,1] for a chamber:
   * c = 0 at that chamber's maximum volume, c = 1 at its minimum volume.
   * This is the scalar the volume -> deformation calibration is expressed in.
   */
  contraction (volKey, tMs) {
    const r = this.vRange[volKey];
    return (r.max - this.sample(volKey, tMs)) / r.span;
  }

  beatIndex (tMs) {
    const onsets = this.meta.beat_onset_ms;
    let t = tMs % this.duration; if (t < 0) t += this.duration;
    let lo = 0; let hi = onsets.length - 1; let res = -1;
    while (lo <= hi) {
      const mid = (lo + hi) >> 1;
      if (onsets[mid] <= t) { res = mid; lo = mid + 1; } else { hi = mid - 1; }
    }
    return res;
  }

  provenanceLine () {
    const m = this.meta;
    return `${m.data_source} | state=${m.state} | HR=${m.hr_bpm_mean.toFixed(1)} bpm | ` +
           `is_user_personal_data=${m.is_user_personal_data}`;
  }
}

export async function loadTimeseries (url) {
  const r = await fetch(url);
  if (!r.ok) throw new Error(`cannot fetch ${url}: ${r.status}`);
  return new HeartTimeseries(await r.json());
}
