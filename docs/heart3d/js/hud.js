/**
 * hud.js — 2-D overlay: scrolling action potential, ventricular pressure /
 * volume, valve state lamps, pressure-volume loop and live frame rate.
 * Everything plotted is read straight out of the contract so the panel can be
 * compared against the cardiac-model track without reinterpretation.
 */

export class HUD {
  constructor (canvas) {
    this.c = canvas;
    this.ctx = canvas.getContext('2d');
    this.frames = [];
    this.fps = 0;
    this._lastT = performance.now();
    this._acc = 0;
    this._n = 0;
  }

  resize (w, h, dpr) {
    this.c.width = Math.floor(w * dpr);
    this.c.height = Math.floor(h * dpr);
    this.c.style.width = w + 'px';
    this.c.style.height = h + 'px';
    this.ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    this.w = w; this.h = h;
  }

  tickFps (now) {
    const dt = now - this._lastT;
    this._lastT = now;
    this._acc += dt; this._n++;
    if (this._acc > 400) {
      this.fps = (this._n * 1000) / this._acc;
      this._acc = 0; this._n = 0;
    }
  }

  _trace (x, y, w, h, arr, i0, span, lo, hi, color, label, unit) {
    const ctx = this.ctx;
    ctx.save();
    ctx.strokeStyle = 'rgba(255,255,255,0.13)';
    ctx.lineWidth = 1;
    ctx.strokeRect(x, y, w, h);
    ctx.beginPath();
    for (let k = 0; k < w; k++) {
      const idx = Math.floor(i0 + (k / w) * span) % arr.length;
      const v = arr[idx < 0 ? idx + arr.length : idx];
      const yy = y + h - ((v - lo) / (hi - lo)) * h;
      if (k === 0) ctx.moveTo(x + k, yy); else ctx.lineTo(x + k, yy);
    }
    ctx.strokeStyle = color;
    ctx.lineWidth = 1.6;
    ctx.stroke();
    ctx.fillStyle = 'rgba(230,236,245,0.82)';
    ctx.font = '11px ui-monospace, Menlo, monospace';
    ctx.fillText(`${label} [${lo.toFixed(0)}..${hi.toFixed(0)} ${unit}]`, x + 4, y + 12);
    ctx.restore();
  }

  draw (ts, tMs, field, flow, extra) {
    const ctx = this.ctx;
    const W = this.w; const H = this.h;
    ctx.clearRect(0, 0, W, H);

    const panelW = Math.min(330, W * 0.34);
    const x0 = W - panelW - 12;
    let y = 12;

    ctx.fillStyle = 'rgba(8,11,18,0.62)';
    ctx.strokeStyle = 'rgba(120,140,170,0.25)';
    ctx.lineWidth = 1;
    const panelH = 410;
    ctx.beginPath();
    ctx.roundRect(x0 - 10, y - 6, panelW + 20, panelH, 10);
    ctx.fill(); ctx.stroke();

    ctx.fillStyle = '#e9eef6';
    ctx.font = '600 13px ui-sans-serif, system-ui, sans-serif';
    ctx.fillText('Cardiac contract signals', x0, y + 10);
    y += 22;

    ctx.font = '10px ui-monospace, Menlo, monospace';
    ctx.fillStyle = 'rgba(200,212,230,0.72)';
    const prov = ts.provenanceLine();
    for (const line of wrap(prov, 46)) { ctx.fillText(line, x0, y); y += 12; }
    y += 4;

    const windowMs = 3000;
    const span = windowMs / ts.dt;
    const i0 = ts.index(tMs) - span;

    const tw = panelW; const th = 52;
    this._trace(x0, y, tw, th, ts.sig.Vm_mV, i0, span, -95, 45, '#56d4ff', 'Vm', 'mV');
    y += th + 8;
    this._trace(x0, y, tw, th, ts.sig.P_LV_mmHg, i0, span, 0, 140, '#ff6b6b', 'P_LV', 'mmHg');
    y += th + 8;
    this._trace(x0, y, tw, th, ts.sig.V_LV_mL, i0, span,
      ts.vRange.V_LV_mL.min - 5, ts.vRange.V_LV_mL.max + 5, '#ffd166', 'V_LV', 'mL');
    y += th + 8;
    this._trace(x0, y, tw, th, ts.sig.Q_coronary_mL_s, i0, span,
      0, ts.qPeak.Q_coronary_mL_s, '#b388ff', 'Q_coronary', 'mL/s');
    y += th + 10;

    // valve lamps
    ctx.font = '10px ui-monospace, Menlo, monospace';
    const lamps = [['mitral', '#ffd166'], ['aortic', '#ff6b6b'],
      ['tricuspid', '#7ee787'], ['pulmonary', '#56d4ff']];
    let lx = x0;
    for (const [name, col] of lamps) {
      const on = ts.valveOpen(name, tMs);
      ctx.beginPath();
      ctx.arc(lx + 5, y + 5, 5, 0, Math.PI * 2);
      ctx.fillStyle = on ? col : 'rgba(255,255,255,0.14)';
      ctx.fill();
      ctx.fillStyle = on ? '#e9eef6' : 'rgba(200,212,230,0.45)';
      ctx.fillText(name.slice(0, 4), lx + 13, y + 9);
      lx += panelW / 4;
    }
    y += 22;

    const beat = ts.beatIndex(tMs);
    const rr = beat >= 0 ? ts.perBeat.rr_ms[beat] : NaN;
    const sv = beat >= 0 ? ts.perBeat.sv_mL[beat] : NaN;
    const ef = beat >= 0 ? ts.perBeat.ef_pct[beat] : NaN;
    const qt = beat >= 0 ? ts.perBeat.qt_ms[beat] : NaN;
    ctx.fillStyle = '#e9eef6';
    ctx.font = '11px ui-monospace, Menlo, monospace';
    const rows = [
      `beat ${beat + 1}/${ts.perBeat.rr_ms.length}  RR ${fmt(rr, 0)} ms`,
      `SV ${fmt(sv, 1)} mL   EF ${fmt(ef, 1)} %   QT ${fmt(qt, 0)} ms`,
      `c_LV ${field.state.cLV.toFixed(2)}  c_RV ${field.state.cRV.toFixed(2)}`,
      `c_LA ${field.state.cLA.toFixed(2)}  c_RA ${field.state.cRA.toFixed(2)}`,
      /* Two different numbers, deliberately both shown. The geometry's excursion is
         fixed by the shape calibration (field.calib.provenance.inputs), which was solved
         against an earlier development contract; the contract now driving the volumes
         states its own strain. Printing only one of them would misreport the scene. */
      `long.short geom ${fmt(field.calib?.provenance?.inputs?.longitudinal_shortening_pct, 1)} %` +
        ` / contract ${ts.derived.longitudinal_shortening_pct.toFixed(1)} %`,
      `wall+      geom ${fmt(field.calib?.provenance?.inputs?.wall_thickening_pct, 1)} %` +
        ` / contract ${ts.derived.wall_thickening_pct.toFixed(1)} %`,
      `fps ${this.fps.toFixed(1)}   tri ${(extra.triangles / 1000).toFixed(0)}k   ` +
        `part ${(extra.particles / 1000).toFixed(0)}k`
    ];
    for (const r of rows) { ctx.fillText(r, x0, y); y += 14; }

    // pressure-volume loop, bottom-left
    this._pvLoop(ts, tMs, 16, H - 186, 170, 170);
  }

  _pvLoop (ts, tMs, x, y, w, h) {
    const ctx = this.ctx;
    ctx.save();
    ctx.fillStyle = 'rgba(8,11,18,0.62)';
    ctx.strokeStyle = 'rgba(120,140,170,0.25)';
    ctx.beginPath(); ctx.roundRect(x - 8, y - 8, w + 16, h + 26, 10);
    ctx.fill(); ctx.stroke();

    const V = ts.sig.V_LV_mL; const P = ts.sig.P_LV_mmHg;
    const vr = ts.vRange.V_LV_mL;
    const vlo = vr.min - 8; const vhi = vr.max + 8;
    let plo = 0; let phi = 0;
    for (let i = 0; i < P.length; i++) phi = Math.max(phi, P[i]);
    phi *= 1.08;

    ctx.beginPath();
    const step = Math.max(1, Math.floor(P.length / 1400));
    for (let i = 0; i < P.length; i += step) {
      const xx = x + ((V[i] - vlo) / (vhi - vlo)) * w;
      const yy = y + h - ((P[i] - plo) / (phi - plo)) * h;
      if (i === 0) ctx.moveTo(xx, yy); else ctx.lineTo(xx, yy);
    }
    ctx.strokeStyle = 'rgba(255,209,102,0.45)';
    ctx.lineWidth = 1.1; ctx.stroke();

    const vNow = ts.sample('V_LV_mL', tMs); const pNow = ts.sample('P_LV_mmHg', tMs);
    const cx = x + ((vNow - vlo) / (vhi - vlo)) * w;
    const cy = y + h - ((pNow - plo) / (phi - plo)) * h;
    ctx.beginPath(); ctx.arc(cx, cy, 4, 0, Math.PI * 2);
    ctx.fillStyle = '#ff6b6b'; ctx.fill();

    ctx.fillStyle = 'rgba(230,236,245,0.82)';
    ctx.font = '10px ui-monospace, Menlo, monospace';
    ctx.fillText('LV pressure-volume loop', x, y + h + 14);
    ctx.restore();
  }
}

function fmt (v, d) { return Number.isFinite(v) ? v.toFixed(d) : 'n/a'; }

function wrap (s, n) {
  const out = []; let cur = '';
  for (const word of s.split(' ')) {
    if ((cur + ' ' + word).trim().length > n) { out.push(cur.trim()); cur = word; } else cur += ' ' + word;
  }
  if (cur.trim()) out.push(cur.trim());
  return out;
}
