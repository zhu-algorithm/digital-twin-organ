/**
 * hud.js — in-canvas overlay.
 *
 * The heart scene shipped a HUD laid out for a desktop viewport: on a phone its
 * right-hand half was cut off outside the viewport, and because it is painted
 * inside a canvas, DOM bounding-box checks reported zero overflow. Nothing short
 * of looking at a screenshot caught it.
 *
 * So this HUD does two things differently:
 *   1. Layout is a function of the viewport: below NARROW_PX it switches to a
 *      single stacked column, drops the per-tree table to the two headline rows
 *      and shrinks the type.
 *   2. Every element it draws is recorded in `this.boxes` as a rectangle in CSS
 *      pixels, so an automated check can assert that nothing lies outside the
 *      viewport, which is the failure that actually happened. `overflow()`
 *      returns the offending boxes.
 */

export const NARROW_PX = 520;

/** Round a length in mm down to 1/2/5 x 10^n so the bar reads cleanly. */
function niceLength (mm) {
  const e = Math.floor(Math.log10(Math.max(mm, 1e-9)));
  const base = Math.pow(10, e);
  const f = mm / base;
  const m = f >= 5 ? 5 : (f >= 2 ? 2 : 1);
  return m * base;
}

/** mm -> the unit a reader expects at that scale. */
function formatLength (mm) {
  if (mm >= 10) return (mm / 10).toPrecision(2).replace(/\.0+$/, '') + ' cm';
  if (mm >= 1) return mm.toPrecision(2).replace(/\.0+$/, '') + ' mm';
  if (mm >= 0.001) return (mm * 1000).toPrecision(2).replace(/\.0+$/, '') + ' \u00b5m';
  return (mm * 1e6).toPrecision(2) + ' nm';
}

export class HUD {
  constructor (canvas) {
    this.c = canvas;
    this.ctx = canvas.getContext('2d');
    this.w = 0; this.h = 0; this.dpr = 1;
    this.fps = 0;
    this._frames = 0; this._t0 = performance.now();
    this.boxes = [];
    this.visible = true;
  }

  resize (w, h, dpr) {
    this.w = w; this.h = h; this.dpr = dpr;
    this.c.width = Math.max(1, Math.floor(w * dpr));
    this.c.height = Math.max(1, Math.floor(h * dpr));
    this.c.style.width = w + 'px';
    this.c.style.height = h + 'px';
    this.ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    this.narrow = w < NARROW_PX;
  }

  tickFps (now) {
    this._frames++;
    if (now - this._t0 >= 500) {
      this.fps = this._frames * 1000 / (now - this._t0);
      this._frames = 0; this._t0 = now;
    }
  }

  _panel (x, y, w, h, label) {
    const g = this.ctx;
    g.fillStyle = 'rgba(8,11,18,0.70)';
    g.strokeStyle = 'rgba(120,140,170,0.28)';
    g.lineWidth = 1;
    g.beginPath();
    const r = 8;
    g.moveTo(x + r, y); g.arcTo(x + w, y, x + w, y + h, r);
    g.arcTo(x + w, y + h, x, y + h, r); g.arcTo(x, y + h, x, y, r);
    g.arcTo(x, y, x + w, y, r); g.closePath();
    g.fill(); g.stroke();
    this.boxes.push({ label, x, y, w, h });
  }

  _text (s, x, y, color, size, label, align = 'left') {
    const g = this.ctx;
    g.font = `${size}px ui-monospace, Menlo, monospace`;
    g.fillStyle = color;
    g.textAlign = align;
    g.textBaseline = 'alphabetic';
    g.fillText(s, x, y);
    const m = g.measureText(s);
    const bx = align === 'right' ? x - m.width : (align === 'center' ? x - m.width / 2 : x);
    this.boxes.push({ label: label || s.slice(0, 18), x: bx, y: y - size, w: m.width, h: size * 1.25 });
  }

  /**
   * @param {object} st  scene state supplied by main.js
   */
  draw (st) {
    const g = this.ctx;
    g.clearRect(0, 0, this.w, this.h);
    this.boxes = [];
    if (!this.visible) return;

    const narrow = this.narrow;
    const fs = narrow ? 10 : 11.5;
    const lh = fs * 1.55;
    const pad = narrow ? 8 : 12;
    const margin = narrow ? 8 : 14;

    // ------------------------------------------------- respiration strip (top)
    const stripH = narrow ? 40 : 46;
    const stripW = Math.min(this.w - margin * 2, narrow ? this.w - margin * 2 : 330);
    this._panel(margin, margin, stripW, stripH, 'respiration-panel');
    const u = st.respiration.phase01;
    const inspFrac = 1 / (1 + 1.5);
    const bx = margin + pad; const bw = stripW - pad * 2;
    const barY = margin + stripH - pad - 4;
    this._text(
      `${st.respiration.breaths_per_min.toFixed(2)}/min  ${u < inspFrac ? 'insp' : 'exp'}  ` +
      `${st.respiration.displacement_mm.toFixed(1)} mm`,
      bx, margin + pad + fs, '#d8e2f0', fs, 'respiration-text');
    g.fillStyle = 'rgba(120,150,190,0.22)';
    g.fillRect(bx, barY, bw, 4);
    g.fillStyle = u < inspFrac ? '#6fd3a0' : '#8fb6e6';
    g.fillRect(bx, barY, bw * u, 4);
    this.boxes.push({ label: 'respiration-bar', x: bx, y: barY, w: bw, h: 4 });

    // ------------------------------------------------- flow panel (bottom-left)
    const rows = narrow
      ? [['portal vein', st.flow.portal], ['hepatic vein', st.flow.hepatic]]
      : [['portal vein', st.flow.portal], ['hepatic artery', st.flow.arterial],
         ['hepatic vein', st.flow.hepatic], ['bile duct', st.flow.bile]];
    const panelH = pad * 2 + lh * (rows.length + 1.6);
    const panelW = Math.min(this.w - margin * 2, narrow ? this.w - margin * 2 : 290);
    const py = this.h - margin - panelH;
    this._panel(margin, py, panelW, panelH, 'flow-panel');
    this._text(st.view === 'lobule' ? 'LOBULE — velocities' : 'FLOW — v = Q/A',
               margin + pad, py + pad + fs, '#9fb4cf', fs, 'flow-title');
    rows.forEach((r, i) => {
      const y = py + pad + fs + lh * (i + 1.3);
      this._text(r[0], margin + pad, y, r[1].color, fs, 'flow-row-' + i);
      this._text(r[1].text, margin + panelW - pad, y, '#e6edf7', fs, 'flow-val-' + i, 'right');
    });

    // ------------------------------------------------- stats (top-right / stacked)
    const lines = [
      `${st.fps.toFixed(0)} fps`,
      `${(st.triangles / 1000).toFixed(0)}k tri`,
      `${st.particlesDrawn} part`
    ];
    if (st.view === 'lobule') lines.push(`${st.magnification}x`);
    if (narrow) {
      const sw = Math.min(this.w - margin * 2, 200);
      const sh = pad * 2 + fs * 1.3;
      const sy = margin + stripH + 6;
      this._panel(margin, sy, sw, sh, 'stats-panel');
      this._text(lines.join('  '), margin + pad, sy + pad + fs, '#9fb4cf', fs, 'stats-text');
    } else {
      const sw = 128;
      const sh = pad * 2 + lh * lines.length;
      const sx = this.w - margin - sw;
      this._panel(sx, margin, sw, sh, 'stats-panel');
      lines.forEach((l, i) => this._text(
        l, sx + sw - pad, margin + pad + fs + lh * i, '#9fb4cf', fs, 'stats-' + i, 'right'));
    }

    // ------------------------------------------------- scale bar + tier
    // The bar length is computed from the camera every frame (mm per CSS pixel
    // at the orbit target), so it is an independent readout of the real size on
    // screen. The magnification printed beside it is the tier's DECLARED scene
    // units per millimetre of tissue. If the declared magnification were wrong
    // the two would disagree, which is the point of showing both.
    if (st.mmPerPixel && isFinite(st.mmPerPixel)) {
      const targetPx = narrow ? 86 : 130;
      const rawMm = st.mmPerPixel * targetPx;
      const nice = niceLength(rawMm);
      const barPx = nice / st.mmPerPixel;
      const bx2 = margin;
      const panelTopY = this.h - margin - (pad * 2 + lh * (rows.length + 1.6));
      const by2 = panelTopY - (narrow ? 30 : 34);
      const bw2 = Math.min(barPx, this.w - margin * 2);
      g.strokeStyle = '#dce6f3'; g.lineWidth = 2;
      g.beginPath();
      g.moveTo(bx2, by2 + 10); g.lineTo(bx2 + bw2, by2 + 10);
      g.moveTo(bx2, by2 + 5); g.lineTo(bx2, by2 + 15);
      g.moveTo(bx2 + bw2, by2 + 5); g.lineTo(bx2 + bw2, by2 + 15);
      g.stroke();
      this.boxes.push({ label: 'scalebar', x: bx2, y: by2 + 5, w: bw2, h: 10 });
      this._text(formatLength(nice), bx2, by2 + 2, '#dce6f3', fs, 'scalebar-label');
      const tierTxt = (st.tier ? st.tier.zh + ' · ' : '')
        + (st.magnification >= 1000 ? (st.magnification / 1000).toFixed(0) + '000×'
                                    : st.magnification.toFixed(0) + '×');
      this._text(tierTxt, bx2 + bw2 + 8, by2 + 14, '#9fb4cf', fs, 'tier-label');
    }

    // ------------------------------------------------- provenance footer
    if (!narrow && st.note) {
      this._text(st.note, this.w - margin, this.h - margin - 2,
                 'rgba(150,166,188,0.85)', 10, 'note', 'right');
    }
  }

  /** Boxes that stick out of the viewport — the failure this HUD exists to avoid. */
  overflow (tol = 0.5) {
    return this.boxes.filter(b =>
      b.x < -tol || b.y < -tol || b.x + b.w > this.w + tol || b.y + b.h > this.h + tol);
  }
}
