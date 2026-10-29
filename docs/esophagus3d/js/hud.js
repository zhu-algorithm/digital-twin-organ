/**
 * hud.js — in-canvas overlay for the oesophagus scene.
 *
 * Carried over from this project's liver HUD, which exists because the heart
 * page shipped a HUD laid out for a desktop viewport: on a phone its right-hand
 * half fell outside the viewport, and because it is painted inside a canvas,
 * DOM bounding-box checks reported zero overflow. Only a screenshot caught it.
 *
 * So, as there:
 *   1. layout is a function of the viewport — below NARROW_PX it switches to a
 *      stacked single-column form, drops the table to its headline rows and
 *      shrinks the type;
 *   2. EVERY element drawn is recorded in `this.boxes` as a rectangle in CSS
 *      pixels, so an automated check can assert that nothing lies outside the
 *      viewport. `overflow()` returns the offending boxes.
 *
 * Projected 3-D labels, carried over from the kidney page: the anatomy is only
 * legible if the three measured constrictions, both sphincters, the GOJ and
 * each vascular territory are NAMED in place, so main.js hands the HUD a list
 * of world-space anchors and it projects them. Labels are recorded in
 * `this.boxes` like everything else, and are culled when they would collide or
 * fall off screen rather than being allowed to overflow.
 *
 * Changed for the oesophagus: the kidney's respiration strip is replaced by a
 * PERISTALSIS strip showing the wave front's position along the measured
 * centreline, because that is this organ's signature motion. The cardiac strip
 * stays — the four arterial territories really are pulsatile.
 */

export const NARROW_PX = 520;

export class HUD {
  constructor (canvas) {
    this.c = canvas;
    this.ctx = canvas.getContext('2d');
    this.w = 0; this.h = 0; this.dpr = 1;
    this.fps = 0;
    this._frames = 0; this._t0 = performance.now();
    this.boxes = [];
    this.visible = true;
    this.labelsDrawn = 0;
    this.labelsCulled = 0;
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
    g.fillStyle = 'rgba(8,11,18,0.72)';
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
    return m.width;
  }

  /**
   * Project world-space anatomical labels onto the canvas.
   * @param {Array} labels [{ pos: THREE.Vector3, text, color, pri }]
   * @param {THREE.Camera} camera
   * @param {object} reserved rectangles already occupied by panels
   */
  drawLabels (labels, camera, THREE) {
    if (!labels || !labels.length) return;
    const g = this.ctx;
    const size = this.narrow ? 9.5 : 11;
    const padx = 5; const pady = 3;
    const placed = this.boxes.slice();
    const v = new THREE.Vector3();
    const cands = [];
    for (const L of labels) {
      v.copy(L.pos).project(camera);
      if (v.z < -1 || v.z > 1) continue;
      const x = (v.x * 0.5 + 0.5) * this.w;
      const y = (-v.y * 0.5 + 0.5) * this.h;
      cands.push({ ...L, x, y, depth: v.z });
    }
    cands.sort((a, b) => (b.pri || 0) - (a.pri || 0) || a.depth - b.depth);
    this.labelsDrawn = 0; this.labelsCulled = 0;
    const limit = this.narrow ? 5 : 14;
    for (const c of cands) {
      if (this.labelsDrawn >= limit) { this.labelsCulled++; continue; }
      g.font = `${size}px ui-monospace, Menlo, monospace`;
      const tw = g.measureText(c.text).width;
      const bw = tw + padx * 2; const bh = size + pady * 2;
      // keep the whole chip inside the viewport: this is the failure the HUD
      // exists to prevent, so a label that cannot fit is dropped, not clipped
      let bx = c.x + 9; let by = c.y - bh / 2;
      if (bx + bw > this.w - 4) bx = c.x - 9 - bw;
      bx = Math.max(4, Math.min(bx, this.w - bw - 4));
      by = Math.max(4, Math.min(by, this.h - bh - 4));
      if (bw > this.w - 8 || bh > this.h - 8) { this.labelsCulled++; continue; }
      const hit = placed.some(p => !(bx + bw < p.x - 2 || bx > p.x + p.w + 2 ||
                                     by + bh < p.y - 2 || by > p.y + p.h + 2));
      if (hit) { this.labelsCulled++; continue; }
      g.fillStyle = 'rgba(6,9,15,0.80)';
      g.strokeStyle = c.color || 'rgba(140,160,190,0.65)';
      g.lineWidth = 1;
      g.beginPath();
      const r = 4;
      g.moveTo(bx + r, by); g.arcTo(bx + bw, by, bx + bw, by + bh, r);
      g.arcTo(bx + bw, by + bh, bx, by + bh, r); g.arcTo(bx, by + bh, bx, by, r);
      g.arcTo(bx, by, bx + bw, by, r); g.closePath();
      g.fill(); g.stroke();
      // leader line to the anchor
      g.beginPath();
      g.moveTo(c.x, c.y);
      g.lineTo(bx < c.x ? bx + bw : bx, by + bh / 2);
      g.strokeStyle = 'rgba(150,170,200,0.45)';
      g.stroke();
      g.fillStyle = c.color || '#d8e2f0';
      g.textAlign = 'left'; g.textBaseline = 'alphabetic';
      g.fillText(c.text, bx + padx, by + bh - pady - 1);
      const box = { label: 'anat-' + c.text.slice(0, 16), x: bx, y: by, w: bw, h: bh };
      this.boxes.push(box); placed.push(box);
      this.labelsDrawn++;
    }
  }

  /** @param {object} st scene state supplied by main.js */
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

    // -------------------------------------- cardiac + peristalsis strip (top)
    const stripH = narrow ? 52 : 60;
    const stripW = Math.min(this.w - margin * 2, narrow ? this.w - margin * 2 : 350);
    this._panel(margin, margin, stripW, stripH, 'rhythm-panel');
    const bx = margin + pad; const bw = stripW - pad * 2;
    this._text(
      `${st.heart.bpm.toFixed(1)} bpm   pulse x${st.heart.pulse.toFixed(2)}`,
      bx, margin + pad + fs, '#ffb4a0', fs, 'cardiac-text');
    const cy = margin + pad + fs + 5;
    g.fillStyle = 'rgba(190,120,120,0.22)'; g.fillRect(bx, cy, bw, 3);
    g.fillStyle = '#e4705f'; g.fillRect(bx, cy, bw * st.heart.phase01, 3);
    this.boxes.push({ label: 'cardiac-bar', x: bx, y: cy, w: bw, h: 3 });
    const ry = cy + lh;
    // peristalsis strip: the front position along the MEASURED centreline.
    // The bar is the wave front's arc fraction, so a viewer can see directly
    // that the contraction is at one place at a time rather than everywhere.
    const wv = st.wave;
    this._text(
      wv.mode === 'idle'
        ? `蠕动 静息   v=${wv.declared_speed_mm_s} mm/s`
        : `${wv.mode === 'primary' ? '原发性' : '继发性'}蠕动  ` +
          `s=${Math.max(0, wv.front_now).toFixed(0)}/${wv.L.toFixed(0)} mm  ` +
          `v=${wv.declared_speed_mm_s} mm/s`,
      bx, ry, '#a9d8c0', fs, 'peristalsis-text');
    const by2 = ry + 5;
    const frac = wv.mode === 'idle' ? 0 : Math.max(0, Math.min(1, wv.front_now / wv.L));
    g.fillStyle = 'rgba(120,150,190,0.22)'; g.fillRect(bx, by2, bw, 3);
    g.fillStyle = '#6fd3a0'; g.fillRect(bx, by2, bw * frac, 3);
    this.boxes.push({ label: 'peristalsis-bar', x: bx, y: by2, w: bw, h: 3 });

    // -------------------------------------------- flow panel (bottom-left)
    const rows = narrow ? st.flowRows.slice(0, 3) : st.flowRows;
    const panelH = pad * 2 + lh * (rows.length + 1.6);
    const panelW = Math.min(this.w - margin * 2, narrow ? this.w - margin * 2 : 320);
    const py = this.h - margin - panelH;
    this._panel(margin, py, panelW, panelH, 'flow-panel');
    this._text(st.flowTitle, margin + pad, py + pad + fs, '#9fb4cf', fs, 'flow-title');
    rows.forEach((r, i) => {
      const y = py + pad + fs + lh * (i + 1.3);
      const nameW = this._text(r.name, margin + pad, y, r.color, fs, 'flow-row-' + i);
      // The right-aligned value must not run back past the label. A long value
      // (an error string, a long cell-type name) otherwise extends LEFT out of
      // the panel and out of the viewport — the first headless audit caught
      // exactly that in the molecular view, where the value was a fetch error
      // message 400 px wide in a 390 px viewport. Values are ellipsised to the
      // space actually available.
      g.font = `${fs}px ui-monospace, Menlo, monospace`;
      const avail = panelW - pad * 2 - nameW - 10;
      let val = String(r.value);
      if (g.measureText(val).width > avail) {
        while (val.length > 1 && g.measureText(val + '…').width > avail) val = val.slice(0, -1);
        val += '…';
      }
      this._text(val, margin + panelW - pad, y, '#e6edf7', fs, 'flow-val-' + i, 'right');
    });

    // ------------------------------- stats (top-right on wide, stacked narrow)
    const lines = [`${st.fps.toFixed(0)} fps`, `${(st.triangles / 1000).toFixed(0)}k tri`,
                   `${st.particlesDrawn} part`];
    if (st.magnification) lines.push(`${st.magnification}x`);
    if (narrow) {
      const sw = Math.min(this.w - margin * 2, 230);
      const sh = pad * 2 + fs * 1.3;
      const sy = margin + stripH + 6;
      this._panel(margin, sy, sw, sh, 'stats-panel');
      this._text(lines.join('  '), margin + pad, sy + pad + fs, '#9fb4cf', fs, 'stats-text');
    } else {
      const sw = 132;
      const sh = pad * 2 + lh * lines.length;
      const sx = this.w - margin - sw;
      this._panel(sx, margin, sw, sh, 'stats-panel');
      lines.forEach((l, i) => this._text(
        l, sx + sw - pad, margin + pad + fs + lh * i, '#9fb4cf', fs, 'stats-' + i, 'right'));
    }

    // -------------------------------------------------- provenance footer
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
