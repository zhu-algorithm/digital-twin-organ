/**
 * section.js — an interactive anatomical cut through the segmented liver.
 *
 * What a section view has to get right, and what this module does about it:
 *
 *  1. SOLID CUT FACE (capping). A clipped triangle mesh is hollow: without a
 *     cap you look straight through the cut into the inside of the shell and
 *     see the back of the far wall. The cut face here is built with the
 *     two-pass stencil technique: the clipped solid's BACK faces increment the
 *     stencil, its FRONT faces decrement it, and a quad lying on the plane is
 *     then painted only where the stencil is non-zero, i.e. exactly where the
 *     plane is inside the solid. Nothing is painted where the plane is outside
 *     the liver, so the cap is the true cross-section, not a disc.
 *
 *  2. VESSEL LUMENS AT TRUE CALIBRE. The four tubular trees are open-ended
 *     swept rings (zero wall thickness), so the stencil parity trick produces
 *     no cap for them — a ray through a lumen crosses one front face and one
 *     back face and nets to zero. Their openings are therefore computed
 *     ANALYTICALLY: every tree segment that straddles the plane contributes one
 *     opening, placed at the interpolated crossing point with the radius
 *     interpolated between the two node radii of liver_vessels.json. The
 *     rendered opening diameter is that radius x 2, unscaled, so the calibres on
 *     the cut face are the calibres in the Murray audit by construction.
 *     Each opening is drawn as an opaque wall annulus plus a translucent lumen
 *     disc, so the flow particles inside the kept half of the vessel remain
 *     visible through the opening, and carries a direction mark taken from the
 *     sign of (segment direction . plane normal) — into-page or out-of-page is
 *     read from the tree's own geometry, not chosen.
 *
 *  3. COUINAUD ON THE CUT FACE. The cap is tinted by segment using the SAME
 *     classifier that baked _field0.y onto the parenchyma (nearest portal
 *     terminal, i.e. the Selle-style Voronoi partition), evaluated on a
 *     256 x 256 grid over the cut face on the CPU whenever the plane moves, and
 *     uploaded as a texture. Segment boundaries are drawn where neighbouring
 *     grid cells disagree.
 *
 *  4. INSIDE/OUTSIDE AND AREA. The same pass slices the parenchyma triangle
 *     soup with the plane, assembles the intersection loops and gets the
 *     cross-section area by the shoelace formula. That area is an independent
 *     CPU computation; the cap audit measures the rendered cap in pixels with a
 *     known orthographic scale and compares the two. A missing or wrongly sized
 *     cap cannot pass that comparison.
 *
 *  5. RESPIRATION. The plane is defined in the breathing body's local frame and
 *     the cap geometry is parented to it, so the cut face travels with the
 *     organ exactly as the rest of the scene does. The audit records the cut
 *     face centroid's world displacement over a cycle.
 */

import * as THREE from 'three';

export const SECTION_PRESETS = {
  sagittal: { normal: [1, 0, 0], label: '矢状面' },
  coronal: { normal: [0, 0, 1], label: '冠状面' },
  axial: { normal: [0, 1, 0], label: '轴位面' }
};

const CAP_OFFSET_MM = 0.0;       // cap sits exactly on the plane; z-fight handled
                                 // by polygonOffset, not by displacing geometry
const GRID = 256;                // cut-face classification grid

/** Lumen tints: same hues as the tubes, lifted so they read on the cut face. */
const LUMEN_COLOR = {
  portal_vein: 0x4f8fe0,
  hepatic_vein: 0x8a68d0,
  hepatic_artery: 0xe8453f,
  bile_duct: 0xaacb45
};

/* ------------------------------------------------------------------ slicing */

/**
 * Intersect an indexed triangle mesh with a plane and return the cut loops in
 * the plane's 2-D frame.
 * @returns {{segments: Float32Array, area_mm2: number, loops: number}}
 */
export function sliceGeometry (geometry, planeNormal, planeConst, u, v, origin) {
  const pos = geometry.getAttribute('position');
  const index = geometry.getIndex();
  const n = planeNormal;
  const nTri = index ? index.count / 3 : pos.count / 3;
  const segs = [];
  const a = new THREE.Vector3(); const b = new THREE.Vector3(); const c = new THREE.Vector3();
  const p = new THREE.Vector3();
  const tri = [a, b, c];
  const d = [0, 0, 0];
  for (let t = 0; t < nTri; t++) {
    for (let k = 0; k < 3; k++) {
      const i = index ? index.getX(t * 3 + k) : t * 3 + k;
      tri[k].fromBufferAttribute(pos, i);
      d[k] = n.dot(tri[k]) + planeConst;
    }
    // triangle must straddle the plane
    const npos = (d[0] > 0) + (d[1] > 0) + (d[2] > 0);
    if (npos === 0 || npos === 3) continue;
    // Each triangle contributes one segment, ORIENTED: it runs from the edge
    // crossing where the signed distance goes positive -> negative to the one
    // where it goes negative -> positive. Without that the segments come out
    // with arbitrary direction and the shoelace sum below partly cancels (on
    // this mesh it under-reports the sagittal cut as 1676 mm2 instead of 7775).
    const pts = [];
    for (let k = 0; k < 3; k++) {
      const k2 = (k + 1) % 3;
      if ((d[k] > 0) === (d[k2] > 0)) continue;
      const tt = d[k] / (d[k] - d[k2]);
      p.copy(tri[k]).lerp(tri[k2], tt).sub(origin);
      pts.push([u.dot(p), v.dot(p), d[k] > 0]);
    }
    if (pts.length === 2) {
      const a0 = pts[0][2] ? pts[0] : pts[1];
      const b0 = pts[0][2] ? pts[1] : pts[0];
      segs.push(a0[0], a0[1], b0[0], b0[1]);
    }
  }
  // Shoelace over the oriented boundary segments. Because every segment runs the
  // same way round the solid, the cross products accumulate instead of
  // cancelling, and the loops never have to be stitched into rings.
  let a2 = 0;
  for (let i = 0; i < segs.length; i += 4) {
    a2 += segs[i] * segs[i + 3] - segs[i + 2] * segs[i + 1];
  }
  return {
    segments: new Float32Array(segs),
    area_mm2: Math.abs(a2) * 0.5,
    loops: segs.length / 4
  };
}

/** Even-odd point-in-polygon test against a soup of unordered cut segments. */
function insideSegments (segs, x, y) {
  let cross = 0;
  for (let i = 0; i < segs.length; i += 4) {
    const x1 = segs[i]; const y1 = segs[i + 1];
    const x2 = segs[i + 2]; const y2 = segs[i + 3];
    if ((y1 > y) !== (y2 > y)) {
      const xi = x1 + (y - y1) / (y2 - y1) * (x2 - x1);
      if (xi > x) cross++;
    }
  }
  return (cross & 1) === 1;
}

/* ------------------------------------------------------------------- class */

export class SectionPlane {
  /**
   * @param {object} o
   * @param {THREE.Group} o.liverBody   the breathing group
   * @param {THREE.Mesh}  o.parenchyma  closed parenchymal shell
   * @param {object} o.solids           name -> closed Mesh to cap
   * @param {object} o.vessels          parsed liver_vessels.json
   * @param {object} o.couinaud         parsed couinaud.json
   * @param {Array}  o.clipMaterials    materials to attach the clipping plane to
   */
  constructor (o) {
    this.liverBody = o.liverBody;
    this.parenchyma = o.parenchyma;
    this.solids = o.solids || {};
    this.vessels = o.vessels;
    this.couinaud = o.couinaud;
    this.terminalsAsset = o.terminals || {};
    this.clipMaterials = o.clipMaterials || [];
    this.flowSystems = o.flowSystems || [];
    this.enabled = false;
    this.flip = 1;

    const box = new THREE.Box3().setFromBufferAttribute(
      this.parenchyma.geometry.getAttribute('position'));
    this.box = box;
    this.centre = box.getCenter(new THREE.Vector3());
    this.extent = box.getSize(new THREE.Vector3());
    this.diag = this.extent.length();

    // local-frame plane: normal + constant (keep dot(p,n)+c <= 0)
    this.normal = new THREE.Vector3(1, 0, 0);
    this.offset01 = 0.5;
    this.plane = new THREE.Plane(this.normal.clone(), 0);
    this.worldPlane = new THREE.Plane();
    this.u = new THREE.Vector3();
    this.v = new THREE.Vector3();
    this.origin = new THREE.Vector3();

    this.group = new THREE.Group();          // everything that lives on the cut
    this.group.name = 'section_cut';
    this.liverBody.add(this.group);

    this._buildStencilGroups();
    this._buildSegmentTexture();
    this._buildCapQuads();
    this._buildLumens();
    this.setPreset('sagittal');
    this.setEnabled(false);
  }

  /* ------------------------------------------------------------- stencils */

  /**
   * Stencil-pass material. The clip test is written into the shader against the
   * BODY-LOCAL position, the same test the parenchyma and vessel materials
   * already use, rather than relying on three.js local clipping planes: the two
   * must agree exactly or the stencil region and the visible cut would describe
   * different solids, and a single shared uniform makes that impossible.
   */
  _stencilMaterial (side, op) {
    const m = new THREE.ShaderMaterial({
      uniforms: { uClipPlane: { value: new THREE.Vector4(1, 0, 0, 0) },
                  uClipOn: { value: 1.0 } },
      vertexShader: `
        varying vec3 vLocal;
        void main () {
          vLocal = position;
          gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
        }`,
      fragmentShader: `
        precision highp float;
        varying vec3 vLocal;
        uniform vec4 uClipPlane;
        uniform float uClipOn;
        void main () {
          if (uClipOn > 0.5 && dot(vLocal, uClipPlane.xyz) > uClipPlane.w) discard;
          gl_FragColor = vec4(1.0);
        }`,
      side,
      colorWrite: false,
      depthWrite: false,
      depthTest: false,
      stencilWrite: true,
      stencilFunc: THREE.AlwaysStencilFunc,
      stencilFail: op,
      stencilZFail: op,
      stencilZPass: op
    });
    this.stencilMats.push(m);
    return m;
  }

  _buildStencilGroups () {
    this.stencilGroups = {};
    this.capOrder = [];
    this.stencilMats = [];
    let order = 1;
    for (const [name, mesh] of Object.entries(this.solids)) {
      const g = new THREE.Group();
      g.name = 'stencil_' + name;
      const back = this._stencilMaterial(THREE.BackSide, THREE.IncrementWrapStencilOp);
      const front = this._stencilMaterial(THREE.FrontSide, THREE.DecrementWrapStencilOp);
      const mb = new THREE.Mesh(mesh.geometry, back);
      const mf = new THREE.Mesh(mesh.geometry, front);
      mb.renderOrder = order; mf.renderOrder = order + 0.1;
      mb.frustumCulled = false; mf.frustumCulled = false;
      g.add(mb, mf);
      this.group.add(g);
      this.stencilGroups[name] = { group: g, order };
      this.capOrder.push(name);
      order += 1;
    }
  }

  /* ------------------------------------------------- Couinaud on the face */

  _buildSegmentTexture () {
    // Portal terminals carrying the geometric Couinaud label, and the hepatic
    // venous terminals: exactly the two point sets build_liver_fields.py used to
    // bake _field0.y (segment id) and _field0.z (zone coordinate). Supplied by
    // tools/build_section_assets.py as couinaud_terminals.json so the cut face
    // and the capsule cannot disagree.
    const T = this.terminalsAsset || {};
    this.terminals = T.portal_terminals || [];
    this.hvTerminals = T.hepatic_vein_terminals || [];
    this.hasSegmentIds = this.terminals.length > 0;
    const data = new Uint8Array(GRID * GRID * 4);
    this.segTex = new THREE.DataTexture(data, GRID, GRID, THREE.RGBAFormat);
    this.segTex.minFilter = THREE.LinearFilter;
    this.segTex.magFilter = THREE.LinearFilter;
    this.segTex.needsUpdate = true;
    // nearest-neighbour copy used for the segment id (must not be interpolated)
    this.idTex = new THREE.DataTexture(data, GRID, GRID, THREE.RGBAFormat);
    this.idTex.minFilter = THREE.NearestFilter;
    this.idTex.magFilter = THREE.NearestFilter;
    this.idTex.needsUpdate = true;
  }

  /**
   * Recompute cut-face mask, segment id, zone coordinate and the analytic area.
   * R = segment id, G = zone01 * 255, B = distance-to-nearest-portal (mm, /2),
   * A = inside mask.
   */
  _recomputeFace () {
    const t0 = (typeof performance !== 'undefined') ? performance.now() : 0;
    const slice = this.slice || sliceGeometry(
      this.parenchyma.geometry, this.plane.normal, this.plane.constant,
      this.u, this.v, this.origin);
    const half = this.diag * 0.55;
    const data = this.segTex.image.data;
    const segCount = new Float64Array(9);
    const zoneHist = new Float64Array(10);
    let inside = 0;
    const P = new THREE.Vector3();
    const terms = this.terminals; const hvt = this.hvTerminals;
    for (let j = 0; j < GRID; j++) {
      const fy = (j + 0.5) / GRID * 2 - 1;
      for (let i = 0; i < GRID; i++) {
        const fx = (i + 0.5) / GRID * 2 - 1;
        const o = (j * GRID + i) * 4;
        const x = fx * half; const y = fy * half;
        if (!insideSegments(slice.segments, x, y)) {
          data[o] = 0; data[o + 1] = 0; data[o + 2] = 0; data[o + 3] = 0;
          continue;
        }
        inside++;
        P.copy(this.origin).addScaledVector(this.u, x).addScaledVector(this.v, y);
        let sid = 0; let bd = Infinity;
        for (let k = 0; k < terms.length; k++) {
          const t = terms[k];
          const dx = t[0] - P.x; const dy = t[1] - P.y; const dz = t[2] - P.z;
          const d2 = dx * dx + dy * dy + dz * dz;
          if (d2 < bd) { bd = d2; sid = t[3]; }
        }
        let bh = Infinity;
        for (let k = 0; k < hvt.length; k++) {
          const t = hvt[k];
          const dx = t[0] - P.x; const dy = t[1] - P.y; const dz = t[2] - P.z;
          const d2 = dx * dx + dy * dy + dz * dz;
          if (d2 < bh) bh = d2;
        }
        const dp = Math.sqrt(bd); const dh = Math.sqrt(bh);
        const zone = dp / Math.max(dp + dh, 1e-6);
        segCount[Math.max(0, Math.min(8, sid))]++;
        zoneHist[Math.min(9, Math.floor(zone * 10))]++;
        data[o] = sid;
        data[o + 1] = Math.round(Math.min(zone, 1) * 255);
        data[o + 2] = Math.round(Math.min(dp / 2, 255));
        data[o + 3] = 255;
      }
    }
    this.segTex.needsUpdate = true;
    this.idTex.image.data = data;
    this.idTex.needsUpdate = true;
    const cell = (2 * half / GRID) ** 2;
    this.faceStats = {
      area_mm2_mesh_slice: Number(slice.area_mm2.toFixed(2)),
      area_mm2_grid: Number((inside * cell).toFixed(2)),
      cut_segments: Array.from(segCount.slice(1))
        .map((cnt, i) => ({ segment: i + 1, area_mm2: Number((cnt * cell).toFixed(1)) }))
        .filter(s => s.area_mm2 > 0),
      zone01_decile_area_mm2: Array.from(zoneHist).map(cnt => Number((cnt * cell).toFixed(1))),
      grid: GRID,
      grid_cell_mm: Number((2 * half / GRID).toFixed(4)),
      slice_edge_count: slice.loops,
      cpu_ms: Number((((typeof performance !== 'undefined') ? performance.now() : 0) - t0).toFixed(1))
    };
    this.faceDirty = false;
  }

  /* -------------------------------------------------------------- cap quad */

  _buildCapQuads () {
    const uniforms = {
      uSegTex: { value: this.segTex },
      uIdTex: { value: this.idTex },
      uHalf: { value: this.diag * 0.55 },
      uSegMix: { value: 0.0 },
      uSegColors: { value: CAP_SEG_COLORS.map(c => new THREE.Color(c)) },
      uBoundary: { value: 1.0 },
      uZoneMix: { value: 0.0 },
      uZoneRatio: { value: 1.0 },
      uFuncMix: { value: 0.0 },
      uFuncMode: { value: 0 },
      uExtraction: { value: 0.5 },
      uPO2: { value: new THREE.Vector2(62.5, 32.5) },
      uTime: { value: 0.0 }
    };
    this.capUniforms = uniforms;
    const mat = new THREE.ShaderMaterial({
      uniforms,
      side: THREE.DoubleSide,
      stencilWrite: true,
      stencilRef: 0,
      stencilFunc: THREE.NotEqualStencilFunc,
      stencilFail: THREE.ReplaceStencilOp,
      stencilZFail: THREE.ReplaceStencilOp,
      stencilZPass: THREE.ReplaceStencilOp,
      polygonOffset: true,
      polygonOffsetFactor: -4,
      polygonOffsetUnits: -4,
      vertexShader: /* glsl */`
        varying vec2 vFace;
        void main () {
          vFace = position.xy;
          gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
        }`,
      fragmentShader: /* glsl */`
        precision highp float;
        uniform sampler2D uSegTex;
        uniform sampler2D uIdTex;
        uniform float uHalf, uSegMix, uBoundary, uZoneMix, uZoneRatio, uFuncMix, uTime;
        uniform float uExtraction;
        uniform vec2 uPO2;
        uniform int uFuncMode;
        uniform vec3 uSegColors[9];
        varying vec2 vFace;
        float h21 (vec2 p) {
          p = fract(p * vec2(123.34, 456.21)); p += dot(p, p + 45.32);
          return fract(p.x * p.y);
        }
        vec2 h22 (vec2 p) {
          return fract(sin(vec2(dot(p, vec2(127.1, 311.7)), dot(p, vec2(269.5, 183.3)))) * 43758.5453);
        }
        float cells (vec2 p) {
          vec2 n = floor(p), f = fract(p); float f1 = 8.0, f2 = 8.0;
          for (int j = -1; j <= 1; j++) for (int i = -1; i <= 1; i++) {
            vec2 g = vec2(float(i), float(j));
            float d = length(g + h22(n + g) - f);
            if (d < f1) { f2 = f1; f1 = d; } else if (d < f2) { f2 = d; }
          }
          return f2 - f1;
        }
        void main () {
          vec2 uv = vFace / uHalf * 0.5 + 0.5;
          vec4 s = texture2D(uIdTex, uv);
          // Outside the measured cross-section nothing is painted at all, so the
          // cut face can never be a disc that overshoots the organ.
          if (s.a < 0.5) discard;
          float sid = s.r * 255.0;
          float zone = texture2D(uSegTex, uv).g;      // 0 periportal -> 1 pericentral
          // cut parenchyma: a darker, wetter red than the capsule, with the
          // ~1 mm lobular stipple that a fresh cut surface actually shows
          float stip = smoothstep(0.02, 0.28, cells(vFace * 1.05));
          vec3 base = vec3(0.42, 0.13, 0.11) * (0.78 + 0.42 * stip);
          base += vec3(0.05, 0.015, 0.01) * h21(vFace * 7.0);

          int si = int(clamp(sid, 0.0, 8.0) + 0.5);
          vec3 segc = uSegColors[si];
          // segment boundary: neighbouring grid cells disagree
          float px = 1.0 / float(${GRID});
          float nb = 0.0;
          nb += abs(texture2D(uIdTex, uv + vec2(px, 0.0)).r - s.r);
          nb += abs(texture2D(uIdTex, uv - vec2(px, 0.0)).r - s.r);
          nb += abs(texture2D(uIdTex, uv + vec2(0.0, px)).r - s.r);
          nb += abs(texture2D(uIdTex, uv - vec2(0.0, px)).r - s.r);
          float edge = step(0.0001, nb);

          // Couinaud tint first, the active function channel ON TOP of it: with
          // the order the other way round a segment tint at full strength hides
          // whichever channel the reader just switched on.
          base = mix(base, mix(base, segc, 0.68), uSegMix);

          // ---------------------------------------------- function channels
          // Every channel is a function of the SAME zone coordinate that the
          // capsule uses: zone01 = d(nearest portal terminal) / (d_portal +
          // d_hepatic-vein terminal), from this subject's own two vascular trees.
          if (uFuncMix > 0.001) {
            vec3 fc = base;
            if (uFuncMode == 0) {
              // oxygen: pO2 falls from the periportal to the pericentral end.
              // Endpoints come from the literature (uPO2), the SHAPE is linear in
              // the zone coordinate on purpose — nothing here measures the profile.
              float po2 = mix(uPO2.x, uPO2.y, zone);
              float t = clamp((po2 - uPO2.y) / max(uPO2.x - uPO2.y, 1e-3), 0.0, 1.0);
              fc = mix(vec3(0.16, 0.20, 0.52), vec3(0.95, 0.30, 0.25), 1.0 - t);
              fc = mix(fc, vec3(1.0), 0.10 * sin(6.2831 * (zone * 3.0 - uTime * 0.35)) + 0.10);
            } else if (uFuncMode == 1) {
              // metabolic zonation: zone 1 (periportal) vs zone 3 (pericentral)
              float z1 = 1.0 - smoothstep(0.25, 0.55, zone);
              float z3 = smoothstep(0.45, 0.75, zone);
              fc = vec3(0.20, 0.52, 0.86) * z1 + vec3(0.92, 0.55, 0.14) * z3
                 + vec3(0.30, 0.30, 0.34) * (1.0 - z1 - z3);
            } else if (uFuncMode == 2) {
              // per-gene measured zonal ramp; direction follows the measured
              // central/periportal ratio, so a reversed gene ramps the other way
              float z = uZoneRatio >= 1.0 ? zone : 1.0 - zone;
              float amp = clamp(abs(log(max(uZoneRatio, 1e-3)) / log(2.21)), 0.0, 1.0);
              fc = mix(vec3(0.13, 0.22, 0.46), vec3(0.95, 0.68, 0.16), z);
              fc = mix(vec3(0.42, 0.33, 0.30), fc, 0.25 + 0.75 * amp);
            } else if (uFuncMode == 3) {
              // substrate concentration falling along the sinusoid. SCHEMATIC:
              // the decay constant is the extraction-ratio slider, an INPUT. It
              // is not computed from this project's expression data.
              float cconc = exp(-zone * (-log(max(1.0 - uExtraction, 0.02))));
              fc = mix(vec3(0.10, 0.12, 0.18), vec3(0.45, 0.95, 0.70), cconc);
              float wave = 0.5 + 0.5 * sin(6.2831 * (zone * 4.0 - uTime * 0.5));
              fc *= 0.85 + 0.3 * wave;
            } else {
              // bile: canalicular flow runs the OTHER way, parenchyma -> porta
              float b = 1.0 - zone;
              fc = mix(vec3(0.12, 0.14, 0.10), vec3(0.78, 0.86, 0.26), b);
              float wave = 0.5 + 0.5 * sin(6.2831 * (b * 5.0 - uTime * 0.25));
              fc *= 0.8 + 0.35 * wave;
            }
            base = mix(base, mix(base * 0.45, fc, 0.9), uFuncMix);
          }

          base = mix(base, vec3(0.97, 0.95, 0.86), edge * uBoundary * 0.85);

          gl_FragColor = vec4(base, 1.0);
          #include <tonemapping_fragment>
          #include <colorspace_fragment>
        }`
    });
    this.capMat = mat;
    // Built at its true size, NOT a unit quad scaled by the model matrix: the
    // fragment shader reads position.xy as the face coordinate in millimetres, so
    // a unit quad would collapse the whole cut face onto one texel.
    const capSize = this.diag * 1.1;
    this.capSize = capSize;
    const geo = new THREE.PlaneGeometry(capSize, capSize);
    this.caps = {};
    for (const name of this.capOrder) {
      const m = name === 'parenchyma' ? mat : mat.clone();
      if (name !== 'parenchyma') {
        m.uniforms = THREE.UniformsUtils.clone(uniforms);
        m.uniforms.uSegTex.value = this.segTex;
        m.polygonOffsetFactor = -4 - this.stencilGroups[name].order;
      }
      const q = new THREE.Mesh(geo, m);
      q.name = 'cap_' + name;
      q.renderOrder = this.stencilGroups[name].order + 0.2;
      q.frustumCulled = false;
      q.onAfterRender = (renderer) => { renderer.clearStencil(); };
      this.group.add(q);
      this.caps[name] = q;
    }
  }

  /* ---------------------------------------------------------- vessel lumens */

  _buildLumens () {
    this.lumenGroup = new THREE.Group();
    this.lumenGroup.name = 'section_lumens';
    this.group.add(this.lumenGroup);
    this.lumenMeshes = {};
    const ringGeo = new THREE.RingGeometry(0.80, 1.0, 20);
    const discGeo = new THREE.CircleGeometry(1.0, 20);
    const MAXN = 600;
    for (const tree of this.vessels.trees) {
      const col = LUMEN_COLOR[tree.name] || 0x888888;
      const wallMat = new THREE.MeshBasicMaterial({
        color: col, side: THREE.DoubleSide, depthTest: false, transparent: true, opacity: 0.98
      });
      const lumMat = new THREE.MeshBasicMaterial({
        color: col, side: THREE.DoubleSide, depthTest: false,
        transparent: true, opacity: 0.30, blending: THREE.AdditiveBlending
      });
      const wall = new THREE.InstancedMesh(ringGeo, wallMat, MAXN);
      const lum = new THREE.InstancedMesh(discGeo, lumMat, MAXN);
      wall.renderOrder = 40; lum.renderOrder = 39;
      wall.frustumCulled = false; lum.frustumCulled = false;
      wall.count = 0; lum.count = 0;
      this.lumenGroup.add(lum, wall);
      this.lumenMeshes[tree.name] = { wall, lum, max: MAXN };
    }
  }

  /**
   * Recompute the openings every tree makes on the current plane.
   * Radius at the crossing is linearly interpolated between the two node radii
   * straight out of liver_vessels.json — no scaling anywhere in this path.
   */
  _recomputeLumens () {
    const n = this.plane.normal; const cst = this.plane.constant;
    const rep = {};
    const m = new THREE.Matrix4();
    const q = new THREE.Quaternion();
    const Z = new THREE.Vector3(0, 0, 1);
    const P = new THREE.Vector3();
    const A = new THREE.Vector3(); const B = new THREE.Vector3();
    q.setFromUnitVectors(Z, n);
    for (const tree of this.vessels.trees) {
      const holder = this.lumenMeshes[tree.name];
      const nodes = tree.nodes; const parent = tree.parent; const rad = tree.radius_mm;
      const diam = []; let k = 0; let into = 0; let out = 0;
      for (let i = 0; i < nodes.length; i++) {
        const p = parent[i];
        if (p < 0) continue;
        A.fromArray(nodes[p]); B.fromArray(nodes[i]);
        const da = n.dot(A) + cst; const db = n.dot(B) + cst;
        if ((da > 0) === (db > 0)) continue;
        const t = da / (da - db);
        P.copy(A).lerp(B, t);
        const r = rad[p] + (rad[i] - rad[p]) * t;
        if (k < holder.max) {
          m.compose(P, q, new THREE.Vector3(r, r, r));
          holder.wall.setMatrixAt(k, m);
          holder.lum.setMatrixAt(k, m);
        }
        // direction of travel relative to the cut face, from the segment itself
        const dirDot = B.clone().sub(A).normalize().dot(n);
        if (dirDot > 0) out++; else into++;
        diam.push(2 * r);
        k++;
      }
      holder.wall.count = Math.min(k, holder.max);
      holder.lum.count = Math.min(k, holder.max);
      holder.wall.instanceMatrix.needsUpdate = true;
      holder.lum.instanceMatrix.needsUpdate = true;
      diam.sort((a, b) => a - b);
      rep[tree.name] = {
        openings: k,
        rendered: holder.wall.count,
        diameter_mm_min: diam.length ? Number(diam[0].toFixed(4)) : null,
        diameter_mm_median: diam.length ? Number(diam[diam.length >> 1].toFixed(4)) : null,
        diameter_mm_max: diam.length ? Number(diam[diam.length - 1].toFixed(4)) : null,
        crossing_with_normal: out,
        crossing_against_normal: into,
        tree_segments_total: parent.filter(x => x >= 0).length,
        radius_scale_applied: 1.0
      };
    }
    this.lumenReportData = rep;
  }

  /* ------------------------------------------------------------- plane set */

  setPreset (name) {
    const pr = SECTION_PRESETS[name];
    if (!pr) return;
    this.presetName = name;
    this.setNormal(new THREE.Vector3().fromArray(pr.normal));
  }

  setNormal (nrm) {
    this.normal.copy(nrm).normalize();
    this.presetName = this.presetName || 'free';
    this._apply();
  }

  /** Free rotation in two angles (radians), used by the drag handle. */
  setFreeAngles (az, el) {
    this.presetName = 'free';
    const ce = Math.cos(el);
    this.setNormal(new THREE.Vector3(ce * Math.cos(az), Math.sin(el), ce * Math.sin(az)));
  }

  setOffset (t01) {
    this.offset01 = Math.min(Math.max(t01, 0.02), 0.98);
    this._apply();
  }

  setFlip (f) { this.flip = f < 0 ? -1 : 1; this._apply(); }

  _apply () {
    const n = this.normal.clone().multiplyScalar(this.flip);
    // span of the mesh along the normal
    const half = new THREE.Vector3(
      Math.abs(n.x) * this.extent.x, Math.abs(n.y) * this.extent.y,
      Math.abs(n.z) * this.extent.z).length() * 0.5;
    const d0 = n.dot(this.centre);
    const d = d0 + (this.offset01 * 2 - 1) * half;
    this.plane.set(n, -d);
    // The face frame's origin must be the ORGAN centre projected onto the plane,
    // not the plane's closest point to the world origin. The segmentation sits at
    // roughly (-49, 458, 184) mm in the scan frame, so taking n*d would place the
    // 2-D face coordinates several hundred millimetres away from the liver and the
    // whole cross-section would fall outside the classification grid.
    this.origin.copy(this.centre).addScaledVector(n, d - n.dot(this.centre));
    // face frame
    const ref = Math.abs(n.y) < 0.9 ? new THREE.Vector3(0, 1, 0) : new THREE.Vector3(1, 0, 0);
    this.u.crossVectors(ref, n).normalize();
    this.v.crossVectors(n, this.u).normalize();

    const basis = new THREE.Matrix4().makeBasis(this.u, this.v, n);
    for (const name of this.capOrder) {
      const q = this.caps[name];
      q.matrixAutoUpdate = false;
      q.matrix.copy(basis);
      q.matrix.setPosition(this.origin.clone().addScaledVector(n, -CAP_OFFSET_MM));
      q.matrixWorldNeedsUpdate = true;
      if (q.material.uniforms) q.material.uniforms.uHalf.value = this.diag * 0.55;
    }
    // Cheap every time: slice the parenchyma (137k triangles, a few ms) and
    // recompute the vessel openings. Expensive (the 313 x 263 nearest-terminal
    // classification over the grid) is deferred so dragging the plane stays
    // interactive; the face texture catches up when the drag settles.
    this.slice = sliceGeometry(this.parenchyma.geometry, this.plane.normal,
                               this.plane.constant, this.u, this.v, this.origin);
    this._recomputeLumens();
    this.faceDirty = true;
    this._faceDeadline = (typeof performance !== 'undefined' ? performance.now() : 0) + 90;
    if (!this.faceStats) this._recomputeFace();
    this._syncMaterials();
  }

  /**
   * One clip test, one place. Everything that can be cut — the parenchyma and
   * gallbladder surfaces, the four vessel trees, the flow particles and the
   * stencil pass that builds the cut face — reads the SAME body-local plane out
   * of the same four numbers. three.js local clipping planes are deliberately
   * NOT used: they live in world space and would need a second, independently
   * updated copy of the plane, which is exactly how a cut face ends up
   * describing a different solid from the one that was cut.
   */
  _syncMaterials () {
    const n = this.plane.normal;
    // the shader test is `dot(vLocal, n) > w  ->  discard`, and the plane keeps
    // n.p + constant < 0, so w = -constant
    const v4 = new THREE.Vector4(n.x, n.y, n.z, -this.plane.constant);
    const on = this.enabled ? 1.0 : 0.0;
    for (const m of this.clipMaterials) {
      const u = m.userData && m.userData.uniforms;
      if (!u || !u.uClipPlane) continue;
      u.uClipPlane.value.copy(v4);
      u.uClipOn.value = on;
    }
    for (const m of this.stencilMats) {
      m.uniforms.uClipPlane.value.copy(v4);
      m.uniforms.uClipOn.value = on;
    }
    for (const f of (this.flowSystems || [])) {
      f.material.uniforms.uClipPlane.value.copy(v4);
      f.material.uniforms.uClipOn.value = on;
    }
    this.clipVec4 = v4;
  }

  setEnabled (on) {
    this.enabled = !!on;
    this.group.visible = this.enabled;
    this._syncMaterials();
  }

  setSegmentMix (x) { this._setU('uSegMix', x); }

  _setU (key, value) {
    for (const name of this.capOrder) {
      const u = this.caps[name].material.uniforms;
      if (u && u[key]) u[key].value = value;
    }
  }

  /** @param {number} mode 0 oxygen, 1 zonation, 2 per-gene, 3 clearance, 4 bile */
  setFunctionChannel (mode, mix) {
    this._setU('uFuncMode', mode | 0);
    this._setU('uFuncMix', mix === undefined ? 1.0 : mix);
    this.funcMode = mode | 0;
  }

  setZoneRatio (r) { this._setU('uZoneRatio', r); }
  setExtraction (e) { this._setU('uExtraction', Math.min(Math.max(e, 0.02), 0.98)); }

  /** Keep the world-space clipping plane in step with the breathing body. */
  update (simMs) {
    this.liverBody.updateMatrixWorld();
    this.worldPlane.copy(this.plane).applyMatrix4(this.liverBody.matrixWorld);
    for (const name of this.capOrder) {
      const u = this.caps[name].material.uniforms;
      if (u) u.uTime.value = simMs * 1e-3;
    }
    if (this.faceDirty && this.enabled) {
      const now = (typeof performance !== 'undefined') ? performance.now() : 0;
      if (now >= this._faceDeadline) this._recomputeFace();
    }
  }

  /** Force the deferred classification to finish (used by the audit probe). */
  flush () { if (this.faceDirty) this._recomputeFace(); }

  /** World position of the cut-face centroid (for the respiration audit). */
  faceCentreWorld () {
    return this.origin.clone().applyMatrix4(this.liverBody.matrixWorld);
  }

  report () {
    return {
      preset: this.presetName,
      normal_local: this.plane.normal.toArray().map(x => Number(x.toFixed(6))),
      plane_constant_local_mm: Number(this.plane.constant.toFixed(4)),
      offset01: Number(this.offset01.toFixed(4)),
      capping: {
        method: 'two-pass stencil (back faces increment, front faces decrement), '
              + 'cap quad painted only where stencil != 0',
        solids_capped: this.capOrder.slice(),
        cap_quad_geometric_offset_mm: CAP_OFFSET_MM,
        depth_fight_handling: 'polygonOffset per cap, no geometric displacement'
      },
      face: this.faceStats,
      lumens: this.lumenReportData
    };
  }
}

export const CAP_SEG_COLORS = [
  '#000000', '#c8a24a', '#7fbf7f', '#5fa8d3', '#b07fc0',
  '#e08a5a', '#d4606a', '#8a8fd0', '#52b0a0'
];
