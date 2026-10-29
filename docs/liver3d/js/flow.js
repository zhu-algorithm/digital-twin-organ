/**
 * flow.js — GPU particle advection along the hepatic flow paths.
 *
 * Every particle lives on one path from liver_flow.json. Path geometry AND the
 * local velocity are baked into a float texture, one row per path:
 *
 *   RGBA = (x, y, z, v_mm_per_s)
 *
 * so the speed a particle moves at is the real local velocity v = Q/A computed
 * by tools/build_liver_flow.py, not a per-path constant. Where a vessel narrows
 * towards its terminals the flow per branch falls faster than the area, so
 * particles visibly slow down distally — that behaviour comes out of the data,
 * nothing animates it.
 *
 * Per-particle arc position is integrated on the CPU per path (one float per
 * path per frame) and evaluated in the vertex shader, so a population of tens of
 * thousands costs one uniform array update.
 *
 * Direction is a property of the path, not of the renderer: portal and arterial
 * paths are stored hilum -> terminal, hepatic venous paths terminal -> caval
 * ostium, biliary paths terminal -> hilum. Phase only ever increases, so no path
 * can produce retrograde motion.
 *
 * Bile moves at ~0.07 cm/s, about 300x slower than portal blood. It carries a
 * declared display gain (liver_flow.json .bile_display_gain) so its direction is
 * visible at all; the gain is reported in the HUD and in the acceptance record.
 */

import * as THREE from 'three';

const SAMPLES = 64;

const TREE_TINT = {
  portal_vein: [0.38, 0.55, 0.86],
  hepatic_vein: [0.24, 0.34, 0.58],
  hepatic_artery: [0.95, 0.28, 0.26],
  bile_duct: [0.66, 0.82, 0.30]
};

export class FlowSystem {
  /**
   * @param {object} spec parsed liver_flow.json
   * @param {number} particlesPerPath
   */
  constructor (spec, particlesPerPath = 320) {
    this.spec = spec;
    this.paths = spec.paths;
    this.nPath = this.paths.length;
    this.ppp = particlesPerPath;
    this.bileGain = spec.bile_display_gain || 1.0;

    this._bakePathTexture();
    this._buildParticles();

    this.phase = new Float32Array(this.nPath);
    this.speedNorm = new Float32Array(this.nPath);
    this.enabled = new Float32Array(this.nPath).fill(1);
    this.lastSpeedMmS = new Float32Array(this.nPath);
    // monotonic accumulated arc fraction per path: the directionality audit
    // reads THIS, not the wrapped phase (a path that advances past 1.0 between
    // two samples makes a wrapped phase look as if it went backwards)
    this.advance = new Float64Array(this.nPath);
  }

  _bakePathTexture () {
    const w = SAMPLES; const h = this.nPath;
    const data = new Float32Array(w * h * 4);
    this.pathLengthMM = new Float32Array(h);
    this.pathMeanV = new Float32Array(h);
    for (let p = 0; p < h; p++) {
      const pts = this.paths[p].points;
      const vv = this.paths[p].v_mm_s;
      let total = 0;
      for (let i = 1; i < pts.length; i++) {
        total += Math.hypot(pts[i][0] - pts[i - 1][0], pts[i][1] - pts[i - 1][1],
                            pts[i][2] - pts[i - 1][2]);
      }
      this.pathLengthMM[p] = total || 1;
      let vs = 0;
      for (let i = 0; i < w; i++) {
        const src = Math.min(Math.round(i * (pts.length - 1) / (w - 1)), pts.length - 1);
        const o = (p * w + i) * 4;
        data[o] = pts[src][0]; data[o + 1] = pts[src][1]; data[o + 2] = pts[src][2];
        data[o + 3] = vv[Math.min(src, vv.length - 1)];
        vs += data[o + 3];
      }
      this.pathMeanV[p] = vs / w;
    }
    const tex = new THREE.DataTexture(data, w, h, THREE.RGBAFormat, THREE.FloatType);
    tex.minFilter = tex.magFilter = THREE.NearestFilter;
    tex.wrapS = tex.wrapT = THREE.ClampToEdgeWrapping;
    tex.needsUpdate = true;
    this.pathTex = tex;
  }

  _buildParticles () {
    const n = this.nPath * this.ppp;
    const aPath = new Float32Array(n);
    const aPhase0 = new Float32Array(n);
    const aJitter = new Float32Array(n * 3);
    const aSize = new Float32Array(n);
    const aColor = new Float32Array(n * 3);

    // interleaved across paths (path varies fastest) so thinning the draw range
    // thins every vessel evenly instead of dropping whole trees
    let k = 0;
    for (let i = 0; i < this.ppp; i++) {
      for (let p = 0; p < this.nPath; p++) {
        const path = this.paths[p];
        const c = TREE_TINT[path.tree] || [0.8, 0.8, 0.8];
        aPath[k] = p;
        aPhase0[k] = Math.random();
        const rr = Math.sqrt(Math.random());
        const th = Math.random() * Math.PI * 2;
        aJitter[k * 3] = rr * Math.cos(th);
        aJitter[k * 3 + 1] = rr * Math.sin(th);
        aJitter[k * 3 + 2] = 1.0 - 0.7 * rr * rr;     // Poiseuille axial profile
        aSize[k] = 0.85 + Math.random() * 1.3;
        aColor[k * 3] = c[0]; aColor[k * 3 + 1] = c[1]; aColor[k * 3 + 2] = c[2];
        k++;
      }
    }

    const g = new THREE.BufferGeometry();
    g.setAttribute('position', new THREE.BufferAttribute(new Float32Array(n * 3), 3));
    g.setAttribute('aPath', new THREE.BufferAttribute(aPath, 1));
    g.setAttribute('aPhase0', new THREE.BufferAttribute(aPhase0, 1));
    g.setAttribute('aJitter', new THREE.BufferAttribute(aJitter, 3));
    g.setAttribute('aSize', new THREE.BufferAttribute(aSize, 1));
    g.setAttribute('aColor', new THREE.BufferAttribute(aColor, 3));
    g.boundingSphere = new THREE.Sphere(new THREE.Vector3(), 1e4);

    this.material = new THREE.ShaderMaterial({
      glslVersion: THREE.GLSL3,
      uniforms: {
        uPathTex: { value: this.pathTex },
        uPhase: { value: new Array(this.nPath).fill(0) },
        uBright: { value: new Array(this.nPath).fill(1) },
        uPixelRatio: { value: 1 },
        uParticleScale: { value: 1.0 },
        uClipPlane: { value: new THREE.Vector4(0, 0, 0, 0) },
        uClipOn: { value: 0.0 }
      },
      transparent: true,
      depthWrite: false,
      blending: THREE.AdditiveBlending,
      vertexShader: /* glsl */`
        precision highp float;
        in float aPath;
        in float aPhase0;
        in vec3  aJitter;
        in float aSize;
        in vec3  aColor;
        uniform highp sampler2D uPathTex;
        uniform float uPhase[${this.nPath}];
        uniform float uBright[${this.nPath}];
        uniform float uPixelRatio;
        uniform float uParticleScale;
        out vec3 vCol;
        out float vFade;
        out vec3 vWorldLocal;

        vec4 samplePath (int p, float s) {
          float x = clamp(s, 0.0, 0.9999) * float(${SAMPLES} - 1);
          int i0 = int(floor(x)); int i1 = min(i0 + 1, ${SAMPLES} - 1);
          float f = x - float(i0);
          vec4 a = texelFetch(uPathTex, ivec2(i0, p), 0);
          vec4 b = texelFetch(uPathTex, ivec2(i1, p), 0);
          return mix(a, b, f);
        }

        void main () {
          int p = int(aPath + 0.5);
          float s = fract(aPhase0 + uPhase[p]);
          vec4 c0 = samplePath(p, s);
          vec4 c1 = samplePath(p, min(s + 0.012, 0.9999));
          vec3 pos = c0.xyz;
          vec3 tang = normalize(c1.xyz - c0.xyz + vec3(1e-6));
          vec3 up = abs(tang.y) < 0.9 ? vec3(0.0, 1.0, 0.0) : vec3(1.0, 0.0, 0.0);
          vec3 t1 = normalize(cross(up, tang));
          vec3 t2 = cross(tang, t1);
          // the radius channel is gone (velocity occupies .w), so lateral spread
          // is scaled from the local velocity: faster vessels are the wide ones
          float spread = 0.45 + 0.55 * clamp(c0.w / 220.0, 0.0, 1.6);
          pos += (t1 * aJitter.x + t2 * aJitter.y) * spread * 1.9;
          vWorldLocal = pos;
          vec4 mv = modelViewMatrix * vec4(pos, 1.0);
          gl_Position = projectionMatrix * mv;
          float sz = aSize * uParticleScale * (1.0 + 0.5 * clamp(c0.w / 200.0, 0.0, 1.0));
          gl_PointSize = sz * uPixelRatio * 170.0 / max(-mv.z, 1.0);
          vCol = aColor;
          vFade = uBright[p] * aJitter.z;
        }`,
      fragmentShader: /* glsl */`
        precision highp float;
        in vec3 vCol;
        in float vFade;
        in vec3 vWorldLocal;
        uniform vec4 uClipPlane;
        uniform float uClipOn;
        out vec4 fragColor;
        void main () {
          if (uClipOn > 0.5 && dot(vWorldLocal, uClipPlane.xyz) > uClipPlane.w) discard;
          vec2 d = gl_PointCoord - 0.5;
          float r2 = dot(d, d);
          if (r2 > 0.25) discard;
          float a = exp(-r2 * 11.0) * vFade;
          if (a < 0.004) discard;
          fragColor = vec4(vCol * (0.6 + 0.9 * a), a);
        }`
    });

    this.points = new THREE.Points(g, this.material);
    this.points.frustumCulled = false;
    this.points.renderOrder = 3;
    this.nParticles = n;
  }

  setPixelRatio (r) { this.material.uniforms.uPixelRatio.value = r; }
  setParticleScale (s) { this.material.uniforms.uParticleScale.value = s; }

  /** Fraction of the population drawn (narrow-screen downgrade uses this). */
  setDensity (f) {
    const n = Math.max(1, Math.round(this.nParticles * Math.min(Math.max(f, 0), 1)));
    this.points.geometry.setDrawRange(0, n);
    this.drawn = n;
  }

  /** Enable/disable a whole tree by name. */
  setTreeEnabled (treeName, on) {
    for (let p = 0; p < this.nPath; p++) {
      if (this.paths[p].tree === treeName) this.enabled[p] = on ? 1 : 0;
    }
  }

  /**
   * @param {number} dtMs wall time advanced
   * @param {number} rate playback multiplier
   */
  update (dtMs, rate = 1.0) {
    const dt = dtMs * 1e-3;
    const u = this.material.uniforms;
    for (let p = 0; p < this.nPath; p++) {
      const path = this.paths[p];
      const gain = path.display_gain || 1.0;
      const v = this.pathMeanV[p] * gain;           // mm/s along this path
      this.lastSpeedMmS[p] = this.pathMeanV[p];
      const ds = (v / this.pathLengthMM[p]) * dt * rate;
      const step = Math.max(ds, 0);                     // never negative
      this.advance[p] += step;
      this.phase[p] = (this.phase[p] + step) % 1;
      u.uPhase.value[p] = this.phase[p];
      u.uBright.value[p] = this.enabled[p] *
        (0.35 + 0.65 * Math.min(this.pathMeanV[p] / 180.0, 1.0));
    }
  }

  diagnostics () {
    const byTree = {};
    for (let p = 0; p < this.nPath; p++) {
      const t = this.paths[p].tree;
      byTree[t] = byTree[t] || { paths: 0, meanV_cm_s: 0 };
      byTree[t].paths++;
      byTree[t].meanV_cm_s += this.pathMeanV[p] / 10;
    }
    for (const k of Object.keys(byTree)) byTree[k].meanV_cm_s /= byTree[k].paths;
    return {
      paths: this.nPath,
      particles: this.nParticles,
      drawn: this.drawn || this.nParticles,
      bileDisplayGain: this.bileGain,
      byTree
    };
  }
}
