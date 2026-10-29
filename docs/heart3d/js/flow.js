/**
 * flow.js — GPU blood-particle advection along chamber and vessel centrelines.
 *
 * Every particle lives on one named flow path (a resampled centreline through a
 * chamber or vessel). All per-particle position evaluation happens in the
 * vertex shader: the path polylines are baked into a float texture and the
 * particle's arc-length coordinate is `fract(phase0 + pathPhase[path])`, so a
 * single uniform per path moves the whole population.
 *
 * Physiological driving (all from the `heart_timeseries.json` contract, no
 * cosmetic animation anywhere):
 *
 *   path phase rate  d(s)/dt = v(t) / L ,  v(t) = Q(t) / A
 *
 * with Q read from the contract flow signal assigned to the path, A the path's
 * cross-sectional area and L its length. Q is clamped at zero, and each path is
 * additionally gated by its valve's `valve_open` flag, so:
 *   - systole, aortic valve open  -> LV outflow + aortic particles accelerate
 *   - diastole, mitral valve open -> LA -> LV filling particles move
 *   - valve closed                -> that path's particles stop advancing
 *   - no path ever advances backwards, so no regurgitant jet can appear
 *
 * Colouring: right heart (venae cavae -> RA -> RV -> pulmonary artery)
 * deoxygenated venous blue-purple; left heart (pulmonary veins -> LA -> LV ->
 * aorta) oxygenated arterial scarlet. The coronary paths are driven by
 * `Q_coronary_mL_s`, which the contract requires to be diastole-dominant, so
 * coronary filling is visibly out of phase with ejection.
 */

import * as THREE from 'three';
import { DEFORM_GLSL } from './deform.js';

const PATH_SAMPLES = 128;

export class FlowSystem {
  /**
   * @param {object} spec parsed flow_paths.json
   * @param {DeformField} field
   * @param {number} particlesPerPath
   */
  constructor (spec, field, particlesPerPath = 900) {
    this.spec = spec;
    this.field = field;
    this.paths = spec.paths;
    this.nPath = this.paths.length;
    this.ppp = particlesPerPath;

    this._bakePathTexture();
    this._buildParticles();

    this.phase = new Float32Array(this.nPath);
    this.speedNorm = new Float32Array(this.nPath);
    this._prepared = null;
  }

  /**
   * Precompute, per path, the cycle-mean of its bound flow signal (for paths
   * declared `flow_mode: "mean"`, i.e. the near-continuous caval and
   * pulmonary-venous returns) and the peak EFFECTIVE flow, used to normalise
   * particle brightness so a low-flow path is still visible.
   * @param {HeartTimeseries} ts
   */
  prepare (ts) {
    const meanQ = new Float32Array(this.nPath);
    const peakQ = new Float32Array(this.nPath);
    for (let p = 0; p < this.nPath; p++) {
      const path = this.paths[p];
      const a = ts.sig[path.flow_signal];
      const gate = path.valve_gate ? ts.valve[path.valve_gate] : null;
      const frac = path.flow_fraction !== undefined ? path.flow_fraction : 1.0;
      let sum = 0; let peak = 0;
      for (let i = 0; i < a.length; i++) {
        const q = Math.max(a[i], 0) * (gate ? gate[i] : 1);
        sum += q;
        if (q > peak) peak = q;
      }
      meanQ[p] = (sum / a.length) * frac;
      peakQ[p] = (path.flow_mode === 'mean' ? meanQ[p] : peak * frac) || 1e-6;
    }
    this._prepared = { meanQ, peakQ, src: ts };
  }

  _bakePathTexture () {
    // RGBA32F, one row per path: xyz = heart-space point, w = lumen radius (mm)
    const w = PATH_SAMPLES; const h = this.nPath;
    const data = new Float32Array(w * h * 4);
    this.pathLengthMM = new Float32Array(h);
    for (let p = 0; p < h; p++) {
      const pts = this.paths[p].points;
      const rad = this.paths[p].radius_mm;
      // cumulative arc length of the source polyline
      const cum = [0];
      for (let i = 1; i < pts.length; i++) {
        const dx = pts[i][0] - pts[i - 1][0];
        const dy = pts[i][1] - pts[i - 1][1];
        const dz = pts[i][2] - pts[i - 1][2];
        cum.push(cum[i - 1] + Math.sqrt(dx * dx + dy * dy + dz * dz));
      }
      const total = cum[cum.length - 1] || 1;
      this.pathLengthMM[p] = total;
      // uniform arc-length resampling
      let j = 0;
      for (let i = 0; i < w; i++) {
        const target = (i / (w - 1)) * total;
        while (j < cum.length - 2 && cum[j + 1] < target) j++;
        const seg = Math.max(cum[j + 1] - cum[j], 1e-9);
        const u = Math.min(Math.max((target - cum[j]) / seg, 0), 1);
        const a = pts[j]; const b = pts[j + 1];
        const o = (p * w + i) * 4;
        data[o] = a[0] + (b[0] - a[0]) * u;
        data[o + 1] = a[1] + (b[1] - a[1]) * u;
        data[o + 2] = a[2] + (b[2] - a[2]) * u;
        const ra = Array.isArray(rad) ? rad[Math.min(j, rad.length - 1)] : rad;
        const rb = Array.isArray(rad) ? rad[Math.min(j + 1, rad.length - 1)] : rad;
        data[o + 3] = ra + (rb - ra) * u;
      }
    }
    const tex = new THREE.DataTexture(data, w, h, THREE.RGBAFormat, THREE.FloatType);
    tex.minFilter = tex.magFilter = THREE.NearestFilter;
    tex.wrapS = THREE.ClampToEdgeWrapping;
    tex.wrapT = THREE.ClampToEdgeWrapping;
    tex.needsUpdate = true;
    this.pathTex = tex;
  }

  _buildParticles () {
    const n = this.nPath * this.ppp;
    const aPath = new Float32Array(n);
    const aPhase0 = new Float32Array(n);
    const aJitter = new Float32Array(n * 3);
    const aSize = new Float32Array(n);
    const aChamber = new Float32Array(n);
    const aTrans = new Float32Array(n);
    const aColor = new Float32Array(n * 3);

    // Particles are interleaved across paths (path varies fastest) so that
    // thinning the draw range reduces density evenly on every path instead of
    // dropping whole vessels.
    let k = 0;
    const col = new THREE.Color();
    for (let i = 0; i < this.ppp; i++) {
      for (let p = 0; p < this.nPath; p++) {
        const path = this.paths[p];
        col.set(path.color);
        aPath[k] = p;
        aPhase0[k] = Math.random();
        // Parabolic (Poiseuille-like) radial seeding: denser near the axis.
        const rr = Math.sqrt(Math.random());
        const th = Math.random() * Math.PI * 2;
        aJitter[k * 3] = rr * Math.cos(th);
        aJitter[k * 3 + 1] = rr * Math.sin(th);
        aJitter[k * 3 + 2] = 1.0 - 0.75 * rr * rr; // axial velocity profile factor
        aSize[k] = 0.9 + Math.random() * 1.5;
        aChamber[k] = path.chamber_id;
        aTrans[k] = path.transmural !== undefined ? path.transmural : 0.0;
        aColor[k * 3] = col.r; aColor[k * 3 + 1] = col.g; aColor[k * 3 + 2] = col.b;
        k++;
      }
    }

    const g = new THREE.BufferGeometry();
    g.setAttribute('position', new THREE.BufferAttribute(new Float32Array(n * 3), 3));
    g.setAttribute('aPath', new THREE.BufferAttribute(aPath, 1));
    g.setAttribute('aPhase0', new THREE.BufferAttribute(aPhase0, 1));
    g.setAttribute('aJitter', new THREE.BufferAttribute(aJitter, 3));
    g.setAttribute('aSize', new THREE.BufferAttribute(aSize, 1));
    g.setAttribute('aChamber', new THREE.BufferAttribute(aChamber, 1));
    g.setAttribute('aTrans', new THREE.BufferAttribute(aTrans, 1));
    g.setAttribute('aColor', new THREE.BufferAttribute(aColor, 3));
    g.boundingSphere = new THREE.Sphere(new THREE.Vector3(), 1e4);

    const u = Object.assign({}, this.field.uniforms, {
      uPathTex: { value: this.pathTex },
      uNPath: { value: this.nPath },
      uPhase: { value: new Array(this.nPath).fill(0) },
      uSpeedNorm: { value: new Array(this.nPath).fill(0) },
      uPixelRatio: { value: 1 },
      uParticleScale: { value: 1.0 }
    });

    this.material = new THREE.ShaderMaterial({
      glslVersion: THREE.GLSL3,
      uniforms: u,
      transparent: true,
      depthWrite: false,
      blending: THREE.AdditiveBlending,
      vertexShader: /* glsl */`
        precision highp float;
        in float aPath;
        in float aPhase0;
        in vec3  aJitter;
        in float aSize;
        in float aChamber;
        in float aTrans;
        in vec3  aColor;

        uniform highp sampler2D uPathTex;
        uniform int   uNPath;
        uniform float uPhase[${this.nPath}];
        uniform float uSpeedNorm[${this.nPath}];
        uniform float uPixelRatio;
        uniform float uParticleScale;

        out vec3  vCol;
        out float vAlpha;

        ${DEFORM_GLSL}

        void main() {
          int pi = int(aPath + 0.5);
          float s = fract(aPhase0 + uPhase[pi]);
          float v = float(pi) + 0.5;
          float fs = s * float(${PATH_SAMPLES} - 1);
          int i0 = int(floor(fs));
          int i1 = min(i0 + 1, ${PATH_SAMPLES} - 1);
          float w = fs - float(i0);
          vec4 samp = mix(texelFetch(uPathTex, ivec2(i0, pi), 0),
                          texelFetch(uPathTex, ivec2(i1, pi), 0), w);
          vec3 pH = samp.xyz;
          float lumen = samp.w;

          // local tangent for building the cross-section frame
          vec3 pA = texelFetch(uPathTex, ivec2(max(i0 - 1, 0), pi), 0).xyz;
          vec3 pB = texelFetch(uPathTex, ivec2(min(i1 + 1, ${PATH_SAMPLES} - 1), pi), 0).xyz;
          vec3 tang = normalize(pB - pA + vec3(1e-6));
          vec3 ref = abs(tang.z) < 0.9 ? vec3(0.0, 0.0, 1.0) : vec3(1.0, 0.0, 0.0);
          vec3 e1 = normalize(cross(tang, ref));
          vec3 e2 = cross(tang, e1);
          pH += (e1 * aJitter.x + e2 * aJitter.y) * lumen * 0.82;

          int ci = int(aChamber + 0.5);
          vec3 J; float phi;
          vec3 dH = deformHeart(pH, aTrans, ci, J, phi);
          // modelViewMatrix, not viewMatrix: the heart sits under a rotated
          // group, and skipping modelMatrix leaves the particles in a different
          // frame from the myocardium.
          vec4 mv = modelViewMatrix * uHeartToWorld * vec4(dH, 1.0);

          float spd = uSpeedNorm[pi];
          vCol = aColor * (0.55 + 0.75 * spd);
          // Particles fade out where flow has stopped, so a shut valve reads as
          // a dark, still column rather than a frozen bright streak.
          vAlpha = (0.10 + 0.90 * smoothstep(0.0, 0.25, spd)) * (0.55 + 0.45 * aJitter.z);

          gl_Position = projectionMatrix * mv;
          gl_PointSize = aSize * uParticleScale * uPixelRatio * (260.0 / max(-mv.z, 1.0))
                         * (0.75 + 0.5 * spd);
        }
      `,
      fragmentShader: /* glsl */`
        precision highp float;
        in vec3  vCol;
        in float vAlpha;
        out vec4 fragColour;
        void main() {
          vec2 d = gl_PointCoord - 0.5;
          float r2 = dot(d, d);
          if (r2 > 0.25) discard;
          float a = exp(-r2 * 11.0) * vAlpha;
          fragColour = vec4(vCol, a);
        }
      `
    });

    this.points = new THREE.Points(g, this.material);
    this.points.frustumCulled = false;
    this.points.renderOrder = 10;
  }

  /**
   * Advance every path phase using the contract flows.
   * @param {HeartTimeseries} ts
   * @param {number} tMs simulated time
   * @param {number} dtMs simulated time advanced this frame
   */
  update (ts, tMs, dtMs) {
    if (!this._prepared || this._prepared.src !== ts) this.prepare(ts);
    const { meanQ, peakQ } = this._prepared;
    const uPhase = this.material.uniforms.uPhase.value;
    const uSpd = this.material.uniforms.uSpeedNorm.value;
    for (let p = 0; p < this.nPath; p++) {
      const path = this.paths[p];
      const frac = path.flow_fraction !== undefined ? path.flow_fraction : 1.0;
      let qEff;
      if (path.flow_mode === 'mean') {
        // caval / pulmonary-venous return: near-continuous, so the cycle-mean
        // flow drives a steady column rather than a pulsatile one
        qEff = meanQ[p];
      } else {
        const q = Math.max(ts.sample(path.flow_signal, tMs), 0);   // no backflow
        const gate = path.valve_gate ? ts.valveOpen(path.valve_gate, tMs) : 1;
        qEff = q * gate * frac;
      }
      const vMMs = (qEff * 1000.0) / path.area_mm2;      // mm/s  (1 mL = 1000 mm^3)
      const dPhase = (vMMs * (dtMs / 1000.0)) / Math.max(this.pathLengthMM[p], 1e-6);
      uPhase[p] = (uPhase[p] + dPhase * path.phase_gain) % 1.0;
      uSpd[p] = Math.min(qEff / peakQ[p], 1.0);
      this.speedNorm[p] = uSpd[p];
      this.phase[p] = uPhase[p];
      this.velocityCmS = this.velocityCmS || new Float32Array(this.nPath);
      this.velocityCmS[p] = vMMs / 10.0;
    }
  }

  /** Diagnostic: per-path instantaneous normalised speed, by name. */
  speeds () {
    const o = {};
    for (let p = 0; p < this.nPath; p++) o[this.paths[p].name] = this.speedNorm[p];
    return o;
  }

  setPixelRatio (r) { this.material.uniforms.uPixelRatio.value = r; }
  setDensity (frac) {
    const n = this.nPath * this.ppp;
    this.points.geometry.setDrawRange(0, Math.floor(n * frac));
  }
}
