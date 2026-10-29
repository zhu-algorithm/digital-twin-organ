/**
 * materials.js — kidney surface appearance, entirely procedural.
 *
 * No imported texture assets of unverified provenance. All surface detail is
 * generated in the fragment shader from fields BAKED ON THE MESH by
 * tools/build_kidney_fields.py, so the detail belongs to this particular kidney
 * rather than sliding over it:
 *
 *   _field0 = (capsule_depth_01, pyramid_id, medulla_01, concavity)
 *   _field1 = (column_of_bertin, cmj_distance_norm, corticomedullary_zone, side)
 *
 *  - Renal capsule: a thin, tough, fibrous coat. Rendered as a tight clearcoat
 *    over a brown-red cortex. It is much less glassy than the liver's serous
 *    Glisson capsule, and the difference is deliberate.
 *  - Cortical granularity: the ~0.2 mm renal corpuscles make the cortical
 *    surface and cut face finely granular. Built from a Worley cell pattern at
 *    that real length scale, not from generic fractal noise.
 *  - Medullary striation: a pyramid is a bundle of parallel collecting ducts
 *    and loops of Henle running apex-ward, so the medulla is STRIATED along the
 *    pyramid axis while the cortex is granular. The shader switches between the
 *    two using the baked medulla_01 field, which is why the corticomedullary
 *    junction is visible on the cut face without any painted line.
 *  - Columns of Bertin: cortical tissue diving between pyramids. Tinted as
 *    cortex (because that is what it is) wherever the baked column field is set,
 *    so the columns read correctly even deep inside the organ.
 *  - Hilar groove: driven by the measured concavity field, not painted.
 */

import * as THREE from 'three';

const VARYINGS = /* glsl */`
varying vec4 vF0;
varying vec4 vF1;
varying vec3 vLocal;
varying vec3 vNrm;
`;

const FS_LIB = /* glsl */`
float h21 (vec2 p) {
  p = fract(p * vec2(123.34, 456.21));
  p += dot(p, p + 45.32);
  return fract(p.x * p.y);
}
vec2 h22 (vec2 p) {
  return fract(sin(vec2(dot(p, vec2(127.1, 311.7)), dot(p, vec2(269.5, 183.3)))) * 43758.5453);
}
float vnoise (vec2 p) {
  vec2 i = floor(p), f = fract(p);
  f = f * f * (3.0 - 2.0 * f);
  return mix(mix(h21(i), h21(i + vec2(1, 0)), f.x),
             mix(h21(i + vec2(0, 1)), h21(i + vec2(1, 1)), f.x), f.y);
}
float fbm (vec2 p) {
  float s = 0.0, a = 0.5;
  for (int i = 0; i < 5; i++) { s += a * vnoise(p); p *= 2.03; a *= 0.5; }
  return s;
}
// Worley F2-F1: polygonal cell edges. At cell size ~0.2 mm this is the renal
// corpuscle granularity visible on a cut cortex.
float cellEdges (vec2 p) {
  vec2 n = floor(p), f = fract(p);
  float f1 = 8.0, f2 = 8.0;
  for (int j = -1; j <= 1; j++) {
    for (int i = -1; i <= 1; i++) {
      vec2 g = vec2(float(i), float(j));
      float d = length(g + h22(n + g) - f);
      if (d < f1) { f2 = f1; f1 = d; } else if (d < f2) { f2 = d; }
    }
  }
  return f2 - f1;
}
vec3 tangentU (vec3 nrm) {
  return abs(nrm.y) < 0.9 ? normalize(cross(vec3(0, 1, 0), nrm))
                          : normalize(cross(vec3(1, 0, 0), nrm));
}
`;

export const CORTEX_COLOR = 0x9d4a3c;
export const MEDULLA_COLOR = 0x6e2f33;
export const COLUMN_COLOR = 0xa9594a;

const COMMON_UNIFORMS = () => ({
  uCortexMedulla: { value: 1.0 },   // show the corticomedullary distinction
  uColumns: { value: 1.0 },         // highlight the columns of Bertin
  uPyramidTint: { value: 0.0 },     // colour each pyramid separately
  uZone: { value: 0.0 },            // corticomedullary axis ramp
  uClipPlane: { value: new THREE.Vector4(0, 0, 0, 0) },
  uClipOn: { value: 0.0 },
  uCortexCol: { value: new THREE.Color(CORTEX_COLOR) },
  uMedullaCol: { value: new THREE.Color(MEDULLA_COLOR) },
  uColumnCol: { value: new THREE.Color(COLUMN_COLOR) },
  uTime: { value: 0.0 }
});

function patch (shader, uniforms) {
  Object.assign(shader.uniforms, uniforms);
  shader.vertexShader = shader.vertexShader
    .replace('#include <common>', `#include <common>
      attribute vec4 _field0;
      attribute vec4 _field1;
      ${VARYINGS}`)
    .replace('#include <begin_vertex>', `#include <begin_vertex>
      vF0 = _field0; vF1 = _field1; vLocal = position; vNrm = normal;`);
  shader.fragmentShader = shader.fragmentShader
    .replace('#include <common>', `#include <common>
      ${VARYINGS}
      uniform float uCortexMedulla;
      uniform float uColumns;
      uniform float uPyramidTint;
      uniform float uZone;
      uniform vec4  uClipPlane;
      uniform float uClipOn;
      uniform vec3  uCortexCol;
      uniform vec3  uMedullaCol;
      uniform vec3  uColumnCol;
      uniform float uTime;
      ${FS_LIB}`);
  return shader;
}

/** Renal parenchyma: cortex, medulla, columns, capsule. */
export function makeParenchymaMaterial (opts = {}) {
  const mat = new THREE.MeshPhysicalMaterial({
    color: new THREE.Color(opts.color || CORTEX_COLOR),
    roughness: 0.58,
    metalness: 0.0,
    clearcoat: 0.55,             // the fibrous renal capsule: tight, not glassy
    clearcoatRoughness: 0.28,
    // No sheen layer and single-sided: the capsule and cut meshes are closed,
    // so a back face is never visible, and both features cost fill rate that
    // the phone viewports cannot spare.
    side: THREE.FrontSide
  });
  const u = COMMON_UNIFORMS();
  mat.userData.uniforms = u;
  mat.onBeforeCompile = (shader) => {
    patch(shader, u);
    shader.fragmentShader = shader.fragmentShader.replace(
      '#include <map_fragment>', /* glsl */`
      #include <map_fragment>
      if (uClipOn > 0.5 && dot(vLocal, uClipPlane.xyz) > uClipPlane.w) discard;

      float depth01 = vF0.x;
      float pyrId   = vF0.y;
      float med     = clamp(vF0.z, 0.0, 1.0);
      float conc    = vF0.w;
      float column  = clamp(vF1.x, 0.0, 1.0);
      float zone    = clamp(vF1.z, 0.0, 1.0);

      vec3 nrm = normalize(vNrm);
      vec3 t = tangentU(nrm);
      vec3 b = cross(nrm, t);
      vec2 uvmm = vec2(dot(vLocal, t), dot(vLocal, b));

      // --- cortex: granular at the renal-corpuscle scale (~0.2 mm) ----------
      float gran = smoothstep(0.02, 0.30, cellEdges(uvmm * 5.0));
      vec3 cortex = uCortexCol * (0.84 + 0.30 * gran);
      cortex = mix(cortex, cortex * vec3(1.06, 0.92, 0.90),
                   fbm(uvmm * 0.35) * 0.45);

      // --- medulla: STRIATED along the pyramid axis, not granular ----------
      // a pyramid is a parallel bundle of collecting ducts and loops of Henle,
      // so its cut face is streaked rather than speckled
      float stripe = 0.5 + 0.5 * sin(dot(vLocal, b) * 7.5 + fbm(uvmm * 0.5) * 2.2);
      stripe = smoothstep(0.25, 0.85, stripe);
      vec3 medulla = uMedullaCol * (0.78 + 0.36 * stripe);

      vec3 base = mix(cortex, medulla, med * uCortexMedulla);
      // columns of Bertin are cortex, and are tinted as cortex even deep in
      base = mix(base, uColumnCol * (0.86 + 0.28 * gran), column * uColumns * 0.85);

      // optional: a distinct hue per pyramid so the lobar pattern is countable
      if (uPyramidTint > 0.001 && pyrId > 0.5) {
        float h = fract(pyrId * 0.1371);
        vec3 pc = 0.5 + 0.5 * cos(6.2831853 * (h + vec3(0.0, 0.33, 0.67)));
        base = mix(base, mix(base, pc, 0.68), uPyramidTint);
      }
      // optional: the corticomedullary axis, capsule -> papilla
      vec3 zc = mix(vec3(0.92, 0.72, 0.30), vec3(0.18, 0.28, 0.58), zone);
      base = mix(base, mix(base * 0.55, zc, 0.85), uZone);

      // --- capsule and the measured hilar groove ---------------------------
      base *= 0.90 + 0.20 * smoothstep(0.0, 0.35, 1.0 - depth01);
      base = mix(base, base * vec3(0.62, 0.52, 0.50),
                 clamp(-conc / 0.45, 0.0, 1.0) * 0.75);

      diffuseColor.rgb = base;
      `);
    shader.fragmentShader = shader.fragmentShader.replace(
      '#include <roughnessmap_fragment>', /* glsl */`
      #include <roughnessmap_fragment>
      // the capsule is tauter and glossier where the surface is convex; a cut
      // face is matt
      roughnessFactor = clamp(roughnessFactor - 0.10 * smoothstep(0.0, 0.5, vF0.w)
                              + 0.26 * clamp(vF0.z, 0.0, 1.0), 0.06, 1.0);
      `);
  };
  mat.customProgramCacheKey = () => 'kidney-parenchyma-v1';
  return mat;
}

/** Vessel wall: thin, tinted by contents. */
export function makeVesselMaterial (color, opts = {}) {
  const mat = new THREE.MeshPhysicalMaterial({
    color: new THREE.Color(color),
    roughness: opts.roughness !== undefined ? opts.roughness : 0.54,
    metalness: 0.0,
    clearcoat: 0.2,
    clearcoatRoughness: 0.4,
    transparent: opts.opacity !== undefined && opts.opacity < 1.0,
    opacity: opts.opacity !== undefined ? opts.opacity : 1.0,
    vertexColors: true,
    side: THREE.FrontSide
  });
  const u = COMMON_UNIFORMS();
  mat.userData.uniforms = u;
  mat.onBeforeCompile = (shader) => {
    patch(shader, u);
    shader.fragmentShader = shader.fragmentShader.replace(
      '#include <map_fragment>', /* glsl */`
      #include <map_fragment>
      if (uClipOn > 0.5 && dot(vLocal, uClipPlane.xyz) > uClipPlane.w) discard;
      `);
  };
  mat.customProgramCacheKey = () => 'kidney-vessel-' + color;
  return mat;
}

/** Great vessels (aorta, cava): seen from outside, slightly translucent. */
export function makeGreatVesselMaterial (color) {
  const mat = new THREE.MeshPhysicalMaterial({
    color: new THREE.Color(color),
    roughness: 0.42, metalness: 0.0, clearcoat: 0.0,
    transparent: true, opacity: 0.55, side: THREE.FrontSide, depthWrite: false
  });
  const u = COMMON_UNIFORMS();
  mat.userData.uniforms = u;
  mat.onBeforeCompile = (shader) => { patch(shader, u); };
  mat.customProgramCacheKey = () => 'kidney-greatvessel-' + color;
  return mat;
}
