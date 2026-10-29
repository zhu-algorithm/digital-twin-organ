/**
 * materials.js — liver surface appearance, entirely procedural.
 *
 * No imported texture assets of unverified provenance. All surface detail is
 * generated in the fragment shader from fields BAKED ON THE MESH by
 * tools/build_liver_fields.py, so the detail belongs to this particular liver
 * rather than sliding over it:
 *
 *   _field0 = (capsule_depth_01, couinaud_segment_id, zone_01, concavity)
 *   _field1 = (falciform_mask, hemiliver, umbilical_fissure_distance, lobule_phase)
 *
 *  - Glisson's capsule: a thin, wet, specular serous layer. Rendered as a
 *    clearcoat sheen over a deep red-brown parenchyma; the capsule is what gives
 *    a fresh liver its characteristic glassy highlight.
 *  - Parenchyma grain: the ~1 mm lobular mottling of the liver surface, built
 *    from a hexagonal Worley-style cell pattern in the tangent plane at the real
 *    lobular length scale (lobule diameter ~1 mm), not from fractal noise alone.
 *  - Falciform ligament and ligamentum teres: drawn ONLY where the baked
 *    falciform mask is non-zero, i.e. where the umbilical-fissure plane meets the
 *    capsule AND the segmented surface is genuinely concave. It is a measured
 *    groove, not a painted line — a hand-drawn ligament was the main complaint
 *    about the earlier 2-D page.
 *  - Couinaud tint: optional, driven by the baked segment id.
 *  - Zonal tint: optional, driven by the baked zone_01 field (the normalised
 *    portal-to-central position computed from the two vascular trees). The
 *    colour ramp is bound to one measured gene at a time by main.js, because the
 *    zonal gradient is NOT the same for every gene.
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
// Worley F2-F1: polygonal cell edges. At cell size ~1 mm this is the lobular
// mottling actually visible on a fresh liver capsule.
float lobuleCells (vec2 p, out vec2 cellId) {
  vec2 n = floor(p), f = fract(p);
  float f1 = 8.0, f2 = 8.0; vec2 best = n;
  for (int j = -1; j <= 1; j++) {
    for (int i = -1; i <= 1; i++) {
      vec2 g = vec2(float(i), float(j));
      vec2 o = h22(n + g);
      float d = length(g + o - f);
      if (d < f1) { f2 = f1; f1 = d; best = n + g; }
      else if (d < f2) { f2 = d; }
    }
  }
  cellId = best;
  return f2 - f1;
}
vec3 tangentFrameU (vec3 nrm, vec3 pos) {
  vec3 t = abs(nrm.y) < 0.9 ? normalize(cross(vec3(0, 1, 0), nrm))
                            : normalize(cross(vec3(1, 0, 0), nrm));
  return t;
}
`;

const COMMON_UNIFORMS = () => ({
  uSegmentMix: { value: 0.0 },
  uZoneMix: { value: 0.0 },
  uZoneRatio: { value: 1.0 },       // central/periportal ratio of the bound gene
  uSegColors: { value: SEGMENT_COLORS.map(c => new THREE.Color(c)) },
  uLigament: { value: 1.0 },
  uClipPlane: { value: new THREE.Vector4(0, 0, 0, 0) },   // xyz normal, w offset
  uClipOn: { value: 0.0 },
  uTime: { value: 0.0 }
});

/** Couinaud segment palette, index 0 unused (id 1..8). */
export const SEGMENT_COLORS = [
  '#000000', '#c8a24a', '#7fbf7f', '#5fa8d3', '#b07fc0',
  '#e08a5a', '#d4606a', '#8a8fd0', '#52b0a0'
];

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
      uniform float uSegmentMix;
      uniform float uZoneMix;
      uniform float uZoneRatio;
      uniform vec3  uSegColors[9];
      uniform float uLigament;
      uniform vec4  uClipPlane;
      uniform float uClipOn;
      uniform float uTime;
      ${FS_LIB}`);
  return shader;
}

/**
 * Parenchyma + Glisson capsule.
 * @param {object} opts
 */
export function makeParenchymaMaterial (opts = {}) {
  const mat = new THREE.MeshPhysicalMaterial({
    color: new THREE.Color(opts.color || 0x7c2f28),
    roughness: 0.52,
    metalness: 0.0,
    clearcoat: 0.85,            // the serous capsule
    clearcoatRoughness: 0.12,
    sheen: 0.25,
    sheenColor: new THREE.Color(0xd08878),
    side: THREE.DoubleSide
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
      float segId   = vF0.y;
      float zone    = clamp(vF0.z, 0.0, 1.0);
      float conc    = vF0.w;
      float lig     = vF1.x;

      // --- parenchymal grain at the real lobular scale (~1 mm cells) -------
      vec3 t = tangentFrameU(normalize(vNrm), vLocal);
      vec3 b = cross(normalize(vNrm), t);
      vec2 uvmm = vec2(dot(vLocal, t), dot(vLocal, b));
      vec2 cid;
      float edge = lobuleCells(uvmm * 1.05, cid);        // 1 unit ~= 0.95 mm
      float grain = smoothstep(0.02, 0.26, edge);
      float mottle = fbm(uvmm * 0.22) * 0.5 + fbm(uvmm * 0.9) * 0.3;

      vec3 base = diffuseColor.rgb;
      base *= 0.82 + 0.28 * grain;
      base = mix(base, base * vec3(1.08, 0.86, 0.82), mottle * 0.5);
      // a subtly paler, yellower capsule where the surface is convex and taut
      base = mix(base, base * vec3(1.12, 1.02, 0.92), smoothstep(0.0, 0.6, conc));

      // --- measured falciform / round-ligament groove ----------------------
      float ligm = clamp(lig, 0.0, 1.0) * uLigament;
      // the ligament itself: a pale fibrous band inside the measured groove
      float band = smoothstep(0.18, 0.55, ligm);
      base = mix(base, vec3(0.86, 0.80, 0.66), band * 0.9);

      // --- optional overlays ----------------------------------------------
      int si = int(clamp(segId, 0.0, 8.0) + 0.5);
      vec3 segc = uSegColors[si];
      base = mix(base, mix(base, segc, 0.72), uSegmentMix);

      // Zonal ramp: periportal -> pericentral, with the direction taken from the
      // MEASURED central/periportal ratio (uZoneRatio < 1 means the gene is higher
      // periportally, so the ramp runs the other way). The ramp is linear in the
      // zone coordinate on purpose — shaping it with an exponent would make a
      // shallow gradient look steep, which is exactly the claim this scene must
      // not make. The magnitude is stated numerically instead.
      float z = uZoneRatio >= 1.0 ? zone : 1.0 - zone;
      vec3 zonec = mix(vec3(0.13, 0.22, 0.46), vec3(0.88, 0.60, 0.12), clamp(z, 0.0, 1.0));
      base = mix(base, mix(base * 0.5, zonec, 0.88), uZoneMix);

      diffuseColor.rgb = base;
      `);
    shader.fragmentShader = shader.fragmentShader.replace(
      '#include <roughnessmap_fragment>', /* glsl */`
      #include <roughnessmap_fragment>
      // the capsule is wettest on convex, exposed surface
      roughnessFactor = clamp(roughnessFactor - 0.10 * smoothstep(0.0, 0.5, vF0.w)
                              + 0.22 * clamp(vF1.x, 0.0, 1.0), 0.05, 1.0);
      `);
  };
  mat.customProgramCacheKey = () => 'liver-parenchyma-v2';
  return mat;
}

/** Vessel wall: thin, translucent, tinted by contents. */
export function makeVesselMaterial (color, opts = {}) {
  const mat = new THREE.MeshPhysicalMaterial({
    color: new THREE.Color(color),
    roughness: opts.roughness !== undefined ? opts.roughness : 0.58,
    metalness: 0.0,
    clearcoat: 0.18,
    clearcoatRoughness: 0.4,
    transparent: opts.opacity !== undefined && opts.opacity < 1.0,
    opacity: opts.opacity !== undefined ? opts.opacity : 1.0,
    side: THREE.DoubleSide
  });
  const u = COMMON_UNIFORMS();
  mat.userData.uniforms = u;
  mat.onBeforeCompile = (shader) => {
    patch(shader, u);
    shader.fragmentShader = shader.fragmentShader.replace(
      '#include <map_fragment>', /* glsl */`
      #include <map_fragment>
      if (uClipOn > 0.5 && dot(vLocal, uClipPlane.xyz) > uClipPlane.w) discard;
      float v = fbm(vec2(vLocal.y * 0.6, atan(vLocal.z, vLocal.x) * 2.0));
      diffuseColor.rgb *= 0.88 + 0.24 * v;
      `);
  };
  mat.customProgramCacheKey = () => 'liver-vessel-' + color;
  return mat;
}

/** Gallbladder: bile-filled, translucent green-tinged. */
export function makeGallbladderMaterial () {
  const mat = new THREE.MeshPhysicalMaterial({
    color: new THREE.Color(0x5f7a3a),
    roughness: 0.22,
    metalness: 0.0,
    clearcoat: 1.0,
    clearcoatRoughness: 0.06,
    transmission: 0.55,
    thickness: 8.0,
    ior: 1.36,
    side: THREE.DoubleSide
  });
  const u = COMMON_UNIFORMS();
  mat.userData.uniforms = u;
  mat.onBeforeCompile = (shader) => { patch(shader, u); };
  mat.customProgramCacheKey = () => 'liver-gallbladder-v1';
  return mat;
}

export const ALL_MATERIAL_UNIFORMS = (mats) =>
  mats.filter(Boolean).map(m => m.userData.uniforms).filter(Boolean);
