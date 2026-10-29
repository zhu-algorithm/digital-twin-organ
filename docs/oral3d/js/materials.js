/**
 * materials.js — oral surface appearance, entirely procedural.
 *
 * Every mesh in oral.glb carries two custom vertex attributes baked by
 * tools/build_oral_fields.py:
 *   _field0 = (region_code, keratinised, epithelium_mm, permeability_rank/8)
 *   _field1 = (u_along_tongue, v_height_in_part, part_id, kinematic_group)
 *
 * The shader reads them, so what you see encodes what was measured or declared:
 *  - keratinised mucosa (hard palate, gingiva, tongue dorsum) is paler, drier
 *    and less translucent than non-keratinised mucosa (buccal, floor of mouth,
 *    ventral tongue, soft palate), which is the actual optical difference and
 *    also the anatomical basis for sublingual drug delivery
 *  - epithelial thickness drives subsurface depth, so the thin sublingual
 *    epithelium looks thin
 *  - the permeability ORDERING can be switched on as a false-colour overlay.
 *    It is an ordering, not a coefficient; the HUD says so whenever it is on.
 *
 * No texture files, no image maps: the bundle ships geometry and shaders only.
 */

import * as THREE from 'three';

export const COLOR = {
  bone: 0xd9cdb4,
  tooth: 0xf2ece0,
  tongue: 0xb4615c,
  mucosa_kerat: 0xe2a79c,
  mucosa_nonkerat: 0xc9736e,
  gland_serous: 0xc9a87c,
  gland_mixed: 0xc29a73,
  gland_mucous: 0xb8946e,
  muscle: 0x9c4a44,
  airway: 0x7f8c9e,
  vein: 0x3b64ab,
  artery: 0xcf3328,
  papilla: 0xc06a62
};

const VARYINGS = /* glsl */`
  varying vec4 vF0;
  varying vec4 vF1;
  varying vec3 vObj;
`;

const FS_LIB = /* glsl */`
  float hash31(vec3 p){
    p = fract(p * 0.3183099 + vec3(0.1, 0.2, 0.3));
    p *= 17.0;
    return fract(p.x * p.y * p.z * (p.x + p.y + p.z));
  }
  float noise3(vec3 x){
    vec3 i = floor(x); vec3 f = fract(x);
    f = f * f * (3.0 - 2.0 * f);
    float n000 = hash31(i + vec3(0,0,0)), n100 = hash31(i + vec3(1,0,0));
    float n010 = hash31(i + vec3(0,1,0)), n110 = hash31(i + vec3(1,1,0));
    float n001 = hash31(i + vec3(0,0,1)), n101 = hash31(i + vec3(1,0,1));
    float n011 = hash31(i + vec3(0,1,1)), n111 = hash31(i + vec3(1,1,1));
    return mix(mix(mix(n000,n100,f.x), mix(n010,n110,f.x), f.y),
               mix(mix(n001,n101,f.x), mix(n011,n111,f.x), f.y), f.z);
  }
  float fbm(vec3 p){
    return 0.55 * noise3(p) + 0.28 * noise3(p * 2.07) + 0.17 * noise3(p * 4.11);
  }
`;

const COMMON_UNIFORMS = () => ({
  uPermOverlay: { value: 0.0 },
  uKeratOverlay: { value: 0.0 },
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
      vF0 = _field0; vF1 = _field1; vObj = position;`);
  shader.fragmentShader = shader.fragmentShader
    .replace('#include <common>', `#include <common>
      ${VARYINGS}
      uniform float uPermOverlay;
      uniform float uKeratOverlay;
      uniform float uTime;
      ${FS_LIB}`);
  return shader;
}

/**
 * Soft tissue: mucosa, tongue, glands, muscle. `kind` picks the base colour and
 * the grain scale; the shader then modulates by the BAKED keratinisation and
 * epithelial thickness, so two sheets with the same base colour still look
 * different if their measured epithelium differs.
 */
export function makeTissueMaterial (kind, opts = {}) {
  const base = new THREE.Color(opts.color ?? COLOR.mucosa_nonkerat);
  const mat = new THREE.MeshPhysicalMaterial({
    color: base,
    roughness: opts.roughness ?? 0.62,
    metalness: 0.0,
    clearcoat: opts.clearcoat ?? 0.25,
    clearcoatRoughness: 0.5,
    sheen: opts.sheen ?? 0.5,
    sheenColor: new THREE.Color(0xffd9d2),
    transmission: 0.0,
    side: opts.side ?? THREE.DoubleSide,
    transparent: opts.opacity !== undefined,
    opacity: opts.opacity ?? 1.0,
    flatShading: false
  });
  const uniforms = COMMON_UNIFORMS();
  uniforms.uGrain = { value: opts.grain ?? 1.0 };
  mat.userData.uniforms = uniforms;
  mat.onBeforeCompile = (shader) => {
    patch(shader, uniforms);
    shader.fragmentShader = shader.fragmentShader
      .replace('#include <common>', `#include <common>
        uniform float uGrain;`)
      .replace('#include <color_fragment>', `#include <color_fragment>
        float kerat = vF0.y;
        float epi   = vF0.z;            // mm
        float perm  = vF0.w;

        // keratinised epithelium: paler, drier, matte
        vec3 dry = mix(diffuseColor.rgb, vec3(0.93, 0.80, 0.74), 0.42);
        diffuseColor.rgb = mix(diffuseColor.rgb, dry, kerat);

        // thin epithelium lets the lamina propria's vessels show through:
        // the thinner the measured epithelium, the redder the surface
        float thin = clamp(1.0 - epi / 0.80, 0.0, 1.0);
        diffuseColor.rgb = mix(diffuseColor.rgb,
                               diffuseColor.rgb * vec3(1.18, 0.72, 0.70), thin * 0.55);

        // procedural grain at the papillary scale
        float g = fbm(vObj * (1.6 * uGrain));
        diffuseColor.rgb *= 0.90 + 0.20 * g;

        // false-colour overlays, both off by default
        vec3 permC = mix(vec3(0.12, 0.20, 0.45), vec3(1.0, 0.85, 0.25), perm);
        diffuseColor.rgb = mix(diffuseColor.rgb, permC, uPermOverlay * step(0.01, perm));
        vec3 keratC = kerat > 0.5 ? vec3(0.95, 0.80, 0.35) : vec3(0.30, 0.65, 0.85);
        diffuseColor.rgb = mix(diffuseColor.rgb, keratC,
                               uKeratOverlay * step(0.5, vF0.x));`);
  };
  mat.customProgramCacheKey = () => 'oral-tissue-' + kind + '-' +
    (opts.opacity ?? 1);
  return mat;
}

/** Bone: mandible, midface, cervical vertebrae. */
export function makeBoneMaterial (opts = {}) {
  const mat = new THREE.MeshPhysicalMaterial({
    color: new THREE.Color(opts.color ?? COLOR.bone),
    roughness: 0.78, metalness: 0.0, clearcoat: 0.12,
    transparent: opts.opacity !== undefined, opacity: opts.opacity ?? 1.0,
    side: THREE.FrontSide
  });
  const uniforms = COMMON_UNIFORMS();
  mat.userData.uniforms = uniforms;
  mat.onBeforeCompile = (shader) => {
    patch(shader, uniforms);
    shader.fragmentShader = shader.fragmentShader
      .replace('#include <color_fragment>', `#include <color_fragment>
        float g = fbm(vObj * 0.75);
        diffuseColor.rgb *= 0.88 + 0.22 * g;
        diffuseColor.rgb = mix(diffuseColor.rgb,
                               diffuseColor.rgb * vec3(0.92, 0.95, 1.02),
                               0.35 * smoothstep(0.3, 0.9, g));`);
  };
  mat.customProgramCacheKey = () => 'oral-bone-' + (opts.opacity ?? 1);
  return mat;
}

/** Enamel: the crowns. Bright, slightly translucent, low roughness. */
export function makeToothMaterial () {
  const mat = new THREE.MeshPhysicalMaterial({
    color: new THREE.Color(COLOR.tooth), roughness: 0.22, metalness: 0.0,
    clearcoat: 0.85, clearcoatRoughness: 0.12,
    // NO transmission: three.js renders the whole scene into a backbuffer once
    // per transmissive mesh, which costs more than enamel translucency is worth
    ior: 1.63, side: THREE.FrontSide
  });
  const uniforms = COMMON_UNIFORMS();
  mat.userData.uniforms = uniforms;
  mat.onBeforeCompile = (shader) => {
    patch(shader, uniforms);
    shader.fragmentShader = shader.fragmentShader
      .replace('#include <color_fragment>', `#include <color_fragment>
        // cervical third is slightly more yellow (thinner enamel over dentine)
        diffuseColor.rgb = mix(diffuseColor.rgb * vec3(1.0, 0.95, 0.84),
                               diffuseColor.rgb, smoothstep(0.0, 0.6, vF1.y));`);
  };
  mat.customProgramCacheKey = () => 'oral-tooth-v1';
  return mat;
}

export function makeVesselMaterial (color, opts = {}) {
  const mat = new THREE.MeshPhysicalMaterial({
    color: new THREE.Color(color), vertexColors: true,
    roughness: 0.42, metalness: 0.0,
    clearcoat: 0.5, clearcoatRoughness: 0.3,
    transparent: opts.opacity !== undefined, opacity: opts.opacity ?? 1.0,
    side: THREE.FrontSide
  });
  const uniforms = COMMON_UNIFORMS();
  mat.userData.uniforms = uniforms;
  mat.onBeforeCompile = (shader) => {
    patch(shader, uniforms);
  };
  mat.customProgramCacheKey = () => 'oral-vessel-' + color + '-' +
    (opts.opacity ?? 1);
  return mat;
}

export function makeGreatVesselMaterial (color) {
  const mat = new THREE.MeshPhysicalMaterial({
    color: new THREE.Color(color), roughness: 0.5, metalness: 0.0,
    clearcoat: 0.3, transparent: true, opacity: 0.88, side: THREE.DoubleSide
  });
  const uniforms = COMMON_UNIFORMS();
  mat.userData.uniforms = uniforms;
  mat.onBeforeCompile = (shader) => { patch(shader, uniforms); };
  mat.customProgramCacheKey = () => 'oral-greatvessel-' + color;
  return mat;
}
