/**
 * materials.js — myocardial surface appearance.
 *
 * Everything here is procedural: no imported texture assets of unverified
 * provenance. Surface detail is generated in the fragment shader from the baked
 * cardiac coordinate system (transmural depth, apex->base coordinate zeta,
 * circumferential angle theta) which travels with the mesh as vertex
 * attributes, so the detail deforms with the muscle instead of sliding over it.
 *
 *  - Myocardial fibre anisotropy: helix angle rotates from about +60 deg on the
 *    endocardium to about -60 deg on the epicardium (Streeter's transmural
 *    fibre rotation). Striation runs along the local fibre direction in the
 *    (theta, zeta) tangent frame, so the apical vortex appears where theta
 *    collapses at the apex and the free walls show oblique fibres.
 *  - Epicardial fat: mask baked per-vertex from proximity to the
 *    atrioventricular groove and the anterior/posterior interventricular
 *    grooves, measured on the segmentation rather than painted. Fat is paler,
 *    yellower, lobulated and matte.
 *  - Wet epicardium: clearcoat serous sheen; endocardial surfaces darker.
 */

import * as THREE from 'three';
import { DEFORM_GLSL } from './deform.js';

const VARYINGS = /* glsl */`
varying vec3 vField0;   // transmural, zeta, radial
varying vec3 vField1;   // chamberId, fat, region
varying vec3 vHeartPos;
varying float vTheta;
`;

/**
 * Vertex injection. three.js evaluates <defaultnormal_vertex> BEFORE
 * <begin_vertex>, so the deformation is computed at <beginnormal_vertex> (where
 * `position` is already available and `objectNormal` has just been set) and the
 * resulting position is handed to <begin_vertex> through a global.
 */
const VS_HEAD = /* glsl */`
attribute vec4 _field0;
attribute vec4 _field1;
${VARYINGS}
vec3 gDeformed;
`;

const VS_NORMAL = /* glsl */`
{
  vec4 f0 = _field0;
  vec4 f1 = _field1;
  int  ci = int(f1.x + 0.5);
  vec3 pH = (uWorldToHeart * vec4(position, 1.0)).xyz;
  vec3 J; float phi;
  vec3 dH = deformHeart(pH, f0.x, ci, J, phi);
  gDeformed = (uHeartToWorld * vec4(dH, 1.0)).xyz;

  vec3 nH = normalize(mat3(uWorldToHeart) * objectNormal);
  nH = deformNormalHeart(nH, J, phi);
  objectNormal = normalize(mat3(uHeartToWorld) * nH);

  vField0 = f0.xyz;
  vField1 = f1.xyz;
  vHeartPos = dH;
  vTheta = f0.w;
}
`;

const FS_LIB = /* glsl */`
float h21 (vec2 p) {
  p = fract(p * vec2(123.34, 456.21));
  p += dot(p, p + 45.32);
  return fract(p.x * p.y);
}
float vnoise (vec2 p) {
  vec2 i = floor(p), f = fract(p);
  f = f * f * (3.0 - 2.0 * f);
  return mix(mix(h21(i), h21(i + vec2(1.0, 0.0)), f.x),
             mix(h21(i + vec2(0.0, 1.0)), h21(i + vec2(1.0, 1.0)), f.x), f.y);
}
float fbm (vec2 p) {
  float s = 0.0, a = 0.5;
  for (int i = 0; i < 5; i++) { s += a * vnoise(p); p *= 2.03; a *= 0.5; }
  return s;
}
vec2 fibreDir (float transmural) {
  float alpha = radians(mix(60.0, -60.0, clamp(transmural, 0.0, 1.0)));
  return vec2(cos(alpha), sin(alpha));
}
`;

function patchVertex (shader) {
  shader.vertexShader = shader.vertexShader
    .replace('#include <common>', '#include <common>\n' + VS_HEAD + DEFORM_GLSL)
    .replace('#include <beginnormal_vertex>', '#include <beginnormal_vertex>\n' + VS_NORMAL)
    .replace('#include <begin_vertex>', '#include <begin_vertex>\n  transformed = gDeformed;');
}

export function makeMyocardiumMaterial (field) {
  const mat = new THREE.MeshPhysicalMaterial({
    color: 0xffffff,
    roughness: 0.60,
    metalness: 0.0,
    clearcoat: 0.9,
    clearcoatRoughness: 0.12,
    side: THREE.DoubleSide
  });

  mat.onBeforeCompile = (shader) => {
    Object.assign(shader.uniforms, field.uniforms);
    shader.uniforms.uFibreScale = { value: 42.0 };
    patchVertex(shader);

    shader.fragmentShader = shader.fragmentShader
      .replace('#include <common>', `#include <common>
${VARYINGS}
uniform float uFibreScale;
${FS_LIB}
float gFat = 0.0;
float gMuscle = 0.5;`)
      .replace('#include <color_fragment>', `#include <color_fragment>
{
  float transmural = clamp(vField0.x, 0.0, 1.0);
  float zeta = vField0.y;
  float fat = clamp(vField1.y, 0.0, 1.0);
  float theta = vTheta;

  vec2 uvF = vec2(theta * 1.9, zeta * 3.1) * uFibreScale;
  vec2 fd = fibreDir(transmural);
  float along  = dot(uvF, fd);
  float across = dot(uvF, vec2(-fd.y, fd.x));
  float stria = fbm(vec2(along * 0.18, across * 1.35));
  float bands = 0.5 + 0.5 * sin(across * 2.6 + fbm(vec2(along * 0.25, across * 0.4)) * 5.0);
  float muscle = mix(stria, bands, 0.45);

  vec3 myoDeep = vec3(0.330, 0.072, 0.070);
  vec3 myoLit  = vec3(0.690, 0.215, 0.160);
  vec3 myo = mix(myoDeep, myoLit, muscle * 0.85 + 0.15);
  myo *= mix(0.74, 1.0, transmural);

  float lob = fbm(vec2(theta * 11.0, zeta * 9.0) * 1.7);
  vec3 fatCol = mix(vec3(0.80, 0.67, 0.36), vec3(0.90, 0.80, 0.52), lob);
  float fatA = smoothstep(0.45, 0.88, fat);

  diffuseColor.rgb *= mix(myo, fatCol, fatA);
  gFat = fatA;
  gMuscle = muscle;
}`)
      .replace('#include <roughnessmap_fragment>', `#include <roughnessmap_fragment>
  roughnessFactor = clamp(roughnessFactor + 0.30 * gFat - 0.12 * (gMuscle - 0.5), 0.05, 1.0);`)
      .replace('#include <lights_physical_fragment>', `#include <lights_physical_fragment>
  #ifdef USE_CLEARCOAT
    material.clearcoat = mix(material.clearcoat, 0.18, gFat);
    material.clearcoatRoughness = mix(material.clearcoatRoughness, 0.60, gFat);
  #endif`);
  };
  mat.customProgramCacheKey = () => 'myo-v2';
  return mat;
}

/** Translucent endocardial shell for the chamber interiors (seen from inside). */
export function makeBloodPoolMaterial (field, color, opacity = 0.16) {
  const mat = new THREE.MeshPhysicalMaterial({
    color, roughness: 0.3, metalness: 0.0,
    transparent: true, opacity, depthWrite: false, side: THREE.BackSide
  });
  mat.onBeforeCompile = (shader) => {
    Object.assign(shader.uniforms, field.uniforms);
    patchVertex(shader);
    shader.fragmentShader = shader.fragmentShader
      .replace('#include <common>', '#include <common>\n' + VARYINGS);
  };
  mat.customProgramCacheKey = () => 'pool-v2';
  return mat;
}

/** Great-vessel wall: pale, smooth, faint longitudinal adventitial texture. */
export function makeVesselWallMaterial (field, color) {
  const mat = new THREE.MeshPhysicalMaterial({
    color, roughness: 0.45, metalness: 0.0,
    clearcoat: 0.55, clearcoatRoughness: 0.28, side: THREE.DoubleSide
  });
  mat.onBeforeCompile = (shader) => {
    Object.assign(shader.uniforms, field.uniforms);
    patchVertex(shader);
    shader.fragmentShader = shader.fragmentShader
      .replace('#include <common>', '#include <common>\n' + VARYINGS + FS_LIB)
      .replace('#include <color_fragment>', `#include <color_fragment>
{
  float f = fbm(vec2(vTheta * 24.0, vField0.y * 40.0));
  diffuseColor.rgb *= 0.86 + 0.24 * f;
}`);
  };
  mat.customProgramCacheKey = () => 'vwall-v2';
  return mat;
}

/** Coronary artery wall: deep arterial red, glossy, deforms with the epicardium. */
export function makeCoronaryMaterial (field) {
  const mat = new THREE.MeshPhysicalMaterial({
    color: 0xd6433a,
    roughness: 0.48, metalness: 0.0,
    clearcoat: 0.45, clearcoatRoughness: 0.22,
    emissive: 0x2a0806, emissiveIntensity: 0.6
  });
  mat.onBeforeCompile = (shader) => {
    Object.assign(shader.uniforms, field.uniforms);
    patchVertex(shader);
    shader.fragmentShader = shader.fragmentShader
      .replace('#include <common>', '#include <common>\n' + VARYINGS);
  };
  mat.customProgramCacheKey = () => 'coro-v2';
  return mat;
}
