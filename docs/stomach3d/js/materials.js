// materials.js — shaders for the gastric surfaces.
//
// The baked vertex attributes do the work (see build_stomach_fields.py):
//   _field0 = (region_id, s_norm, wall_depth_01, concavity)
//   _field1 = (luminal_pH, oxyntic_density_01, theta, surface_class)
//
// region_id 1 cardia · 2 fundus · 3 body · 4 antrum · 5 pylorus · 0 neighbour
import * as THREE from 'three';

// Deepened and spread apart deliberately: at the first render the five regions
// were all pale pink and, under the environment map plus a bright key light,
// read as one undifferentiated white surface.  Regional partition is a headline
// result of this page, so the five tints have to be separable on screen.
export const REGION_COLOUR = [
  [0.34, 0.37, 0.42],   // 0 neighbour organ — neutral grey
  [0.78, 0.60, 0.34],   // 1 cardia         — amber
  [0.80, 0.33, 0.30],   // 2 fundus         — deep red
  [0.62, 0.26, 0.30],   // 3 body           — darker red
  [0.52, 0.34, 0.56],   // 4 antrum         — violet
  [0.30, 0.40, 0.62],   // 5 pylorus        — blue
];

// pH colour ramp: acid (red) -> neutral (blue-green).  Used by the section and
// gland views; the legend on the page states the endpoints.
export function phColour(ph) {
  const t = Math.min(Math.max((ph - 1.0) / 6.4, 0), 1);
  return new THREE.Color().setHSL(0.02 + 0.52 * t, 0.78 - 0.22 * t, 0.42 + 0.14 * t);
}

// `attribute` is ILLEGAL in a fragment shader (and a reserved word under
// GLSL ES 3.0), so the two stages get different declaration blocks.  Injecting
// one shared block into both silently kills the fragment shader: three.js then
// renders with an unlinked program, the surface falls back to flat shading, and
// the lost fragment work makes the frame budget look better than it is.  The
// headless audit asserts LINK_STATUS on every program for exactly this reason.
const VERT_DEFS = /* glsl */`
  attribute vec4 _field0;
  attribute vec4 _field1;
  varying vec4 vF0;
  varying vec4 vF1;
  varying vec3 vN;
  varying vec3 vWorld;
`;
const FRAG_DEFS = /* glsl */`
  varying vec4 vF0;
  varying vec4 vF1;
  varying vec3 vN;
  varying vec3 vWorld;
`;

const COMMON_VERT = /* glsl */`
  vF0 = _field0; vF1 = _field1;
  vN = normalize(normalMatrix * normal);
`;

/** Wall / serosa / mucosa material.
 *  mode 0 = anatomical region tint, 1 = luminal pH, 2 = oxyntic gland density */
export function wallMaterial(opts = {}) {
  const u = {
    uMode:       { value: 0 },
    uOpacity:    { value: opts.opacity ?? 1.0 },
    uRegionCol:  { value: REGION_COLOUR.map(c => new THREE.Vector3(...c)) },
    uHighlight:  { value: -1 },
    uTime:       { value: 0 },
    uCut:        { value: 0 },
    uPhLo:       { value: 1.5 },
    uPhSurf:     { value: 7.0 },
  };
  const m = new THREE.MeshStandardMaterial({
    color: 0xffffff, roughness: 0.78, metalness: 0.0,
    transparent: (opts.opacity ?? 1.0) < 1.0,
    opacity: opts.opacity ?? 1.0,
    side: opts.side ?? THREE.FrontSide,
    flatShading: false,
  });
  m.userData.u = u;
  m.onBeforeCompile = (sh) => {
    Object.assign(sh.uniforms, u);
    sh.vertexShader = sh.vertexShader
      .replace('#include <common>', '#include <common>\n' + VERT_DEFS)
      .replace('#include <begin_vertex>',
        '#include <begin_vertex>\n' + COMMON_VERT + '\n  vWorld = (modelMatrix * vec4(transformed,1.0)).xyz;');
    sh.fragmentShader = sh.fragmentShader
      .replace('#include <common>', '#include <common>\n' + FRAG_DEFS + `
        uniform int  uMode;
        uniform vec3 uRegionCol[6];
        uniform int  uHighlight;
        uniform float uOpacity, uPhLo, uPhSurf;
        vec3 phRamp(float ph){
          float t = clamp((ph - 1.0) / 6.4, 0.0, 1.0);
          return mix(vec3(0.86,0.20,0.16), vec3(0.22,0.62,0.72), t);
        }`)
      .replace('#include <color_fragment>', /* glsl */`
        #include <color_fragment>
        int rid = int(vF0.x + 0.5);
        vec3 base = uRegionCol[rid < 0 ? 0 : (rid > 5 ? 0 : rid)];
        if (uMode == 1) {
          base = phRamp(vF1.x);
        } else if (uMode == 2) {
          base = mix(vec3(0.30,0.32,0.36), vec3(0.95,0.72,0.30), vF1.y);
        }
        // the measured concavity makes the incisura read as a real groove
        float groove = clamp(-vF0.w, 0.0, 1.0);
        base *= (1.0 - 0.38 * groove);
        // mucosal surfaces are wetter and lighter than serosa
        if (vF1.w > 0.5 && vF1.w < 1.5) base = mix(base, base * 1.10 + 0.03, 0.45);
        if (uHighlight >= 0 && rid == uHighlight) base = mix(base, vec3(1.0,0.92,0.55), 0.42);
        diffuseColor.rgb = base;
        diffuseColor.a *= uOpacity;
      `);
  };
  return m;
}

export function vesselMaterial(colour, emissive = 0.18) {
  return new THREE.MeshStandardMaterial({
    color: colour, roughness: 0.42, metalness: 0.08,
    emissive: new THREE.Color(colour).multiplyScalar(emissive),
  });
}

export function particleMaterial(colour, size) {
  return new THREE.PointsMaterial({
    color: colour, size, sizeAttenuation: true,
    transparent: true, opacity: 0.95, depthWrite: false,
    blending: THREE.AdditiveBlending,
  });
}
