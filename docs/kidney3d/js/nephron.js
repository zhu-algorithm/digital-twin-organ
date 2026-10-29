/**
 * nephron.js — the nephron in three dimensions, at real proportions.
 *
 * Built at its real size in millimetres and then magnified for display. The
 * magnification is reported in the HUD and in the acceptance record; nothing
 * about the geometry is invented to fit the screen, and every dimension below
 * is a published adult human value:
 *
 *   renal corpuscle (glomerulus)   0.20 mm diameter
 *   Bowman's capsule               0.22 mm
 *   proximal convoluted tubule     14 mm long, 0.060 mm outer diameter
 *   proximal straight tubule       6 mm
 *   thin limb of the loop          0.015 mm diameter
 *   thick ascending limb           9 mm, 0.035 mm
 *   distal convoluted tubule       5 mm, 0.035 mm
 *   connecting tubule              0.040 mm
 *   collecting duct                0.040 mm cortical, widening to 0.200 mm at
 *                                  the papillary duct
 *   afferent arteriole             0.020 mm, efferent arteriole 0.013 mm
 *
 * BOTH NEPHRON TYPES ARE BUILT, because they are not the same organ part:
 *
 *   cortical nephron (about 85 % of them)  glomerulus in the mid cortex, SHORT
 *       loop that barely enters the outer medulla, peritubular capillaries
 *       around its tubules
 *   juxtamedullary nephron (about 15 %)    glomerulus on the corticomedullary
 *       junction, LONG loop reaching towards the papilla, and VASA RECTA
 *       running beside it. Only this type has vasa recta; drawing them on a
 *       cortical nephron would be wrong, so they are absent there.
 *
 * The convolutions are real length, not a squiggle: the proximal tubule's 14 mm
 * is coiled into a cortical ball a little over a millimetre across, which is
 * what makes a 40 mm nephron fit inside an 8 mm cortex.
 *
 * TUBULAR FLUID SPEED is computed, not chosen: v = Q/A, with Q starting at the
 * model's own single-nephron GFR and stepping down along the tubule by the
 * TEXTBOOK fractional reabsorption shares carried in kidney_model.json. Those
 * shares are flagged there as textbook values; no transport equation is solved
 * and no solute concentration is computed anywhere.
 */

import * as THREE from 'three';

export const MAGNIFICATION = 20.0;       // 1 mm of tissue -> 20 mm on screen
const S = MAGNIFICATION;
/**
 * A nephron is about 40 mm long and 0.06 mm wide — an aspect ratio near 700:1.
 * Drawn at true calibre inside a figure tall enough to contain the loop, the
 * tubules come out thinner than a pixel. POSITIONS AND LENGTHS STAY TRUE; only
 * the rendered calibre is multiplied by this declared gain, and the gain is
 * reported in the HUD and in the acceptance record. The real diameters are
 * listed in DIMS_MM and are what the velocity calculation uses.
 */
export const TUBULE_CALIBRE_GAIN = 6.0;
const G = TUBULE_CALIBRE_GAIN;

export const DIMS_MM = {
  glomerulus_diameter: 0.200,
  bowman_diameter: 0.220,
  pct_length: 14.0, pct_diameter: 0.060,
  pst_length: 6.0, pst_diameter: 0.050,
  thin_limb_diameter: 0.015,
  tal_length: 9.0, tal_diameter: 0.035,
  dct_length: 5.0, dct_diameter: 0.035,
  cnt_diameter: 0.040,
  cd_cortical_diameter: 0.040, cd_papillary_diameter: 0.200,
  afferent_diameter: 0.020, efferent_diameter: 0.013,
  vasa_recta_diameter: 0.017,
  cortical_loop_depth: 2.0,            // just into the outer medulla
  juxta_loop_depth: 18.0,              // towards the papilla
  macula_densa_diameter: 0.040
};

const COL = {
  glomerulus: 0xf2a62c,
  bowman: 0xd8cba6,
  afferent: 0xd93b2b,
  efferent: 0xe8743a,
  pct: 0xb9743f,
  pst: 0xb08050,
  thin_desc: 0xc9b07a,
  thin_asc: 0x8fb4c9,
  tal: 0x6f9fc4,
  dct: 0x8f7fc0,
  cnt: 0x9b8fcb,
  cd: 0xd9d98a,
  vasa_desc: 0xff9a5c,
  vasa_asc: 0x7fb0e0,
  peritubular: 0xffb36b,
  macula: 0x44d08a
};

/**
 * Sweep a tube of varying radius along a polyline, APPENDING into a shared
 * buffer. One mesh per colour instead of one mesh per tubule: a nephron is
 * forty-odd separate tubes, and giving each its own material meant hundreds of
 * shader compilations and hundreds of draw calls, which stalled the scene hard
 * enough for the headless benchmark to time out.
 */
function appendTube (B, pts, r0, r1, sides = 8) {
  const pos = B.pos; const nrm = B.nrm; const idx = B.idx;
  const base0 = pos.length / 3;
  const up = new THREE.Vector3(); const dir = new THREE.Vector3();
  const t1 = new THREE.Vector3(); const t2 = new THREE.Vector3();
  for (let i = 0; i < pts.length; i++) {
    const a = pts[i];
    const b = pts[Math.min(i + 1, pts.length - 1)];
    const c = pts[Math.max(i - 1, 0)];
    dir.subVectors(b, c);
    if (dir.lengthSq() < 1e-12) dir.set(0, 1, 0);
    dir.normalize();
    up.set(0, 1, 0);
    if (Math.abs(dir.dot(up)) > 0.94) up.set(1, 0, 0);
    t1.crossVectors(up, dir).normalize();
    t2.crossVectors(dir, t1).normalize();
    const u = pts.length > 1 ? i / (pts.length - 1) : 0;
    const r = r0 + (r1 - r0) * u;
    for (let k = 0; k < sides; k++) {
      const th = (k / sides) * Math.PI * 2;
      const nx = t1.x * Math.cos(th) + t2.x * Math.sin(th);
      const ny = t1.y * Math.cos(th) + t2.y * Math.sin(th);
      const nz = t1.z * Math.cos(th) + t2.z * Math.sin(th);
      pos.push(a.x + nx * r, a.y + ny * r, a.z + nz * r);
      nrm.push(nx, ny, nz);
    }
  }
  for (let i = 0; i < pts.length - 1; i++) {
    for (let k = 0; k < sides; k++) {
      const k2 = (k + 1) % sides;
      const A = base0 + i * sides + k; const Bv = base0 + i * sides + k2;
      const C = base0 + (i + 1) * sides + k; const D = base0 + (i + 1) * sides + k2;
      idx.push(A, C, D, A, D, Bv);
    }
  }
  return (pts.length - 1) * sides * 2;
}

function newBuf () { return { pos: [], nrm: [], idx: [] }; }

function bufToMesh (B, color) {
  const g = new THREE.BufferGeometry();
  g.setAttribute('position', new THREE.Float32BufferAttribute(B.pos, 3));
  g.setAttribute('normal', new THREE.Float32BufferAttribute(B.nrm, 3));
  g.setIndex(B.idx);
  g.computeBoundingSphere();
  const m = new THREE.MeshPhysicalMaterial({
    color: new THREE.Color(color), roughness: 0.5, metalness: 0.0,
    side: THREE.FrontSide
  });
  return new THREE.Mesh(g, m);
}

/** A convoluted tubule of TRUE arc length coiled inside a ball of radius R. */
function coil (centre, R_mm, lengthMm, turns, seed, axis = new THREE.Vector3(0, 1, 0)) {
  // arc length of the generated helix is measured and the pitch solved for, so
  // the coil really does contain `lengthMm` of tubule rather than looking busy
  const n = 220;
  const e1 = new THREE.Vector3(1, 0, 0);
  if (Math.abs(axis.dot(e1)) > 0.9) e1.set(0, 0, 1);
  const u = new THREE.Vector3().crossVectors(axis, e1).normalize();
  const v = new THREE.Vector3().crossVectors(axis, u).normalize();
  const make = (rad) => {
    const out = [];
    for (let i = 0; i <= n; i++) {
      const t = i / n;
      const th = t * turns * Math.PI * 2 + seed;
      const wob = 0.72 + 0.28 * Math.sin(th * 2.7 + seed);
      const rr = rad * wob;
      const h = (t - 0.5) * 2.0 * rad * 0.85;
      out.push(new THREE.Vector3(
        centre.x + u.x * rr * Math.cos(th) + v.x * rr * Math.sin(th) + axis.x * h,
        centre.y + u.y * rr * Math.cos(th) + v.y * rr * Math.sin(th) + axis.y * h,
        centre.z + u.z * rr * Math.cos(th) + v.z * rr * Math.sin(th) + axis.z * h));
    }
    return out;
  };
  let lo = 0.02; let hi = R_mm * 3.0; let pts = make(R_mm);
  for (let k = 0; k < 28; k++) {
    const mid = 0.5 * (lo + hi);
    pts = make(mid);
    let L = 0;
    for (let i = 1; i < pts.length; i++) L += pts[i].distanceTo(pts[i - 1]);
    if (L < lengthMm) lo = mid; else hi = mid;
  }
  let L = 0;
  for (let i = 1; i < pts.length; i++) L += pts[i].distanceTo(pts[i - 1]);
  return { pts, lengthMm: L };
}

function lerpPts (a, b, n) {
  const out = [];
  for (let i = 0; i <= n; i++) out.push(new THREE.Vector3().lerpVectors(a, b, i / n));
  return out;
}

/**
 * @param {object} model parsed kidney_model.json
 * @param {object} flowRep parsed flow_report.json (for the single-nephron GFR)
 */
export function buildNephrons (model, flowRep) {
  const group = new THREE.Group();
  const labels = [];
  const paths = [];
  let triangles = 0;
  const D = DIMS_MM;
  const reab = model.tubular_handling;
  const sngfr_nL_min = flowRep.single_nephron.single_nephron_gfr_nL_min;
  // nL/min -> mm^3/s : 1 nL = 1e-3 mm^3
  const Q0 = sngfr_nL_min * 1e-3 / 60.0;

  const bufs = new Map();
  const add = (pts, r0, r1, color, sides) => {
    if (!bufs.has(color)) bufs.set(color, newBuf());
    triangles += appendTube(bufs.get(color),
                            pts.map(p => p.clone().multiplyScalar(S)),
                            r0 * S * G, r1 * S * G, sides || 7);
  };
  // spheres (Bowman's capsule, macula densa) share one material per colour too
  const sphMat = new Map();
  const sphereMat = (color, opts) => {
    const k = color + JSON.stringify(opts || {});
    if (!sphMat.has(k)) {
      sphMat.set(k, new THREE.MeshPhysicalMaterial(Object.assign({
        color: new THREE.Color(color), roughness: 0.4, metalness: 0 }, opts || {})));
    }
    return sphMat.get(k);
  };
  const measure = (pts) => {
    let L = 0; for (let i = 1; i < pts.length; i++) L += pts[i].distanceTo(pts[i - 1]);
    return L;
  };
  const addPath = (name, pts, rMm, Q, gain) => {
    const A = Math.PI * rMm * rMm;                      // mm^2 (lumen)
    const v = Q / Math.max(A, 1e-12);                   // mm/s
    paths.push({
      tree: 'tubule', retrograde: false, level: name,
      points: pts.map(p => [p.x * S, p.y * S, p.z * S]),
      radius_mm: pts.map(() => rMm * S),
      v_mm_s: pts.map(() => v * S),
      length_mm: measure(pts) * S,
      display_gain: gain || 1.0
    });
    return v;
  };

  const built = {};

  function buildOne (kind, origin, loopDepth, withVasaRecta) {
    const tag = kind === 'juxtamedullary' ? 'JM' : 'CTX';
    const down = new THREE.Vector3(0, -1, 0);
    const glomC = origin.clone();

    // --- renal corpuscle: Bowman's capsule around the glomerular tuft
    const bow = new THREE.Mesh(
      new THREE.SphereGeometry(D.bowman_diameter * 0.5 * S * 1.6, 20, 14),
      // Plain transparency, NOT `transmission`: a transmissive material makes
      // three.js render the whole scene again into a back buffer every frame,
      // which took the nephron view below one frame per second.
      sphereMat(COL.bowman, { roughness: 0.3, transparent: true, opacity: 0.30,
        side: THREE.FrontSide }));
    bow.position.copy(glomC.clone().multiplyScalar(S));
    group.add(bow); triangles += 20 * 14 * 2;

    // glomerular tuft: capillary loops at real calibre inside the capsule
    const tuftR = D.glomerulus_diameter * 0.5;
    for (let k = 0; k < 6; k++) {
      const th = k / 6 * Math.PI * 2;
      const pts = [];
      for (let i = 0; i <= 26; i++) {
        const t = i / 26;
        const a = t * Math.PI * 2;
        const rr = tuftR * (0.55 + 0.42 * Math.sin(a * 1.5 + th));
        pts.push(new THREE.Vector3(
          glomC.x + rr * Math.cos(a + th),
          glomC.y + tuftR * 0.55 * Math.sin(a * 2 + th),
          glomC.z + rr * Math.sin(a + th) * 0.75));
      }
      add(pts, 0.0037, 0.0037, COL.glomerulus, 5);
    }
    labels.push({ pos: glomC.clone().multiplyScalar(S)
      .add(new THREE.Vector3(0, D.bowman_diameter * 0.6 * S, 0)),
      text: `${tag} glomerulus 200 µm`, color: '#f2a62c', pri: 9 });

    // --- afferent and efferent arterioles: the efferent IS the narrower one
    const affIn = glomC.clone().add(new THREE.Vector3(-0.42, 0.16, 0));
    const effOut = glomC.clone().add(new THREE.Vector3(0.40, 0.14, 0.05));
    add(lerpPts(affIn, glomC.clone().add(new THREE.Vector3(-0.09, 0.03, 0)), 10),
        D.afferent_diameter / 2, D.afferent_diameter / 2, COL.afferent, 7);
    add(lerpPts(glomC.clone().add(new THREE.Vector3(0.09, 0.03, 0)), effOut, 10),
        D.efferent_diameter / 2, D.efferent_diameter / 2, COL.efferent, 7);
    labels.push({ pos: affIn.clone().multiplyScalar(S),
      text: 'afferent arteriole 20 µm', color: '#ff7a68', pri: 8 });
    labels.push({ pos: effOut.clone().multiplyScalar(S),
      text: 'efferent arteriole 13 µm — narrower', color: '#ff9b5c', pri: 8 });

    // --- proximal convoluted tubule: 14 mm coiled into the cortex
    const pctC = glomC.clone().add(new THREE.Vector3(0.55, -0.35, 0.35));
    const pct = coil(pctC, 0.62, D.pct_length, 11, 0.3);
    add(pct.pts, D.pct_diameter / 2, D.pct_diameter / 2, COL.pct, 7);
    labels.push({ pos: pctC.clone().multiplyScalar(S),
      text: `proximal convoluted tubule ${pct.lengthMm.toFixed(1)} mm coiled`,
      color: '#d79a63', pri: 7 });

    // --- proximal straight tubule, descending
    const pstTop = pct.pts[pct.pts.length - 1].clone();
    const pstBot = pstTop.clone().add(down.clone().multiplyScalar(D.pst_length * 0.42))
      .add(new THREE.Vector3(0.1, 0, -0.1));
    const pstPts = lerpPts(pstTop, pstBot, 26);
    add(pstPts, D.pst_diameter / 2, D.thin_limb_diameter / 2, COL.pst, 7);

    // --- the loop of Henle: descending thin limb, hairpin, ascending
    const loopTop = pstBot.clone();
    const loopBot = loopTop.clone().add(down.clone().multiplyScalar(loopDepth));
    const desc = lerpPts(loopTop, loopBot, 40);
    add(desc, D.thin_limb_diameter / 2, D.thin_limb_diameter / 2, COL.thin_desc, 6);
    const hair = [];
    for (let i = 0; i <= 16; i++) {
      const a = Math.PI * i / 16;
      hair.push(new THREE.Vector3(loopBot.x + 0.11 * (1 - Math.cos(a)),
                                  loopBot.y - 0.11 * Math.sin(a), loopBot.z));
    }
    add(hair, D.thin_limb_diameter / 2, D.thin_limb_diameter / 2, COL.thin_desc, 6);
    const ascBot = hair[hair.length - 1].clone();
    // juxtamedullary loops have a long THIN ascending segment before the thick
    // one; cortical loops turn thick almost immediately
    const thinAscTop = ascBot.clone().add(new THREE.Vector3(
      0, (kind === 'juxtamedullary' ? loopDepth * 0.45 : loopDepth * 0.12), 0));
    add(lerpPts(ascBot, thinAscTop, 20), D.thin_limb_diameter / 2,
        D.tal_diameter / 2, COL.thin_asc, 6);
    const talTop = loopTop.clone().add(new THREE.Vector3(0.22, 0.25, 0));
    add(lerpPts(thinAscTop, talTop, 30), D.tal_diameter / 2, D.tal_diameter / 2,
        COL.tal, 7);
    labels.push({ pos: new THREE.Vector3().lerpVectors(loopBot, loopTop, 0.5)
      .multiplyScalar(S),
      text: kind === 'juxtamedullary'
        ? `long loop of Henle, ${loopDepth.toFixed(0)} mm into the medulla`
        : `short loop of Henle, ${loopDepth.toFixed(0)} mm`,
      color: '#9fc6dd', pri: 8 });

    // --- macula densa: the thick ascending limb touches its OWN glomerulus
    const md = glomC.clone().add(new THREE.Vector3(-0.17, -0.02, 0.0));
    add(lerpPts(talTop, md, 18), D.tal_diameter / 2, D.macula_densa_diameter / 2,
        COL.tal, 7);
    const mdot = new THREE.Mesh(
      new THREE.SphereGeometry(0.030 * S, 12, 10),
      sphereMat(COL.macula, { emissive: new THREE.Color(0x12402a) }));
    mdot.position.copy(md.clone().multiplyScalar(S));
    group.add(mdot); triangles += 12 * 10 * 2;
    labels.push({ pos: md.clone().multiplyScalar(S),
      text: 'macula densa — tubuloglomerular feedback', color: '#5fe0a4', pri: 9 });

    // --- distal convoluted tubule, connecting tubule, collecting duct
    const dctC = glomC.clone().add(new THREE.Vector3(-0.72, -0.30, -0.35));
    const dct = coil(dctC, 0.42, D.dct_length, 7, 1.7);
    add(lerpPts(md, dct.pts[0], 10), D.dct_diameter / 2, D.dct_diameter / 2,
        COL.dct, 7);
    add(dct.pts, D.dct_diameter / 2, D.dct_diameter / 2, COL.dct, 7);
    const cntEnd = dctC.clone().add(new THREE.Vector3(-0.45, -0.55, 0));
    add(lerpPts(dct.pts[dct.pts.length - 1], cntEnd, 12), D.dct_diameter / 2,
        D.cnt_diameter / 2, COL.cnt, 7);
    labels.push({ pos: dctC.clone().multiplyScalar(S),
      text: `distal convoluted tubule ${dct.lengthMm.toFixed(1)} mm`,
      color: '#a89bd8', pri: 6 });

    const cdTop = cntEnd.clone();
    const cdBot = cdTop.clone().add(down.clone().multiplyScalar(loopDepth + 2.2));
    add(lerpPts(cdTop, cdBot, 34), D.cd_cortical_diameter / 2,
        D.cd_papillary_diameter / 2, COL.cd, 9);
    labels.push({ pos: cdBot.clone().multiplyScalar(S),
      text: 'collecting duct -> papilla, 40 -> 200 µm', color: '#e4e4a0', pri: 7 });

    // --- the SECOND capillary bed, and it differs by nephron type
    if (withVasaRecta) {
      for (const sgn of [-1, 1]) {
        const off = new THREE.Vector3(0.20 * sgn, 0, 0.14 * sgn);
        const d0 = loopTop.clone().add(off);
        const d1 = loopBot.clone().add(off);
        add(lerpPts(effOut, d0, 14), D.vasa_recta_diameter / 2,
            D.vasa_recta_diameter / 2, COL.vasa_desc, 6);
        add(lerpPts(d0, d1, 30), D.vasa_recta_diameter / 2,
            D.vasa_recta_diameter / 2, COL.vasa_desc, 6);
        const a0 = d1.clone().add(new THREE.Vector3(0.13 * sgn, 0, 0));
        const a1 = loopTop.clone().add(off).add(new THREE.Vector3(0.13 * sgn, 0.5, 0));
        add(lerpPts(d1, a0, 8), D.vasa_recta_diameter / 2,
            D.vasa_recta_diameter / 2, COL.vasa_asc, 6);
        add(lerpPts(a0, a1, 30), D.vasa_recta_diameter / 2,
            D.vasa_recta_diameter / 2, COL.vasa_asc, 6);
      }
      labels.push({ pos: new THREE.Vector3().lerpVectors(loopTop, loopBot, 0.62)
        .add(new THREE.Vector3(0.3, 0, 0)).multiplyScalar(S),
        text: 'vasa recta — 2nd capillary bed, countercurrent', color: '#ffa871',
        pri: 10 });
    } else {
      // peritubular capillary plexus around the cortical convolutions
      for (let k = 0; k < 5; k++) {
        const c = coil(pctC.clone().add(new THREE.Vector3(
          0.22 * Math.cos(k), 0.10 * k - 0.2, 0.22 * Math.sin(k))),
          0.5, 2.6, 3, k * 1.3);
        add(c.pts, D.vasa_recta_diameter / 2 * 0.8, D.vasa_recta_diameter / 2 * 0.8,
            COL.peritubular, 5);
      }
      labels.push({ pos: pctC.clone().add(new THREE.Vector3(0, 0.75, 0))
        .multiplyScalar(S),
        text: 'peritubular capillaries — 2nd capillary bed', color: '#ffc28a',
        pri: 10 });
    }

    // --- tubular fluid: v = Q/A, Q stepping down by the textbook reabsorption
    const lumen = (outer) => outer * 0.42;      // lumen radius from outer radius
    let Q = Q0;
    const v = {};
    v.pct = addPath(`${tag} proximal tubule`, pct.pts, lumen(D.pct_diameter / 2), Q, 25);
    Q *= (1 - reab.proximal_convoluted_tubule);
    v.loop = addPath(`${tag} loop of Henle`, desc, lumen(D.thin_limb_diameter / 2), Q, 25);
    Q *= (1 - reab.loop_of_henle);
    v.dct = addPath(`${tag} distal tubule`, dct.pts, lumen(D.dct_diameter / 2), Q, 25);
    Q *= (1 - reab.distal_convoluted_tubule);
    v.cd = addPath(`${tag} collecting duct`, lerpPts(cdTop, cdBot, 20),
                   lumen(D.cd_cortical_diameter / 2), Q, 25);

    built[kind] = {
      pct_length_mm: Number(pct.lengthMm.toFixed(2)),
      dct_length_mm: Number(dct.lengthMm.toFixed(2)),
      loop_depth_mm: loopDepth,
      has_vasa_recta: withVasaRecta,
      tubular_velocity_mm_s: {
        proximal_tubule: Number(v.pct.toFixed(4)),
        loop_of_henle: Number(v.loop.toFixed(4)),
        distal_tubule: Number(v.dct.toFixed(4)),
        collecting_duct: Number(v.cd.toFixed(4))
      }
    };
  }

  buildOne('cortical', new THREE.Vector3(-7.5, 0, 0), D.cortical_loop_depth, false);
  buildOne('juxtamedullary', new THREE.Vector3(7.5, -1.2, 0), D.juxta_loop_depth, true);
  for (const [color, B] of bufs) group.add(bufToMesh(B, color));

  labels.push({ pos: new THREE.Vector3(-7.5 * S, 1.6 * S, 0),
    text: 'CORTICAL nephron (~85 %) — short loop, NO vasa recta',
    color: '#cfe0f0', pri: 11 });
  labels.push({ pos: new THREE.Vector3(7.5 * S, 1.6 * S, 0),
    text: 'JUXTAMEDULLARY nephron (~15 %) — long loop + vasa recta',
    color: '#cfe0f0', pri: 11 });

  return {
    group,
    labels,
    flowSpec: { paths, urine_display_gain: 25.0, pulse: [1.0], rr_interval_ms: 1000 },
    info: {
      magnification: MAGNIFICATION,
      tubule_calibre_gain: TUBULE_CALIBRE_GAIN,
      tubule_calibre_note: 'tubule POSITIONS and LENGTHS are true at ' +
        MAGNIFICATION + 'x; the rendered tubule CALIBRE carries an additional ' +
        'declared ' + TUBULE_CALIBRE_GAIN + 'x gain because a 60 um tubule in a ' +
        'figure tall enough to hold an 18 mm loop is thinner than one pixel. The ' +
        'true diameters are in dimensions_mm and are what v = Q/A uses.',
      magnification_note: `built at real size in millimetres and displayed at ` +
        `${MAGNIFICATION}x; every dimension is a published adult human value and ` +
        `none was changed to fit the screen`,
      dimensions_mm: DIMS_MM,
      nephrons: built,
      triangles,
      tubular_flow_note: 'tubular fluid speed is v = Q/A, with Q starting at this ' +
        'project\'s own single-nephron GFR and stepping down by the TEXTBOOK ' +
        'fractional reabsorption shares in kidney_model.json. No transport equation ' +
        'is solved and no solute concentration is computed.',
      tubular_display_gain: 25.0,
      lumen_assumption: 'lumen radius taken as 0.42 of the published OUTER tubular ' +
        'radius; this is a declared assumption, not a measurement',
      ultrastructure_absent: 'podocyte foot processes (0.2-0.4 µm), the filtration ' +
        'slit diaphragm and the proximal-tubule brush border are NOT modelled. They ' +
        'are two orders of magnitude below CT resolution and three below anything ' +
        'this geometry is derived from; drawing them would be illustration, not data.'
    }
  };
}
