// eclipses.ts — the nearest solar or lunar eclipse before or after a moment.
//
// Method, for both kinds:
//   1. find the next new (or full) Moon exactly, from the apparent longitudes;
//   2. skip it unless the Moon is near enough to a node for an eclipse;
//   3. find GREATEST ECLIPSE: for the Sun, the instant the shadow axis passes
//      closest to the Earth's centre; for the Moon, the instant the Moon's
//      centre passes closest to the axis of the Earth's shadow;
//   4. classify it from the cone geometry at that instant.
//
// Greatest SOLAR eclipse is the instant the shadow axis passes nearest the
// Earth's centre, in true distance. Whether the axis or the shadow cones then
// meet the Earth is judged on the Earth's true figure: with the polar axis
// stretched by 1/(1 − f), which turns the ellipsoid into a sphere of the
// equatorial radius. (Timing the minimum in the stretched space instead moved
// high-latitude eclipses by up to half a minute.)
//
// LUNAR eclipses use Danjon's rule for the Earth's shadow: the Earth's radius
// enlarged by 1/100 for its atmosphere when the Moon's parallax is taken. It
// classifies all 1,446 lunar eclipses of 1800–2400 exactly as Swiss Ephemeris
// does; Chauvenet's older 1/50 enlargement of the whole shadow misclassifies
// six. The rule decides which faint penumbral grazes count at all.
//
// The Moon's radius is the eclipse-prediction standard, k = 0.2725076 Earth
// radii (IAU), which also decides a few borderline hybrids.
//
// Swiss Ephemeris's return-flag bits are reproduced, so callers that decode
// them keep working.

import { norm360 } from "./angles.ts";

export const SE_ECL_CENTRAL = 1, SE_ECL_NONCENTRAL = 2, SE_ECL_TOTAL = 4, SE_ECL_ANNULAR = 8,
  SE_ECL_PARTIAL = 16, SE_ECL_ANNULAR_TOTAL = 32, SE_ECL_PENUMBRAL = 64;

export type EclipseCategory = "total" | "annular" | "hybrid" | "partial" | "penumbral";

export interface Eclipse {
  /** Greatest eclipse, as a UT Julian Day (and as TT). */
  jd: number;
  jdTT: number;
  category: EclipseCategory;
  /** Swiss Ephemeris return-flag bits for the same eclipse. */
  flags: number;
  /** Solar: the shadow axis met the Earth. */
  central: boolean;
  /**
   * Solar: the axis's closest approach to the Earth's centre, in Earth radii
   * (gamma). Lunar: the Moon's distance from the shadow axis, in umbral radii.
   */
  gamma: number;
}

type Vec = [number, number, number];

/** What the eclipse search needs from the engine. */
export interface SkySource {
  /** Apparent geocentric Sun and Moon at TT: km, true equator and equinox of date. */
  sunMoon(jdTT: number): { sun: Vec; moon: Vec; sunLon: number; moonLon: number; moonLat: number };
  deltaT(jdUT: number): number;
}

const R_EARTH = 6378.137;
const FLATTENING = 1 / 298.257223563;
const R_SUN = 696000;
const R_MOON = 0.2725076 * R_EARTH;
const SYNODIC = 29.530588853;
const ELONGATION_RATE = 360 / SYNODIC; // °/day, on average
const D2R = Math.PI / 180;

const sub = (a: Vec, b: Vec): Vec => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const dot = (a: Vec, b: Vec) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const len = (a: Vec) => Math.sqrt(dot(a, a));
const signed = (d: number) => {
  const r = norm360(d);
  return r > 180 ? r - 360 : r;
};

/** The syzygy (new Moon at 0°, full Moon at 180°) nearest after or before `t`, in TT. */
function syzygy(sky: SkySource, t: number, target: 0 | 180, backward: boolean): number {
  const elong = (tt: number) => {
    const p = sky.sunMoon(tt);
    return norm360(p.moonLon - p.sunLon);
  };
  const e0 = elong(t);
  let x = backward ? t - norm360(e0 - target) / ELONGATION_RATE : t + norm360(target - e0) / ELONGATION_RATE;
  for (let i = 0; i < 12; i++) {
    const f = signed(elong(x) - target);
    const step = f / ELONGATION_RATE;
    x -= step;
    if (Math.abs(step) < 1e-9) break;
  }
  return x;
}

/** Minimise a unimodal function on [a, b] (golden section, to ~0.05 s). */
function minimise(f: (t: number) => number, a: number, b: number): number {
  const g = (Math.sqrt(5) - 1) / 2;
  let c = b - g * (b - a), d = a + g * (b - a);
  let fc = f(c), fd = f(d);
  while (b - a > 5e-7) {
    if (fc < fd) {
      b = d; d = c; fd = fc; c = b - g * (b - a); fc = f(c);
    } else {
      a = c; c = d; fc = fd; d = a + g * (b - a); fd = f(d);
    }
  }
  return (a + b) / 2;
}

/**
 * The Moon's shadow at an instant. `stretch` turns the Earth's ellipsoid into a
 * sphere of the equatorial radius (for deciding what meets the Earth); without
 * it the distances are true (for timing greatest eclipse).
 */
function moonShadow(sky: SkySource, tt: number, stretch: boolean) {
  const { sun, moon } = sky.sunMoon(tt);
  const k = stretch ? 1 / (1 - FLATTENING) : 1;
  const S: Vec = [sun[0], sun[1], sun[2] * k];
  const M: Vec = [moon[0], moon[1], moon[2] * k];
  const sm = sub(M, S);
  const L = len(sub(moon, sun));
  const d: Vec = [sm[0] / len(sm), sm[1] / len(sm), sm[2] / len(sm)]; // Sun → Moon → onward
  const s = -dot(M, d); // from the Moon to the point of the axis nearest the Earth's centre
  const near: Vec = [M[0] + s * d[0], M[1] + s * d[1], M[2] + s * d[2]];
  const D = len(near);
  const tanUmbra = (R_SUN - R_MOON) / L;
  const tanPenumbra = (R_SUN + R_MOON) / L;
  return { D, s, l1: R_MOON + s * tanPenumbra, l2: R_MOON - s * tanUmbra, tanUmbra };
}

function solarAt(sky: SkySource, tt: number): Eclipse | null {
  const best = minimise((x) => moonShadow(sky, x, false).D, tt - 0.25, tt + 0.25);
  // Whether any eclipse happens at all is decided on a sphere of the
  // equatorial radius, as NASA's eclipse canon and Swiss Ephemeris decide it:
  // the polar grazes this admits (1935 January 5, magnitude 0.001) are in
  // every catalogue. Centrality and type are decided on the true figure.
  const sphere = moonShadow(sky, best, false);
  if (sphere.D > R_EARTH + sphere.l1) return null; // even the penumbra misses
  const { D, s, l2, tanUmbra } = moonShadow(sky, best, true);
  const gamma = D / R_EARTH;
  let category: EclipseCategory, flags: number;
  if (D < R_EARTH) {
    // Where the axis meets the surface, nearer the Moon than the fundamental plane.
    const sSurface = s - Math.sqrt(R_EARTH * R_EARTH - D * D);
    const umbraThere = R_MOON - sSurface * tanUmbra;
    if (umbraThere > 0 && l2 > 0) [category, flags] = ["total", SE_ECL_TOTAL | SE_ECL_CENTRAL];
    else if (umbraThere > 0) [category, flags] = ["hybrid", SE_ECL_ANNULAR_TOTAL | SE_ECL_CENTRAL];
    else [category, flags] = ["annular", SE_ECL_ANNULAR | SE_ECL_CENTRAL];
  } else if (D < R_EARTH + Math.abs(l2)) {
    // The umbra or antumbra grazes the Earth while the axis misses it.
    [category, flags] = l2 > 0 ? ["total", SE_ECL_TOTAL | SE_ECL_NONCENTRAL] : ["annular", SE_ECL_ANNULAR | SE_ECL_NONCENTRAL];
  } else {
    [category, flags] = ["partial", SE_ECL_PARTIAL];
  }
  return { jd: best, jdTT: best, category, flags, central: D < R_EARTH, gamma };
}

/** The Moon against the Earth's shadow at an instant: angles in radians. */
function earthShadow(sky: SkySource, tt: number) {
  const { sun, moon } = sky.sunMoon(tt);
  const ds = len(sun), dm = len(moon);
  const anti: Vec = [-sun[0] / ds, -sun[1] / ds, -sun[2] / ds];
  const theta = Math.acos(Math.min(1, dot(anti, moon) / dm));
  const parallaxMoon = Math.asin(R_EARTH / dm);
  const parallaxSun = Math.asin(R_EARTH / ds);
  const sunRadius = Math.asin(R_SUN / ds);
  const moonRadius = Math.asin(R_MOON / dm);
  // Danjon: the Earth's radius enlarged 1/100 for its atmosphere, in the parallax.
  const umbra = 1.01 * parallaxMoon + parallaxSun - sunRadius;
  const penumbra = 1.01 * parallaxMoon + parallaxSun + sunRadius;
  return { theta, umbra, penumbra, moonRadius };
}

function lunarAt(sky: SkySource, tt: number): Eclipse | null {
  const best = minimise((x) => earthShadow(sky, x).theta, tt - 0.25, tt + 0.25);
  const { theta, umbra, penumbra, moonRadius } = earthShadow(sky, best);
  let category: EclipseCategory, flags: number;
  if (theta + moonRadius < umbra) [category, flags] = ["total", SE_ECL_TOTAL];
  else if (theta - moonRadius < umbra) [category, flags] = ["partial", SE_ECL_PARTIAL];
  else if (theta - moonRadius < penumbra) [category, flags] = ["penumbral", SE_ECL_PENUMBRAL];
  else return null;
  return { jd: best, jdTT: best, category, flags, central: false, gamma: theta / umbra };
}

/** Moon's ecliptic latitude limits beyond which no eclipse is possible at a syzygy. */
const SOLAR_LIMIT = 1.8;
const LUNAR_LIMIT = 1.8;

function search(sky: SkySource, jdUT: number, backward: boolean, kind: "solar" | "lunar"): Eclipse | null {
  const target = kind === "solar" ? 0 : 180;
  // Start half a day back (or on): greatest eclipse can fall an hour or two
  // either side of the syzygy, so the syzygy just behind the start time may
  // still own the next eclipse. The order check below drops it if not.
  let t = jdUT + sky.deltaT(jdUT) / 86400 + (backward ? 0.5 : -0.5);
  for (let i = 0; i < 40; i++) {
    const s = syzygy(sky, t, target, backward);
    const lat = Math.abs(sky.sunMoon(s).moonLat);
    if (lat < (kind === "solar" ? SOLAR_LIMIT : LUNAR_LIMIT)) {
      const e = kind === "solar" ? solarAt(sky, s) : lunarAt(sky, s);
      if (e) {
        // Report in UT; ΔT changes by microseconds across the search.
        const ut = e.jdTT - sky.deltaT(e.jdTT) / 86400;
        const inOrder = backward ? ut < jdUT : ut > jdUT;
        if (inOrder) return { ...e, jd: ut };
      }
    }
    t = s + (backward ? -1 : 1) * (SYNODIC / 2); // well clear of this syzygy, short of the next
  }
  return null;
}

/** Nearest global solar eclipse after (or, with `backward`, before) a UT Julian Day. */
export function solarEclipse(sky: SkySource, jdUT: number, backward: boolean): Eclipse | null {
  return search(sky, jdUT, backward, "solar");
}

/** Nearest lunar eclipse after (or, with `backward`, before) a UT Julian Day. */
export function lunarEclipse(sky: SkySource, jdUT: number, backward: boolean): Eclipse | null {
  return search(sky, jdUT, backward, "lunar");
}

/** The one word a reader wants, from Swiss-style flags. */
export function categoryOf(flags: number): EclipseCategory | "unknown" {
  if (flags & SE_ECL_ANNULAR_TOTAL) return "hybrid";
  if (flags & SE_ECL_TOTAL) return "total";
  if (flags & SE_ECL_ANNULAR) return "annular";
  if (flags & SE_ECL_PENUMBRAL) return "penumbral";
  if (flags & SE_ECL_PARTIAL) return "partial";
  return "unknown";
}

/** Exposed for tests: the elongation-based syzygy finder. */
export const _internal = { syzygy, D2R };
