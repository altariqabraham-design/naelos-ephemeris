// houses.ts — sidereal time, the angles, and the house systems.
//
// Everything here is spherical geometry on the sky of date: the local meridian
// and horizon, the true equator and the true ecliptic. Where a textbook gives a
// tangent formula with a quadrant to fix up, the code intersects the two great
// circles instead and chooses the intersection by which side of the sky it is
// on, which cannot pick the wrong quadrant. Placidus and Koch, which are
// defined by time rather than by circles, are solved from their definitions.

import { frameOfDate, centuries, delaunay, iauWeight, meanEquationOfOrigins } from "./frames.ts";
import { norm360 } from "./angles.ts";

const D2R = Math.PI / 180;
const R2D = 180 / Math.PI;
const AS2R = Math.PI / 648000;
const TWO_PI = 2 * Math.PI;

type Vec = [number, number, number];

/** Swiss Ephemeris's letters, which Naelos already stores. */
export type HouseSystem = "P" | "K" | "W" | "A" | "E" | "C" | "R" | "O" | "X" | "M";

export interface Houses {
  /** Twelve cusp longitudes, cusp 1 first, on the true ecliptic of date. */
  cusps: number[];
  ascendant: number;
  mc: number;
  /** Right ascension of the meridian (local apparent sidereal time), degrees. */
  armc: number;
  vertex: number;
  /** The ecliptic point rising at the equator: the "East Point". */
  equatorialAscendant: number;
  /** The system the cusps were drawn in. */
  system: HouseSystem;
  /**
   * True when Placidus or Koch was asked for where it cannot exist (inside the
   * polar circles some ecliptic degrees never rise or set) and Porphyry was
   * drawn instead, as Swiss Ephemeris does.
   */
  fallback: boolean;
}

/** Earth rotation angle (IAU 2000), radians, from a UT1 Julian Day. */
export function earthRotationAngle(jdUT: number): number {
  const t = jdUT - 2451545.0;
  const f = jdUT % 1;
  const theta = TWO_PI * (f + 0.779057273264 + 0.00273781191135448 * t);
  return ((theta % TWO_PI) + TWO_PI) % TWO_PI;
}

/**
 * Greenwich apparent sidereal time, radians: mean sidereal time from the Earth
 * rotation angle, plus the equation of the equinoxes (nutation in longitude
 * projected on the equator, and its leading complementary terms).
 *
 * Mean sidereal time is IAU 2006's polynomial inside the DE440 tables. That
 * polynomial is a fit to the present and diverges beyond it, so outside them
 * it is ERA minus the equation of the origins integrated through the
 * long-term precession model (frames.ts), blended in over the same years as
 * the frame itself.
 */
export function apparentSiderealTime(jdUT: number, jdTT: number): number {
  const t = centuries(jdTT);
  const era = earthRotationAngle(jdUT);
  const w = iauWeight(jdTT);
  let gmst = era + (0.014506 + (4612.156534 + (1.3915817 + (-0.00000044 + (-0.000029956 - 0.0000000368 * t) * t) * t) * t) * t) * AS2R;
  if (w < 1) {
    const long = era - meanEquationOfOrigins(jdTT);
    // Move from one to the other by the shortest way round.
    let d = (long - gmst) % TWO_PI;
    if (d > Math.PI) d -= TWO_PI;
    else if (d < -Math.PI) d += TWO_PI;
    gmst += (1 - w) * d;
  }
  const { dpsi, epsA } = frameOfDate(jdTT);
  const om = delaunay(t)[4]; // the Moon's mean node
  const eect = (0.00264096 * Math.sin(om) + 0.00006352 * Math.sin(2 * om)) * AS2R;
  const gast = gmst + dpsi * Math.cos(epsA) + eect;
  return ((gast % TWO_PI) + TWO_PI) % TWO_PI;
}

const cross = (a: Vec, b: Vec): Vec => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const dot = (a: Vec, b: Vec) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];

/** The local sky of date, as unit vectors in the true equatorial frame. */
interface Sky {
  theta: number; // ARMC, radians
  phi: number;
  eps: number;
  zenith: Vec;
  north: Vec; // the north point of the horizon
  east: Vec; // the east point of the horizon
  pole: Vec; // the ecliptic pole
}

function sky(theta: number, phi: number, eps: number): Sky {
  const cp = Math.cos(phi), sp = Math.sin(phi), ct = Math.cos(theta), st = Math.sin(theta);
  const zenith: Vec = [cp * ct, cp * st, sp];
  const north: Vec = [-sp * ct, -sp * st, cp];
  return { theta, phi, eps, zenith, north, east: cross(north, zenith), pole: [0, -Math.sin(eps), Math.cos(eps)] };
}

/** Ecliptic longitude (degrees) of an equatorial-of-date direction. */
function longitudeOf(v: Vec, eps: number): number {
  const y = v[1] * Math.cos(eps) + v[2] * Math.sin(eps);
  return norm360(Math.atan2(y, v[0]) * R2D);
}

/**
 * Where the ecliptic crosses the great circle with normal `n`, on the side of
 * the sky `side` points into.
 */
function eclipticCrossing(s: Sky, n: Vec, side: Vec): number {
  const d = cross(n, s.pole);
  return longitudeOf(dot(d, side) >= 0 ? d : [-d[0], -d[1], -d[2]], s.eps);
}

/** Ecliptic point with a given right ascension (radians). */
function longitudeAtRA(ra: number, eps: number): number {
  return norm360(Math.atan2(Math.sin(ra), Math.cos(ra) * Math.cos(eps)) * R2D);
}

/** Declination (radians) of an ecliptic longitude (degrees). */
const declinationOf = (lon: number, eps: number) => Math.asin(Math.sin(eps) * Math.sin(lon * D2R));
/** Right ascension (radians) of an ecliptic longitude (degrees). */
const raOf = (lon: number, eps: number) => Math.atan2(Math.sin(lon * D2R) * Math.cos(eps), Math.cos(lon * D2R));

/** Ascendant for a meridian RA `theta` and latitude `phi` (radians): the eastern horizon crossing. */
function ascendantAt(theta: number, phi: number, eps: number): number {
  const s = sky(theta, phi, eps);
  return eclipticCrossing(s, s.zenith, s.east);
}

/** Trisect each quadrant in longitude. Porphyry, and the polar fallback. */
function porphyry(asc: number, mc: number): number[] {
  const c = new Array<number>(12);
  const q1 = norm360(asc - mc); // MC → Asc, through the 11th and 12th
  const q2 = norm360(mc + 180 - asc); // Asc → IC, through the 2nd and 3rd
  c[0] = asc;
  c[9] = mc;
  c[10] = norm360(mc + q1 / 3);
  c[11] = norm360(mc + (2 * q1) / 3);
  c[1] = norm360(asc + q2 / 3);
  c[2] = norm360(asc + (2 * q2) / 3);
  return opposite(c);
}

/** Fill cusps 4–9 from 10–3: each house opposes the house six away. */
function opposite(c: number[]): number[] {
  c[3] = norm360(c[9] + 180);
  c[4] = norm360(c[10] + 180);
  c[5] = norm360(c[11] + 180);
  c[6] = norm360(c[0] + 180);
  c[7] = norm360(c[1] + 180);
  c[8] = norm360(c[2] + 180);
  return c;
}

/**
 * Placidus: a cusp is the ecliptic degree that has covered a third (or two
 * thirds) of its own diurnal or nocturnal semi-arc. Solved by iterating on the
 * degree's declination. Returns null where a degree never rises or sets.
 */
function placidus(s: Sky, asc: number, mc: number): number[] | null {
  const c = new Array<number>(12);
  c[0] = asc;
  c[9] = mc;
  const tp = Math.tan(s.phi);
  // [cusp index, fraction, above the horizon?]
  const specs: Array<[number, number, boolean]> = [[10, 1 / 3, true], [11, 2 / 3, true], [1, 2 / 3, false], [2, 1 / 3, false]];
  for (const [idx, f, above] of specs) {
    let lon = norm360(mc + (above ? 30 : 120) + (idx === 11 || idx === 2 ? 30 : 0));
    for (let i = 0; i < 100; i++) {
      const x = -tp * Math.tan(declinationOf(lon, s.eps));
      if (Math.abs(x) > 1) return null;
      const sa = Math.acos(x); // diurnal semi-arc
      const ra = above ? s.theta + f * sa : s.theta + Math.PI - f * (Math.PI - sa);
      const next = longitudeAtRA(ra, s.eps);
      const step = Math.abs(((next - lon + 540) % 360) - 180);
      lon = next;
      if (step < 1e-11) break;
    }
    c[idx] = lon;
  }
  return opposite(c);
}

/**
 * Koch: the Ascendants at the moments when the Midheaven degree had covered
 * thirds of its diurnal semi-arc. Returns null where that degree never rises.
 */
function koch(s: Sky, asc: number, mc: number): number[] | null {
  const x = Math.tan(s.phi) * Math.tan(declinationOf(mc, s.eps));
  if (Math.abs(x) > 1) return null;
  const dsa = Math.PI / 2 + Math.asin(x); // the MC degree's diurnal semi-arc
  const c = new Array<number>(12);
  c[0] = asc;
  c[9] = mc;
  c[10] = ascendantAt(s.theta - (2 * dsa) / 3, s.phi, s.eps);
  c[11] = ascendantAt(s.theta - dsa / 3, s.phi, s.eps);
  c[1] = ascendantAt(s.theta + dsa / 3, s.phi, s.eps);
  c[2] = ascendantAt(s.theta + (2 * dsa) / 3, s.phi, s.eps);
  return opposite(c);
}

/**
 * Systems whose house circles all pass through the north and south points of
 * the horizon: Campanus divides the prime vertical into equal parts,
 * Regiomontanus the celestial equator. Each cusp is where the ecliptic crosses
 * the circle through the poles of the horizon and a dividing point.
 */
function axial(s: Sky, asc: number, kind: "C" | "R"): number[] {
  const c = new Array<number>(12);
  c[0] = asc;
  // The 10th cusp is where the meridian crosses the ecliptic ABOVE the
  // horizon. Everywhere outside the polar circles that is the Midheaven;
  // inside them the Midheaven can sit below the horizon, and the cusp is the
  // opposite point (as Swiss Ephemeris draws Campanus and Regiomontanus there).
  // Every other cusp is likewise taken on its own circle's side of the sky.
  // Swiss differs inside the polar circles for some Regiomontanus cusps; the
  // geometry here is the definition.
  c[9] = eclipticCrossing(s, s.east, s.zenith);
  const point = (a: number): Vec => {
    if (kind === "C") {
      // a from the zenith through the east point, along the prime vertical
      return [
        Math.cos(a) * s.zenith[0] + Math.sin(a) * s.east[0],
        Math.cos(a) * s.zenith[1] + Math.sin(a) * s.east[1],
        Math.cos(a) * s.zenith[2] + Math.sin(a) * s.east[2],
      ];
    }
    const ra = s.theta + a; // a east of the meridian, along the equator
    return [Math.cos(ra), Math.sin(ra), 0];
  };
  for (const [idx, deg] of [[10, 30], [11, 60], [1, 120], [2, 150]] as const) {
    const p = point(deg * D2R);
    c[idx] = eclipticCrossing(s, cross(s.north, p), p);
  }
  return opposite(c);
}

/**
 * The house cusps, angles and sidereal time for a moment and place.
 *
 * `jdUT` is the moment; `deltaTSeconds` converts it to TT for precession and
 * nutation. Latitude and longitude are geographic degrees, east positive.
 */
export function houses(jdUT: number, deltaTSeconds: number, lat: number, lon: number, system: HouseSystem): Houses {
  const jdTT = jdUT + deltaTSeconds / 86400;
  const eps = frameOfDate(jdTT).eps;
  const theta = apparentSiderealTime(jdUT, jdTT) + lon * D2R;
  const phi = lat * D2R;
  const s = sky(theta, phi, eps);

  const mc = longitudeAtRA(theta, eps);
  const asc = eclipticCrossing(s, s.zenith, s.east);
  // The Vertex: where the ecliptic crosses the prime vertical, in the west.
  const vertex = eclipticCrossing(s, s.north, [-s.east[0], -s.east[1], -s.east[2]]);
  const equatorialAscendant = ascendantAt(theta, 0, eps);

  let cusps: number[] | null;
  let used: HouseSystem = system;
  switch (system) {
    case "W": {
      const first = Math.floor(asc / 30) * 30;
      cusps = Array.from({ length: 12 }, (_, i) => norm360(first + 30 * i));
      break;
    }
    case "A":
    case "E":
      cusps = Array.from({ length: 12 }, (_, i) => norm360(asc + 30 * i));
      break;
    case "O":
      cusps = porphyry(asc, mc);
      break;
    case "C":
    case "R":
      cusps = axial(s, asc, system);
      break;
    case "X":
      // Meridian (axial rotation): equal divisions of the equator from the
      // meridian, carried to the ecliptic along hour circles.
      cusps = Array.from({ length: 12 }, (_, i) => longitudeAtRA(theta + ((i + 3) % 12) * 30 * D2R, eps));
      break;
    case "M":
      // Morinus: equal divisions of the equator, carried to the ecliptic along
      // circles through the ecliptic poles.
      cusps = Array.from({ length: 12 }, (_, i) => {
        const ra = theta + ((i + 3) % 12) * 30 * D2R;
        return norm360(Math.atan2(Math.sin(ra) * Math.cos(eps), Math.cos(ra)) * R2D);
      });
      break;
    // Inside the polar circles some ecliptic degrees never rise or set, and
    // Placidus and Koch are undefined there. The rule is the latitude, not
    // whether this moment's cusps happen to compute: a chart must not change
    // system as the sky turns. Porphyry stands in, as in Swiss Ephemeris.
    case "P":
      cusps = Math.abs(phi) > Math.PI / 2 - eps ? null : placidus(s, asc, mc);
      break;
    case "K":
      cusps = Math.abs(phi) > Math.PI / 2 - eps ? null : koch(s, asc, mc);
      break;
    default:
      throw new Error(`unknown house system ${system as string}`);
  }
  const fallback = cusps === null;
  if (fallback) {
    cusps = porphyry(asc, mc);
    used = "O";
  }
  return { cusps: cusps as number[], ascendant: asc, mc, armc: norm360(theta * R2D), vertex, equatorialAscendant, system: used, fallback };
}

/** Right ascension of an ecliptic longitude, degrees — exposed for callers that need it. */
export const rightAscensionOf = (lon: number, eps: number) => norm360(raOf(lon, eps) * R2D);
