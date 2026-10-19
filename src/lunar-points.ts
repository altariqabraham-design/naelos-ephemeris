// lunar-points.ts — the Moon's nodes and apogee: the North Node and Lilith.
//
// MEAN points come from the mean elements of the lunar theory (ELP2000-82B as
// published by Meeus, chapter 47): smooth, slowly turning, the Node always
// retrograde. Those are polynomials fitted to the present, and they drift
// away from it — 2.4° from the true node's average by 13,000 BC — so beyond
// the DE440 tables they take a correction measured on DE441 itself (the true
// points averaged over a century, minus the polynomial; see
// tools/build_lunar_mean.mts), blended in exactly as the long-term frame is.
// From 1550 to 2650 they are the polynomials alone, as they always were. TRUE points come from the Moon's OSCULATING orbit — the ellipse
// it would follow from this instant if the Sun and planets let go — built from
// the engine's own lunar position and velocity:
//
//   angular momentum h = r × v           the orbit's plane
//   eccentricity e = (v × h)/μ − r/|r|   points to perigee
//   True Node = where that plane crosses the ecliptic of date, northbound
//   True Lilith = the apogee, −e, as a point in space (it has a latitude)
//
// Both are computed in the inertial frame and only then turned into the frame
// of date: rotating r and v first would add the frame's own slow turning to
// the velocity. Positions are geometric (no light-time or aberration) — these
// are points of an orbit, not bodies anyone sees — on the true ecliptic and
// equinox of date, as Swiss Ephemeris gives them.

import { apply, iauWeight, type Frame } from "./frames.ts";
import {
  LUNAR_MEAN_FROM_EPOCH, LUNAR_MEAN_STEP_YEARS, MEAN_NODE_CORRECTION, MEAN_APOGEE_CORRECTION,
} from "./generated/lunar-mean.ts";
import type { Tables } from "./apparent.ts";

/** G(M_Earth + M_Moon) from DE440, km³/day². */
const GM = (398600.435507 + 4902.800118) * 86400 * 86400;
const AU_KM = 149597870.7;
const D2R = Math.PI / 180;
/** Mean distance and eccentricity of the lunar orbit, and its mean inclination. */
const MEAN_DISTANCE_KM = 384400;
const MEAN_ECCENTRICITY = 0.054900489;
const MEAN_INCLINATION = 5.1453964 * D2R;

type Vec = [number, number, number];
const cross = (a: number[], b: number[]): Vec => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const dot = (a: number[], b: number[]) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const norm = (a: number[]) => Math.sqrt(dot(a, a));

/** A point given as a unit vector on the true equator of date, and a distance in AU. */
export interface LunarPoint {
  u: Vec;
  dist: number;
}

/** The ecliptic pole of date, in true-equator-of-date coordinates. */
const eclipticPole = (eps: number): Vec => [0, -Math.sin(eps), Math.cos(eps)];

/** An ecliptic-of-date longitude and latitude (radians) as a true-equator-of-date unit vector. */
function fromEcliptic(lon: number, lat: number, eps: number): Vec {
  const x = Math.cos(lat) * Math.cos(lon), y = Math.cos(lat) * Math.sin(lon), z = Math.sin(lat);
  const c = Math.cos(eps), s = Math.sin(eps);
  return [x, c * y - s * z, s * y + c * z];
}

function osculating(tables: Tables, tdb: number, frame: Frame) {
  const r: Vec = [0, 0, 0], v: Vec = [0, 0, 0];
  tables.moonGeocentric(tdb, r, v);
  const h = cross(r, v);
  const rv = cross(v, h);
  const rl = norm(r);
  const e: Vec = [rv[0] / GM - r[0] / rl, rv[1] / GM - r[1] / rl, rv[2] / GM - r[2] / rl];
  const p = dot(h, h) / GM; // semi-latus rectum, km
  return { hDate: apply(frame.npb, h), eDate: apply(frame.npb, e), p, ecc: norm(e) };
}

/** True (osculating) North Node. Its distance is where the osculating ellipse crosses the node. */
export function trueNode(tables: Tables, tdb: number, frame: Frame): LunarPoint {
  const { hDate, eDate, p } = osculating(tables, tdb, frame);
  const n = cross(eclipticPole(frame.eps), hDate);
  const nl = norm(n);
  const u: Vec = [n[0] / nl, n[1] / nl, n[2] / nl];
  return { u, dist: p / (1 + dot(eDate, u)) / AU_KM };
}

/** True Lilith: the apogee of the osculating orbit, with its latitude. */
export function trueApogee(tables: Tables, tdb: number, frame: Frame): LunarPoint {
  const { eDate, p, ecc } = osculating(tables, tdb, frame);
  const u: Vec = [-eDate[0] / ecc, -eDate[1] / ecc, -eDate[2] / ecc];
  return { u, dist: p / (1 - ecc) / AU_KM }; // a(1 + e) = p / (1 − e)
}

/** A correction table at a TT moment, by cubic Hermite, degrees. */
function correction(table: readonly number[], jdTT: number): number {
  const last = table.length - 1;
  const x = (2000 + (jdTT - 2451545.0) / 365.25 - LUNAR_MEAN_FROM_EPOCH) / LUNAR_MEAN_STEP_YEARS;
  let i = Math.floor(x);
  if (i < 1) i = 1;
  else if (i > last - 2) i = last - 2;
  const u = Math.min(Math.max(x - i, -1), 2); // held at the ends beyond the table
  const p0 = table[i - 1], p1 = table[i], p2 = table[i + 1], p3 = table[i + 2];
  const m1 = (p2 - p0) / 2, m2 = (p3 - p1) / 2;
  const u2 = u * u, u3 = u2 * u;
  return (2 * u3 - 3 * u2 + 1) * p1 + (u3 - 2 * u2 + u) * m1 + (-2 * u3 + 3 * u2) * p2 + (u3 - u2) * m2;
}

/**
 * Mean elements of the lunar orbit (Meeus, Astronomical Algorithms, ch. 47),
 * degrees, with DE441's correction beyond the DE440 tables.
 */
function meanElements(jdTT: number) {
  const T = (jdTT - 2451545.0) / 36525;
  let node = 125.0445479 - 1934.1362891 * T + 0.0020754 * T * T + (T * T * T) / 467441 - (T * T * T * T) / 60616000;
  let perigee = 83.3532465 + 4069.0137287 * T - 0.01032 * T * T - (T * T * T) / 80053 + (T * T * T * T) / 18999000;
  const c = meanCorrection(jdTT);
  node += c.node;
  perigee += c.apogee; // the apogee is the perigee turned half round
  return { node, perigee };
}

/** The correction applied to the mean node and mean apogee at a TT moment, degrees (zero 1550–2650). */
export function meanCorrection(jdTT: number): { node: number; apogee: number } {
  const far = 1 - iauWeight(jdTT);
  if (far === 0) return { node: 0, apogee: 0 };
  return { node: far * correction(MEAN_NODE_CORRECTION, jdTT), apogee: far * correction(MEAN_APOGEE_CORRECTION, jdTT) };
}

/** Mean North Node, on the true equinox of date (mean node plus nutation in longitude). */
export function meanNode(jdTT: number, frame: Frame): LunarPoint {
  const lon = meanElements(jdTT).node * D2R + frame.dpsi;
  return { u: fromEcliptic(lon, 0, frame.eps), dist: MEAN_DISTANCE_KM / AU_KM };
}

/**
 * Mean Lilith: the apogee of the mean orbit, as a point on that orbit — so it
 * has a latitude, up to the orbit's 5.1° inclination — projected onto the
 * ecliptic of date.
 */
export function meanApogee(jdTT: number, frame: Frame): LunarPoint {
  const { node, perigee } = meanElements(jdTT);
  const om = node * D2R;
  const u = (perigee + 180 - node) * D2R; // the apogee's argument of latitude
  const lon = om + Math.atan2(Math.sin(u) * Math.cos(MEAN_INCLINATION), Math.cos(u)) + frame.dpsi;
  const lat = Math.asin(Math.sin(u) * Math.sin(MEAN_INCLINATION));
  return { u: fromEcliptic(lon, lat, frame.eps), dist: (MEAN_DISTANCE_KM * (1 + MEAN_ECCENTRICITY)) / AU_KM };
}
