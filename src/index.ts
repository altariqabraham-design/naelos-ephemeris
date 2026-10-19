// Naelos Ephemeris — planetary positions from NASA JPL's ephemerides: DE440,
// the most accurate for modern times, across 1550–2649, and DE441, JPL's
// long-span companion to it, before and after — 13,200 BC to AD 17,191 in all.
// Chiron and the asteroids 1600–2497, as far as JPL publishes their orbits.
//
// The interface deliberately answers the way Swiss Ephemeris's calc_ut does
// (same body numbers, same flags, same result fields), so a chart engine
// written against Swiss can switch by changing what it constructs.

import { parseTable, type TableFile } from "./tables.ts";
import { frameOfDate } from "./frames.ts";
import { julday, tdbMinusTT } from "./time.ts";
import { Tables, geocentric, heliocentric, spherical, EARTH } from "./apparent.ts";
import { deltaT as builtInDeltaT } from "./deltat.ts";
import { norm360 } from "./angles.ts";
import { houses, type Houses, type HouseSystem } from "./houses.ts";
import { meanNode, trueNode, meanApogee, trueApogee } from "./lunar-points.ts";
import { solarEclipse, lunarEclipse, type Eclipse, type SkySource } from "./eclipses.ts";
import { stations as findStations, crossings as findCrossings, combined, MOTION, type AngleFn, type Motion } from "./events.ts";
import { chartUncertainty, type UncertaintyInput, type Uncertainty } from "./uncertainty.ts";
import { rectify, type RectifyInput, type Rectification } from "./rectification.ts";
import { provenanceOf, type Provenance } from "./provenance.ts";
import { FILES } from "./generated/files.ts";

export { horizonsUrl, ENGINE, type Provenance } from "./provenance.ts";

export { CERTAIN, type UncertaintyInput, type Uncertainty, type Placement, type Share, type AspectShare } from "./uncertainty.ts";
export { type RectifyInput, type Rectification, type LifeEvent, type Candidate, type Evidence, type Technique } from "./rectification.ts";

/** A station: the moment a body stands still before turning. Times in UT and TT. */
export interface StationEvent { jd: number; jdTT: number; longitude: number; turning: "retrograde" | "direct" }
/** The moment a longitude reaches a degree; `direction` −1 means while retrograde. */
export interface CrossingEvent { jd: number; jdTT: number; target: number; direction: 1 | -1 }
export const LUNATION_PHASES = ["new", "first quarter", "full", "last quarter"] as const;
export type LunationPhase = (typeof LUNATION_PHASES)[number];

export { categoryOf, type Eclipse, type EclipseCategory } from "./eclipses.ts";
export {
  SE_ECL_CENTRAL, SE_ECL_NONCENTRAL, SE_ECL_TOTAL, SE_ECL_ANNULAR, SE_ECL_PARTIAL, SE_ECL_ANNULAR_TOTAL, SE_ECL_PENUMBRAL,
} from "./eclipses.ts";

export { julday, juldayJulian, juldayHistorical, revjul, formatYear, GREGORIAN_REFORM_JD, type CalendarDate } from "./time.ts";
export { FILES as TABLE_FILES, NEAR_SOURCE, FAR_SOURCE } from "./generated/files.ts";
export { norm360 } from "./angles.ts";
export { parseTable } from "./tables.ts";
export { deltaT, deltaTUncertainty, MEASURED_UNTIL, PREDICTED_UNTIL } from "./deltat.ts";
export { apparentSiderealTime, earthRotationAngle, type Houses, type HouseSystem } from "./houses.ts";

export const SE_SUN = 0, SE_MOON = 1, SE_MERCURY = 2, SE_VENUS = 3, SE_MARS = 4, SE_JUPITER = 5,
  SE_SATURN = 6, SE_URANUS = 7, SE_NEPTUNE = 8, SE_PLUTO = 9, SE_EARTH = EARTH;
/** The Moon's nodes and apogee: mean and true North Node, mean and true Lilith. */
export const SE_MEAN_NODE = 10, SE_TRUE_NODE = 11, SE_MEAN_APOG = 12, SE_OSCU_APOG = 13;
/** Chiron and the four largest main-belt asteroids, from JPL Horizons. */
export const SE_CHIRON = 15, SE_CERES = 17, SE_PALLAS = 18, SE_JUNO = 19, SE_VESTA = 20;
/**
 * Not a body: the obliquity and nutation of date, laid out as Swiss Ephemeris
 * lays them out — `longitude` the true obliquity, `latitude` the mean
 * obliquity, `distance` the nutation in longitude and `longitudeSpeed` the
 * nutation in obliquity, all in degrees. With the true obliquity an ecliptic
 * place converts to right ascension and declination exactly.
 */
export const SE_ECL_NUT = -1;
const LUNAR_POINTS = new Set([10, 11, 12, 13]);

/** Where the tables begin and end (TDB Julian Days): 13,200 BC to AD 17,191. */
export const FIRST_JD = FILES[0][1];
export const LAST_JD = FILES[FILES.length - 1][2];

/** Whether a file is one of the far tables (DE441), held outside the repository. */
export function isFarTable(name: string): boolean {
  return FILES.some(([n, , , far]) => n === name && far === 1);
}

/**
 * The century files needed for every moment from `from` to `to` (Julian Days,
 * UT or TT), from each file's exact span (src/generated/files.ts). A margin
 * either side covers light-time and speeds — the first hours of a century can
 * need the century before, because light left the body before it began — and
 * ΔT, which grows to 8.5 days by 13,200 BC: a UT moment there is TT a week
 * later, so the margin grows with ΔT's size. Moments outside the tables get
 * the nearest file, which then refuses them.
 */
export function tablesFor(from: number, to: number): string[] {
  const margin = (jd: number) => 1 + Math.abs(builtInDeltaT(jd)) / 86400;
  const lo = Math.min(from, to), hi = Math.max(from, to);
  const a = lo - margin(lo), b = hi + margin(hi);
  const names = FILES.filter(([, start, end], i) =>
    (i === 0 ? -Infinity : start) <= b && (i === FILES.length - 1 ? Infinity : end) > a).map(([n]) => n);
  return names;
}

export const SEFLG_SWIEPH = 2;
export const SEFLG_HELCTR = 8;
export const SEFLG_SPEED = 256;
export const SEFLG_EQUATORIAL = 2048;

export interface Position {
  /** Ecliptic longitude, or right ascension with SEFLG_EQUATORIAL (degrees). */
  longitude: number;
  /** Ecliptic latitude, or declination with SEFLG_EQUATORIAL (degrees). */
  latitude: number;
  /** AU. */
  distance: number;
  /** Per day; zero unless SEFLG_SPEED was asked for. */
  longitudeSpeed: number;
  latitudeSpeed: number;
  distanceSpeed: number;
  flags: number;
}

export interface PositionError {
  error: string;
}

export interface Options {
  /**
   * ΔT = TT − UT in seconds, for a UT Julian Day. Defaults to the built-in
   * model: IERS measurements, the 2021 historical splines and a published
   * long-term curve (see deltat.ts).
   */
  deltaT?: (jdUT: number) => number;
}

const BODIES = new Set([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, EARTH, 15, 17, 18, 19, 20]);
/** Half-interval for speeds by central difference, in days. */
const H = 1 / 1440;

export class NaelosEphemeris {
  private readonly tables: Tables;
  private readonly deltaTFn: (jdUT: number) => number;

  constructor(files: TableFile[], options: Options = {}) {
    this.tables = new Tables(files);
    this.deltaTFn = options.deltaT ?? builtInDeltaT;
  }

  /** Build from the raw century files (any subset, in any order). */
  static fromBuffers(buffers: Array<ArrayBuffer | Uint8Array>, options: Options = {}): NaelosEphemeris {
    return new NaelosEphemeris(buffers.map(parseTable), options);
  }

  /**
   * Load more century files, as they are needed: a browser fetches only the
   * centuries a chart falls in (see `tablesFor`). Files already loaded are
   * ignored.
   */
  addTables(buffers: Array<ArrayBuffer | Uint8Array>): void {
    this.tables.add(buffers.map(parseTable));
  }

  /** The TDB span from the first loaded file to the last. */
  get range(): { start: number; end: number } {
    return { start: this.tables.start, end: this.tables.end };
  }

  julday(year: number, month: number, day: number, hour = 0): number {
    return julday(year, month, day, hour);
  }

  /** ΔT in days, as Swiss Ephemeris's swe_deltat returns it. */
  deltat(jdUT: number): number {
    return this.deltaTFn(jdUT) / 86400;
  }

  /** Position at a UT Julian Day. */
  calc_ut(jdUT: number, body: number, flags: number): Position | PositionError {
    return this.calc(jdUT + this.deltaTFn(jdUT) / 86400, body, flags);
  }

  /**
   * House cusps and angles for a UT moment and a place (degrees, east
   * positive). Placidus and Koch fall back to Porphyry where they cannot
   * exist; `fallback` says when that happened.
   */
  houses(jdUT: number, lat: number, lon: number, system: HouseSystem): Houses {
    return houses(jdUT, this.deltaTFn(jdUT), lat, lon, system);
  }

  /** Nearest global solar eclipse after (or, `backward`, before) a UT Julian Day. */
  solarEclipse(jdUT: number, backward = false): Eclipse | null {
    return solarEclipse(this.sky, jdUT, backward);
  }

  /** Nearest lunar eclipse after (or, `backward`, before) a UT Julian Day. */
  lunarEclipse(jdUT: number, backward = false): Eclipse | null {
    return lunarEclipse(this.sky, jdUT, backward);
  }

  // ── exact moments ─────────────────────────────────────────────────────────

  private longitudeOf(body: number): AngleFn {
    return (tt) => {
      const r = this.calc(tt, body, 0);
      if ("error" in r) throw new Error(r.error);
      return r.longitude;
    };
  }

  private motionOf(body: number): Motion {
    const m = MOTION[body];
    if (!m) throw new Error(`no motion profile for body ${body}`);
    return m;
  }

  private toTT(jdUT: number) {
    return jdUT + this.deltaTFn(jdUT) / 86400;
  }

  private toUT(jdTT: number) {
    return jdTT - this.deltaTFn(jdTT - this.deltaTFn(jdTT) / 86400) / 86400;
  }

  /** Every station of a body between two UT moments: when it turns retrograde or direct. */
  stations(body: number, startUT: number, endUT: number): StationEvent[] {
    return findStations(this.longitudeOf(body), this.toTT(startUT), this.toTT(endUT), this.motionOf(body))
      .map((s) => ({ jd: this.toUT(s.jdTT), jdTT: s.jdTT, longitude: norm360(s.longitude), turning: s.turning }));
  }

  /**
   * Every moment a body's longitude reaches any of `targets` (degrees) between
   * two UT moments — all passes, including the three of a retrograde loop.
   * A transit to a natal point at an aspect is a target of natal + angle.
   */
  crossings(body: number, targets: number[], startUT: number, endUT: number): CrossingEvent[] {
    return findCrossings(this.longitudeOf(body), targets, this.toTT(startUT), this.toTT(endUT), this.motionOf(body))
      .map((c) => ({ jd: this.toUT(c.jdTT), jdTT: c.jdTT, target: c.target, direction: c.direction }));
  }

  /** Sign ingresses (crossings of each 30° boundary), with the sign entered (0 = Aries). */
  ingresses(body: number, startUT: number, endUT: number): Array<CrossingEvent & { sign: number }> {
    return this.crossings(body, Array.from({ length: 12 }, (_, i) => i * 30), startUT, endUT).map((c) => ({
      ...c,
      // Moving forward the body enters the sign that starts at the boundary;
      // retrograde it re-enters the one before.
      sign: c.direction === 1 ? Math.round(c.target / 30) % 12 : (Math.round(c.target / 30) + 11) % 12,
    }));
  }

  /** Returns: every moment a body comes back to a longitude (a solar return is the Sun's). */
  returns(body: number, longitude: number, startUT: number, endUT: number): CrossingEvent[] {
    return this.crossings(body, [longitude], startUT, endUT);
  }

  /** New, first-quarter, full and last-quarter Moons between two UT moments. */
  lunations(startUT: number, endUT: number): Array<{ jd: number; jdTT: number; phase: LunationPhase }> {
    const moon = this.longitudeOf(SE_MOON), sun = this.longitudeOf(SE_SUN);
    const elongation: AngleFn = (tt) => norm360(moon(tt) - sun(tt));
    return findCrossings(elongation, [0, 90, 180, 270], this.toTT(startUT), this.toTT(endUT), { step: 2, maxSpeed: 16, monotone: true })
      .map((c) => ({ jd: this.toUT(c.jdTT), jdTT: c.jdTT, phase: LUNATION_PHASES[Math.round(c.target / 90) % 4] }));
  }

  /**
   * Exact aspects between two moving bodies (transit to transit): every moment
   * their separation reaches an angle, either way round, between two UT moments.
   */
  aspectsBetween(bodyA: number, bodyB: number, angles: number[], startUT: number, endUT: number): Array<CrossingEvent & { angle: number }> {
    const a = this.longitudeOf(bodyA), b = this.longitudeOf(bodyB);
    const separation: AngleFn = (tt) => norm360(a(tt) - b(tt));
    const targets = [...new Set(angles.flatMap((x) => [norm360(x), norm360(-x)]))];
    return findCrossings(separation, targets, this.toTT(startUT), this.toTT(endUT), combined(this.motionOf(bodyA), this.motionOf(bodyB)))
      .map((c) => ({ jd: this.toUT(c.jdTT), jdTT: c.jdTT, target: c.target, direction: c.direction, angle: c.target > 180 ? 360 - c.target : c.target }));
  }

  // ── showing the working ───────────────────────────────────────────────────

  /**
   * Where a position came from: the JPL data and table segment, every
   * correction in order, how certain the clock was, and a JPL Horizons link
   * that computes the same quantity independently.
   */
  provenance(jdUT: number, body: number, flags = 0): Provenance | PositionError {
    const r = this.calc_ut(jdUT, body, flags);
    if ("error" in r) return r;
    const tdb = this.toTT(jdUT) + tdbMinusTT(this.toTT(jdUT));
    // The Earth rests on the Earth-Moon barycentre's series; the Moon and her
    // nodes and apogee on the Moon's.
    const seg = this.tables.locate(body === EARTH ? 200 : body >= 10 && body <= 13 ? SE_MOON : body, tdb);
    return provenanceOf({
      body, jdUT, deltaT: this.deltaTFn(jdUT), longitude: r.longitude, latitude: r.latitude, distanceAU: r.distance,
      heliocentric: (flags & SEFLG_HELCTR) !== 0, segment: seg ?? undefined,
    });
  }

  // ── knowing how sure ──────────────────────────────────────────────────────

  /**
   * What a birth window allows: every body's and angle's possible signs and
   * houses with the share of the window each holds, which aspects hold and for
   * how much of it, and how uncertain the clock itself is.
   */
  chartUncertainty(input: UncertaintyInput): Uncertainty {
    return chartUncertainty({
      longitude: (b, t) => this.longitudeOf(b)(this.toTT(t)),
      houses: (t, lat, lon, sys) => this.houses(t, lat, lon, sys),
      deltaT: (t) => this.deltaTFn(t),
    }, input);
  }

  /**
   * Rectification: support for each candidate birth minute from dated life
   * events, with the contacts behind it (transits, solar arcs and the
   * progressed Moon to and from the angles).
   */
  rectify(input: RectifyInput): Rectification {
    return rectify({
      longitude: (b, t) => this.longitudeOf(b)(this.toTT(t)),
      angles: (t, lat, lon, sys) => this.houses(t, lat, lon, sys),
    }, input);
  }

  /** The Sun and Moon as the eclipse search needs them. */
  private readonly sky: SkySource = {
    deltaT: (jdUT) => this.deltaTFn(jdUT),
    sunMoon: (jdTT) => {
      const frame = frameOfDate(jdTT);
      const tdb = jdTT + tdbMinusTT(jdTT);
      const s = geocentric(this.tables, SE_SUN, tdb, frame);
      const m = geocentric(this.tables, SE_MOON, tdb, frame);
      const km = 149597870.7;
      const [sunLon] = spherical(s, frame, false);
      const [moonLon, moonLat] = spherical(m, frame, false);
      return {
        sun: [s.u[0] * s.dist * km, s.u[1] * s.dist * km, s.u[2] * s.dist * km],
        moon: [m.u[0] * m.dist * km, m.u[1] * m.dist * km, m.u[2] * m.dist * km],
        sunLon, moonLon, moonLat,
      };
    },
  };

  /** Position at a TT Julian Day. */
  calc(jdTT: number, body: number, flags: number): Position | PositionError {
    if (body === SE_ECL_NUT) {
      const f = frameOfDate(jdTT);
      const deg = 180 / Math.PI;
      return {
        longitude: f.eps * deg, latitude: f.epsA * deg, distance: f.dpsi * deg,
        longitudeSpeed: f.deps * deg, latitudeSpeed: 0, distanceSpeed: 0, flags: 0,
      };
    }
    if (!BODIES.has(body)) return { error: `body ${body} is not available yet` };
    const helio = (flags & SEFLG_HELCTR) !== 0;
    if (helio && body === SE_SUN) return { error: "the Sun has no heliocentric position" };
    if (!helio && body === EARTH) return { error: "the Earth has no geocentric position" };
    if (helio && LUNAR_POINTS.has(body)) return { error: "the lunar nodes and apogee are geocentric points" };
    const equatorial = (flags & SEFLG_EQUATORIAL) !== 0;
    try {
      const at = (tt: number) => {
        const frame = frameOfDate(tt);
        const tdb = tt + tdbMinusTT(tt);
        const ap = body === SE_MEAN_NODE ? meanNode(tt, frame)
          : body === SE_TRUE_NODE ? trueNode(this.tables, tdb, frame)
          : body === SE_MEAN_APOG ? meanApogee(tt, frame)
          : body === SE_OSCU_APOG ? trueApogee(this.tables, tdb, frame)
          : helio ? heliocentric(this.tables, body, tdb, frame) : geocentric(this.tables, body, tdb, frame);
        const [lon, lat] = spherical(ap, frame, equatorial);
        return { lon, lat, dist: ap.dist };
      };
      const p = at(jdTT);
      let lonSpeed = 0, latSpeed = 0, distSpeed = 0;
      if (flags & SEFLG_SPEED) {
        const a = at(jdTT - H), b = at(jdTT + H);
        let dl = b.lon - a.lon;
        if (dl > 180) dl -= 360;
        else if (dl < -180) dl += 360;
        lonSpeed = dl / (2 * H);
        latSpeed = (b.lat - a.lat) / (2 * H);
        distSpeed = (b.dist - a.dist) / (2 * H);
      }
      return {
        longitude: p.lon, latitude: p.lat, distance: p.dist,
        longitudeSpeed: lonSpeed, latitudeSpeed: latSpeed, distanceSpeed: distSpeed,
        flags: flags & (SEFLG_HELCTR | SEFLG_SPEED | SEFLG_EQUATORIAL),
      };
    } catch (e) {
      return { error: e instanceof Error ? e.message : String(e) };
    }
  }
}
