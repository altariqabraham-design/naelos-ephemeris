// apparent.ts — where a body appears, from where it is.
//
// The tables give where bodies ARE (DE440: barycentric, ICRF, TDB). An
// astrologer's chart gives where they APPEAR from the centre of the Earth, on
// the true ecliptic of date. In between, in order:
//   1. light-time  — the body is seen where it was when its light left;
//   2. deflection  — the Sun's gravity bends that light (relativistic, PPN);
//   3. aberration  — the Earth's own motion tilts it (relativistic);
//   4. precession and nutation — ICRF to the true equator of date;
//   5. the true obliquity — equator of date to ecliptic of date.
// These are the steps of IERS Conventions / ERFA / NOVAS, in the same order.

import { covers, evaluate, type Series, type TableFile } from "./tables.ts";
import { apply, toEcliptic, type Frame } from "./frames.ts";
import { FILES } from "./generated/files.ts";
import { norm360 } from "./angles.ts";

export const C_KM_PER_DAY = 299792.458 * 86400;
export const AU_KM = 149597870.7;
/** DE440's Earth/Moon mass ratio. */
export const EMRAT = 81.3005682214972;
/** 2GM☉/c² in AU (the Sun's Schwarzschild radius). */
const SRS_AU = 1.97412574336e-8;
/** The Sun's radius (IAU 2015 nominal), km. */
const R_SUN_KM = 695700;

/**
 * Fraction of the Sun's mass within a radius (in solar radii), from a standard
 * solar model, coarsely: about half within 0.25 R☉, nine-tenths within 0.5.
 * Used only to taper light deflection for bodies hidden behind the Sun.
 */
const MASS_R = [0, 0.1, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1];
const MASS_F = [0, 0.07, 0.34, 0.5, 0.62, 0.8, 0.91, 0.96, 0.985, 0.996, 1];
function enclosedSolarMass(x: number): number {
  if (x >= 1) return 1;
  let i = 1;
  while (MASS_R[i] < x) i++;
  const u = (x - MASS_R[i - 1]) / (MASS_R[i] - MASS_R[i - 1]);
  return MASS_F[i - 1] + u * (MASS_F[i] - MASS_F[i - 1]);
}

// Series ids inside the tables (Swiss numbering where a body has one).
const SUN = 0, MOON = 1, EMB = 200;
/** Body ids the engine answers for directly from the tables. */
export const EARTH = 14;
/** A series stored relative to the Sun (frame 2): the asteroids. */
const HELIOCENTRIC = 2;

type Vec = [number, number, number];

export class Tables {
  private files: TableFile[] = [];
  private a: Vec = [0, 0, 0];
  private b: Vec = [0, 0, 0];
  private av: Vec = [0, 0, 0];
  private bv: Vec = [0, 0, 0];

  constructor(files: TableFile[]) {
    if (!files.length) throw new Error("no tables loaded");
    this.add(files);
  }

  /** Load more century files. One already loaded (the same start) is ignored. */
  add(files: TableFile[]): void {
    for (const f of files) if (!this.files.some((g) => g.start === f.start)) this.files.push(f);
    this.files.sort((x, y) => x.start - y.start);
  }

  /** The TDB span from the first loaded file to the last (gaps between them are refused). */
  get start(): number {
    return this.files[0].start;
  }
  get end(): number {
    return this.files[this.files.length - 1].end;
  }

  /**
   * The file and series that answer for a body at a TDB moment: the file for
   * that moment's century, or — at the edge of what is loaded — any loaded
   * file whose copy of the series reaches it. Each file carries every segment
   * that overlaps its century, so its first and last reach a little past its
   * ends, which is what light-time needs in the first hours of a century.
   */
  private find(id: number, t: number): { file: TableFile; series: Series } {
    const home = this.files.find((f) => t >= f.start && (t < f.end || t === this.end));
    const own = home?.series.get(id);
    if (home && own && covers(own, t)) return { file: home, series: own };
    for (const file of this.files) {
      const series = file.series.get(id);
      if (series && covers(series, t)) return { file, series };
    }
    if (!home) throw new RangeError(`outside the loaded tables (TDB JD ${t})`);
    if (!own) throw new RangeError(`no series ${id} in the tables for TDB JD ${t}`);
    throw new RangeError(`series ${id} covers TDB JD ${own.t0} to ${own.t0 + own.segments * own.segmentDays}, not ${t}`);
  }

  private raw(id: number, t: number, pos: number[], vel?: number[]): number {
    const { series } = this.find(id, t);
    evaluate(series, t, pos, vel);
    return series.frame;
  }

  /** Which century file and segment serve a series at a TDB moment (for provenance). */
  locate(id: number, t: number): { file: string; index: number; from: number; to: number } | null {
    try {
      const { file: f, series: s } = this.find(id, t);
      const index = Math.min(s.segments - 1, Math.max(0, Math.floor((t - s.t0) / s.segmentDays)));
      // The file's name, from the generated index of every file's span.
      const entry = FILES.find(([, start]) => start === f.start);
      const file = entry ? entry[0] : `a table starting at TDB JD ${f.start}`;
      return { file, index, from: s.t0 + index * s.segmentDays, to: s.t0 + (index + 1) * s.segmentDays };
    } catch {
      return null;
    }
  }

  /** The Moon relative to the Earth: position (km) and velocity (km/day), geometric. */
  moonGeocentric(t: number, pos: number[], vel: number[]): void {
    this.raw(MOON, t, pos, vel);
  }

  /** Barycentric position (km) of a body; velocity (km/day) too if asked. */
  bary(body: number, t: number, pos: number[], vel?: number[]): void {
    if (body !== MOON && body !== EARTH) {
      // The asteroids are stored around the Sun, as Horizons gives them.
      if (this.raw(body, t, pos, vel) === HELIOCENTRIC) {
        this.raw(SUN, t, this.a, vel ? this.av : undefined);
        for (let i = 0; i < 3; i++) {
          pos[i] += this.a[i];
          if (vel) vel[i] += this.av[i];
        }
      }
      return;
    }
    this.raw(EMB, t, this.a, vel ? this.av : undefined);
    this.raw(MOON, t, this.b, vel ? this.bv : undefined);
    // Earth and Moon on either side of their barycentre, by mass.
    const f = body === MOON ? EMRAT / (1 + EMRAT) : -1 / (1 + EMRAT);
    for (let i = 0; i < 3; i++) {
      pos[i] = this.a[i] + f * this.b[i];
      if (vel) vel[i] = this.av[i] + f * this.bv[i];
    }
  }
}

const dot = (a: number[], b: number[]) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const len = (a: number[]) => Math.sqrt(dot(a, a));

export interface Apparent {
  /** Unit vector, true equator and equinox of date. */
  u: Vec;
  /** Distance in AU at the moment the light left. */
  dist: number;
}

/**
 * Light-time: the body's position when the light now arriving at `from` left
 * it. Returns the vector from observer to body, in km.
 */
function retarded(tables: Tables, body: number, t: number, from: number[], out: number[]): number {
  let tau = 0;
  for (let i = 0; i < 4; i++) {
    tables.bary(body, t - tau, out);
    const next = Math.hypot(out[0] - from[0], out[1] - from[1], out[2] - from[2]) / C_KM_PER_DAY;
    if (Math.abs(next - tau) < 1e-12) {
      tau = next;
      break;
    }
    tau = next;
  }
  tables.bary(body, t - tau, out);
  out[0] -= from[0];
  out[1] -= from[1];
  out[2] -= from[2];
  return tau;
}

/** Geocentric apparent place of a body at a TDB Julian Day. */
export function geocentric(tables: Tables, body: number, tdb: number, frame: Frame): Apparent {
  const E: Vec = [0, 0, 0], V: Vec = [0, 0, 0], S: Vec = [0, 0, 0], p: Vec = [0, 0, 0];
  tables.bary(EARTH, tdb, E, V);
  tables.bary(SUN, tdb, S);
  retarded(tables, body, tdb, E, p);
  const dist = len(p);
  let u: Vec = [p[0] / dist, p[1] / dist, p[2] / dist];

  // Deflection by the Sun. Not for the Sun itself, and not along a line that
  // runs straight through it, where the formula has no meaning.
  //
  // BEHIND THE SOLAR DISC the formula's deflection grows without bound as the
  // line of sight nears the Sun's centre: Uranus a minute of arc from it in
  // 2280 was thrown 9″ and its speed turned backward. A body there cannot be
  // seen, but its number must still be sane and continuous. Light passing
  // inside the Sun is bent only by the mass it passes outside of, so the
  // deflection is scaled by the Sun's enclosed mass at the impact parameter —
  // which stays near 1 until well inside the limb, because the Sun's mass is
  // concentrated toward its centre (half of it within a quarter of the
  // radius). Swiss Ephemeris tapers it the same way.
  if (body !== SUN) {
    const e: Vec = [E[0] - S[0], E[1] - S[1], E[2] - S[2]];
    const emag = len(e);
    const q: Vec = [E[0] + p[0] - S[0], E[1] + p[1] - S[1], E[2] + p[2] - S[2]];
    const qmag = len(q);
    const eh: Vec = [e[0] / emag, e[1] / emag, e[2] / emag];
    const qh: Vec = [q[0] / qmag, q[1] / qmag, q[2] / qmag];
    const edotp = dot(eh, u);
    if (Math.abs(edotp) < 0.99999999999) {
      let w = (SRS_AU / (emag / AU_KM)) / (1 + dot(qh, eh));
      // Where does the ray pass the Sun? Along it from the observer, at
      // -e·u (km); that point must lie between observer and body to matter.
      const along = -dot(e, u);
      if (along > 0 && along < dist) {
        const b = Math.sqrt(Math.max(0, emag * emag - along * along)); // impact parameter, km
        if (b < R_SUN_KM) w *= enclosedSolarMass(b / R_SUN_KM);
      }
      const pdotq = dot(u, qh);
      u = [u[0] + w * (pdotq * eh[0] - edotp * qh[0]),
        u[1] + w * (pdotq * eh[1] - edotp * qh[1]),
        u[2] + w * (pdotq * eh[2] - edotp * qh[2])];
    }
  }

  // Aberration, relativistic: the observer's barycentric velocity in units of c.
  const v: Vec = [V[0] / C_KM_PER_DAY, V[1] / C_KM_PER_DAY, V[2] / C_KM_PER_DAY];
  const beta2 = dot(v, v);
  const gi = Math.sqrt(1 - beta2);
  const pv = dot(u, v);
  const qf = 1 + pv / (1 + gi);
  const r = 1 + pv;
  const a: Vec = [(gi * u[0] + qf * v[0]) / r, (gi * u[1] + qf * v[1]) / r, (gi * u[2] + qf * v[2]) / r];
  const am = len(a);

  return { u: apply(frame.npb, [a[0] / am, a[1] / am, a[2] / am]), dist: dist / AU_KM };
}

/**
 * Heliocentric place, Swiss Ephemeris's convention: as seen from the Sun, with
 * light-time from the body to the Sun, on the true ecliptic of date — and no
 * aberration or deflection, which belong to an observer on the Earth.
 */
export function heliocentric(tables: Tables, body: number, tdb: number, frame: Frame): Apparent {
  const S: Vec = [0, 0, 0], p: Vec = [0, 0, 0];
  tables.bary(SUN, tdb, S);
  retarded(tables, body, tdb, S, p);
  const dist = len(p);
  return { u: apply(frame.npb, [p[0] / dist, p[1] / dist, p[2] / dist]), dist: dist / AU_KM };
}

/** Longitude and latitude (degrees) on the true ecliptic of date, or RA/Dec. */
export function spherical(ap: Apparent, frame: Frame, equatorial: boolean): [number, number] {
  const w = equatorial ? ap.u : toEcliptic(ap.u, frame.eps);
  const lon = norm360((Math.atan2(w[1], w[0]) * 180) / Math.PI);
  const lat = (Math.atan2(w[2], Math.hypot(w[0], w[1])) * 180) / Math.PI;
  return [lon, lat];
}
