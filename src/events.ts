// events.ts — exact moments, found rather than sampled.
//
// A transit, an ingress or a return is the moment a longitude reaches a given
// degree. Software usually steps through time and looks for the degree being
// passed, which misses any pass shorter than its step — and a planet turning
// retrograde can cross the same degree three times in a few weeks.
//
// Here every search is COMPLETE by construction:
//   1. Sample the longitude at under half the shortest time between two
//      stations that body ever has (measured across the whole span of the
//      tables). Then at least two samples fall inside every retrograde arc, so
//      the sampled track turns wherever the body does, and the two-step window
//      each turn is refined in holds no other turn.
//   2. Refine each turning point (station) exactly. A "reversal" smaller than
//      an arcsecond is numerical, not a station, and is dropped.
//   3. Between consecutive stations the longitude moves one way only, so each
//      target degree is reached at most once per lap there; solve for it.
// The result is every pass, each to a few milliseconds.

import { norm360 } from "./angles.ts";

/** A longitude-like angle in degrees as a function of TT. */
export type AngleFn = (jdTT: number) => number;

export interface Motion {
  /** Sampling step, days: under half the shortest gap between stations. */
  step: number;
  /** Fastest speed, °/day: bounds how far one sample can move. */
  maxSpeed: number;
  /** The body never turns (Sun, Moon, mean node, mean Lilith). */
  monotone?: boolean;
}

/**
 * Per body: Swiss numbering. The shortest gaps between stations, measured
 * across 1550–2649 (the asteroids 1600–2497): Mercury 19.75 days, Venus 40.5,
 * Mars 59.5, Jupiter 117, Saturn 133, Uranus 148, Neptune and Pluto 156,
 * Chiron 128, Ceres 91, Pallas 74.5, Juno 72, Vesta 84.5 — each at least ten
 * steps. The true node and true Lilith turn far more often, and most of their
 * turns are sub-arcsecond wobbles that step 2 drops; between turns that are
 * not (a reversal of 1″ or more), the shortest gaps found by brute force in
 * three windows across the span were 0.525 days for the node and 0.45 for
 * Lilith: 3.5 and 4.5 steps. (The node's step was 0.25, which left it 2.1 —
 * too near the two that step 1 needs, from a sample of a ninth of the span.)
 */
export const MOTION: Record<number, Motion> = {
  0: { step: 5, maxSpeed: 1.03, monotone: true },
  1: { step: 4, maxSpeed: 15.5, monotone: true },
  2: { step: 2, maxSpeed: 2.5 },
  3: { step: 4, maxSpeed: 1.3 },
  4: { step: 6, maxSpeed: 0.8 },
  5: { step: 10, maxSpeed: 0.25 },
  6: { step: 10, maxSpeed: 0.17 },
  7: { step: 12, maxSpeed: 0.07 },
  8: { step: 12, maxSpeed: 0.045 },
  9: { step: 12, maxSpeed: 0.045 },
  10: { step: 10, maxSpeed: 0.06, monotone: true },
  11: { step: 0.15, maxSpeed: 0.3 },
  12: { step: 10, maxSpeed: 0.12, monotone: true },
  13: { step: 0.1, maxSpeed: 6.5 },
  15: { step: 10, maxSpeed: 0.15 },
  17: { step: 6, maxSpeed: 0.47 },
  18: { step: 6, maxSpeed: 0.62 },
  19: { step: 6, maxSpeed: 0.6 },
  20: { step: 6, maxSpeed: 0.55 },
};

/** A reversal smaller than this (degrees, 1″) is numerical, not a station. */
const MIN_REVERSAL = 1 / 3600;
/** Solve times to this (days): about 9 ms. */
const TIME_TOL = 1e-7;

const signed = (d: number) => {
  const r = norm360(d);
  return r > 180 ? r - 360 : r;
};

export interface Station {
  jdTT: number;
  /** Unwrapped longitude at the station, degrees (norm360 it to display). */
  longitude: number;
  /** "retrograde": turning backward here; "direct": turning forward. */
  turning: "retrograde" | "direct";
}

export interface Crossing {
  jdTT: number;
  /** The degree reached, 0–360. */
  target: number;
  /** +1 moving forward through it, −1 moving backward (retrograde). */
  direction: 1 | -1;
}

/** Golden-section search for a minimum of `g` on [a, b]. */
function golden(g: (t: number) => number, a: number, b: number): number {
  const r = (Math.sqrt(5) - 1) / 2;
  let c = b - r * (b - a), d = a + r * (b - a), gc = g(c), gd = g(d);
  while (b - a > TIME_TOL) {
    if (gc < gd) { b = d; d = c; gd = gc; c = b - r * (b - a); gc = g(c); }
    else { a = c; c = d; gc = gd; d = a + r * (b - a); gd = g(d); }
  }
  return (a + b) / 2;
}

/** Sampled, unwrapped track of an angle. */
function track(f: AngleFn, start: number, end: number, step: number) {
  const n = Math.max(2, Math.ceil((end - start) / step) + 1);
  const h = (end - start) / (n - 1);
  const t = new Float64Array(n), u = new Float64Array(n);
  let prev = f(start);
  t[0] = start;
  u[0] = prev;
  for (let i = 1; i < n; i++) {
    t[i] = start + i * h;
    const v = f(t[i]);
    u[i] = u[i - 1] + signed(v - prev);
    prev = v;
  }
  return { t, u };
}

/** Every station of an angle between `start` and `end` (TT). */
export function stations(f: AngleFn, start: number, end: number, motion: Motion): Station[] {
  if (motion.monotone) return [];
  const { t, u } = track(f, start, end, motion.step);
  const found: Station[] = [];
  for (let i = 1; i < t.length - 1; i++) {
    const before = u[i] - u[i - 1], after = u[i + 1] - u[i];
    if (before === 0 || after === 0 || Math.sign(before) === Math.sign(after)) continue;
    const isMax = before > 0; // forward then backward: turning retrograde
    const ref = u[i], anchor = f(t[i]);
    const unwrap = (x: number) => ref + signed(f(x) - anchor);
    const at = golden((x) => (isMax ? -unwrap(x) : unwrap(x)), t[i - 1], t[i + 1]);
    found.push({ jdTT: at, longitude: unwrap(at), turning: isMax ? "retrograde" : "direct" });
  }
  // Drop numerical wobbles: neighbouring turns that undo each other by less
  // than an arcsecond were never a real retrograde loop. Smallest first, so a
  // wobble beside a real station takes its own partner with it, not the
  // station: the True Node steps by ~0.03″ where two Moon segments join, and
  // swept left to right, a join just after a station would pair the station
  // with the step and leave the step standing in its place.
  for (;;) {
    let k = -1, least = MIN_REVERSAL;
    for (let i = 0; i + 1 < found.length; i++) {
      const d = Math.abs(found[i + 1].longitude - found[i].longitude);
      if (d < least) { least = d; k = i; }
    }
    if (k < 0) break;
    found.splice(k, 2);
  }
  return found;
}

/** Illinois (regula falsi) root of g on [a, b], given g(a) and g(b) of opposite sign. */
function root(g: (t: number) => number, a: number, b: number, ga: number, gb: number): number {
  let side = 0;
  for (let i = 0; i < 100 && b - a > TIME_TOL; i++) {
    const c = (a * gb - b * ga) / (gb - ga);
    const gc = g(c);
    if (gc === 0) return c;
    if (Math.sign(gc) === Math.sign(gb)) {
      b = c; gb = gc;
      if (side === -1) ga /= 2;
      side = -1;
    } else {
      a = c; ga = gc;
      if (side === 1) gb /= 2;
      side = 1;
    }
  }
  return (a + b) / 2;
}

/**
 * Every moment between `start` and `end` (TT) when the angle reaches any of
 * `targets` (degrees), in time order — every pass, retrograde ones included.
 */
export function crossings(f: AngleFn, targets: number[], start: number, end: number, motion: Motion, known?: Station[]): Crossing[] {
  const turns = (known ?? stations(f, start, end, motion)).map((s) => s.jdTT).filter((x) => x > start && x < end);
  const edges = [start, ...turns, end];
  // Chunks short enough that the angle moves under 90° in each.
  const chunk = Math.min(motion.step, 90 / motion.maxSpeed);
  const wanted = targets.map(norm360);
  const out: Crossing[] = [];
  for (let k = 0; k + 1 < edges.length; k++) {
    const a0 = edges[k], b0 = edges[k + 1];
    const pieces = Math.max(1, Math.ceil((b0 - a0) / chunk));
    let a = a0, fa = f(a);
    for (let p = 1; p <= pieces; p++) {
      const b = p === pieces ? b0 : a0 + ((b0 - a0) * p) / pieces;
      const fb = f(b);
      const move = signed(fb - fa); // monotone and under 90°: the swept arc
      if (move !== 0) {
        const dir: 1 | -1 = move > 0 ? 1 : -1;
        for (const target of wanted) {
          // Is the target strictly inside the arc swept from fa to fb?
          const from = signed(target - fa) * dir; // how far ahead, in the direction of motion
          const reach = Math.abs(move);
          if (from > 0 && from <= reach) {
            const g = (x: number) => signed(f(x) - target) * dir;
            const ga = signed(fa - target) * dir, gb = signed(fb - target) * dir;
            const at = gb === 0 ? b : root(g, a, b, ga, gb);
            out.push({ jdTT: at, target, direction: dir });
          }
        }
      }
      a = b;
      fa = fb;
    }
  }
  return out.sort((x, y) => x.jdTT - y.jdTT);
}

/** The step and speed bound for the difference of two bodies' longitudes. */
export function combined(a: Motion, b: Motion): Motion {
  return { step: Math.min(a.step, b.step), maxSpeed: a.maxSpeed + b.maxSpeed, monotone: false };
}
