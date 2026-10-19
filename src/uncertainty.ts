// uncertainty.ts — a chart that knows how sure it is.
//
// A birth time is rarely exact: "about half past two", "in the morning",
// "unknown". The ephemeris is precise to a hundredth of an arcsecond; a birth
// time known to ±30 minutes blurs the Ascendant by about ±7°. Software has
// treated the chart as one exact object anyway.
//
// Given a WINDOW for the birth instead of an instant, this says what is
// certain and what is not: for every body and angle, the signs it can be in
// and how much of the window each holds; its house the same way; how much of
// the window each aspect is within orb. With a uniform prior over the window
// — every minute in it equally likely — a probability is the share of the
// window. Every boundary is found exactly (each change refined by bisection to
// about a second), so the shares are exact, not sampled estimates.
//
// Far from the present the clock itself is uncertain (ΔT), which moves the
// planets but not the houses. Its two-sigma spread widens the window the
// planets' SIGNS are judged over, and is reported.

import { norm360 } from "./angles.ts";
import { deltaTUncertainty } from "./deltat.ts";
import type { HouseSystem } from "./houses.ts";

export interface UncertaintyInput {
  /** The birth window, UT Julian Days. A zero-width window is an exact time. */
  from: number;
  to: number;
  lat: number;
  lon: number;
  houseSystem?: HouseSystem;
  /** Bodies to report, Swiss numbering. Defaults to Sun–Pluto, True Node, True Lilith, Chiron. */
  bodies?: number[];
  /** Aspects to test between every pair, with orbs (degrees). */
  aspects?: Array<{ angle: number; orb: number }>;
  /**
   * Instead of `aspects` between every pair: exactly these pairs, each with
   * its own orb, named as the placements are ("body:4", "asc", "mc"). For a
   * chart whose orbs depend on the bodies — a luminary's wider than an outer
   * planet's — so the shares describe the aspects that chart draws.
   */
  pairs?: Array<{ a: string; b: string; angle: number; orb: number }>;
}

export interface Share<T> {
  value: T;
  /** Share of the window, 0–1. */
  probability: number;
}

export interface Placement {
  /** "body:4" for Mars, "asc", "mc". */
  key: string;
  body?: number;
  /** Sign indices (0 = Aries) by share, largest first. */
  signs: Share<number>[];
  /** House numbers (1–12) by share; absent for the angles. */
  houses?: Share<number>[];
  /** Longitude at the start and end of the window, unwrapped: end − start is the motion. */
  range: [number, number];
  /** Sign (and house, for a body) hold for the whole window. */
  certain: boolean;
}

export interface AspectShare {
  a: string;
  b: string;
  angle: number;
  /** Share of the window during which the aspect is within orb. */
  probability: number;
}

export interface Uncertainty {
  from: number;
  to: number;
  minutes: number;
  placements: Placement[];
  aspects: AspectShare[];
  /** The clock: ΔT at the window and its one-sigma uncertainty, seconds. */
  deltaT: { seconds: number; sigma: number };
}

/** What the uncertainty pass needs from the engine. */
export interface ChartSource {
  /** Apparent ecliptic longitude of a body at a UT moment. */
  longitude(body: number, jdUT: number): number;
  houses(jdUT: number, lat: number, lon: number, system: HouseSystem): { cusps: number[]; ascendant: number; mc: number };
  deltaT(jdUT: number): number;
}

const DEFAULT_BODIES = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 13, 15];
const DEFAULT_ASPECTS = [
  { angle: 0, orb: 8 }, { angle: 60, orb: 6 }, { angle: 90, orb: 8 }, { angle: 120, orb: 8 }, { angle: 180, orb: 8 },
];
/** A share at or above this is reported as certain. */
export const CERTAIN = 0.999;
const SECOND = 1 / 86400;

const signOf = (lon: number) => Math.floor(norm360(lon) / 30) % 12;
const sep = (a: number, b: number) => Math.abs(((a - b + 540) % 360) - 180);

function houseOf(lon: number, cusps: number[]): number {
  const x = norm360(lon);
  for (let i = 0; i < 12; i++) {
    const a = cusps[i], b = cusps[(i + 1) % 12];
    const inside = b > a ? x >= a && x < b : x >= a || x < b;
    if (inside) return i + 1;
  }
  return 1;
}

/**
 * How long a discrete property holds each value over a grid `t` (whose
 * samples' values are `cats`). Each change between neighbouring samples is
 * refined by bisection with `cat`, to about a second.
 */
function shares(t: Float64Array, cats: Int32Array, cat: (x: number) => number): Share<number>[] {
  const a = t[0], b = t[t.length - 1];
  const total = new Map<number, number>();
  if (b <= a) return [{ value: cats[0], probability: 1 }];
  let segStart = a, current = cats[0];
  for (let i = 1; i < t.length; i++) {
    if (cats[i] === current) continue;
    let lo = t[i - 1], hi = t[i];
    while (hi - lo > SECOND) {
      const mid = (lo + hi) / 2;
      if (cat(mid) === current) lo = mid; else hi = mid;
    }
    total.set(current, (total.get(current) ?? 0) + (hi - segStart));
    segStart = hi;
    current = cats[i];
  }
  total.set(current, (total.get(current) ?? 0) + (b - segStart));
  return [...total]
    .map(([value, d]) => ({ value, probability: d / (b - a) }))
    .sort((x, y) => y.probability - x.probability);
}

function grid(a: number, b: number, n: number): Float64Array {
  const t = new Float64Array(n);
  for (let i = 0; i < n; i++) t[i] = n === 1 ? a : a + ((b - a) * i) / (n - 1);
  return t;
}

const unwrapRange = (lons: Float64Array): [number, number] => {
  let end = lons[0];
  for (let i = 1; i < lons.length; i++) end += ((lons[i] - lons[i - 1] + 540) % 360) - 180;
  return [lons[0], end];
};

export function chartUncertainty(src: ChartSource, input: UncertaintyInput): Uncertainty {
  const { from, to, lat, lon } = input;
  const system = input.houseSystem ?? "P";
  const bodies = input.bodies ?? DEFAULT_BODIES;
  const aspects = input.aspects ?? DEFAULT_ASPECTS;
  const mid = (from + to) / 2;
  const sigma = deltaTUncertainty(mid);
  const exact = to <= from;
  // About one sample per 20 s, at most 2,000; the bisection does the rest.
  const n = exact ? 1 : Math.max(3, Math.min(2000, Math.ceil(((to - from) * 86400) / 20) + 1));
  const t = grid(from, to, n);

  // Everything at every sample, once.
  const houses = Array.from(t, (x) => src.houses(x, lat, lon, system));
  const lons = new Map<string, Float64Array>();
  const at = new Map<string, (x: number) => number>();
  for (const b of bodies) {
    const key = `body:${b}`;
    lons.set(key, Float64Array.from(t, (x) => src.longitude(b, x)));
    at.set(key, (x) => src.longitude(b, x));
  }
  lons.set("asc", Float64Array.from(houses, (h) => h.ascendant));
  at.set("asc", (x) => src.houses(x, lat, lon, system).ascendant);
  lons.set("mc", Float64Array.from(houses, (h) => h.mc));
  at.set("mc", (x) => src.houses(x, lat, lon, system).mc);

  const placements: Placement[] = [];
  const spread = (2 * sigma) * SECOND; // ΔT's two sigma, as days
  for (const b of bodies) {
    const key = `body:${b}`;
    const L = lons.get(key)!, f = at.get(key)!;
    // Signs over the window widened by ΔT's uncertainty (when it matters).
    let signs: Share<number>[];
    if (spread > SECOND && !exact) {
      const tw = grid(from - spread, to + spread, n);
      signs = shares(tw, Int32Array.from(tw, (x) => signOf(f(x))), (x) => signOf(f(x)));
    } else {
      signs = exact ? [{ value: signOf(L[0]), probability: 1 }] : shares(t, Int32Array.from(L, signOf), (x) => signOf(f(x)));
    }
    const houseShares = exact
      ? [{ value: houseOf(L[0], houses[0].cusps), probability: 1 }]
      : shares(t, Int32Array.from(L, (v, i) => houseOf(v, houses[i].cusps)), (x) => houseOf(f(x), src.houses(x, lat, lon, system).cusps));
    placements.push({
      key, body: b, signs, houses: houseShares, range: unwrapRange(L),
      certain: signs[0].probability >= CERTAIN && houseShares[0].probability >= CERTAIN,
    });
  }
  for (const key of ["asc", "mc"]) {
    const L = lons.get(key)!, f = at.get(key)!;
    const signs = exact ? [{ value: signOf(L[0]), probability: 1 }] : shares(t, Int32Array.from(L, signOf), (x) => signOf(f(x)));
    placements.push({ key, signs, range: unwrapRange(L), certain: signs[0].probability >= CERTAIN });
  }

  const keys = [...lons.keys()];
  const tests = input.pairs ?? keys.flatMap((a, i) =>
    keys.slice(i + 1).flatMap((b) => aspects.map(({ angle, orb }) => ({ a, b, angle, orb }))));
  const out: AspectShare[] = [];
  for (const { a, b, angle, orb } of tests) {
    const A = lons.get(a), B = lons.get(b), fa = at.get(a), fb = at.get(b);
    if (!A || !B || !fa || !fb) continue; // not among the bodies asked for
    const inOrb = (x: number, y: number) => (Math.abs(sep(x, y) - angle) <= orb ? 1 : 0);
    const cats = Int32Array.from(A, (v, k) => inOrb(v, B[k]));
    if (!cats.some((c) => c === 1)) continue;
    const s = exact ? [{ value: cats[0], probability: 1 }] : shares(t, cats, (x) => inOrb(fa(x), fb(x)));
    const p = s.find((e) => e.value === 1)?.probability ?? 0;
    if (p > 0) out.push({ a, b, angle, probability: p });
  }
  return { from, to, minutes: (to - from) * 1440, placements, aspects: out, deltaT: { seconds: src.deltaT(mid), sigma } };
}
