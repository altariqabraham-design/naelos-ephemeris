// rectification.ts — narrowing an uncertain birth time from the life lived.
//
// The rectification astrologers do by hand, made explicit and checkable: for
// each candidate birth minute in the window, how strongly do the person's
// dated life events coincide with contacts to that minute's ANGLES — the
// Ascendant and Midheaven, which move a degree every four minutes and so are
// what a birth time decides? Three classical techniques are weighed:
//
//   transits    Jupiter–Pluto (and Mars, lightly) on the event date, to the
//               natal angles, by conjunction, opposition or square;
//   solar arc   every natal planet carried forward by the solar arc to the
//               angles, and the angles carried forward to the planets;
//   progressed  the secondary-progressed Moon to the angles.
//
// Each contact contributes by how exact it is (a Gaussian in the orb, wider
// for transits than for the slow-moving directions). The result is SUPPORT —
// relative, summing to 1 over the window — with every contact behind it
// listed. It is not a measured probability: it says which minutes the chosen
// techniques favour, and why, so the reasoning can be checked and argued with.

import { norm360 } from "./angles.ts";
import type { HouseSystem } from "./houses.ts";

export interface LifeEvent {
  /** When it happened, as a UT Julian Day (noon if only the day is known). */
  jd: number;
  /** How much it should count (default 1): a marriage more than a short trip. */
  weight?: number;
  label?: string;
}

export interface RectifyInput {
  /** The window the birth lies in, UT Julian Days. */
  from: number;
  to: number;
  lat: number;
  lon: number;
  events: LifeEvent[];
  /** Candidate spacing, minutes (default 1). */
  stepMinutes?: number;
  houseSystem?: HouseSystem;
}

export type Technique = "transit" | "solar arc" | "progressed Moon";

export interface Evidence {
  event: number;
  technique: Technique;
  /** What moved: a body number, or "Asc"/"MC" when a directed angle reached a planet. */
  mover: number | "Asc" | "MC";
  /** What it reached: "Asc"/"MC", or a natal body number. */
  point: number | "Asc" | "MC";
  aspect: 0 | 90 | 180;
  /** How far from exact, degrees. */
  orb: number;
  contribution: number;
}

export interface Candidate {
  jd: number;
  ascendant: number;
  mc: number;
  score: number;
  /** Share of the window's total support, 0–1. */
  support: number;
  evidence: Evidence[];
}

export interface Rectification {
  candidates: Candidate[];
  /** The strongest few, best first. */
  best: Candidate[];
  /** The narrowest span of candidates holding 68% of the support. */
  window68: { from: number; to: number };
  techniques: string;
}

/** What rectification needs from the engine. */
export interface RectifySource {
  longitude(body: number, jdUT: number): number;
  angles(jdUT: number, lat: number, lon: number, system: HouseSystem): { ascendant: number; mc: number };
}

const TRANSITING: Array<[number, number]> = [[9, 1.0], [8, 0.9], [7, 1.0], [6, 0.9], [5, 0.6], [4, 0.3]];
const NATAL = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9];
const ASPECTS: Array<[0 | 90 | 180, number]> = [[0, 1.0], [180, 0.8], [90, 0.7]];
/** Gaussian widths, degrees: transits spread wider than the slow directions. */
const SIGMA: Record<Technique, number> = { transit: 1.0, "solar arc": 0.33, "progressed Moon": 0.5 };
const TROPICAL_YEAR = 365.2422;

const sep = (a: number, b: number) => Math.abs(((a - b + 540) % 360) - 180);

/** The closest of the three aspects between two longitudes, and its orb. */
function closest(a: number, b: number): { aspect: 0 | 90 | 180; orb: number; weight: number } {
  const s = sep(a, b);
  let best = { aspect: 0 as 0 | 90 | 180, orb: Infinity, weight: 0 };
  for (const [aspect, weight] of ASPECTS) {
    const orb = Math.abs(s - aspect);
    if (orb < best.orb) best = { aspect, orb, weight };
  }
  return best;
}

export function rectify(src: RectifySource, input: RectifyInput): Rectification {
  const { from, to, lat, lon, events } = input;
  const system = input.houseSystem ?? "P";
  const step = (input.stepMinutes ?? 1) / 1440;
  const n = Math.max(1, Math.floor((to - from) / step + 1e-9) + 1);

  // The sky on each event date does not depend on the birth minute: once each.
  const transits = events.map((e) => TRANSITING.map(([b, w]) => ({ body: b, weight: w, lon: src.longitude(b, e.jd) })));

  const candidates: Candidate[] = [];
  for (let i = 0; i < n; i++) {
    const tb = Math.min(to, from + i * step);
    const { ascendant, mc } = src.angles(tb, lat, lon, system);
    const natal = NATAL.map((b) => ({ body: b, lon: src.longitude(b, tb) }));
    const sunAtBirth = natal[0].lon;
    const evidence: Evidence[] = [];
    let score = 0;
    const add = (ev: Omit<Evidence, "contribution">, weight: number) => {
      const sigma = SIGMA[ev.technique];
      if (ev.orb > 3 * sigma) return;
      const contribution = weight * Math.exp(-0.5 * (ev.orb / sigma) ** 2);
      score += contribution;
      evidence.push({ ...ev, contribution });
    };
    events.forEach((e, k) => {
      const w = e.weight ?? 1;
      const angles: Array<["Asc" | "MC", number]> = [["Asc", ascendant], ["MC", mc]];
      // Transits to the angles.
      for (const tr of transits[k]) {
        for (const [name, a] of angles) {
          const c = closest(tr.lon, a);
          add({ event: k, technique: "transit", mover: tr.body, point: name, aspect: c.aspect, orb: c.orb }, w * tr.weight * c.weight);
        }
      }
      // A day after birth for each year of life: the secondary progression.
      const tp = tb + (e.jd - tb) / TROPICAL_YEAR;
      const arc = norm360(src.longitude(0, tp) - sunAtBirth);
      for (const p of natal) {
        for (const [name, a] of angles) {
          const c1 = closest(p.lon + arc, a); // directed planet to natal angle
          add({ event: k, technique: "solar arc", mover: p.body, point: name, aspect: c1.aspect, orb: c1.orb }, w * c1.weight);
          const c2 = closest(a + arc, p.lon); // directed angle to natal planet
          add({ event: k, technique: "solar arc", mover: name, point: p.body, aspect: c2.aspect, orb: c2.orb }, w * c2.weight * 0.8);
        }
      }
      const moon = src.longitude(1, tp);
      for (const [name, a] of angles) {
        const c = closest(moon, a);
        add({ event: k, technique: "progressed Moon", mover: 1, point: name, aspect: c.aspect, orb: c.orb }, w * 0.6 * c.weight);
      }
    });
    candidates.push({ jd: tb, ascendant, mc, score, support: 0, evidence: evidence.sort((x, y) => y.contribution - x.contribution) });
  }

  // Support: a softmax of the scores, so one clearly better minute stands out
  // without the rest vanishing. With no events, it is uniform.
  const top = Math.max(...candidates.map((c) => c.score));
  const weights = candidates.map((c) => Math.exp(c.score - top));
  const sum = weights.reduce((s, x) => s + x, 0);
  candidates.forEach((c, i) => (c.support = weights[i] / sum));

  // The narrowest run of consecutive candidates holding 68% of the support.
  let bestSpan = { from: from, to: to }, width = Infinity;
  for (let a = 0, b = 0, acc = 0; a < candidates.length; a++) {
    while (b < candidates.length && acc < 0.68) acc += candidates[b++].support;
    if (acc >= 0.68 && candidates[b - 1].jd - candidates[a].jd < width) {
      width = candidates[b - 1].jd - candidates[a].jd;
      bestSpan = { from: candidates[a].jd, to: candidates[b - 1].jd };
    }
    acc -= candidates[a].support;
  }

  return {
    candidates,
    best: [...candidates].sort((x, y) => y.score - x.score).slice(0, 5),
    window68: bestSpan,
    techniques: "transits of Mars–Pluto, solar-arc directions and the secondary-progressed Moon, to and from the Ascendant and Midheaven; conjunction, opposition and square",
  };
}

/** Exposed for tests. */
export const _internal = { closest };
