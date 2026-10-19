// Exact moments: stations, crossings, ingresses, returns, lunations, aspects
// between moving bodies. The claim is COMPLETENESS — every pass found, the
// three passes of a retrograde loop included — so each search is held to a
// brute-force count over the same span, and each moment to its definition.

import { describe, it, expect } from "vitest";
import { readFileSync } from "node:fs";
import { resolve } from "node:path";
import { NaelosEphemeris, julday, SEFLG_SPEED, SE_SUN, SE_MOON, SE_MERCURY, SE_MARS, SE_JUPITER, SE_SATURN, SE_TRUE_NODE, SE_OSCU_APOG } from "../src/index.ts";
import { stations } from "../src/events.ts";

const ROOT = resolve(__dirname, "..");
const manifest = JSON.parse(readFileSync(resolve(ROOT, "data/manifest.json"), "utf8"));
const eph = NaelosEphemeris.fromBuffers(manifest.files.map((f: { name: string }) => readFileSync(resolve(ROOT, "data", f.name))));
const wrap = (d: number) => ((d + 540) % 360) - 180;
const lon = (jdUT: number, body: number) => {
  const r = eph.calc_ut(jdUT, body, SEFLG_SPEED);
  if ("error" in r) throw new Error(r.error);
  return r;
};

/**
 * Count the times a body passes a degree by stepping finely: the slow, sure
 * way. A sample landing exactly on the degree is skipped rather than counted
 * twice (arriving at zero, then leaving it).
 */
function brute(body: number, target: number, a: number, b: number, step: number) {
  let n = 0, prev = wrap(lon(a, body).longitude - target);
  for (let t = a + step; t <= b; t += step) {
    const d = wrap(lon(t, body).longitude - target);
    if (d === 0) continue;
    if (Math.sign(d) !== Math.sign(prev) && Math.abs(d - prev) < 90) n++;
    prev = d;
  }
  return n;
}

describe("completeness: the exact search finds what brute force finds", () => {
  const cases: Array<[string, number, number, number, number]> = [
    // name, body, from, to, brute step (days)
    ["Mercury, five years", SE_MERCURY, julday(2015, 1, 1), julday(2020, 1, 1), 1 / 24],
    ["Mars, ten years", SE_MARS, julday(1990, 1, 1), julday(2000, 1, 1), 1 / 8],
    ["the True Node, three years", SE_TRUE_NODE, julday(1995, 1, 1), julday(1998, 1, 1), 1 / 24],
    ["True Lilith, a year", SE_OSCU_APOG, julday(2020, 1, 1), julday(2021, 1, 1), 1 / 96],
  ];
  for (const [name, body, a, b, step] of cases) {
    it(`${name}`, () => {
      // A degree the body certainly reaches: where it stands mid-span.
      const target = lon((a + b) / 2, body).longitude;
      const found = eph.crossings(body, [target], a, b);
      expect(found.length).toBeGreaterThan(0);
      expect(found.length).toBe(brute(body, target, a, b, step));
      for (const c of found) expect(Math.abs(wrap(lon(c.jd, body).longitude - target))).toBeLessThan(0.5 / 3600);
    });
  }

  it("a retrograde loop crosses its degree three times, and says which way", () => {
    // Mercury's 2023 April–May retrograde in Taurus: a degree inside the loop.
    const st = eph.stations(SE_MERCURY, julday(2023, 4, 1), julday(2023, 6, 1));
    expect(st.map((s) => s.turning)).toEqual(["retrograde", "direct"]);
    const inside = (st[0].longitude + st[1].longitude) / 2;
    const passes = eph.crossings(SE_MERCURY, [inside], julday(2023, 3, 1), julday(2023, 7, 1));
    expect(passes.map((p) => p.direction)).toEqual([1, -1, 1]);
  });
});

describe("stations are where the body stands still", () => {
  it("Mercury: speed zero at each, to a ten-thousandth of a degree a day", () => {
    const st = eph.stations(SE_MERCURY, julday(2000, 1, 1), julday(2010, 1, 1));
    expect(st.length).toBeGreaterThan(55); // three retrogrades a year, two stations each
    for (const s of st) expect(Math.abs(lon(s.jd, SE_MERCURY).longitudeSpeed)).toBeLessThan(1e-4);
    for (let i = 1; i < st.length; i++) expect(st[i].turning).not.toBe(st[i - 1].turning);
  });

  it("the Sun and Moon never station", () => {
    expect(eph.stations(SE_SUN, julday(2000, 1, 1), julday(2001, 1, 1))).toEqual([]);
    expect(eph.stations(SE_MOON, julday(2000, 1, 1), julday(2001, 1, 1))).toEqual([]);
  });

  it("drops a sub-arcsecond wobble beside a real station, not the station", () => {
    // A minimum at 10, a 0.1″ wobble at 20–30, a maximum at 40 (arcseconds,
    // piecewise linear). Swept left to right, the minimum and the wobble's
    // first turn undo each other by 0.5″ and were dropped together, leaving
    // the wobble's second turn reported as the station, 20 days late.
    const knots: Array<[number, number]> = [[0, 5], [10, 0], [20, 0.5], [30, 0.4], [40, 3.4], [50, 0]];
    const f = (t: number) => {
      const next = knots.findIndex(([x]) => x > t);
      const i = next < 0 ? knots.length - 2 : Math.max(0, next - 1);
      const [x0, y0] = knots[i], [x1, y1] = knots[i + 1];
      return (y0 + ((t - x0) * (y1 - y0)) / (x1 - x0)) / 3600;
    };
    const found = stations(f, 0, 50, { step: 1, maxSpeed: 1 });
    expect(found.map((s) => [Math.round(s.jdTT * 1000) / 1000, s.turning])).toEqual([[10, "direct"], [40, "retrograde"]]);
  });
});

describe("ingresses, returns and lunations", () => {
  it("the Sun enters each sign once a year, in order", () => {
    const ing = eph.ingresses(SE_SUN, julday(2025, 1, 1), julday(2026, 1, 1));
    expect(ing.length).toBe(12);
    expect(ing.map((i) => i.sign)).toEqual([10, 11, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9]); // Aquarius first
  });

  it("a retrograde re-entry enters the sign behind", () => {
    const back = eph.ingresses(SE_MERCURY, julday(2000, 1, 1), julday(2010, 1, 1)).filter((i) => i.direction === -1);
    expect(back.length).toBeGreaterThan(0);
    for (const i of back) expect(i.sign).toBe((Math.round(i.target / 30) + 11) % 12);
  });

  it("a solar return comes once a year, at the natal degree", () => {
    const natal = 164.0854;
    const ret = eph.returns(SE_SUN, natal, julday(2000, 1, 1), julday(2010, 1, 1));
    expect(ret.length).toBe(10);
    for (const r of ret) expect(Math.abs(wrap(lon(r.jd, SE_SUN).longitude - natal)) * 3600).toBeLessThan(0.05);
  });

  it("lunations come in order, about 7.4 days apart", () => {
    const l = eph.lunations(julday(2024, 1, 1), julday(2025, 1, 1));
    expect(l.length).toBeGreaterThanOrEqual(49);
    const order = ["new", "first quarter", "full", "last quarter"];
    for (let i = 1; i < l.length; i++) {
      expect(order.indexOf(l[i].phase)).toBe((order.indexOf(l[i - 1].phase) + 1) % 4);
      const gap = l[i].jd - l[i - 1].jd;
      expect(gap).toBeGreaterThan(6.5);
      expect(gap).toBeLessThan(8.5);
    }
  });

  it("the new Moon of the 2024-04-08 eclipse falls within an hour of greatest eclipse", () => {
    const e = eph.solarEclipse(julday(2024, 4, 1))!;
    const nm = eph.lunations(julday(2024, 4, 7), julday(2024, 4, 10)).find((x) => x.phase === "new")!;
    expect(Math.abs(nm.jd - e.jd) * 24).toBeLessThan(1);
  });
});

describe("aspects between moving bodies", () => {
  it("the great conjunction of Jupiter and Saturn on 2020-12-21, and only one that year", () => {
    const c = eph.aspectsBetween(SE_JUPITER, SE_SATURN, [0], julday(2020, 1, 1), julday(2021, 1, 1));
    expect(c.length).toBe(1);
    expect(Math.abs(c[0].jd - julday(2020, 12, 21, 12))).toBeLessThan(0.75);
  });
});
