// A chart that knows how sure it is — and rectification from the life lived.

import { describe, it, expect } from "vitest";
import { readFileSync } from "node:fs";
import { resolve } from "node:path";
import { NaelosEphemeris, julday, norm360, CERTAIN } from "../src/index.ts";

const ROOT = resolve(__dirname, "..");
const manifest = JSON.parse(readFileSync(resolve(ROOT, "data/manifest.json"), "utf8"));
const eph = NaelosEphemeris.fromBuffers(manifest.files.map((f: { name: string }) => readFileSync(resolve(ROOT, "data", f.name))));
const CHARLOTTE = { lat: 35.2271, lon: -80.8431 };
const BIRTH = julday(1992, 9, 6, 10 + 41 / 60); // 06:41 EDT
const signOf = (lon: number) => Math.floor(norm360(lon) / 30);

describe("a birth window", () => {
  const u = eph.chartUncertainty({ from: BIRTH - 30 / 1440, to: BIRTH + 30 / 1440, ...CHARLOTTE });
  const get = (key: string) => u.placements.find((p) => p.key === key)!;

  it("gives every placement shares that add up to the whole window", () => {
    for (const p of u.placements) {
      expect(p.signs.reduce((s, x) => s + x.probability, 0)).toBeCloseTo(1, 9);
      if (p.houses) expect(p.houses.reduce((s, x) => s + x.probability, 0)).toBeCloseTo(1, 9);
    }
  });

  it("holds the Ascendant in Virgo for all of ±30 minutes, as the geometry says", () => {
    expect(get("asc").signs).toEqual([{ value: 5, probability: 1 }]);
    expect(get("asc").certain).toBe(true);
  });

  it("knows the Sun's house is not certain when a cusp sweeps past it", () => {
    const sun = get("body:0");
    expect(sun.signs[0].probability).toBe(1);
    expect(sun.houses!.length).toBe(2);
    expect(sun.certain).toBe(false);
  });

  it("agrees with a second-by-second count", () => {
    // Independent check: each MC sign's share among 3,601 samples a second apart.
    const counts = new Map<number, number>();
    const N = 3601;
    for (let i = 0; i < N; i++) {
      const t = BIRTH - 30 / 1440 + (i / (N - 1)) * (60 / 1440);
      const s = signOf(eph.houses(t, CHARLOTTE.lat, CHARLOTTE.lon, "P").mc);
      counts.set(s, (counts.get(s) ?? 0) + 1);
    }
    for (const share of get("mc").signs) expect(Math.abs(share.probability - (counts.get(share.value) ?? 0) / N)).toBeLessThan(2 / N);
  });

  it("an exact time is certain in everything", () => {
    const e = eph.chartUncertainty({ from: BIRTH, to: BIRTH, ...CHARLOTTE });
    for (const p of e.placements) expect(p.certain).toBe(true);
    for (const a of e.aspects) expect(a.probability).toBe(1);
  });

  it("tests exactly the pairs it is given, each with its own orb", () => {
    // The Moon's aspects across a whole unknown day, where a wider or narrower
    // orb genuinely changes the share. Asked pair by pair, each share must be
    // what the every-pair pass gives at that orb.
    // On this day the Moon squares Venus for about half the day at 8° and a
    // sixteenth at 3°; she trines Jupiter for half at 8° and a tenth at 3°.
    const day = { from: julday(1992, 9, 6, 4), to: julday(1992, 9, 7, 4), ...CHARLOTTE };
    const both = [{ angle: 90, orb: 0 }, { angle: 120, orb: 0 }];
    const wide = eph.chartUncertainty({ ...day, aspects: both.map((x) => ({ ...x, orb: 8 })) });
    const tight = eph.chartUncertainty({ ...day, aspects: both.map((x) => ({ ...x, orb: 3 })) });
    const pairs = [
      { a: "body:1", b: "body:3", angle: 90, orb: 8 },
      { a: "body:1", b: "body:5", angle: 120, orb: 3 },
      { a: "body:1", b: "body:99", angle: 90, orb: 8 }, // not cast: ignored
    ];
    const asked = eph.chartUncertainty({ ...day, pairs });
    const share = (u: typeof asked, b: string, angle: number) =>
      u.aspects.find((x) => x.a === "body:1" && x.b === b && x.angle === angle)?.probability ?? 0;
    expect(asked.aspects.map((x) => `${x.b}@${x.angle}`)).toEqual(["body:3@90", "body:5@120"]);
    expect(share(asked, "body:3", 90)).toBe(share(wide, "body:3", 90));
    expect(share(asked, "body:5", 120)).toBe(share(tight, "body:5", 120));
    // and the orb mattered: these are genuine fractions, different at the other orb
    expect(share(wide, "body:3", 90)).toBeGreaterThan(0.4);
    expect(share(tight, "body:3", 90)).toBeLessThan(0.2);
    expect(share(tight, "body:5", 120)).toBeGreaterThan(0.05);
    expect(share(wide, "body:5", 120)).toBeGreaterThan(share(tight, "body:5", 120) * 3);
  });
});

describe("an unknown birth time (the whole day)", () => {
  const u = eph.chartUncertainty({ from: julday(1992, 9, 6, 4), to: julday(1992, 9, 7, 4), ...CHARLOTTE });

  it("can still say the Moon's sign for certain when she stays in one all day", () => {
    const moon = u.placements.find((p) => p.key === "body:1")!;
    expect(moon.signs[0].probability).toBeGreaterThanOrEqual(CERTAIN);
  });

  it("gives the Ascendant every sign, with the long-ascension signs holding longest", () => {
    const asc = u.placements.find((p) => p.key === "asc")!;
    expect(asc.signs.length).toBe(12);
    // At 35° N, Virgo and Libra rise slowly, Pisces and Aries fast.
    const share = (s: number) => asc.signs.find((x) => x.value === s)!.probability;
    expect(share(5)).toBeGreaterThan(share(11));
    expect(share(6)).toBeGreaterThan(share(0));
  });
});

describe("the clock's own uncertainty", () => {
  it("is reported, and grows far from the present", () => {
    const near = eph.chartUncertainty({ from: julday(2000, 1, 1), to: julday(2000, 1, 1, 0.1), ...CHARLOTTE });
    const far = eph.chartUncertainty({ from: julday(2350, 1, 1), to: julday(2350, 1, 1, 0.1), ...CHARLOTTE });
    expect(near.deltaT.sigma).toBeLessThan(0.01);
    expect(far.deltaT.sigma).toBeGreaterThan(100);
  });
});

describe("rectification", () => {
  // Events built on a known birth minute's own angles: the dates transiting
  // Saturn, Uranus, Jupiter and Neptune reached them.
  const h = eph.houses(BIRTH, CHARLOTTE.lat, CHARLOTTE.lon, "P");
  const events: Array<{ jd: number }> = [];
  for (const [b, pts] of [[6, [h.ascendant, h.mc, h.ascendant + 180, h.mc + 180]], [7, [h.ascendant, h.mc, h.mc + 180]], [5, [h.ascendant, h.mc]], [8, [h.ascendant, h.mc, h.ascendant + 90]]] as const) {
    for (const x of eph.crossings(b, [...pts], BIRTH + 365 * 18, BIRTH + 365 * 60).slice(0, 3)) events.push({ jd: x.jd });
  }

  it("finds the minute the events were built on, from a ±2-hour window", () => {
    const r = eph.rectify({ from: BIRTH - 2 / 24, to: BIRTH + 2 / 24, ...CHARLOTTE, events: events.slice(0, 9) });
    expect(Math.abs(r.best[0].jd - BIRTH) * 1440).toBeLessThanOrEqual(2);
    expect((r.window68.to - r.window68.from) * 1440).toBeLessThanOrEqual(10);
    expect(r.candidates.reduce((s, c) => s + c.support, 0)).toBeCloseTo(1, 9);
  });

  it("shows its working: every contact it counted, with its orb", () => {
    const r = eph.rectify({ from: BIRTH - 2 / 24, to: BIRTH + 2 / 24, ...CHARLOTTE, events: events.slice(0, 9) });
    const ev = r.best[0].evidence;
    expect(ev.length).toBeGreaterThan(0);
    expect(ev.some((e) => e.technique === "transit" && e.orb < 0.1)).toBe(true);
    for (const e of ev) expect(e.contribution).toBeGreaterThan(0);
  });

  it("with no events, favours no minute", () => {
    const r = eph.rectify({ from: BIRTH - 1 / 24, to: BIRTH + 1 / 24, ...CHARLOTTE, events: [] });
    const s = r.candidates.map((c) => c.support);
    expect(Math.max(...s) - Math.min(...s)).toBeLessThan(1e-12);
  });
});
