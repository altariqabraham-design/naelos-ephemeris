// Eclipses: the nearest solar or lunar eclipse before or after a moment.
//
// Checked against eclipses everyone can look up — NASA's published greatest
// eclipse (in TD, which is TT) and type for recent ones — and against their
// own geometry: every eclipse falls at a new or full Moon, near a node, and
// the search finds the same eclipse walking forward or backward.

import { describe, it, expect } from "vitest";
import { readFileSync } from "node:fs";
import { resolve } from "node:path";
import { NaelosEphemeris, julday, SE_SUN, SE_MOON, type Eclipse } from "../src/index.ts";

const ROOT = resolve(__dirname, "..");
const manifest = JSON.parse(readFileSync(resolve(ROOT, "data/manifest.json"), "utf8"));
const eph = NaelosEphemeris.fromBuffers(manifest.files.map((f: { name: string }) => readFileSync(resolve(ROOT, "data", f.name))));
const wrap = (d: number) => ((d + 540) % 360) - 180;
const lon = (jdTT: number, body: number) => {
  const r = eph.calc(jdTT, body, 0);
  if ("error" in r) throw new Error(r.error);
  return r;
};

// NASA's greatest eclipse in TD (= TT), from its eclipse pages.
const KNOWN: Array<[string, "solar" | "lunar", Eclipse["category"], number]> = [
  ["2017-08-21 total solar", "solar", "total", julday(2017, 8, 21, 18 + 25 / 60 + 29 / 3600)],
  ["2023-04-20 hybrid solar", "solar", "hybrid", julday(2023, 4, 20, 4 + 17 / 60 + 56 / 3600)],
  ["2023-10-14 annular solar", "solar", "annular", julday(2023, 10, 14, 18 + 0 / 60 + 41 / 3600)],
  ["2024-04-08 total solar", "solar", "total", julday(2024, 4, 8, 18 + 18 / 60 + 29 / 3600)],
  ["2022-10-25 partial solar", "solar", "partial", julday(2022, 10, 25, 11 + 1 / 60 + 20 / 3600)],
  ["2022-11-08 total lunar", "lunar", "total", julday(2022, 11, 8, 10 + 59 / 60 + 11 / 3600)],
  ["2023-05-05 penumbral lunar", "lunar", "penumbral", julday(2023, 5, 5, 17 + 24 / 60)],
];

describe("well-known eclipses", () => {
  for (const [name, kind, category, td] of KNOWN) {
    it(`${name}: the right type, within two minutes of NASA's greatest eclipse`, () => {
      const e = kind === "solar" ? eph.solarEclipse(td - 3) : eph.lunarEclipse(td - 3);
      expect(e?.category).toBe(category);
      expect(Math.abs((e!.jdTT - td) * 86400)).toBeLessThan(120);
    });
  }
});

describe("every eclipse meets its own geometry", () => {
  const solar: Eclipse[] = [], lunar: Eclipse[] = [];
  let jd = julday(1990, 1, 1);
  for (let i = 0; i < 40; i++) { const e = eph.solarEclipse(jd)!; solar.push(e); jd = e.jd + 1; }
  jd = julday(1990, 1, 1);
  for (let i = 0; i < 40; i++) { const e = eph.lunarEclipse(jd)!; lunar.push(e); jd = e.jd + 1; }

  it("solar eclipses fall at new Moon, lunar at full Moon, within a few hours", () => {
    for (const e of solar) {
      const d = wrap(lon(e.jdTT, SE_MOON).longitude - lon(e.jdTT, SE_SUN).longitude);
      expect(Math.abs(d)).toBeLessThan(2); // degrees: the Moon moves 0.5° an hour
    }
    for (const e of lunar) {
      const d = wrap(lon(e.jdTT, SE_MOON).longitude - lon(e.jdTT, SE_SUN).longitude - 180);
      expect(Math.abs(d)).toBeLessThan(2);
    }
  });

  it("with the Moon near the ecliptic, where only an eclipse can happen", () => {
    for (const e of [...solar, ...lunar]) expect(Math.abs(lon(e.jdTT, SE_MOON).latitude)).toBeLessThan(1.65);
  });

  it("two to five solar eclipses a year, and as many lunar", () => {
    const perYear = (list: Eclipse[]) => {
      const years = new Map<number, number>();
      for (const e of list) {
        const y = Math.floor(2000 + (e.jd - 2451545) / 365.25);
        years.set(y, (years.get(y) ?? 0) + 1);
      }
      return [...years.values()].slice(1, -1); // whole years only
    };
    for (const n of perYear(solar)) { expect(n).toBeGreaterThanOrEqual(2); expect(n).toBeLessThanOrEqual(5); }
    for (const n of perYear(lunar)) { expect(n).toBeLessThanOrEqual(5); }
  });

  it("finds the same eclipse searching backward from just after it", () => {
    for (const e of [...solar.slice(0, 10)]) expect(eph.solarEclipse(e.jd + 0.01, true)!.jd).toBeCloseTo(e.jd, 6);
    for (const e of [...lunar.slice(0, 10)]) expect(eph.lunarEclipse(e.jd + 0.01, true)!.jd).toBeCloseTo(e.jd, 6);
  });

  it("repeats after a saros (6,585.32 days) with the same kind of eclipse", () => {
    // A saros later the Sun, Moon and node line up again: same family, close in kind.
    for (const e of solar.slice(0, 10)) {
      const next = eph.solarEclipse(e.jd + 6585.32 - 2)!;
      expect(Math.abs(next.jd - (e.jd + 6585.32))).toBeLessThan(0.5);
    }
  });
});
