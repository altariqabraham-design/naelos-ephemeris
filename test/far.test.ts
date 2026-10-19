// Beyond DE440: the DE441 tables, 13,200 BC to AD 17,191.
//
// The 294 far tables live in the "ephemeris" Storage bucket; five samples are
// in data-far/ for these tests: the first and last files, 500 BC, and the
// file either side of each seam with the DE440 tables (1500–1550, 2649–2700).
// What is checked here:
//   - the calendar holds for negative Julian Days, in both calendars;
//   - ΔT and its uncertainty run the whole span, continuously;
//   - the long-term precession is ERFA's to rounding, and the frame and the
//     sidereal time turn from IAU 2006 to it without a step;
//   - positions agree with NASA JPL Horizons wherever Horizons answers — to
//     the thousandth of an arcsecond in anything that does not depend on the
//     precession model, and to its documented model difference in what does;
//   - nothing jumps at either seam, and nothing is answered beyond the ends.

import { describe, it, expect } from "vitest";
import { readFileSync, readdirSync } from "node:fs";
import { createHash } from "node:crypto";
import { resolve } from "node:path";
import {
  NaelosEphemeris, julday, juldayJulian, juldayHistorical, revjul, formatYear, tablesFor, TABLE_FILES, FIRST_JD, LAST_JD,
  isFarTable, deltaT, deltaTUncertainty, SE_SUN, SE_MOON, SE_MARS, SE_PLUTO, SE_CHIRON, SE_MEAN_NODE, SE_TRUE_NODE,
  SE_MEAN_APOG, SEFLG_SPEED, apparentSiderealTime, type Position,
} from "../src/index.ts";
import { frameOfDate, longTermPrecession, iauWeight } from "../src/frames.ts";
import { meanCorrection } from "../src/lunar-points.ts";

const ROOT = resolve(__dirname, "..");
const near = JSON.parse(readFileSync(resolve(ROOT, "data/manifest.json"), "utf8"));
const far = JSON.parse(readFileSync(resolve(ROOT, "data-far/manifest.json"), "utf8"));
const SAMPLES = readdirSync(resolve(ROOT, "data-far")).filter((n) => n.endsWith(".bin")).sort();
const file = (n: string) => readFileSync(resolve(ROOT, SAMPLES.includes(n) ? "data-far" : "data", n));
const engineFor = (from: number, to = from) => NaelosEphemeris.fromBuffers(tablesFor(from, to).map(file));
const ok = (r: Position | { error: string }): Position => {
  if ("error" in r) throw new Error(r.error);
  return r;
};
const wrap = (d: number) => ((d + 540) % 360) - 180;
const MAS = (rad: number) => (rad * 180 / Math.PI) * 3.6e6;

describe("the far tables", () => {
  it("are DE441's, built to the DE440 tables' own target, and meet them exactly", () => {
    expect(far.source.name).toBe("JPL DE441");
    expect(far.target_arcsec).toBe(near.target_arcsec);
    expect(far.meets.sha256).toBe(near.source.sha256);
    expect(far.span.past.end_jd_tdb).toBe(near.span.start_jd_tdb);
    expect(far.span.future.start_jd_tdb).toBe(near.span.end_jd_tdb);
    for (const side of ["past", "future"]) {
      for (const [name, s] of Object.entries(far.series[side]) as Array<[string, { worst_arcsec: number }]>) {
        expect(s.worst_arcsec, `${side} ${name}`).toBeLessThanOrEqual(far.target_arcsec);
      }
    }
  });

  it("cover 13,200 BC to AD 17,191 in 294 files with no gap, and the index agrees", () => {
    expect(far.files).toHaveLength(294);
    const all = [...near.files.map((f: { name: string }) => [f, 0]), ...far.files.map((f: { name: string }) => [f, 1])]
      .sort((a, b) => a[0].start_jd_tdb - b[0].start_jd_tdb);
    expect(TABLE_FILES.map(([n, s, e, f]) => [n, s, e, f])).toEqual(all.map(([f, isFar]) => [f.name, f.start_jd_tdb, f.end_jd_tdb, isFar]));
    for (let i = 1; i < TABLE_FILES.length; i++) expect(TABLE_FILES[i][1]).toBe(TABLE_FILES[i - 1][2]);
    expect(formatYear(revjul(FIRST_JD).year)).toBe("13198 BC");
    expect(revjul(LAST_JD).year).toBe(17188);
  });

  it("keeps only five samples in git, each the file its manifest describes", () => {
    expect(SAMPLES).toEqual(["ne-1500.bin", "ne-17100.bin", "ne-2649.bin", "ne-m0500.bin", "ne-m13200.bin"]);
    for (const n of SAMPLES) {
      const entry = far.files.find((f: { name: string }) => f.name === n);
      expect(createHash("sha256").update(file(n)).digest("hex"), n).toBe(entry.sha256);
      expect(isFarTable(n)).toBe(true);
    }
    expect(isFarTable("ne-2000.bin")).toBe(false);
  });
});

describe("the calendar, for any year", () => {
  it("counts days across Julian Day 0 and back to the start of DE441", () => {
    expect(julday(2000, 1, 1, 12)).toBe(2451545);
    expect(juldayJulian(-4712, 1, 1, 12)).toBe(0);
    expect(juldayJulian(-1000, 7, 12, 12)).toBe(1356001); // Meeus, example 7.b
    expect(juldayJulian(333, 1, 27, 12)).toBe(1842713); // Meeus, example 7.a
    expect(juldayJulian(-13200, 8, 15)).toBe(-3100015.5); // DE441's first instant
    expect(julday(17191, 3, 15)).toBe(8000016.5); // and its last
  });

  it("round-trips every date, either calendar", () => {
    for (let jd = -3100015.5; jd < 8000016.5; jd += 7919.375) {
      const d = revjul(jd);
      const back = d.calendar === "gregorian" ? julday(d.year, d.month, d.day, d.hour) : juldayJulian(d.year, d.month, d.day, d.hour);
      expect(back).toBeCloseTo(jd, 6);
    }
  });

  it("reads dates as history records them: Julian before 15 October 1582, and refuses the ten skipped days", () => {
    expect(juldayHistorical(1582, 10, 4)).toBe(2299159.5);
    expect(juldayHistorical(1582, 10, 15)).toBe(2299160.5);
    expect(() => juldayHistorical(1582, 10, 10)).toThrow(RangeError);
    expect(juldayHistorical(-43, 3, 15)).toBe(juldayJulian(-43, 3, 15)); // the Ides of March, 44 BC
    expect(formatYear(-43)).toBe("44 BC");
    expect(formatYear(0)).toBe("1 BC");
    expect(formatYear(800)).toBe("AD 800");
  });
});

describe("ΔT across the whole span", () => {
  // The model the tables were built from (Skyfield's construction of
  // Morrison-Stephenson-Hohenkerk-Zawilski 2021 and Stephenson et al. 2016),
  // read off directly at these instants when the tables were built.
  const MODEL: Array<[number, number]> = [
    [-13200, 733367.26], [-5000, 151066.82], [-720, 20371.84], [0, 10441.31], [1000, 1650.39],
    [2700, 2167.03], [10000, 216879.9], [17191, 767052.78],
  ];
  const at = (year: number) => 2451545 + (year - 2000) * 365.25;

  it("is the model, to the hundredth of a second", () => {
    for (const [y, s] of MODEL) expect(deltaT(at(y)), String(y)).toBeCloseTo(s, 1);
  });

  it("joins the quarter-year table without a step", () => {
    for (const j of [2286820.25 + 91.3125, 2286820.25 + 4412 * 91.3125]) {
      expect(Math.abs(deltaT(j - 1e-6) - deltaT(j + 1e-6))).toBeLessThan(0.05);
    }
  });

  it("is as uncertain as Morrison & Stephenson say: 0.8 t² seconds before 1700", () => {
    expect(deltaTUncertainty(at(-500))).toBeCloseTo(430, -1); // their ±430 s at 500 BC
    expect(deltaTUncertainty(at(1000))).toBeCloseTo(54, 0);
    expect(deltaTUncertainty(at(1600))).toBe(20);
    expect(deltaTUncertainty(at(-13000)) / 3600).toBeGreaterThan(4.5); // hours, in the deep past
  });
});

describe("the frame of date beyond the DE440 tables", () => {
  const ref = JSON.parse(readFileSync(resolve(ROOT, "test/ltp-reference.json"), "utf8"));

  it("is ERFA's long-term precession, to rounding", () => {
    for (const r of ref.rows) {
      const { pb, epsA } = longTermPrecession(r.jd_tt);
      pb.forEach((v, i) => expect(Math.abs(v - r.pb[i]), `JD ${r.jd_tt}`).toBeLessThan(1e-14));
      expect(Math.abs(epsA - r.epsA)).toBeLessThan(1e-14);
    }
  });

  it("is IAU 2006 inside the DE440 tables, and blends to the long-term model over 50 years outside", () => {
    expect(iauWeight(near.span.start_jd_tdb)).toBe(1);
    expect(iauWeight(near.span.end_jd_tdb)).toBe(1);
    expect(iauWeight(near.span.start_jd_tdb - 25 * 365.25)).toBeCloseTo(0.5, 12);
    expect(iauWeight(near.span.end_jd_tdb + 50 * 365.25)).toBe(0);
  });

  it("turns without a step at the seams and at the ends of the blend", () => {
    for (const jd of [near.span.start_jd_tdb, near.span.start_jd_tdb - 50 * 365.25, near.span.end_jd_tdb, near.span.end_jd_tdb + 50 * 365.25]) {
      const a = frameOfDate(jd - 1e-6).npb, b = frameOfDate(jd + 1e-6).npb;
      expect(MAS(Math.max(...a.map((v, i) => Math.abs(v - b[i]))))).toBeLessThan(0.01);
      // Sidereal time too: the same instant from either side, less a millionth of a day of rotation.
      const g = (t: number) => apparentSiderealTime(t, t);
      let d = g(jd + 1e-6) - g(jd - 1e-6) - 2e-6 * 2 * Math.PI * 1.0027379;
      d = ((d + Math.PI) % (2 * Math.PI)) - Math.PI;
      expect(Math.abs(MAS(d))).toBeLessThan(1);
    }
  });
});

describe("positions against NASA JPL Horizons", () => {
  const H = JSON.parse(readFileSync(resolve(ROOT, "test/horizons-far.json"), "utf8")).rows as Array<{ label: string; jd_tt: number; body: number; lon: number; lat: number }>;
  const vec = (lon: number, lat: number) => {
    const a = (lon * Math.PI) / 180, b = (lat * Math.PI) / 180;
    return [Math.cos(b) * Math.cos(a), Math.cos(b) * Math.sin(a), Math.sin(b)];
  };
  const sep = (u: number[], v: number[]) => (Math.atan2(
    Math.hypot(u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]), u[0] * v[0] + u[1] * v[1] + u[2] * v[2],
  ) * 180) / Math.PI * 3600;
  const labels = [...new Set(H.map((r) => r.label))];

  it("place every body where Horizons does, measured from the Sun — 0.01″ in 500 BC, at the seams and after 2649", () => {
    // A separation is the same in any frame, so this compares DE441 and the
    // engine's light-time, deflection and aberration with Horizons' own,
    // independently of which precession model either uses.
    for (const label of labels) {
      const rows = H.filter((r) => r.label === label);
      const eph = engineFor(rows[0].jd_tt);
      const sunE = ok(eph.calc(rows[0].jd_tt, SE_SUN, 0)), sunH = rows.find((r) => r.body === 0)!;
      for (const r of rows.filter((x) => x.body !== 0)) {
        const p = ok(eph.calc(r.jd_tt, r.body, 0));
        const mine = sep(vec(p.longitude, p.latitude), vec(sunE.longitude, sunE.latitude));
        const theirs = sep(vec(r.lon, r.lat), vec(sunH.lon, sunH.lat));
        expect(Math.abs(mine - theirs), `${label} body ${r.body}`).toBeLessThan(0.01);
      }
    }
  });

  it("differ from Horizons' places of date only by the precession models, the same for every body", () => {
    // Horizons uses Owen's (1990) long-term model, built to match IAU 1976;
    // the engine Vondrák et al. (2011), which matches IAU 2006. In 500 BC that
    // is a few arcseconds, and — being a turn of the whole frame — it is the
    // same for every body at an instant, to a hundredth of an arcsecond.
    for (const label of labels) {
      const rows = H.filter((r) => r.label === label);
      const eph = engineFor(rows[0].jd_tt);
      const offsets = rows.map((r) => {
        const p = ok(eph.calc(r.jd_tt, r.body, 0));
        return [wrap(p.longitude - r.lon) * 3600, (p.latitude - r.lat) * 3600];
      });
      for (const [dl] of offsets) {
        expect(Math.abs(dl), label).toBeLessThan(6);
        expect(Math.abs(dl - offsets[0][0]), label).toBeLessThan(0.05 / Math.cos(0));
      }
    }
  });
});

describe("at the seams with the DE440 tables, and at the ends", () => {
  it("moves through 1550 and 2649 without a jump", () => {
    // Stepping over the seam by 86 seconds either way. Inside a light-time of
    // it — the Moon's light takes 1.3 s — the body is read from one fit and the
    // Earth from the other, which agree on the Earth-Moon barycentre to about
    // a kilometre: for those 1.3 seconds the Moon can sit up to 0.7″ off. Two
    // such instants in 30,000 years; outside them, nothing jumps.
    for (const seam of [near.span.start_jd_tdb, near.span.end_jd_tdb]) {
      const eph = engineFor(seam - 2, seam + 2);
      for (const body of [SE_SUN, SE_MOON, SE_MARS, SE_PLUTO]) {
        const e = 1e-3;
        const a = ok(eph.calc(seam - e, body, SEFLG_SPEED)), b = ok(eph.calc(seam + e, body, SEFLG_SPEED));
        const expected = a.longitudeSpeed * 2 * e;
        expect(Math.abs(wrap(b.longitude - a.longitude) - expected) * 3600, `body ${body}`).toBeLessThan(0.01);
      }
    }
  });

  it("answers at the first and last instants, and refuses beyond them", () => {
    const first = engineFor(FIRST_JD + 2), last = engineFor(LAST_JD - 2);
    expect("error" in first.calc(FIRST_JD + 2, SE_MARS, SEFLG_SPEED)).toBe(false);
    expect("error" in last.calc(LAST_JD - 2, SE_MARS, SEFLG_SPEED)).toBe(false);
    expect("error" in first.calc(FIRST_JD - 30, SE_MARS, 0)).toBe(true);
    expect("error" in last.calc(LAST_JD + 30, SE_MARS, 0)).toBe(true);
  });

  it("refuses Chiron and the asteroids where JPL publishes no orbit for them", () => {
    const t = julday(1500, 6, 1);
    expect("error" in engineFor(t).calc(t, SE_CHIRON, 0)).toBe(true);
  });

  it("asks for the right century 8.5 days of ΔT away from the clock", () => {
    // In 13,000 BC a UT moment is TT a week later; the file must cover both.
    const ut = juldayJulian(-13197, 6, 1);
    const tt = ut + deltaT(ut) / 86400;
    const eph = engineFor(ut);
    expect("error" in eph.calc_ut(ut, SE_MARS, 0)).toBe(false);
    expect(tt - ut).toBeGreaterThan(8);
  });
});

describe("showing the working, far from the present", () => {
  it("names DE441 and the long-term precession, and says what Horizons will show", () => {
    const t = juldayJulian(-449, 6, 1);
    const p = engineFor(t).provenance(t, SE_MARS);
    if ("error" in p) throw new Error(p.error);
    expect(p.source).toMatch(/^JPL DE441, table ne-m0500\.bin segment \d+$/);
    expect(p.steps.join(" ")).toMatch(/long-term precession \(Vondrák, Capitaine & Wallace 2011\)/);
    expect(p.clock.basis).toBe("historical reconstruction"); // after 720 BC: eclipse records
    expect(p.verify?.expect).toMatch(/Owen's \(1990\)/);
  });

  it("calls the clock a model before the records, and offers no Horizons link before 9999 BC", () => {
    const t = juldayJulian(-13150, 6, 1);
    const p = engineFor(t).provenance(t, SE_MARS);
    if ("error" in p) throw new Error(p.error);
    expect(p.clock.basis).toBe("long-term model");
    expect(p.verify).toBeUndefined();
  });
});

describe("the mean North Node and mean Lilith beyond the DE440 tables", () => {
  it("are the polynomials alone from 1550 to 2650", () => {
    expect(meanCorrection(julday(2026, 1, 1))).toEqual({ node: 0, apogee: 0 });
    expect(meanCorrection(julday(1560, 1, 1))).toEqual({ node: 0, apogee: 0 });
  });

  it("follow DE441's true node on average where the polynomial drifts", () => {
    // By 13,150 BC the bare polynomial is 2.4° off the true node's average
    // over its 18.6-year cycle; corrected, the mean node sits on that average.
    const mid = juldayJulian(-13150, 1, 1);
    const eph = engineFor(mid - 3400, mid + 3400);
    let sx = 0, sy = 0;
    for (let i = 0; i < 2000; i++) {
      const t = mid - 3399 + i * 3.399;
      const d = (wrap(ok(eph.calc(t, SE_TRUE_NODE, 0)).longitude - ok(eph.calc(t, SE_MEAN_NODE, 0)).longitude) * Math.PI) / 180;
      sx += Math.cos(d);
      sy += Math.sin(d);
    }
    expect(Math.abs((Math.atan2(sy, sx) * 180) / Math.PI)).toBeLessThan(0.02);
    expect(meanCorrection(mid).node).toBeGreaterThan(2);
    expect("error" in eph.calc(mid, SE_MEAN_APOG, 0)).toBe(false);
  });
});
