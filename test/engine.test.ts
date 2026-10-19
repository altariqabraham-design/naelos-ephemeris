// Naelos Ephemeris against JPL: DE440 reduced by Skyfield (tools/build_reference.py)
// at 300 instants across 1550–2649. Every body, three kinds of position.
//
// The tables are fitted to 0.01″; the engine's nutation (IAU 2000B's terms)
// and Skyfield's (IAU 2000A) differ by up to 0.0026″ across the span. The
// bounds below leave room for both and nothing more.

import { describe, it, expect } from "vitest";
import { readFileSync } from "node:fs";
import { gunzipSync } from "node:zlib";
import { createHash } from "node:crypto";
import { resolve } from "node:path";
import {
  NaelosEphemeris, julday, tablesFor, TABLE_FILES, SEFLG_SPEED, SEFLG_EQUATORIAL, SEFLG_HELCTR, SE_EARTH, SE_MARS,
  SE_PLUTO, SE_CHIRON, SE_ECL_NUT, type Position,
} from "../src/index.ts";

const ROOT = resolve(__dirname, "..");
const manifest = JSON.parse(readFileSync(resolve(ROOT, "data/manifest.json"), "utf8"));
const buffers = manifest.files.map((f: { name: string }) => readFileSync(resolve(ROOT, "data", f.name)));
const eph = NaelosEphemeris.fromBuffers(buffers);
const ref = JSON.parse(gunzipSync(readFileSync(resolve(ROOT, "test/reference.json.gz"))).toString());
const jds: number[] = ref.jd_tt;

const wrap = (d: number) => ((d + 540) % 360) - 180;
const arcsec = (lon: number, lat: number, rlon: number, rlat: number) =>
  Math.hypot(wrap(lon - rlon) * Math.cos((rlat * Math.PI) / 180), lat - rlat) * 3600;
const ok = (r: ReturnType<typeof eph.calc>): Position => {
  if ("error" in r) throw new Error(r.error);
  return r;
};

const NAMES = ["Sun", "Moon", "Mercury", "Venus", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"];

describe("the tables", () => {
  it("match the checksums the build recorded", () => {
    manifest.files.forEach((f: { name: string; sha256: string }, i: number) => {
      expect(createHash("sha256").update(buffers[i]).digest("hex"), f.name).toBe(f.sha256);
    });
  });

  it("were built from DE440 to 0.01″, and say so", () => {
    expect(manifest.source.name).toBe("JPL DE440");
    expect(manifest.target_arcsec).toBe(0.01);
    for (const s of Object.values(manifest.series) as Array<{ worst_arcsec: number }>) {
      expect(s.worst_arcsec).toBeLessThanOrEqual(0.01);
    }
  });

  it("cover 1550 to 2649 and refuse what lies outside", () => {
    expect(eph.range.start).toBeLessThan(julday(1550, 2, 1));
    expect(eph.range.end).toBeGreaterThan(julday(2649, 1, 1));
    expect("error" in eph.calc(julday(1500, 1, 1), 0, 0)).toBe(true);
    expect("error" in eph.calc(julday(2700, 1, 1), 0, 0)).toBe(true);
  });
});

describe("geocentric apparent places agree with JPL", () => {
  for (let b = 0; b <= 9; b++) {
    it(`${NAMES[b]}: within 0.02″, speed within ${b === 1 ? "0.001" : "0.00001"}°/day`, () => {
      let worst = 0, worstSpeed = 0;
      jds.forEach((jd, i) => {
        const r = ok(eph.calc(jd, b, SEFLG_SPEED));
        worst = Math.max(worst, arcsec(r.longitude, r.latitude, ref.geo[b].lon[i], ref.geo[b].lat[i]));
        worstSpeed = Math.max(worstSpeed, Math.abs(r.longitudeSpeed - ref.geo[b].speed[i]));
      });
      expect(worst).toBeLessThan(0.02);
      // The Moon's speed carries the derivative of a degree-21 fit; 1e-3 °/day
      // is 3.6″ a day, invisible in any chart.
      expect(worstSpeed).toBeLessThan(b === 1 ? 1e-3 : 1e-5);
    });
  }
});

describe("equatorial places agree with JPL", () => {
  for (let b = 0; b <= 9; b++) {
    it(`${NAMES[b]}: RA/Dec within 0.02″`, () => {
      let worst = 0;
      jds.forEach((jd, i) => {
        const r = ok(eph.calc(jd, b, SEFLG_EQUATORIAL));
        worst = Math.max(worst, arcsec(r.longitude, r.latitude, ref.equ[b].ra[i], ref.equ[b].dec[i]));
      });
      expect(worst).toBeLessThan(0.02);
    });
  }
});

describe("heliocentric places agree with JPL", () => {
  for (const b of [1, 2, 3, 4, 5, 6, 7, 8, 9, SE_EARTH]) {
    // Mercury is the hard case: 0.31 AU from the Sun at perihelion, nearer than
    // the Earth ever comes, so a table error looks larger from there.
    const bound = b === 2 ? 0.04 : 0.02;
    it(`${b === SE_EARTH ? "Earth" : NAMES[b]}: within ${bound}″`, () => {
      let worst = 0;
      jds.forEach((jd, i) => {
        const r = ok(eph.calc(jd, b, SEFLG_HELCTR));
        worst = Math.max(worst, arcsec(r.longitude, r.latitude, ref.helio[b].lon[i], ref.helio[b].lat[i]));
      });
      expect(worst).toBeLessThan(bound);
    });
  }

  it("has no heliocentric Sun and no geocentric Earth", () => {
    expect("error" in eph.calc(2451545, 0, SEFLG_HELCTR)).toBe(true);
    expect("error" in eph.calc(2451545, SE_EARTH, 0)).toBe(true);
  });
});

describe("behind the Sun", () => {
  it("keeps a planet's position smooth and moving forward while it passes behind the solar disc", () => {
    // 2280-06-02: Uranus a minute of arc from the Sun's centre. The untapered
    // deflection formula threw it 9″ and turned its speed backward here.
    const T = 2553967.0;
    let prev = -Infinity;
    for (let dt = -0.5; dt <= 0.5; dt += 0.01) {
      const r = eph.calc(T + dt, 7, SEFLG_SPEED);
      if ("error" in r) throw new Error(r.error);
      expect(r.longitude).toBeGreaterThan(prev);
      expect(r.longitudeSpeed).toBeGreaterThan(0.02);
      expect(r.longitudeSpeed).toBeLessThan(0.1);
      prev = r.longitude;
    }
  });
});

describe("loading centuries as they are needed", () => {
  const file = (name: string) => readFileSync(resolve(ROOT, "data", name));

  it("names the files a span needs, and the century before when light left the body in it", () => {
    // The DE440 files, as the index lists them; the far (DE441) files are test/far.test.ts's.
    const near = TABLE_FILES.filter(([, , , far]) => far === 0).map(([n]) => n);
    expect(near).toEqual(manifest.files.map((f: { name: string }) => f.name));
    expect(tablesFor(julday(2026, 10, 18), julday(2026, 10, 18))).toEqual(["ne-2000.bin"]);
    expect(tablesFor(julday(2000, 1, 1, 6), julday(2000, 1, 1, 6))).toEqual(["ne-1900.bin", "ne-2000.bin"]);
    expect(tablesFor(julday(1720, 1, 1), julday(1560, 1, 1))).toEqual(["ne-1550.bin", "ne-1600.bin", "ne-1700.bin"]);
    // Before 1550 the DE441 tables take over.
    expect(tablesFor(julday(1400, 6, 1), julday(1400, 6, 1))).toEqual(["ne-1400.bin"]);
  });

  it("answers from one century, refuses another until it is added, and ignores a repeat", () => {
    const one = NaelosEphemeris.fromBuffers([file("ne-2000.bin")]);
    expect("error" in one.calc(julday(2026, 10, 18), SE_MARS, 0)).toBe(false);
    expect("error" in one.calc(julday(1950, 1, 1), SE_MARS, 0)).toBe(true);
    one.addTables([file("ne-1900.bin"), file("ne-2000.bin")]);
    expect(one.calc(julday(1950, 1, 1), SE_MARS, SEFLG_SPEED)).toEqual(eph.calc(julday(1950, 1, 1), SE_MARS, SEFLG_SPEED));
    expect(one.range).toEqual({ start: julday(1900, 1, 1), end: julday(2100, 1, 1) });
  });

  it("reaches back past the start of its century for light-time", () => {
    // Pluto's light arriving in the first hour of 2000 left it in 1999. The
    // 2000 file carries the segment that straddles the boundary.
    const one = NaelosEphemeris.fromBuffers([file("ne-2000.bin")]);
    const t = julday(2000, 1, 1, 0.5);
    expect(one.calc(t, SE_PLUTO, SEFLG_SPEED)).toEqual(eph.calc(t, SE_PLUTO, SEFLG_SPEED));
  });

  it("asks for the century before where a series starts exactly on the boundary", () => {
    // The asteroids' segment grids happen to meet 2200-01-01 exactly, so the
    // 2200 file holds nothing of them before it — and Chiron's light takes
    // hours to arrive.
    const t = julday(2200, 1, 1, 0.5);
    expect("error" in NaelosEphemeris.fromBuffers([file("ne-2200.bin")]).calc(t, SE_CHIRON, 0)).toBe(true);
    const enough = NaelosEphemeris.fromBuffers(tablesFor(t, t).map(file));
    expect(enough.calc(t, SE_CHIRON, SEFLG_SPEED)).toEqual(eph.calc(t, SE_CHIRON, SEFLG_SPEED));
  });
});

describe("the obliquity and nutation of date (SE_ECL_NUT)", () => {
  const D = Math.PI / 180;
  it("is laid out as Swiss Ephemeris lays it out, and the mean obliquity is IAU 2006's", () => {
    const n = ok(eph.calc(2451545, SE_ECL_NUT, 0));
    expect(n.latitude * 3600).toBeCloseTo(84381.406, 6); // mean obliquity at J2000
    expect(n.longitude - n.latitude).toBeCloseTo(n.longitudeSpeed, 12); // true = mean + Δε
    expect(Math.abs(n.distance * 3600)).toBeLessThan(20); // Δψ is under 20″
  });

  it("turns an ecliptic place into right ascension and declination exactly", () => {
    for (const jd of [julday(1600, 3, 1), julday(2026, 10, 18, 12), julday(2600, 7, 1)]) {
      const e = ok(eph.calc(jd, SE_ECL_NUT, 0)).longitude * D;
      const ecl = ok(eph.calc(jd, SE_MARS, 0)), equ = ok(eph.calc(jd, SE_MARS, SEFLG_EQUATORIAL));
      const l = ecl.longitude * D, b = ecl.latitude * D;
      const dec = Math.asin(Math.sin(b) * Math.cos(e) + Math.cos(b) * Math.sin(e) * Math.sin(l)) / D;
      const ra = Math.atan2(Math.sin(l) * Math.cos(e) - Math.tan(b) * Math.sin(e), Math.cos(l)) / D;
      expect(Math.abs(dec - equ.latitude) * 3600).toBeLessThan(1e-6);
      expect(Math.abs(wrap(ra - equ.longitude)) * 3600).toBeLessThan(1e-6);
    }
  });
});

describe("dates", () => {
  it("numbers days as the Gregorian calendar does", () => {
    expect(julday(2000, 1, 1, 12)).toBe(2451545);
    expect(julday(1800, 1, 1)).toBe(2378496.5);
    expect(julday(2400, 1, 1)).toBe(2597641.5);
    expect(julday(1992, 9, 6, 10 + 41 / 60)).toBeCloseTo(2448871.945138889, 9);
  });
});
