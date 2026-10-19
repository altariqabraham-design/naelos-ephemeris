// ΔT, sidereal time and the house systems.
//
// Sidereal time is checked against ERFA (IAU 2006/2000A). The house systems
// have no external oracle to check against, so each is checked against its
// own DEFINITION: a Placidus cusp has covered the stated fraction of its own
// semi-arc; a Campanus cusp lies on the circle through the horizon's poles and
// the stated point of the prime vertical; and so on. A formula with a wrong
// quadrant or a swapped house fails these, whatever another program says.

import { describe, it, expect } from "vitest";
import { readFileSync } from "node:fs";
import { gunzipSync } from "node:zlib";
import { resolve } from "node:path";
import {
  deltaT, deltaTUncertainty, MEASURED_UNTIL, PREDICTED_UNTIL, julday, apparentSiderealTime, norm360,
} from "../src/index.ts";
import { houses, type HouseSystem } from "../src/houses.ts";
import { frameOfDate } from "../src/frames.ts";

const ROOT = resolve(__dirname, "..");
const ref = JSON.parse(gunzipSync(readFileSync(resolve(ROOT, "test/reference.json.gz"))).toString());
const D2R = Math.PI / 180, R2D = 180 / Math.PI;
const wrap = (d: number) => ((d + 540) % 360) - 180;

describe("ΔT", () => {
  it("is the measured value where the Earth's rotation has been measured", () => {
    expect(deltaT(julday(2000, 1, 1, 12))).toBeCloseTo(63.83, 1);
    expect(deltaT(julday(2026, 10, 1))).toBeCloseTo(69.2, 1);
  });

  it("follows the historical reconstruction before 1973", () => {
    expect(deltaT(julday(1900, 1, 1))).toBeGreaterThan(-3);
    expect(deltaT(julday(1900, 1, 1))).toBeLessThan(-1);
    expect(deltaT(julday(1800, 1, 1))).toBeGreaterThan(13);
    expect(deltaT(julday(1800, 1, 1))).toBeLessThan(20);
  });

  it("knows where measurement ends, and says the future is a guess", () => {
    expect(MEASURED_UNTIL).toBeGreaterThan(julday(2026, 1, 1));
    expect(PREDICTED_UNTIL).toBeGreaterThan(MEASURED_UNTIL);
    expect(deltaTUncertainty(julday(2000, 1, 1))).toBeLessThan(0.01);
    expect(deltaTUncertainty(julday(1850, 1, 1))).toBeGreaterThanOrEqual(0.5);
    const later = [2030, 2100, 2200, 2300, 2400].map((y) => deltaTUncertainty(julday(y, 1, 1)));
    for (let i = 1; i < later.length; i++) expect(later[i]).toBeGreaterThan(later[i - 1]);
    expect(later[4]).toBeGreaterThan(200);
  });

  it("is smooth: no step between neighbouring days anywhere in 1550–2649", () => {
    let worst = 0;
    for (let jd = julday(1550, 2, 1); jd < julday(2649, 1, 1); jd += 1) {
      worst = Math.max(worst, Math.abs(deltaT(jd + 1) - deltaT(jd)));
    }
    expect(worst).toBeLessThan(0.02); // seconds per day
  });
});

describe("sidereal time agrees with ERFA (IAU 2006/2000A)", () => {
  it("within 0.01″ at 300 instants across 1550–2649", () => {
    // The engine's nutation is IAU 2000B, the 77-term abridgement, which
    // tracks 2000A to a milliarcsecond near the present and drifts to a few
    // milliarcseconds by 1550 and 2650. 0.01″ of sidereal time is 0.0007 s.
    let worst = 0;
    ref.gast.jd_ut.forEach((ut: number, i: number) => {
      const ours = apparentSiderealTime(ut, ref.jd_tt[i]) * R2D;
      worst = Math.max(worst, Math.abs(wrap(ours - ref.gast.deg[i])) * 3600);
    });
    expect(worst).toBeLessThan(0.01);
  });
});

// ── the houses, held to their definitions ──────────────────────────────────

const PLACES = [
  { lat: 35.2271, lon: -80.8431 }, // Charlotte
  { lat: 51.5, lon: -0.12 }, // London
  { lat: -33.87, lon: 151.21 }, // Sydney
  { lat: 0.5, lon: 36.8 }, // near the equator
  { lat: 64.13, lon: -21.9 }, // Reykjavík, just short of the polar circle
  { lat: -54.8, lon: -68.3 }, // Ushuaia
];
const MOMENTS = Array.from({ length: 40 }, (_, i) => julday(1800 + i * 15, 1 + (i % 12), 1 + i, (i * 7) % 24));
const DT = 60;

/** Equatorial unit vector of an ecliptic longitude of date. */
function ecl(lon: number, eps: number): [number, number, number] {
  const l = lon * D2R;
  return [Math.cos(l), Math.sin(l) * Math.cos(eps), Math.sin(l) * Math.sin(eps)];
}
/** Local frame: zenith and north point of the horizon, in equatorial coordinates of date. */
function local(armc: number, lat: number) {
  const t = armc * D2R, p = lat * D2R;
  return {
    zenith: [Math.cos(p) * Math.cos(t), Math.cos(p) * Math.sin(t), Math.sin(p)],
    north: [-Math.sin(p) * Math.cos(t), -Math.sin(p) * Math.sin(t), Math.cos(p)],
    east: [-Math.sin(t), Math.cos(t), 0],
  };
}
const dot = (a: number[], b: number[]) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const cross = (a: number[], b: number[]) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];

function each(system: HouseSystem, check: (h: ReturnType<typeof houses>, ctx: { lat: number; eps: number; jd: number }) => void) {
  for (const p of PLACES) for (const jd of MOMENTS) {
    const h = houses(jd, DT, p.lat, p.lon, system);
    check(h, { lat: p.lat, eps: frameOfDate(jd + DT / 86400).eps, jd });
  }
}

describe("the angles", () => {
  it("put the Ascendant on the eastern horizon and the MC on the meridian", () => {
    each("O", (h, { lat, eps }) => {
      const { zenith, east } = local(h.armc, lat);
      const a = ecl(h.ascendant, eps);
      expect(Math.abs(dot(a, zenith))).toBeLessThan(1e-12); // on the horizon
      expect(dot(a, east)).toBeGreaterThan(0); // and in the east
      const m = ecl(h.mc, eps);
      expect(Math.abs(dot(m, east))).toBeLessThan(1e-12); // on the meridian
    });
  });

  it("put the Vertex on the prime vertical, in the west", () => {
    each("O", (h, { lat, eps }) => {
      const { north, east } = local(h.armc, lat);
      const v = ecl(h.vertex, eps);
      expect(Math.abs(dot(v, north))).toBeLessThan(1e-12);
      expect(dot(v, east)).toBeLessThan(0);
    });
  });
});

describe("each house system meets its own definition", () => {
  it("Placidus: every intermediate cusp has covered its third of its own semi-arc", () => {
    each("P", (h, { lat, eps }) => {
      const phi = lat * D2R;
      for (const [i, f, above] of [[10, 1 / 3, true], [11, 2 / 3, true], [1, 2 / 3, false], [2, 1 / 3, false]] as const) {
        const lon = h.cusps[i] * D2R;
        const dec = Math.asin(Math.sin(eps) * Math.sin(lon));
        const ra = Math.atan2(Math.sin(lon) * Math.cos(eps), Math.cos(lon));
        const sa = Math.acos(-Math.tan(phi) * Math.tan(dec));
        const hourAngleEast = wrap((ra * R2D) - h.armc); // degrees east of the meridian
        const expected = above ? f * sa * R2D : 180 - f * (180 - sa * R2D);
        expect(Math.abs(wrap(hourAngleEast - expected))).toBeLessThan(1e-7);
      }
    });
  });

  it("Koch: the intermediate cusps are the Ascendants at thirds of the MC degree's semi-arc", () => {
    // The Ascendant for a sidereal time, worked out here independently: the
    // eastern crossing of the horizon and the ecliptic.
    const ascAt = (armc: number, lat: number, eps: number) => {
      const { zenith, east } = local(armc, lat);
      let d = cross(zenith, [0, -Math.sin(eps), Math.cos(eps)]);
      if (dot(d, east) < 0) d = d.map((x) => -x);
      return norm360(Math.atan2(d[1] * Math.cos(eps) + d[2] * Math.sin(eps), d[0]) * R2D);
    };
    each("K", (h, { lat, eps }) => {
      const dec = Math.asin(Math.sin(eps) * Math.sin(h.mc * D2R));
      const dsa = 90 + Math.asin(Math.tan(lat * D2R) * Math.tan(dec)) * R2D;
      // The MC degree rose dsa ago; the cusps are what was rising at its thirds.
      for (const [i, k] of [[10, -2], [11, -1], [1, 1], [2, 2]] as const) {
        expect(Math.abs(wrap(h.cusps[i] - ascAt(h.armc + (k * dsa) / 3, lat, eps)))).toBeLessThan(1e-9);
      }
    });
  });

  it("Campanus: each cusp lies on the circle through the horizon's poles and its prime-vertical point", () => {
    each("C", (h, { lat, eps }) => {
      const { zenith, north, east } = local(h.armc, lat);
      for (const [i, a] of [[10, 30], [11, 60], [1, 120], [2, 150]] as const) {
        const r = a * D2R;
        const pv = [0, 1, 2].map((k) => Math.cos(r) * zenith[k] + Math.sin(r) * east[k]);
        const normal = cross(north, pv);
        expect(Math.abs(dot(ecl(h.cusps[i], eps), normal))).toBeLessThan(1e-12);
      }
    });
  });

  it("Regiomontanus: each cusp lies on the circle through the horizon's poles and its equator point", () => {
    each("R", (h, { lat, eps }) => {
      const { north } = local(h.armc, lat);
      for (const [i, a] of [[10, 30], [11, 60], [1, 120], [2, 150]] as const) {
        const ra = (h.armc + a) * D2R;
        const normal = cross(north, [Math.cos(ra), Math.sin(ra), 0]);
        expect(Math.abs(dot(ecl(h.cusps[i], eps), normal))).toBeLessThan(1e-12);
      }
    });
  });

  it("Porphyry trisects each quadrant; Equal and Whole Sign step by 30°", () => {
    each("O", (h) => {
      const quadrant = norm360(h.cusps[0] - h.cusps[9]); // MC → Ascendant
      expect(norm360(h.cusps[10] - h.cusps[9])).toBeCloseTo(quadrant / 3, 9);
      expect(norm360(h.cusps[11] - h.cusps[9])).toBeCloseTo((2 * quadrant) / 3, 9);
    });
    each("A", (h) => h.cusps.forEach((c, i) => expect(Math.abs(wrap(c - h.ascendant - 30 * i))).toBeLessThan(1e-9)));
    each("W", (h) => h.cusps.forEach((c, i) => {
      expect(c % 30).toBe(0);
      expect(wrap(c - Math.floor(h.ascendant / 30) * 30 - 30 * i)).toBe(0);
    }));
  });

  it("every quadrant system shares the Ascendant and the MC outside the polar circles", () => {
    for (const s of ["P", "K", "C", "R", "O"] as HouseSystem[]) {
      each(s, (h) => {
        expect(Math.abs(wrap(h.cusps[0] - h.ascendant))).toBeLessThan(1e-9);
        expect(Math.abs(wrap(h.cusps[9] - h.mc))).toBeLessThan(1e-9);
      });
    }
  });

  it("puts houses in opposition: cusp n+6 is cusp n plus 180°", () => {
    for (const s of ["P", "K", "C", "R", "O", "A", "W", "X", "M"] as HouseSystem[]) {
      each(s, (h) => {
        for (let i = 0; i < 6; i++) expect(Math.abs(wrap(h.cusps[i + 6] - h.cusps[i] - 180))).toBeLessThan(1e-9);
      });
    }
  });
});

describe("inside the polar circles", () => {
  it("Placidus and Koch fall back to Porphyry, and say so", () => {
    for (const s of ["P", "K"] as HouseSystem[]) {
      const h = houses(julday(2000, 6, 21, 12), DT, 70, 25, s);
      expect(h.fallback).toBe(true);
      expect(h.system).toBe("O");
      expect(h.cusps).toEqual(houses(julday(2000, 6, 21, 12), DT, 70, 25, "O").cusps);
    }
    expect(houses(julday(2000, 6, 21, 12), DT, 64, 25, "P").fallback).toBe(false);
  });

  it("decides by latitude alone, so a chart never changes system as the sky turns", () => {
    for (let h = 0; h < 24; h += 0.25) {
      expect(houses(julday(2000, 1, 1, h), DT, 66.0, 25, "P").fallback).toBe(false);
      expect(houses(julday(2000, 1, 1, h), DT, 67.0, 25, "P").fallback).toBe(true);
    }
  });

  it("puts the 10th Campanus cusp on the meridian crossing above the horizon", () => {
    for (let hr = 0; hr < 24; hr += 0.5) {
      const h = houses(julday(2000, 12, 21, hr), DT, 75, 15, "C");
      const eps = frameOfDate(julday(2000, 12, 21, hr) + DT / 86400).eps;
      expect(dot(ecl(h.cusps[9], eps), local(h.armc, 75).zenith)).toBeGreaterThan(0);
    }
  });
});
