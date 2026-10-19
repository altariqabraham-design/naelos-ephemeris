// Chiron, Ceres, Pallas, Juno and Vesta against JPL.
//
// Their orbits come from JPL Horizons, and the reference is computed from the
// same Horizons files by a different road: NASA's SPICE toolkit for
// light-time, ERFA for the Sun's light deflection, aberration and the frame of
// date (tools/build_reference.py).

import { describe, it, expect } from "vitest";
import { readFileSync } from "node:fs";
import { gunzipSync } from "node:zlib";
import { resolve } from "node:path";
import {
  NaelosEphemeris, julday, SE_CHIRON, SE_CERES, SE_PALLAS, SE_JUNO, SE_VESTA, SEFLG_SPEED, SEFLG_HELCTR, type Position,
} from "../src/index.ts";

const ROOT = resolve(__dirname, "..");
const manifest = JSON.parse(readFileSync(resolve(ROOT, "data/manifest.json"), "utf8"));
const eph = NaelosEphemeris.fromBuffers(manifest.files.map((f: { name: string }) => readFileSync(resolve(ROOT, "data", f.name))));
const ref = JSON.parse(gunzipSync(readFileSync(resolve(ROOT, "test/reference.json.gz"))).toString());
const jds: number[] = ref.jd_tt;
const wrap = (d: number) => ((d + 540) % 360) - 180;
const ok = (r: ReturnType<typeof eph.calc>): Position => {
  if ("error" in r) throw new Error(r.error);
  return r;
};

const BODIES: Array<[number, string]> = [[SE_CHIRON, "Chiron"], [SE_CERES, "Ceres"], [SE_PALLAS, "Pallas"], [SE_JUNO, "Juno"], [SE_VESTA, "Vesta"]];

describe("the asteroids agree with JPL Horizons", () => {
  for (const [b, name] of BODIES) {
    it(`${name}: geocentric and heliocentric within 0.02″`, () => {
      let geo = 0, helio = 0, n = 0;
      jds.forEach((jd, i) => {
        if (ref.geo[b].lon[i] === null) return; // outside what Horizons publishes
        n++;
        const g = ok(eph.calc(jd, b, SEFLG_SPEED));
        geo = Math.max(geo, Math.hypot(wrap(g.longitude - ref.geo[b].lon[i]) * Math.cos((g.latitude * Math.PI) / 180), g.latitude - ref.geo[b].lat[i]) * 3600);
        const h = ok(eph.calc(jd, b, SEFLG_HELCTR));
        helio = Math.max(helio, Math.hypot(wrap(h.longitude - ref.helio[b].lon[i]) * Math.cos((h.latitude * Math.PI) / 180), h.latitude - ref.helio[b].lat[i]) * 3600);
      });
      expect(n).toBeGreaterThan(200);
      expect(geo).toBeLessThan(0.02);
      expect(helio).toBeLessThan(0.02);
    });
  }

  it("are given only where JPL publishes their orbits, and refused outside", () => {
    const span = manifest.span.asteroids;
    // Horizons publishes 1600–2500. Chiron is stored in 1024-day segments,
    // which stop at the last whole one inside that: March 2497.
    expect(span.start_jd_tdb).toBeLessThan(julday(1600, 2, 1));
    expect(span.end_jd_tdb).toBeGreaterThan(julday(2497, 1, 1));
    for (const [b] of BODIES) {
      expect("error" in eph.calc(julday(1580, 1, 1), b, 0)).toBe(true);
      expect("error" in eph.calc(julday(2550, 1, 1), b, 0)).toBe(true);
    }
    // while the planets carry on
    expect("error" in eph.calc(julday(1580, 1, 1), 4, 0)).toBe(false);
  });

  it("records which Horizons files they came from", () => {
    expect(Object.keys(manifest.asteroids.files).sort()).toEqual(["ceres", "chiron", "juno", "pallas", "vesta"]);
    for (const f of Object.values(manifest.asteroids.files) as Array<{ sha256: string }>) expect(f.sha256).toMatch(/^[0-9a-f]{64}$/);
  });
});
