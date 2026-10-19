// Every position shows its working: its source, its steps, its clock, and an
// independent check. (The live comparison with JPL Horizons needs network and
// lives in tools/check_horizons.mjs; it agrees to about 0.1″.)

import { describe, it, expect } from "vitest";
import { readFileSync } from "node:fs";
import { resolve } from "node:path";
import { NaelosEphemeris, julday, horizonsUrl, SE_MARS, SE_MOON, SE_CERES, SE_TRUE_NODE, SEFLG_HELCTR, type Provenance } from "../src/index.ts";

const ROOT = resolve(__dirname, "..");
const manifest = JSON.parse(readFileSync(resolve(ROOT, "data/manifest.json"), "utf8"));
const eph = NaelosEphemeris.fromBuffers(manifest.files.map((f: { name: string }) => readFileSync(resolve(ROOT, "data", f.name))));
const ok = (p: ReturnType<typeof eph.provenance>): Provenance => {
  if ("error" in p) throw new Error(p.error);
  return p;
};
const FOUNDER = julday(1992, 9, 6, 10 + 41 / 60);

describe("provenance", () => {
  it("names the JPL data and the table segment a planet came from", () => {
    const p = ok(eph.provenance(FOUNDER, SE_MARS));
    expect(p.source).toMatch(/^JPL DE440, table ne-1900\.bin segment \d+$/);
    expect(p.engine).toMatch(/^Naelos Ephemeris /);
  });

  it("lists every correction in order, with the light-time", () => {
    const p = ok(eph.provenance(FOUNDER, SE_MARS));
    expect(p.steps[0]).toMatch(/^light-time \(\d+\.\d min\)/);
    expect(p.steps.join(" | ")).toMatch(/deflection.*aberration.*precession.*ecliptic of date/);
  });

  it("says how certain the clock was", () => {
    expect(ok(eph.provenance(FOUNDER, SE_MARS)).clock.basis).toBe("measured (IERS)");
    expect(ok(eph.provenance(julday(1900, 1, 1), SE_MARS)).clock.basis).toBe("historical reconstruction");
    const far = ok(eph.provenance(julday(2300, 1, 1), SE_MARS)).clock;
    expect(far.basis).toBe("long-term model");
    expect(far.sigma).toBeGreaterThan(60);
  });

  it("gives an independent check against JPL Horizons, asking for the same quantity", () => {
    const p = ok(eph.provenance(FOUNDER, SE_MARS));
    const u = new URL(p.verify!.horizons);
    expect(u.origin + u.pathname).toBe("https://ssd.jpl.nasa.gov/api/horizons.api");
    expect(u.searchParams.get("COMMAND")).toBe("'4'"); // Mars system barycentre, as computed
    expect(u.searchParams.get("CENTER")).toBe("'500@399'"); // from the Earth's centre
    expect(u.searchParams.get("QUANTITIES")).toBe("'31'"); // apparent ecliptic place of date
    expect(Number(u.searchParams.get("TLIST")!.replace(/'/g, ""))).toBeCloseTo(FOUNDER, 7);
  });

  it("asks Horizons for an asteroid by its small-body number", () => {
    expect(new URL(horizonsUrl(SE_CERES, FOUNDER)!).searchParams.get("COMMAND")).toBe("'1;'");
    expect(ok(eph.provenance(FOUNDER, SE_CERES)).source).toMatch(/^JPL Horizons orbit \(SPK 20000001\)/);
  });

  it("explains derived points instead of linking them, and the Moon rests on her own series", () => {
    const node = ok(eph.provenance(FOUNDER, SE_TRUE_NODE));
    expect(node.verify).toBeUndefined();
    expect(node.source).toMatch(/osculating orbit/);
    expect(ok(eph.provenance(FOUNDER, SE_MOON)).source).toMatch(/table ne-1900\.bin segment/);
  });

  it("has no Horizons link for heliocentric places, which Horizons' quantity 31 does not give", () => {
    expect(ok(eph.provenance(FOUNDER, SE_MARS, SEFLG_HELCTR)).verify).toBeUndefined();
  });
});
