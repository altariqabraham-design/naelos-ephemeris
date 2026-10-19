// The Moon's nodes and apogee — mean and true North Node, mean and true Lilith.
//
// The true points are held to DE440's EXACT lunar state (tools/build_reference.py
// uses JPL's own analytic velocity), so this also tests that the engine's Moon
// table carries a good enough velocity: the apogee magnifies a velocity error
// eighteen-fold, which is why the Moon's table rounds its high orders finely.

import { describe, it, expect } from "vitest";
import { readFileSync } from "node:fs";
import { gunzipSync } from "node:zlib";
import { resolve } from "node:path";
import {
  NaelosEphemeris, SE_MEAN_NODE, SE_TRUE_NODE, SE_MEAN_APOG, SE_OSCU_APOG, SEFLG_SPEED, SEFLG_HELCTR,
  SEFLG_EQUATORIAL, type Position,
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

function worst(body: number, lonKey: string, latKey?: string) {
  let lon = 0, lat = 0;
  jds.forEach((jd, i) => {
    const r = ok(eph.calc(jd, body, SEFLG_SPEED));
    lon = Math.max(lon, Math.abs(wrap(r.longitude - ref.lunar[lonKey][i])) * 3600);
    lat = Math.max(lat, Math.abs(r.latitude - (latKey ? ref.lunar[latKey][i] : 0)) * 3600);
  });
  return { lon, lat };
}

describe("the true points, from the Moon's osculating orbit, agree with DE440's exact state", () => {
  it("True Node: within 0.3″, exactly on the ecliptic", () => {
    const w = worst(SE_TRUE_NODE, "true_node");
    expect(w.lon).toBeLessThan(0.3);
    expect(w.lat).toBeLessThan(1e-9);
  });

  it("True Lilith: within 1.5″ in longitude and 0.3″ in latitude", () => {
    const w = worst(SE_OSCU_APOG, "true_apogee_lon", "true_apogee_lat");
    expect(w.lon).toBeLessThan(1.5);
    expect(w.lat).toBeLessThan(0.3);
  });
});

describe("the mean points follow the published lunar theory (ELP2000-82B)", () => {
  // Same formulas on both sides; the only difference is IAU 2000B nutation
  // here and 2000A in the reference, a few milliarcseconds.
  it("Mean Node: within 0.01″", () => {
    expect(worst(SE_MEAN_NODE, "mean_node").lon).toBeLessThan(0.01);
  });

  it("Mean Lilith: within 0.01″, with the latitude of the inclined mean orbit", () => {
    const w = worst(SE_MEAN_APOG, "mean_apogee_lon", "mean_apogee_lat");
    expect(w.lon).toBeLessThan(0.01);
    expect(w.lat).toBeLessThan(0.01);
  });
});

describe("how they move", () => {
  it("the Mean Node always runs backwards, about 0.053° a day", () => {
    for (const jd of jds) {
      const r = ok(eph.calc(jd, SE_MEAN_NODE, SEFLG_SPEED));
      expect(r.longitudeSpeed).toBeLessThan(-0.0525);
      expect(r.longitudeSpeed).toBeGreaterThan(-0.0535);
    }
  });

  it("the True Node wobbles about the mean one but never strays more than about 2°", () => {
    for (const jd of jds) {
      const t = ok(eph.calc(jd, SE_TRUE_NODE, 0)).longitude;
      const m = ok(eph.calc(jd, SE_MEAN_NODE, 0)).longitude;
      expect(Math.abs(wrap(t - m))).toBeLessThan(2.2);
    }
  });

  it("True Lilith can stray up to about 30° from the mean one", () => {
    let far = 0;
    for (const jd of jds) {
      const t = ok(eph.calc(jd, SE_OSCU_APOG, 0)).longitude;
      const m = ok(eph.calc(jd, SE_MEAN_APOG, 0)).longitude;
      far = Math.max(far, Math.abs(wrap(t - m)));
    }
    expect(far).toBeGreaterThan(15);
    expect(far).toBeLessThan(35);
  });

  it("have right ascension and declination too, and no heliocentric place", () => {
    const r = ok(eph.calc(2451545, SE_OSCU_APOG, SEFLG_EQUATORIAL));
    expect(Number.isFinite(r.longitude) && Number.isFinite(r.latitude)).toBe(true);
    for (const b of [SE_MEAN_NODE, SE_TRUE_NODE, SE_MEAN_APOG, SE_OSCU_APOG]) {
      expect("error" in eph.calc(2451545, b, SEFLG_HELCTR)).toBe(true);
    }
  });
});
