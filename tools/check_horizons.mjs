// check_horizons.mjs — ask NASA JPL Horizons, live, whether it agrees.
//
// For a few bodies and moments, takes the verification link provenance()
// gives, asks Horizons for the same quantity (apparent geocentric ecliptic
// longitude and latitude of date), and prints the difference. Needs network.
//
//   node tools/check_horizons.mjs [YYYY-MM-DDTHH:MM (UT) ...]
//
// Expect about 0.05–0.1″: Horizons uses IAU 1976/1980 precession-nutation
// and reads UT as UTC (the engine's UT is UT1, up to 0.9 s different, which
// is up to 0.5″ for the Moon).

import { readFileSync } from "node:fs";
import { fileURLToPath } from "node:url";

const ROOT = new URL("..", import.meta.url);
const { NaelosEphemeris, julday } = await import(new URL("src/index.ts", ROOT).href);
const manifest = JSON.parse(readFileSync(fileURLToPath(new URL("data/manifest.json", ROOT)), "utf8"));
const eph = NaelosEphemeris.fromBuffers(manifest.files.map((f) => readFileSync(fileURLToPath(new URL(`data/${f.name}`, ROOT)))));
const wrap = (d) => ((d + 540) % 360) - 180;

const args = process.argv.slice(2);
const moments = (args.length ? args : ["2000-01-01T12:00", new Date().toISOString().slice(0, 16)]).map((s) => {
  const d = new Date(`${s}Z`);
  return julday(d.getUTCFullYear(), d.getUTCMonth() + 1, d.getUTCDate(), d.getUTCHours() + d.getUTCMinutes() / 60);
});

let worst = 0;
for (const jd of moments) {
  for (const body of [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 15, 17, 18, 19, 20]) {
    const p = eph.provenance(jd, body);
    if ("error" in p || !p.verify) continue;
    const text = await (await fetch(p.verify.horizons)).text();
    const row = text.split("$$SOE")[1]?.split("$$EOE")[0]?.trim();
    if (!row) { console.log(`${p.name}: Horizons gave no data`); continue; }
    const cols = row.split(",").map((s) => s.trim());
    const nums = cols.slice(1).filter((s) => s !== "").map(Number).filter(Number.isFinite);
    const [lon, lat] = nums.slice(-2);
    const d = Math.hypot(wrap(p.longitude - lon) * Math.cos((lat * Math.PI) / 180), p.latitude - lat) * 3600;
    worst = Math.max(worst, d);
    console.log(`${cols[0].padEnd(24)} ${p.name.padEnd(28)} engine ${p.longitude.toFixed(6)}°  JPL ${lon.toFixed(6)}°  ${d.toFixed(3)}″`);
  }
}
console.log(`worst ${worst.toFixed(3)}″`);
process.exitCode = worst < 1 ? 0 : 1;
