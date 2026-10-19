# Naelos Ephemeris

Where the planets appear, from 13,200 BC to AD 17,191, computed from NASA JPL's
ephemerides in plain TypeScript: DE440, JPL's most accurate for modern times,
across 1550–2649, and DE441, its long-span companion, either side. No
WebAssembly, no native code, no runtime dependencies: the same files run in
Node, Deno and the browser.

## How it works

DE440 and DE441 are JPL's planetary ephemerides: where every body *is*, as
Chebyshev series, in kilometres. JPL publishes the two together — DE440 for the
modern era, DE441 integrated across 30,000 years and recommended outside it —
and Naelos Ephemeris uses each where JPL says to. It refits them into compact
tables, one file per century:

- `data/` — DE440, 1550–2649: 12 files, 6.6 MB, with Chiron and the asteroids;
- `data-far/` — DE441, 13,200 BC – 1550 and 2649 – AD 17,191: 294 files,
  162 MB. Too many to put in a package: only their manifest (names and
  sha256) ships; the files are attached to the GitHub release
  [`de441-v1`](https://github.com/altariqabraham-design/naelos-ephemeris/releases/tag/de441-v1)
  (see [The far tables](#the-far-tables)).

Then it does at run time everything that turns *where a body is* into *where
it appears* from the centre of the Earth:

1. light-time: the body is seen where it was when its light left;
2. deflection of that light by the Sun's gravity;
3. aberration from the Earth's own motion (relativistic);
4. precession and nutation, from the ICRF to the true equator of date —
   IAU 2006 precession across 1550–2650, the long-term precession of Vondrák,
   Capitaine & Wallace (2011, valid ±200,000 years) beyond, blended over 50
   years so nothing steps where one hands over; and IAU 2000B's nutation terms,
   on the full fundamental arguments (within 3 mas of IAU 2000A across
   1550–2650);
5. the true obliquity, to the true ecliptic of date.

These are the steps of the IERS Conventions, in the same order as ERFA and
NOVAS. The tables carry no astrological convention at all, so none can be wrong
in them.

## Accuracy

Held to JPL by `test/engine.test.ts`, at 300 instants across 1550–2649, against
DE440 reduced by Skyfield (`tools/build_reference.py`):

| | Agreement with JPL |
|---|---|
| Geocentric longitude/latitude, every body | within 0.02″ (typically 0.01″) |
| Right ascension/declination | within 0.02″ |
| Heliocentric | within 0.02″ (Mercury 0.04″) |
| Speed, Sun to Pluto | within 0.00001°/day |
| Speed, Moon | within 0.001°/day |

For scale, a chart prints positions to the arcminute: 60″.

**Beyond 1550–2649** the far tables are fitted to DE441 to the same 0.01″, and
held to **JPL Horizons** itself (`test/far.test.ts`, from positions captured off
the Horizons API): measured as each body's separation from the Sun — which is
the same in any frame — the engine and Horizons agree within **0.004″** at every
instant tried from 9000 BC to AD 9000, light-time, deflection and aberration
included. Their places *of date* differ by more, and the same for every body
at an instant: Horizons turns DE441 into the sky of date with Owen's (1990)
long-term precession, built to match IAU 1976, and Naelos with Vondrák et
al.'s (2011), which matches IAU 2006 and ERFA. The two models are 1–10″ apart
between 3000 BC and AD 9000, and 3′ by 9000 BC. Horizons computes nothing
before March 9999 BC.

At the two seams with the DE440 tables, DE440 and DE441 agree on every planet
to the millimetre but not on the Moon (740 m in 1550, 98 m in 2649); within 16
days of each seam the far tables carry the Moon smoothly onto DE440, so its
position and velocity — and the True Node and Lilith made from them — meet
exactly.

**What is uncertain that far out is the clock, not the planets.** ΔT, the gap
between the Earth's rotation and uniform time, is known to 20 s in 1600, about
5 minutes at the start of the Christian era and some 5 hours by 13,200 BC
(`deltaTUncertainty`, after Morrison & Stephenson 2004). The planets move
slowly enough not to care; the Moon moves 0.55″ a second and the Ascendant
15″, so a chart of 10,000 BC has sound planets, a Moon good to a few degrees,
and no knowable houses. The engine says so; it does not hide it.

Every table segment is checked against DE440 when it is built, and each matches
DE440's position and velocity exactly at both of its ends, so neighbouring
segments join with no step and no kink.

## Interface

The calls answer the way Swiss Ephemeris's `calc_ut` does (same body numbers,
same flags, same result fields), so a chart engine written against Swiss can
switch by changing what it constructs.

```sh
npm install naelos-ephemeris
```

```ts
import { readFileSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { NaelosEphemeris, julday, tablesFor, SE_MARS, SEFLG_SPEED } from "naelos-ephemeris";

// The DE440 tables (1550–2649) ship in the package.
const table = (name: string) => readFileSync(fileURLToPath(import.meta.resolve(`naelos-ephemeris/data/${name}`)));

const t = julday(2026, 10, 17, 12); // UT
const eph = NaelosEphemeris.fromBuffers(tablesFor(t, t).map(table));
const mars = eph.calc_ut(t, SE_MARS, SEFLG_SPEED);
// { longitude, latitude, distance, longitudeSpeed, latitudeSpeed, distanceSpeed, flags }
```

In a browser, serve the files in `data/` yourself and fetch the ones
`tablesFor` names; each is about 640 KB.

Load only the centuries you need: `tablesFor(from, to)` names them (with a
margin for light-time and for ΔT, which is 8.5 days by 13,200 BC), and a file
answers for its own century and can be loaded on its own. `TABLE_FILES` lists
every file's exact span and `isFarTable` says which come from `data-far/`.
Outside what is loaded, and outside 13,200 BC – AD 17,191, the engine refuses
(`{ error }`) rather than guess; Chiron and the asteroids are refused outside
1600–2497, where JPL publishes no orbit for them.

### The far tables

Before 1550 and after 2649, `tablesFor` names files from `data-far/` (their
names start `ne-m` for BC centuries: `ne-m0500.bin` is 500–400 BC). Fetch them
from the release and check each against the manifest before use:

```ts
import { readFileSync } from "node:fs";
import { createHash } from "node:crypto";
import { fileURLToPath } from "node:url";
import { NaelosEphemeris, juldayHistorical, tablesFor, isFarTable, SE_SUN } from "naelos-ephemeris";
import manifest from "naelos-ephemeris/data-far/manifest.json" with { type: "json" };

const RELEASE = "https://github.com/altariqabraham-design/naelos-ephemeris/releases/download/de441-v1/";

async function table(name: string): Promise<Uint8Array> {
  if (!isFarTable(name)) return readFileSync(fileURLToPath(import.meta.resolve(`naelos-ephemeris/data/${name}`)));
  const entry = manifest.files.find((f) => f.name === name);
  if (!entry) throw new Error(`no far table ${name}`);
  const bytes = new Uint8Array(await (await fetch(RELEASE + name)).arrayBuffer());
  if (createHash("sha256").update(bytes).digest("hex") !== entry.sha256) throw new Error(`${name} failed its checksum`);
  return bytes;
}

const ides = juldayHistorical(-43, 3, 15, 11); // 15 March 44 BC, 11:00 UT, as the Romans dated it
const eph = NaelosEphemeris.fromBuffers(await Promise.all(tablesFor(ides, ides).map(table)));
eph.calc_ut(ides, SE_SUN, 0); // the Sun in Pisces
```

GitHub's release downloads do not allow cross-origin requests, so a browser
cannot fetch them directly: download the release once and serve the files
from your own host, with the same checksum check. Cache what you fetch; a file
never changes under its name within a release, and a rebuilt set of tables
gets a new release tag.

**Dates.** `julday` counts days in the proleptic Gregorian calendar for any
year, `juldayJulian` in the Julian, and `juldayHistorical` reads a date as
history records it (Julian before 15 October 1582, and the ten skipped days
refused); `revjul` goes back. Years are astronomical: year 0 is 1 BC, and
`formatYear(-43)` is "44 BC".

## Exact moments

Transits, stations, ingresses, returns and lunations are **solved, not
sampled**, and every search is complete: no pass is missed, including the three
passes of a retrograde loop.

```ts
eph.crossings(SE_MARS, [natalSun, natalSun + 90, natalSun + 180], startUT, endUT); // transits to a natal point
eph.stations(SE_MERCURY, startUT, endUT);   // { jd, longitude, turning: "retrograde" | "direct" }
eph.ingresses(SE_SATURN, startUT, endUT);   // with the sign entered, re-entries included
eph.returns(SE_SUN, natalSun, startUT, endUT);
eph.lunations(startUT, endUT);              // new, first quarter, full, last quarter
eph.aspectsBetween(SE_JUPITER, SE_SATURN, [0, 90, 180], startUT, endUT);
```

How it is complete: each body is sampled at under half the shortest gap
between its stations across the whole span (measured: Mercury 19.75 days,
Mars 59.5, and so on; the steps are a tenth of that), so at least two samples
fall inside every retrograde arc and the sampled track turns wherever the body
does; each turning point is refined exactly; and between turning points the
longitude moves one way, so each degree is reached once per lap. The tests
hold every search to a brute-force count at fine steps.

## A chart that knows how sure it is

Birth times are rarely exact. Give a window instead of an instant and the
engine says what the window allows:

```ts
const u = eph.chartUncertainty({ from: birth - 30 / 1440, to: birth + 30 / 1440, lat, lon });
// every body and angle: signs and houses with the share of the window each holds,
// `certain` when nothing changes; every aspect with the share it is in orb;
// and ΔT with its own uncertainty (a millisecond today, minutes by 2400)
```

Shares are exact: every boundary — the Ascendant entering a sign, a cusp
passing a planet, an aspect leaving orb — is found to the second. A ±30-minute
window takes about 50 ms; a whole unknown day about 400 ms. Aspects default to
the five majors at fixed orbs; `aspects` changes the set, and `pairs` gives
each pair its own orb, as a chart whose luminaries take wider orbs than its
outer planets needs.

**Rectification** narrows the window from dated life events, the way
astrologers do it by hand — transits of Mars to Pluto, solar-arc directions and
the progressed Moon, to and from the Ascendant and Midheaven — and lists every
contact it counted:

```ts
const r = eph.rectify({ from, to, lat, lon, events: [{ jd: marriage }, { jd: move, weight: 0.5 }] });
// r.best — candidate minutes with their evidence; r.window68 — the narrowed window
```

Its result is *support*, not a measured probability: it says which minutes those
techniques favour, and why. With events built on a known birth minute's own
angles, eight or more of them recover that minute to within two.

## Every position shows its working

```ts
const p = eph.provenance(jdUT, SE_MARS);
// p.source  "JPL DE440, table ne-1900.bin segment 66"
// p.steps   ["light-time (10.8 min): seen where it was when its light left", "deflection …", …]
// p.clock   { deltaT: 58.85, sigma: 0.001, basis: "measured (IERS)" }
// p.verify  { horizons: "https://ssd.jpl.nasa.gov/api/horizons.api?…" }
```

The link asks NASA JPL's public Horizons service for the same quantity (the
apparent geocentric ecliptic place of date), so anyone can check a number
without trusting the engine. `node tools/check_horizons.mjs` does that live for
every body: on 2026-10-18 the worst difference was 0.074″. (Horizons uses the
older IAU 1976/1980 precession-nutation and reads UT as UTC, which accounts for
that.)

## Status

Built:

- the Sun, the Moon and Mercury to Pluto (Jupiter to Pluto as system
  barycentres, as is standard); geocentric, heliocentric and equatorial
  places; speeds;
- **ΔT** (UT to TT): IERS measurements from 1973, the Morrison–Stephenson–
  Hohenkerk–Zawilski (2021) reconstruction before that, and a published
  long-term curve after the IERS predictions end. `deltaTUncertainty` says how
  much of a guess any date's ΔT is: a millisecond today, minutes by 2400;
- **sidereal time** (IAU 2006 across 1550–2650; beyond, the Earth rotation angle
  less the equation of the origins integrated through the long-term
  precession) and **houses**: Placidus, Koch, Campanus,
  Regiomontanus, Porphyry, Equal, Whole Sign, Meridian and Morinus, with the
  Ascendant, MC, Vertex and East Point. Inside the polar circles Placidus and
  Koch are undefined, and Porphyry is drawn instead, flagged as `fallback`.
  Every system is tested against its own definition rather than against
  another program;
- the **Moon's nodes and apogee**: the mean North Node and mean Lilith from the
  published ELP2000-82B mean elements — which drift from the true points'
  average by 2.4° over 15,000 years, so beyond 1550–2650 they take a correction
  measured on DE441 itself (`tools/build_lunar_mean.mts`) — and the True Node
  and True Lilith from
  the Moon's osculating orbit, computed from the engine's own lunar position
  and velocity (within 0.3″ and 1.5″ of DE440's exact state). Lilith has a
  latitude: she is the apogee of an inclined orbit, a point in space;
- **Chiron, Ceres, Pallas, Juno and Vesta**, 1600–2497, from JPL Horizons'
  current orbit solutions, within 0.02″ of those orbits as seen from Earth and
  from the Sun;
- **eclipses**: the nearest solar or lunar eclipse before or after any moment,
  with greatest eclipse and type (total, annular, hybrid, partial, penumbral).
  Solar types are decided on the Earth's true figure; lunar ones by Danjon's
  rule for the Earth's shadow. Over 1800–2400 this finds every eclipse in
  Swiss Ephemeris's list, within 11 s of its time where ΔT is measured.

## Rebuilding

Every number in `data/` and `src/generated/` is produced by a tool. To rebuild:

```sh
python -m venv .venv && .venv/bin/pip install -r tools/requirements.txt
python -I tools/build_tables.py --de440 /path/to/de440.bsp --horizons /path/to/spk-folder
python -I tools/build_frames.py
python -I tools/build_deltat.py --finals /path/to/finals2000A.all
python -I tools/build_reference.py --de440 /path/to/de440.bsp --horizons /path/to/spk-folder
# the far tables, and what depends on them
python -I tools/build_far_tables.py --de441 /path/to/de441.bsp --de440 /path/to/de440.bsp --out /path/to/far
python -I tools/build_deltat_long.py --finals /path/to/finals2000A.all
python -I tools/build_ltp_reference.py
node tools/build_file_index.mjs
npx tsx tools/build_lunar_mean.mts /path/to/far
```

`build_far_tables.py` writes all 294 files; copy its `manifest.json` and the
five samples the tests read into `data-far/`, and attach all 294 to a NEW
release tag (`de441-v2`, …), never replacing the files of an old one: someone
may have pinned it.

- **DE440** (`de440.bsp`, 120 MB) and **DE441** (`de441.bsp`, 3.3 GB):
  <https://ssd.jpl.nasa.gov/ftp/eph/planets/bsp/>.
- **The asteroids**: one SPK per body from the JPL Horizons API, 1600 to 2500,
  the widest span it will write
  (`https://ssd.jpl.nasa.gov/api/horizons.api?format=json&COMMAND='1;'&EPHEM_TYPE=SPK&START_TIME='1600-01-01'&STOP_TIME='2500-01-01'`,
  and `'2;'`, `'3;'`, `'4;'`, `'2060;'`), saved as `ceres.bsp`, `pallas.bsp`,
  `juno.bsp`, `vesta.bsp` and `chiron.bsp`.
- **Earth rotation** (`finals2000A.all`): the IERS rapid series, from the US
  Naval Observatory (<https://maia.usno.navy.mil/ser7/finals2000A.all>).
  Rebuild ΔT from a fresh copy every few months.

The build records every source's SHA-256 in `data/manifest.json` and
`data-far/manifest.json`.

## Sources and licences

- **DE440 and DE441**, NASA/JPL (Park et al. 2021), and the **JPL Horizons**
  orbits of Chiron and the asteroids. Works of the US Government, in the public
  domain. Everything in `data/` and `data-far/` is derived from them.
- **ΔT**: IERS Earth orientation data (freely published); the splines of
  Morrison, Stephenson, Hohenkerk & Zawilski (2021) and the long-term
  parabola of Stephenson, Morrison & Hohenkerk (2016).
- **Long-term precession** (Vondrák, Capitaine & Wallace 2011, corrigendum
  2012), in `src/generated/frames-long.ts`, checked against ERFA's `ltpb`.
- **ΔT uncertainty** before 1700: Morrison & Stephenson (2004), 0.8 t² seconds.
- **IAU 2006 precession** (Capitaine & Wallace 2006), **IAU 2000B nutation**
  (McCarthy & Luzum 2003) and the **fundamental arguments** of Simon et al.
  (1994) as adopted by IERS 2003: published models; the coefficients in
  `src/generated/frames.ts` are checked against ERFA by `tools/build_frames.py`.
- **ERFA** (BSD-3), **Skyfield** and **jplephem** (MIT) and NASA's **SPICE**
  toolkit (public domain, through SpiceyPy, MIT) are used only by the build
  tools, to produce and check numbers. None of them ships with the engine.

## Licence

Copyright 2026 Altariq Abraham, doing business as Apex Node. Licensed under
the [Apache License, Version 2.0](LICENSE); see [NOTICE](NOTICE) for the
attributions above in the form the licence asks to be carried along. The
licence grants no rights to the names "Naelos" or "Apex Node".
