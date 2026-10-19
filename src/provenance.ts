// provenance.ts — every position shows its working.
//
// An astrology reading cites placements; a reader deserves to know where those
// numbers came from and to check them without trusting anyone. For any
// position this says: which JPL data it rests on, which table segment served
// it, every correction applied in order, how certain the clock was — and gives
// a link to NASA JPL's public Horizons service that computes the same quantity
// independently.
//
// Horizons' quantity 31 is the apparent geocentric ecliptic longitude and
// latitude of date, with light-time, deflection and aberration — the same
// quantity as here. Horizons uses the older IAU 1976/1980 precession and
// nutation and its own ΔT, so agreement is to about a tenth of an arcsecond
// near the present: far inside the arcminute a chart prints.

import { MEASURED_UNTIL, PREDICTED_UNTIL, deltaTUncertainty } from "./deltat.ts";
import { FILES, NEAR_SOURCE, FAR_SOURCE } from "./generated/files.ts";
import { iauWeight } from "./frames.ts";

export const ENGINE = "Naelos Ephemeris 0.1.0";

/** What each body is in JPL's numbering, and what it rests on. */
const SOURCE: Record<number, { name: string; horizons: string; data: string }> = {
  0: { name: "Sun", horizons: "10", data: "JPL DE440" },
  1: { name: "Moon", horizons: "301", data: "JPL DE440" },
  2: { name: "Mercury", horizons: "199", data: "JPL DE440" },
  3: { name: "Venus", horizons: "299", data: "JPL DE440" },
  4: { name: "Mars (system barycentre)", horizons: "4", data: "JPL DE440" },
  5: { name: "Jupiter (system barycentre)", horizons: "5", data: "JPL DE440" },
  6: { name: "Saturn (system barycentre)", horizons: "6", data: "JPL DE440" },
  7: { name: "Uranus (system barycentre)", horizons: "7", data: "JPL DE440" },
  8: { name: "Neptune (system barycentre)", horizons: "8", data: "JPL DE440" },
  9: { name: "Pluto (system barycentre)", horizons: "9", data: "JPL DE440" },
  15: { name: "Chiron", horizons: "2060;", data: "JPL Horizons orbit (SPK 20002060)" },
  17: { name: "Ceres", horizons: "1;", data: "JPL Horizons orbit (SPK 20000001)" },
  18: { name: "Pallas", horizons: "2;", data: "JPL Horizons orbit (SPK 20000002)" },
  19: { name: "Juno", horizons: "3;", data: "JPL Horizons orbit (SPK 20000003)" },
  20: { name: "Vesta", horizons: "4;", data: "JPL Horizons orbit (SPK 20000004)" },
};

const DERIVED: Record<number, string> = {
  10: "Mean North Node: ELP2000-82B mean elements (Meeus, Astronomical Algorithms ch. 47), on the true equinox of date",
  11: "True North Node: node of the Moon's osculating orbit, from the Moon's DE440 position and velocity",
  12: "Mean Lilith: apogee of the ELP2000-82B mean lunar orbit, projected on the ecliptic of date",
  13: "True Lilith: apogee of the Moon's osculating orbit, from its DE440 position and velocity",
};

export interface Provenance {
  engine: string;
  body: number;
  name: string;
  jdUT: number;
  jdTT: number;
  /** Where the number rests. */
  source: string;
  /** What was done to it, in order. */
  steps: string[];
  clock: {
    deltaT: number;
    sigma: number;
    basis: "measured (IERS)" | "predicted (IERS)" | "historical reconstruction" | "long-term model";
  };
  longitude: number;
  latitude: number;
  /** An independent check, where one exists. */
  verify?: { horizons: string; expect: string };
}

/** 1 January 720 BC (Julian), where the eclipse records behind ΔT begin. */
const RECORDS_BEGIN = 1458085.5;

function clockBasis(jdUT: number): Provenance["clock"]["basis"] {
  if (jdUT < RECORDS_BEGIN) return "long-term model";
  if (jdUT < 2441317.5) return "historical reconstruction"; // before 1972
  if (jdUT <= MEASURED_UNTIL) return "measured (IERS)";
  if (jdUT <= PREDICTED_UNTIL) return "predicted (IERS)";
  return "long-term model";
}

/** A JPL Horizons API request for the apparent geocentric ecliptic place of a body at a UT instant. */
export function horizonsUrl(body: number, jdUT: number): string | null {
  const s = SOURCE[body];
  if (!s) return null;
  const params: Array<[string, string]> = [
    ["format", "text"],
    ["COMMAND", `'${s.horizons}'`],
    ["OBJ_DATA", "'NO'"],
    ["MAKE_EPHEM", "'YES'"],
    ["EPHEM_TYPE", "'OBSERVER'"],
    ["CENTER", "'500@399'"], // the Earth's centre
    ["TLIST", `'${jdUT.toFixed(8)}'`],
    ["TLIST_TYPE", "'JD'"],
    ["TIME_TYPE", "'UT'"],
    ["QUANTITIES", "'31'"], // apparent ecliptic longitude and latitude of date
    ["ANG_FORMAT", "'DEG'"],
    ["EXTRA_PREC", "'YES'"],
    ["CSV_FORMAT", "'YES'"],
  ];
  return "https://ssd.jpl.nasa.gov/api/horizons.api?" + params.map(([k, v]) => `${k}=${encodeURIComponent(v)}`).join("&");
}

export function provenanceOf(args: {
  body: number; jdUT: number; deltaT: number; longitude: number; latitude: number; distanceAU: number;
  heliocentric: boolean; segment?: { file: string; index: number; from: number; to: number };
}): Provenance {
  const { body, jdUT, deltaT, longitude, latitude, distanceAU, heliocentric } = args;
  const jdTT = jdUT + deltaT / 86400;
  const derived = DERIVED[body];
  const s = SOURCE[body] ?? (body === 14 ? { name: "Earth", horizons: "399", data: "JPL DE440" } : undefined);
  const name = s?.name ?? derived?.split(":")[0] ?? `body ${body}`;
  const minutes = (distanceAU * 499.004784) / 60;
  const steps = derived
    ? ["the Moon's geometric state (no light-time: a point of an orbit, not a body seen)", "IAU 2006 precession and IAU 2000B nutation to the true equinox of date", "true ecliptic of date"]
    : heliocentric
      ? [`light-time from the body to the Sun (${minutes.toFixed(1)} min)`, "IAU 2006 precession and IAU 2000B nutation to the true equinox of date", "true ecliptic of date"]
      : [
          `light-time (${minutes < 1 ? (minutes * 60).toFixed(2) + " s" : minutes.toFixed(1) + " min"}): seen where it was when its light left`,
          "deflection of that light by the Sun's gravity",
          "aberration from the Earth's motion (relativistic)",
          "IAU 2006 precession and IAU 2000B nutation to the true equator of date",
          "true obliquity, to the true ecliptic of date",
        ];
  // The far tables (before 1550, after 2649) were fitted from DE441, not DE440,
  // and beyond the DE440 tables the precession is the long-term model.
  const far = args.segment ? FILES.some(([n, , , f]) => n === args.segment!.file && f === 1) : false;
  const ephemeris = (text: string) => (far ? text.split(NEAR_SOURCE).join(FAR_SOURCE) : text);
  const longTerm = iauWeight(jdTT) < 1;
  const shownSteps = longTerm
    ? steps.map((x) => x.replace("IAU 2006 precession", "long-term precession (Vondrák, Capitaine & Wallace 2011)"))
    : steps;
  const source = ephemeris(derived ?? `${s?.data ?? NEAR_SOURCE}${args.segment ? `, table ${args.segment.file} segment ${args.segment.index}` : ""}`);
  // Horizons answers from 9999 BC; before that there is nothing to check against.
  const url = heliocentric || derived || jdUT < HORIZONS_FIRST ? null : horizonsUrl(body, jdUT);
  const expect = longTerm
    ? "the same separations between bodies to 0.01″; its places of date differ by a few arcseconds within 3000 BC – AD 9000, and by arcminutes before 8000 BC, because Horizons turns DE441 into the sky of date with Owen's (1990) long-term precession, built to match IAU 1976, and Naelos with Vondrák et al.'s (2011), which matches IAU 2006"
    : "within about 0.1″ near the present (JPL uses IAU 1976/1980 precession and nutation, and its own ΔT)";
  return {
    engine: ENGINE, body, name, jdUT, jdTT, source, steps: shownSteps,
    clock: { deltaT, sigma: deltaTUncertainty(jdUT), basis: clockBasis(jdUT) },
    longitude, latitude,
    ...(url ? { verify: { horizons: url, expect } } : {}),
  };
}

/**
 * 16 March 9999 BC (Julian): Horizons computes nothing before 15 March 9999 BC
 * 21:49 UT ("No ephemeris for target … prior to B.C. 9999-MAR-15"), so no
 * link is offered for an earlier moment.
 */
const HORIZONS_FIRST = -1930637.5;
