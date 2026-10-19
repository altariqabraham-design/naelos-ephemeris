// deltat.ts — ΔT = TT − UT: the gap between the Earth's rotation and uniform time.
//
// Ephemerides run on uniform time (TT); birth times are clock times (UT). The
// gap was 13.7 s in 1800 and is about 69 s today; it was 8.5 DAYS in 13,200 BC.
// It is MEASURED from 1973 to the table's MEASURED_UNTIL, reconstructed from
// historical observations before that (eclipse and occultation records reach
// back to 720 BC), a long-term model before those, and PREDICTED after the
// present — and a prediction is a guess that grows worse with every year.
// deltaTUncertainty says how much of a guess.

import { DELTAT, DELTAT_START, DELTAT_STEP, MEASURED_UNTIL, PREDICTED_UNTIL } from "./generated/deltat.ts";
import {
  DELTAT_PAST, DELTAT_PAST_START, DELTAT_PAST_STEP, DELTAT_FUTURE, DELTAT_FUTURE_START, DELTAT_FUTURE_STEP,
} from "./generated/deltat-long.ts";

export { MEASURED_UNTIL, PREDICTED_UNTIL };

/** Cubic Hermite through evenly spaced samples. */
function hermite(values: readonly number[], start: number, step: number, jd: number): number {
  const last = values.length - 1;
  const x = (jd - start) / step;
  let i = Math.floor(x);
  if (i < 1) i = 1;
  else if (i > last - 2) i = last - 2;
  const u = x - i;
  const p0 = values[i - 1], p1 = values[i], p2 = values[i + 1], p3 = values[i + 2];
  const m1 = (p2 - p0) / 2, m2 = (p3 - p1) / 2;
  const u2 = u * u, u3 = u2 * u;
  return (2 * u3 - 3 * u2 + 1) * p1 + (u3 - 2 * u2 + u) * m1 + (-2 * u3 + 3 * u2) * p2 + (u3 - u2) * m2;
}

/** Where the quarter-year table has two neighbours either side. */
const QUARTER_FROM = DELTAT_START + DELTAT_STEP;
const QUARTER_TO = DELTAT_START + (DELTAT.length - 2) * DELTAT_STEP;

/**
 * ΔT in seconds at a UT Julian Day: the quarter-year table across 1549–2651,
 * the 5-year table before it and the 25-year table after (all three from the
 * same model; see tools/build_deltat_long.py).
 */
export function deltaT(jdUT: number): number {
  if (jdUT < QUARTER_FROM) return hermite(DELTAT_PAST, DELTAT_PAST_START, DELTAT_PAST_STEP, jdUT);
  if (jdUT > QUARTER_TO) return hermite(DELTAT_FUTURE, DELTAT_FUTURE_START, DELTAT_FUTURE_STEP, jdUT);
  return hermite(DELTAT, DELTAT_START, DELTAT_STEP, jdUT);
}

const YEAR = 365.25;
const J2400 = 2597641.5;
/** 1820, the epoch Morrison & Stephenson's uncertainty is measured from. */
const J1820 = 2385801.0;

/**
 * One standard deviation of ΔT, in seconds — how far the true value may lie
 * from deltaT(). Approximate by design; its job is to tell a chart for 1650,
 * 2300 or 3000 BC how uncertain its clock is, not to be a model of its own.
 *
 *   atomic era (1955 on, while measured)   0.001 s
 *   1900–1955                              0.1 s
 *   1860–1900                              0.5 s
 *   1800–1860                              1 s
 *   1700–1800 (telescopic timings)         5 s
 *   before 1700                            20 s, or Morrison & Stephenson's
 *                                          (2004) 0.8 t² s with t in centuries
 *                                          from 1820 where that is larger:
 *                                          55 s in AD 1000, 7 minutes in 500
 *                                          BC, an hour in 4000 BC and 5 hours
 *                                          by 13,200 BC
 *   IERS predictions                       grows 0.03 s a year
 *   beyond them                            grows with the square of the years
 *                                          past measurement, reaching 240 s in
 *                                          2400 (the spread between published
 *                                          models there), 11 minutes in 2649
 *                                          and days by AD 17,000
 */
export function deltaTUncertainty(jdUT: number): number {
  if (jdUT <= MEASURED_UNTIL) {
    if (jdUT >= 2435108.5) return 0.001; // 1955-01-01
    if (jdUT >= 2415020.5) return 0.1; // 1900
    if (jdUT >= 2400410.5) return 0.5; // 1860
    if (jdUT >= 2378496.5) return 1; // 1800
    if (jdUT >= 2341972.5) return 5; // 1700
    const t = (jdUT - J1820) / 36525;
    return Math.max(20, 0.8 * t * t);
  }
  const years = (jdUT - MEASURED_UNTIL) / YEAR;
  if (jdUT <= PREDICTED_UNTIL) return 0.001 + 0.03 * years;
  const span = (J2400 - MEASURED_UNTIL) / YEAR;
  return 0.03 * ((PREDICTED_UNTIL - MEASURED_UNTIL) / YEAR) + 240 * (years / span) ** 2;
}
