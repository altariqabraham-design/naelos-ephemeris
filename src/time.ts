// time.ts — calendar dates, Julian Days and the time scales between them.
//
// Years are ASTRONOMICAL: year 0 is 1 BC, year −1 is 2 BC, and so on. The
// tables reach 13,200 BC (year −13199), far before Julian Day 0 (4713 BC), so
// the day count must hold for negative Julian Days too — Meeus's formula as
// printed does not (its floor() assumes y + 4716 > 0). Both calendars repeat
// exactly — the Gregorian every 400 years (146,097 days), the Julian every 4
// (1,461 days) — so a far year is moved into range by whole cycles, counted,
// and moved back.

const GREGORIAN_CYCLE = { years: 400, days: 146097 };
const JULIAN_CYCLE = { years: 4, days: 1461 };

function meeus(year: number, month: number, day: number, hour: number, gregorian: boolean): number {
  let y = year, m = month;
  if (m <= 2) {
    y -= 1;
    m += 12;
  }
  let b = 0;
  if (gregorian) {
    const a = Math.floor(y / 100);
    b = 2 - a + Math.floor(a / 4);
  }
  return Math.floor(365.25 * (y + 4716)) + Math.floor(30.6001 * (m + 1)) + day + b - 1524.5 + hour / 24;
}

function shifted(year: number, month: number, day: number, hour: number, gregorian: boolean): number {
  const cycle = gregorian ? GREGORIAN_CYCLE : JULIAN_CYCLE;
  // Bring the year above −4000, where every floor in Meeus's formula is of a
  // positive number, then take the whole cycles back off.
  const n = year < -4000 ? Math.ceil((-4000 - year) / cycle.years) : 0;
  return meeus(year + n * cycle.years, month, day, hour, gregorian) - n * cycle.days;
}

/**
 * Julian Day of a proleptic GREGORIAN date and decimal hour, for any year.
 * The same day numbering as Swiss Ephemeris's swe_julday with the Gregorian
 * flag, which is how Naelos has always called it.
 */
export function julday(year: number, month: number, day: number, hour = 0): number {
  return shifted(year, month, day, hour, true);
}

/** Julian Day of a proleptic JULIAN-calendar date and decimal hour, for any year. */
export function juldayJulian(year: number, month: number, day: number, hour = 0): number {
  return shifted(year, month, day, hour, false);
}

/** The first Gregorian day, 15 October 1582 (the day after Julian 4 October). */
export const GREGORIAN_REFORM_JD = 2299160.5;

/**
 * Julian Day of a date as history records it: the Julian calendar before 15
 * October 1582, the Gregorian from then on — the convention of astronomers and
 * historians. The ten days between (5–14 October 1582) never happened, and
 * are refused.
 */
export function juldayHistorical(year: number, month: number, day: number, hour = 0): number {
  if (year === 1582 && month === 10 && day >= 5 && day <= 14) {
    throw new RangeError("5–14 October 1582 do not exist: the Gregorian reform skipped them");
  }
  const g = julday(year, month, day, 0);
  return g >= GREGORIAN_REFORM_JD ? g + hour / 24 : juldayJulian(year, month, day, hour);
}

export interface CalendarDate {
  /** Astronomical year: 0 is 1 BC. */
  year: number;
  month: number;
  day: number;
  /** Decimal hours. */
  hour: number;
  calendar: "gregorian" | "julian";
}

/**
 * The calendar date of a Julian Day — the inverse of `juldayHistorical`
 * (Julian before the reform, Gregorian after) unless a calendar is asked for.
 * Meeus ch. 7, made valid for negative Julian Days by whole calendar cycles.
 */
export function revjul(jd: number, calendar?: "gregorian" | "julian"): CalendarDate {
  const cal = calendar ?? (jd >= GREGORIAN_REFORM_JD ? "gregorian" : "julian");
  const cycle = cal === "gregorian" ? GREGORIAN_CYCLE : JULIAN_CYCLE;
  const n = jd < 0 ? Math.ceil(-jd / cycle.days) + 1 : 0;
  const j = jd + 0.5 + n * cycle.days;
  const z = Math.floor(j);
  const f = j - z;
  let a = z;
  if (cal === "gregorian") {
    const alpha = Math.floor((z - 1867216.25) / 36524.25);
    a = z + 1 + alpha - Math.floor(alpha / 4);
  }
  const b = a + 1524;
  const c = Math.floor((b - 122.1) / 365.25);
  const d = Math.floor(365.25 * c);
  const e = Math.floor((b - d) / 30.6001);
  const day = b - d - Math.floor(30.6001 * e);
  const month = e < 14 ? e - 1 : e - 13;
  const year = (month > 2 ? c - 4716 : c - 4715) - n * cycle.years;
  return { year, month, day, hour: f * 24, calendar: cal };
}

/** "13200 BC", "AD 1066", "2026": an astronomical year as people write it. */
export function formatYear(year: number): string {
  if (year <= 0) return `${1 - year} BC`;
  return year < 1000 ? `AD ${year}` : String(year);
}

/**
 * TDB − TT in days. The periodic difference between the two (under 2 ms) from
 * the leading terms of Fairhead & Bretagnon; DE440's argument is TDB.
 */
export function tdbMinusTT(jdTT: number): number {
  const g = ((357.53 + 0.98560028 * (jdTT - 2451545.0)) * Math.PI) / 180;
  return (0.001657 * Math.sin(g) + 0.000014 * Math.sin(2 * g)) / 86400;
}
