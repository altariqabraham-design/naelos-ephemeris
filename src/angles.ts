// angles.ts — the one way to wrap an angle into [0, 360) in Naelos Ephemeris.
//
// The engine is its own package and cannot import the app's helper, so this is
// a copy, held to the same cases as every other copy by the app's
// src/tests/angles.test.ts. The two hand-rolled forms it replaces are each
// wrong at one edge: ((d % 360) + 360) % 360 loses low bits on values already
// in range, and r < 0 ? r + 360 : r returns exactly 360 for a tiny negative.

/**
 * Any angle in degrees, wrapped into [0, 360). A value already in range comes
 * back bit-for-bit unchanged; a negative too small to survive adding 360 comes
 * back as 0. NaN and ±Infinity give NaN.
 */
export function norm360(deg: number): number {
  const r = deg % 360;
  if (r >= 0) return r + 0; // + 0 turns -0 into 0
  const wrapped = r + 360;
  return wrapped === 360 ? 0 : wrapped;
}
