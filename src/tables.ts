// tables.ts — reads the century files written by tools/build_tables.py and
// evaluates their Chebyshev series.
//
// A file holds one series per body: fixed-length segments, each a Chebyshev
// series per axis (ICRF km, TDB). Coefficients arrive as zigzag varints, the
// mean coefficient as a difference from where the previous segment ended. They
// are decoded once, on load, into plain Float64Arrays.
//
// Format 2 lets a series round its high-order coefficients more finely than
// its low ones (the Moon does, so its velocity is good enough to derive the
// True Node and True Lilith from). Coefficient k is in units of step / 2^s_k;
// the continuity prediction runs in the finest unit, so it stays exact.

export interface Series {
  id: number;
  /** 0 = barycentric, 1 = geocentric (the Moon), 2 = heliocentric (the asteroids). */
  frame: number;
  degree: number;
  segmentDays: number;
  /** TDB Julian Day where segment 0 starts. */
  t0: number;
  segments: number;
  /** segments × 3 axes × (degree + 1), km. */
  coeffs: Float64Array;
}

export interface TableFile {
  /** TDB Julian Days this file answers for: [start, end). */
  start: number;
  end: number;
  series: Map<number, Series>;
}

const MAGIC = 0x4850454e; // "NEPH" read little-endian

/** Parse one century file. Throws on anything that is not format 2. */
export function parseTable(input: ArrayBuffer | Uint8Array): TableFile {
  const bytes = input instanceof Uint8Array ? input : new Uint8Array(input);
  const view = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
  if (view.getUint32(0, true) !== MAGIC) throw new Error("not a Naelos Ephemeris table");
  const version = view.getUint8(4);
  if (version !== 2) throw new Error(`unsupported table format ${version}`);
  const count = view.getUint8(5);
  const start = view.getFloat64(6, true);
  const end = view.getFloat64(14, true);
  let pos = 22;
  const series = new Map<number, Series>();
  for (let s = 0; s < count; s++) {
    const id = view.getUint16(pos, true);
    const frame = view.getUint8(pos + 2);
    const degree = view.getUint8(pos + 3);
    const knee = view.getUint8(pos + 4);
    const segmentDays = view.getFloat64(pos + 5, true);
    const t0 = view.getFloat64(pos + 13, true);
    const segments = view.getUint32(pos + 21, true);
    const step = view.getFloat64(pos + 25, true);
    const length = view.getUint32(pos + 33, true);
    pos += 37;
    const coeffs = decode(bytes, pos, pos + length, segments, degree, knee, step);
    pos += length;
    series.set(id, { id, frame, degree, segmentDays, t0, segments, coeffs });
  }
  return { start, end, series };
}

/**
 * s_k for each coefficient: zero at and below the knee, above it the smallest s
 * with 2^s · knee² ≥ k². Integers only, so it matches the builder exactly.
 */
export function shifts(degree: number, knee: number): number[] {
  const out: number[] = [];
  for (let k = 0; k <= degree; k++) {
    let s = 0;
    if (knee && k > knee) while (2 ** s * knee * knee < k * k) s++;
    out.push(s);
  }
  return out;
}

function readZigzag(bytes: Uint8Array, at: { pos: number }): number {
  let value = 0, scale = 1, b: number;
  do {
    b = bytes[at.pos++];
    value += (b & 0x7f) * scale;
    scale *= 128;
  } while (b & 0x80);
  return value % 2 === 0 ? value / 2 : -(value + 1) / 2;
}

function decode(bytes: Uint8Array, pos: number, stop: number, segments: number, degree: number, knee: number, step: number): Float64Array {
  const n = degree + 1;
  const sh = shifts(degree, knee);
  const top = Math.max(...sh);
  const toFinest = sh.map((s) => 2 ** (top - s)); // exact powers of two
  const unit = sh.map((s) => step / 2 ** s);
  const out = new Float64Array(segments * 3 * n);
  const q = new Float64Array(n); // exact integers, well inside 2^53
  const prevEnd = [0, 0, 0]; // finest unit
  const nextEnd = [0, 0, 0];
  const at = { pos };
  for (let s = 0; s < segments; s++) {
    for (let a = 0; a < 3; a++) {
      let alternating = 0; // sum over k >= 1 of Q_k (-1)^k, finest unit
      let tail = 0; // sum over k >= 1 of Q_k, finest unit
      for (let k = 1; k < n; k++) {
        q[k] = readZigzag(bytes, at);
        const Q = q[k] * toFinest[k];
        tail += Q;
        alternating += k % 2 === 0 ? Q : -Q;
      }
      const d = readZigzag(bytes, at);
      // The first segment stores its mean outright; the rest store it as a
      // difference from continuity with the previous segment's end.
      const Q0 = s === 0 ? d * toFinest[0] : prevEnd[a] - alternating + d;
      q[0] = Q0 / toFinest[0];
      nextEnd[a] = Q0 + tail;
      const base = (s * 3 + a) * n;
      for (let k = 0; k < n; k++) out[base + k] = q[k] * unit[k];
    }
    prevEnd[0] = nextEnd[0];
    prevEnd[1] = nextEnd[1];
    prevEnd[2] = nextEnd[2];
  }
  if (at.pos !== stop) throw new Error("table payload length mismatch");
  return out;
}

/** Whether a series has data at a TDB Julian Day (its first and last instants included). */
export function covers(series: Series, t: number): boolean {
  return t >= series.t0 - 1e-9 && t <= series.t0 + series.segments * series.segmentDays + 1e-9;
}

/**
 * Position (km) and, when `vel` is given, velocity (km/day) of a series at a
 * TDB Julian Day. The caller has already chosen the file that covers `t`.
 */
export function evaluate(series: Series, t: number, pos: number[], vel?: number[]): void {
  const { degree, segmentDays, t0, segments, coeffs } = series;
  // A series answers only inside its own data: an asteroid outside the span
  // JPL Horizons publishes is refused, never extrapolated.
  if (!covers(series, t)) {
    throw new RangeError(`series ${series.id} covers TDB JD ${t0} to ${t0 + segments * segmentDays}, not ${t}`);
  }
  let i = Math.floor((t - t0) / segmentDays);
  if (i < 0) i = 0;
  else if (i >= segments) i = segments - 1;
  const x = (2 * (t - (t0 + i * segmentDays))) / segmentDays - 1;
  const n = degree + 1;
  for (let a = 0; a < 3; a++) {
    const base = (i * 3 + a) * n;
    // Clenshaw for the value.
    let b1 = 0, b2 = 0;
    for (let k = degree; k >= 1; k--) {
      const b0 = 2 * x * b1 - b2 + coeffs[base + k];
      b2 = b1;
      b1 = b0;
    }
    pos[a] = x * b1 - b2 + coeffs[base];
    if (vel) {
      // d/dx of sum c_k T_k = sum k c_k U_{k-1}(x).
      let u0 = 1, u1 = 2 * x, d = coeffs[base + 1];
      for (let k = 2; k <= degree; k++) {
        d += k * coeffs[base + k] * u1;
        const u2 = 2 * x * u1 - u0;
        u0 = u1;
        u1 = u2;
      }
      vel[a] = (d * 2) / segmentDays;
    }
  }
}
