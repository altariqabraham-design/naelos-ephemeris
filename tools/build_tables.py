"""build_tables.py — Naelos Ephemeris position tables, fitted from NASA JPL DE440.

Writes one binary file per century (1550–2649, all of DE440 but its last
months) holding a Chebyshev series for every body the engine draws positions
from, plus manifest.json recording what was built from what, and the worst
error found.

    python -I tools/build_tables.py --de440 /path/to/de440.bsp --horizons /path/to/spk-folder

The --horizons folder holds one SPK per asteroid from JPL Horizons
(https://ssd.jpl.nasa.gov/api/horizons.api with EPHEM_TYPE=SPK, 1600 to 2500,
the widest span Horizons will write), named chiron.bsp, ceres.bsp, pallas.bsp,
juno.bsp and vesta.bsp.

WHAT IS STORED. Positions as DE440 itself states them: ICRF, kilometres, TDB.
The Sun, the Earth-Moon barycentre and the planets are barycentric; the Moon is
geocentric. Everything an astrologer sees — light-time, aberration, precession,
nutation, the ecliptic of date — is computed at run time by the engine, so the
tables carry no convention that could be wrong.

HOW SMALL. Each series is cut into fixed-length segments and fitted at
Chebyshev nodes. Coefficients are rounded to a step chosen so that the rounding
cannot move a body by more than the target (default 0.01″) as seen from Earth
at its closest approach, then written as zigzag varints. The segment's mean
coefficient is predicted from where the previous segment ended and only the
difference is stored. Every segment is then checked against DE440 on a dense
grid, and a series that misses the target is rebuilt with a finer step.

THE MOON'S VELOCITY. A velocity is the series' derivative, and the derivative
weights coefficient k by about k², so rounding the high orders as coarsely as
the low ones leaves the Moon's velocity good to only a few parts in a million.
That is invisible in a position and ruinous in the Moon's osculating orbit:
the True Node and True Lilith are computed from position AND velocity, and
the apogee magnifies a velocity error eighteen-fold. So a series may name a
KNEE: coefficients above it are rounded to the step divided by 2^s, the
smallest power of two at least (k / knee)². Only the Moon has one.

FORMAT 2 (little-endian). File: "NEPH", u8 version, u8 series count,
f64 start JD, f64 end JD. Per series: u16 id, u8 frame, u8 degree, u8 knee
(0 = none), f64 segment days, f64 first segment start JD, u32 segment count,
f64 step km, u32 payload bytes, then the payload: per segment, per axis x/y/z,
coefficients 1..degree in their own units, then the mean (coefficient 0) — raw
in a file's first segment, otherwise as a difference from its prediction, in
the finest unit (step / 2^max s), so the prediction is exact integer arithmetic.

DE440 is a work of the US Government (NASA/JPL) and is in the public domain.
"""
import argparse, hashlib, json, math, os, struct, sys, time
import numpy as np
from jplephem.spk import SPK

AU = 149597870.7
RS = 206264.80624709636  # arcseconds per radian

# id, name, frame, segment days, degree, closest approach to Earth (km), knee.
# Frame 0 = barycentric, 1 = geocentric, 2 = heliocentric (the asteroids, as
# Horizons gives them). The closest approach converts the angular target into
# a distance budget: Venus at inferior conjunction is the hardest case for the
# Earth-Moon barycentre too. Each series stops at its last whole segment inside
# its source, so the long segments cost the most at the end: Mars to Pluto (512
# days and longer) end the tables in January 2649, ten months before DE440
# does, and Chiron (1024) ends the asteroids in March 2497.
SERIES = [
    (15, "chiron", 2, 1024, 45, 7.4 * AU, 0),
    (17, "ceres", 2, 512, 39, 1.55 * AU, 0),
    (18, "pallas", 2, 512, 42, 1.2 * AU, 0),
    (19, "juno", 2, 256, 24, 0.98 * AU, 0),
    (20, "vesta", 2, 512, 41, 1.13 * AU, 0),
    (0, "sun", 0, 256, 22, 1.47e8, 0),
    (1, "moon", 1, 16, 21, 356000.0, 2),
    (2, "mercury", 0, 64, 49, 0.53 * AU, 0),
    (3, "venus", 0, 256, 32, 0.26 * AU, 0),
    (4, "mars", 0, 512, 41, 0.37 * AU, 0),
    (5, "jupiter", 0, 1024, 16, 3.9 * AU, 0),
    (6, "saturn", 0, 1024, 12, 8.0 * AU, 0),
    (7, "uranus", 0, 2048, 11, 17.3 * AU, 0),
    (8, "neptune", 0, 2048, 9, 28.8 * AU, 0),
    (9, "pluto", 0, 2048, 7, 28.6 * AU, 0),
    (200, "emb", 0, 256, 26, 0.26 * AU, 0),
]


def shifts(deg, knee):
    """s_k: coefficient k is rounded to step / 2^s_k. Zero at and below the
    knee; above it the smallest s with 2^s · knee² ≥ k² (integers, so the
    engine derives the same numbers with no floating point)."""
    out = []
    for k in range(deg + 1):
        s = 0
        if knee and k > knee:
            while (2 ** s) * knee * knee < k * k:
                s += 1
        out.append(s)
    return np.array(out, dtype=np.int64)

# The span is DE440's own, 1550–2650, less a margin at each end for light-time
# and speed differencing; Chiron and the asteroids cover what JPL Horizons
# publishes for them, 1600–2500. Century boundaries (Gregorian, 1 January):
BOUNDARIES = {1600: 2305447.5, 1700: 2341972.5, 1800: 2378496.5, 1900: 2415020.5, 2000: 2451544.5,
              2100: 2488069.5, 2200: 2524593.5, 2300: 2561117.5, 2400: 2597641.5, 2500: 2634166.5,
              2600: 2670690.5}
MARGIN = 10.0  # days inside each source's coverage


def centuries(start, end):
    """(label, from, to) files: 1550 runs from the start to 1600, 2600 from 2600 to the end."""
    edges = [start] + [jd for jd in BOUNDARIES.values() if start < jd < end] + [end]
    labels = [1550] + [y for y, jd in BOUNDARIES.items() if start < jd < end]
    return [(labels[i], edges[i], edges[i + 1]) for i in range(len(labels))]


# Chiron and the four largest main-belt asteroids are not in DE440. Their
# orbits come from JPL Horizons as SPK files (one per body, heliocentric,
# ICRF), read with NASA's SPICE toolkit; see --horizons.
ASTEROIDS = {"chiron": 20002060, "ceres": 20000001, "pallas": 20000002, "juno": 20000003, "vesta": 20000004}
_spice_loaded = False


def asteroid_state(name, t):
    """Heliocentric position (3, N) km and velocity (3, N) km/day from the Horizons SPK."""
    import spiceypy as sp
    t = np.atleast_1d(t)
    pos = np.empty((3, t.size))
    vel = np.empty((3, t.size))
    for i, jd in enumerate(t):
        s, _ = sp.spkgeo(ASTEROIDS[name], (jd - 2451545.0) * 86400.0, "J2000", 10)
        pos[:, i] = s[:3]
        vel[:, i] = s[3:] * 86400.0
    return pos, vel


def body_position(K, name, t):
    """DE440 position (3, N) in km, ICRF, at TDB `t`."""
    if name in ASTEROIDS:
        return asteroid_state(name, t)[0]
    if name == "sun":
        return K[0, 10].compute(t)
    if name == "emb":
        return K[0, 3].compute(t)
    if name == "moon":  # geocentric
        return K[3, 301].compute(t) - K[3, 399].compute(t)
    if name == "mercury":
        return K[0, 1].compute(t) + K[1, 199].compute(t)
    if name == "venus":
        return K[0, 2].compute(t) + K[2, 299].compute(t)
    code = {"mars": 4, "jupiter": 5, "saturn": 6, "uranus": 7, "neptune": 8, "pluto": 9}[name]
    return K[0, code].compute(t)  # system barycentre, as Swiss Ephemeris reports by default


def body_velocity(K, name, t):
    """DE440 velocity (3, N) in km/day."""
    if name in ASTEROIDS:
        return asteroid_state(name, t)[1]
    v = lambda a, b: K[a, b].compute_and_differentiate(t)[1]
    if name == "sun":
        return v(0, 10)
    if name == "emb":
        return v(0, 3)
    if name == "moon":
        return v(3, 301) - v(3, 399)
    if name == "mercury":
        return v(0, 1) + v(1, 199)
    if name == "venus":
        return v(0, 2) + v(2, 299)
    return v(0, {"mars": 4, "jupiter": 5, "saturn": 6, "uranus": 7, "neptune": 8, "pluto": 9}[name])


def cheb_T(deg, x):
    return np.cos(np.outer(np.arange(deg + 1), np.arccos(np.clip(x, -1.0, 1.0))))


def hermite_operator(deg):
    """The linear map from samples to coefficients for a least-squares Chebyshev
    fit that matches value AND derivative exactly at both ends of the segment.

    Minimises the squared error at 2(deg+1) Chebyshev nodes subject to four
    equality constraints, solved once per degree as a KKT system. Returns the
    nodes and the two halves of the operator: c = Pf @ f(nodes) + Pg @ g, where
    g = [f(+1), f(-1), f'(+1), f'(-1)] with derivatives taken in x."""
    n = deg + 1
    m = 2 * n
    x = np.cos(np.pi * (np.arange(m) + 0.5) / m)
    A = cheb_T(deg, x).T  # (m, n)
    k = np.arange(n)
    B = np.vstack([np.ones(n), (-1.0) ** k, k ** 2.0, (-1.0) ** (k + 1) * k ** 2.0])  # (4, n)
    KKT = np.block([[A.T @ A, B.T], [B, np.zeros((4, 4))]])
    inv = np.linalg.inv(KKT)
    return x, inv[:n, :n] @ A.T, inv[:n, n:]


def zigzag(v):
    return (v << 1) ^ (v >> 63)


def put_varint(out, v):
    z = zigzag(int(v))
    while True:
        b = z & 0x7F
        z >>= 7
        if z:
            out.append(b | 0x80)
        else:
            out.append(b)
            return


def fit_series(K, name, deg, seg, t0, nseg):
    """Chebyshev coefficients, shape (nseg, 3, deg+1).

    Each segment equals DE440 in position AND velocity at both of its ends, as
    JPL's own fits do, so neighbours meet with no jump and no kink. Fitting at
    interior nodes only left joins up to the whole fitting error apart, and a
    speed taken across a join read the jump as motion; matching positions
    alone still left the slopes disagreeing there."""
    x, Pf, Pg = hermite_operator(deg)
    starts = t0 + seg * np.arange(nseg)
    times = (starts[:, None] + (x[None, :] + 1.0) * seg / 2.0).ravel()
    f = body_position(K, name, times).reshape(3, nseg, len(x))
    ends = np.concatenate([starts + seg, starts])  # x = +1, then x = -1
    pe = body_position(K, name, ends).reshape(3, 2, nseg)
    ve = body_velocity(K, name, ends).reshape(3, 2, nseg) * (seg / 2.0)  # d/dx = d/dt · dt/dx
    g = np.stack([pe[:, 0], pe[:, 1], ve[:, 0], ve[:, 1]], axis=-1)  # (3, nseg, 4)
    return np.einsum("km,anm->nak", Pf, f) + np.einsum("kq,anq->nak", Pg, g)


def encode(q, sh):
    """q: int64 (nseg, 3, deg+1), coefficient k in units of step / 2^sh[k].
    The mean coefficient goes last, predicted from continuity in the finest
    unit, where every coefficient is an exact integer."""
    out = bytearray()
    nseg, _, n = q.shape
    top = int(sh.max())
    scale = np.array([2 ** (top - int(s)) for s in sh], dtype=np.int64)  # to the finest unit
    alt = np.array([(-1) ** k for k in range(n)], dtype=np.int64)
    prev_end = None
    for s in range(nseg):
        Q = q[s] * scale  # (3, n), finest unit
        for a in range(3):
            row = q[s, a]
            for k in range(1, n):
                put_varint(out, row[k])
            if prev_end is None:
                put_varint(out, row[0])
            else:
                # start of this segment = sum c_k (-1)^k must equal the previous end
                pred = int(prev_end[a]) - int((Q[a, 1:] * alt[1:]).sum())
                put_varint(out, int(Q[a, 0]) - pred)
        prev_end = Q.sum(axis=1)  # value at x = +1, per axis, finest unit
    return bytes(out)


def worst_velocity_rel(K, name, deg, seg, t0, nseg, cq, chunk=2000):
    """Worst relative velocity error of the dequantised series (diagnostic)."""
    m = 2 * (deg + 1)
    xs = np.linspace(-0.999, 0.999, m)
    th = np.arccos(xs)
    U = np.array([k * np.sin(k * th) / np.sin(th) for k in range(deg + 1)])  # d T_k / dx
    worst = 0.0
    for lo in range(0, nseg, chunk):
        hi = min(nseg, lo + chunk)
        starts = t0 + seg * np.arange(lo, hi)
        times = (starts[:, None] + (xs[None, :] + 1.0) * seg / 2.0).ravel()
        truth = body_velocity(K, name, times).reshape(3, hi - lo, m)
        approx = np.einsum("nak,km->anm", cq[lo:hi], U) * (2.0 / seg)
        rel = np.sqrt(((approx - truth) ** 2).sum(axis=0)) / np.sqrt((truth ** 2).sum(axis=0))
        worst = max(worst, float(rel.max()))
    return worst


def worst_error_km(K, name, deg, seg, t0, nseg, cq, chunk=2000):
    """Max distance between the dequantised series and DE440 on a dense grid."""
    m = 4 * (deg + 1)
    xs = np.linspace(-1.0, 1.0, m)
    Tm = cheb_T(deg, xs)  # (n, m)
    worst = 0.0
    for lo in range(0, nseg, chunk):
        hi = min(nseg, lo + chunk)
        starts = t0 + seg * np.arange(lo, hi)
        times = (starts[:, None] + (xs[None, :] + 1.0) * seg / 2.0).ravel()
        truth = body_position(K, name, times).reshape(3, hi - lo, m)
        approx = np.einsum("nak,km->anm", cq[lo:hi], Tm)
        err = np.sqrt(((approx - truth) ** 2).sum(axis=0)).max()
        worst = max(worst, float(err))
    return worst


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--de440", required=True)
    ap.add_argument("--horizons", required=True,
                    help="folder holding chiron.bsp, ceres.bsp, pallas.bsp, juno.bsp, vesta.bsp from JPL Horizons")
    ap.add_argument("--out", default=os.path.join(os.path.dirname(__file__), "..", "data"))
    ap.add_argument("--target", type=float, default=0.01, help="arcseconds as seen from Earth")
    args = ap.parse_args()

    import spiceypy as sp
    asteroid_sources = {}
    for name, spkid in ASTEROIDS.items():
        path = os.path.join(args.horizons, f"{name}.bsp")
        sp.furnsh(path)
        asteroid_sources[name] = {"spk_id": spkid, "sha256": hashlib.sha256(open(path, "rb").read()).hexdigest()}

    sha = hashlib.sha256(open(args.de440, "rb").read()).hexdigest()
    K = SPK.open(args.de440)
    os.makedirs(args.out, exist_ok=True)
    manifest = {
        "engine": "Naelos Ephemeris",
        "format": 2,
        "source": {"name": "JPL DE440", "file": os.path.basename(args.de440), "sha256": sha,
                   "licence": "public domain (US Government work, NASA/JPL)"},
        "asteroids": {"source": "JPL Horizons small-body SPK files, one per body, heliocentric ICRF",
                      "licence": "public domain (US Government work, NASA/JPL)", "files": asteroid_sources},
        "built": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "generator": "packages/naelos-ephemeris/tools/build_tables.py",
        "target_arcsec": args.target,
        "frame": "ICRF, km, TDB; barycentric, except the Moon (geocentric) and the asteroids (heliocentric)",
        "files": [],
        "series": {},
    }

    # ONE SEGMENT GRID PER SERIES FOR ITS WHOLE SPAN. Each series is fitted and
    # quantised once, on a grid from its source's start; a century file then
    # carries the segments that overlap it, so a segment straddling a boundary
    # is byte-identical in both files. (Per-file grids put two different fits
    # either side of every boundary, and a speed taken across one read the jump
    # between them as motion.) Every segment lies inside its source's coverage;
    # the engine's span ends where the shortest of the planetary series does.
    de_start = max(seg.start_jd for seg in K.segments)
    de_end = min(seg.end_jd for seg in K.segments)
    ast = [sp.spkcov(os.path.join(args.horizons, f"{n}.bsp"), i) for n, i in ASTEROIDS.items()]
    ast_start = max(sp.wnfetd(c, 0)[0] for c in ast) / 86400 + 2451545.0
    ast_end = min(sp.wnfetd(c, sp.wncard(c) - 1)[1] for c in ast) / 86400 + 2451545.0
    G = de_start + MARGIN
    origin = {sid: (ast_start + MARGIN if name in ASTEROIDS else G) for (sid, name, *_) in SERIES}
    count = {}
    for (sid, name, frame, seg, deg, dmin, knee) in SERIES:
        stop = (ast_end if name in ASTEROIDS else de_end) - MARGIN
        count[sid] = int(math.floor((stop - origin[sid]) / seg))
    END = min(origin[sid] + count[sid] * seg for (sid, name, frame, seg, *_) in SERIES if name not in ASTEROIDS)
    CENTURIES = centuries(G, END)
    manifest["span"] = {"start_jd_tdb": G, "end_jd_tdb": END,
                        "asteroids": {"start_jd_tdb": origin[17], "end_jd_tdb": min(origin[s] + count[s] * seg
                                      for (s, n, f, seg, *_) in SERIES if n in ASTEROIDS)}}
    print(f"span TDB {G} .. {END}; asteroids {ast_start:.1f} .. {ast_end:.1f}", flush=True)
    fitted = {}
    for (sid, name, frame, seg, deg, dmin, knee) in SERIES:
        nseg = count[sid]
        c = fit_series(K, name, deg, seg, origin[sid], nseg)
        sh = shifts(deg, knee)
        step = args.target / RS * dmin * 1.5 / math.sqrt(deg + 1)
        for attempt in range(6):
            unit = step / (2.0 ** sh)  # per-coefficient step
            q = np.round(c / unit).astype(np.int64)
            err_as = worst_error_km(K, name, deg, seg, origin[sid], nseg, q * unit) / dmin * RS
            if err_as <= args.target:
                break
            step *= 0.8
        else:
            sys.exit(f"{name}: {err_as:.4f}″ misses {args.target}″ even at the finest step")
        vel = worst_velocity_rel(K, name, deg, seg, origin[sid], nseg, q * unit)
        fitted[sid] = (q, step, sh)
        manifest["series"][name] = {"id": sid, "frame": frame, "segment_days": seg, "degree": deg, "knee": knee,
                                    "step_km": step, "bytes": 0, "worst_arcsec": round(err_as, 5),
                                    "worst_velocity_rel": float(f"{vel:.2e}")}
        print(f"{name:8s} {nseg:6d} seg  step {step*1000:9.3f} m  worst {err_as:.4f}″  velocity {vel:.1e}", flush=True)

    for (label, start, end) in CENTURIES:
        blobs = []
        for (sid, name, frame, seg, deg, dmin, knee) in SERIES:
            q, step, sh = fitted[sid]
            o = origin[sid]
            i0 = max(0, int(math.floor((start - o) / seg + 1e-9)))
            i1 = min(len(q), int(math.ceil((end - o) / seg - 1e-9)))
            if i1 <= i0:
                continue  # this series has no data in this file (asteroids before 1600 or after 2500)
            payload = encode(q[i0:i1], sh)
            head = struct.pack("<HBBBddIdI", sid, frame, deg, knee, float(seg), float(o + i0 * seg), i1 - i0, step, len(payload))
            blobs.append(head + payload)
            manifest["series"][name]["bytes"] += len(payload)
        body = b"NEPH" + struct.pack("<BBdd", 2, len(blobs), float(start), float(end)) + b"".join(blobs)
        print(f"ne-{label}.bin  {len(body)/1024:7.1f} KB", flush=True)
        fname = f"ne-{label}.bin"
        with open(os.path.join(args.out, fname), "wb") as f:
            f.write(body)
        manifest["files"].append({"name": fname, "century": label, "start_jd_tdb": start, "end_jd_tdb": end,
                                  "bytes": len(body), "sha256": hashlib.sha256(body).hexdigest()})

    with open(os.path.join(args.out, "manifest.json"), "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)
    total = sum(f["bytes"] for f in manifest["files"])
    print(f"wrote {len(manifest['files'])} files, {total/1024:.0f} KB")


if __name__ == "__main__":
    main()
