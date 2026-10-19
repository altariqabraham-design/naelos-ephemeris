"""build_far_tables.py — Naelos Ephemeris tables beyond DE440, fitted from DE441.

JPL publishes two ephemerides together: DE440, the most accurate for modern
times (1550–2650), and DE441, integrated for a long span (13,200 BC to AD
17,191) and recommended outside DE440's. data/ holds the DE440 tables; this
writes the rest of DE441's span in the same format, to data-far/ — one file
per century before 1550 and after 2649, and a manifest.

    python -I tools/build_far_tables.py --de441 /path/to/de441.bsp

THE SEAMS. Each body's segment grid is anchored on the DE440 tables' own
first and last instants (data/manifest.json, span): the far-past grid ends
exactly where DE440's begins, the far-future grid begins exactly where it
ends, so there is no gap between the sources and no overlap to choose from.
Only the source changes at the seam (DE441 and DE440 agree there to far
better than the 0.01″ the fit is held to).

NO ASTEROIDS. Chiron and the asteroids come from JPL Horizons orbits, which
reach only 1600–2500; the far tables carry the Sun, the Moon, the planets and
the Earth-Moon barycentre.

The fitting, quantising and encoding are build_tables.py's own, imported, so
the far tables are the same format, held to the same target, checked the same
way against their source on a dense grid.

DE441 is a work of the US Government (NASA/JPL) and is in the public domain.
"""
import argparse, hashlib, json, math, os, struct, sys, time
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # python -I leaves the script's folder off the path
from build_tables import (SERIES, ASTEROIDS, MARGIN, RS, shifts, fit_series, encode,  # noqa: E402
                          worst_error_km, worst_velocity_rel)
from jplephem.spk import SPK  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")  # the arcsecond sign, on a Windows console

PLANETARY = [s for s in SERIES if s[1] not in ASTEROIDS]
CHUNK = 20000  # segments fitted at once, to bound memory


def gregorian_jd(year, month=1, day=1):
    """Proleptic Gregorian Julian Day for any year (astronomical numbering),
    as the engine's julday(): Meeus, moved into range by 400-year cycles."""
    n = math.ceil((-4000 - year) / 400) if year < -4000 else 0
    y, m = year + 400 * n, month
    if m <= 2:
        y, m = y - 1, m + 12
    a = math.floor(y / 100)
    b = 2 - a + math.floor(a / 4)
    return math.floor(365.25 * (y + 4716)) + math.floor(30.6001 * (m + 1)) + day + b - 1524.5 - 146097 * n


def name_for(label):
    """ne-1500.bin, ne-0.bin, ne-m0100.bin (100 BC is year −99; 'm' marks the minus)."""
    return f"ne-{label}.bin" if label >= 0 else f"ne-m{-label:04d}.bin"


class Half:
    """One half of DE441, as build_tables.py's K[center, target] expects.

    DE441 is published as two sets of segments, 13,200 BC–1969 and 1969–AD 17,191;
    jplephem keeps only the last segment it reads for each pair, so the half a
    side needs is chosen here: every segment that covers the moment `at`."""

    def __init__(self, K, at):
        self.map = {(s.center, s.target): s for s in K.segments if s.start_jd <= at <= s.end_jd}
        self.start = max(s.start_jd for s in self.map.values())
        self.end = min(s.end_jd for s in self.map.values())

    def __getitem__(self, key):
        return self.map[key]


# THE MOON AT THE SEAMS. DE440 and DE441 share their planets to the millimetre
# where they meet, but not the Moon: DE441's lunar solution differs by 740 m
# at the 1550 seam (0.43″ at perigee) and 98 m at 2649. Left alone, the Moon
# would jump by that at one instant — under anything a chart shows, but more
# than the tables are fitted to, and True Node and True Lilith, which come
# from the Moon's velocity, would jump with it. So within MOON_BLEND days of
# each seam the Moon is moved from DE441 onto DE440: DE441 plus the difference
# at the seam, faded out with a cosine (and the velocity difference over a
# shorter MOON_VBLEND, which keeps the extra displacement under 300 m). At the
# seam position and velocity are DE440's exactly; nowhere does the Moon stray
# from DE441 by more than the two already disagree.
MOON_BLEND = 16.0   # days: one of the Moon's segments
MOON_VBLEND = 4.0   # days


class _Blended:
    """A DE441 segment (Earth-Moon barycentre → Moon) with the seam blend added."""

    def __init__(self, seg, seam, dp, dv):
        self.seg, self.seam, self.dp, self.dv = seg, seam, dp, dv

    def _terms(self, t):
        x = np.abs(t - self.seam)
        sgn = np.sign(t - self.seam)
        inp, inv = x < MOON_BLEND, x < MOON_VBLEND
        c = np.where(inp, 0.5 * (1 + np.cos(np.pi * x / MOON_BLEND)), 0.0)
        dc = np.where(inp, -0.5 * np.pi / MOON_BLEND * np.sin(np.pi * x / MOON_BLEND) * sgn, 0.0)
        q = np.where(inv, 0.5 * (1 + np.cos(np.pi * x / MOON_VBLEND)), 0.0)
        dq = np.where(inv, -0.5 * np.pi / MOON_VBLEND * np.sin(np.pi * x / MOON_VBLEND) * sgn, 0.0)
        g = (t - self.seam) * q
        dg = q + (t - self.seam) * dq
        pos = np.outer(self.dp, c) + np.outer(self.dv, g)
        vel = np.outer(self.dp, dc) + np.outer(self.dv, dg)
        return pos, vel

    def compute(self, t):
        t = np.asarray(t, dtype=float)
        pos, _ = self._terms(np.atleast_1d(t))
        return self.seg.compute(t) + pos.reshape(self.seg.compute(t).shape)

    def compute_and_differentiate(self, t):
        t = np.asarray(t, dtype=float)
        p, v = self.seg.compute_and_differentiate(t)
        dp, dv = self._terms(np.atleast_1d(t))
        return p + dp.reshape(p.shape), v + dv.reshape(v.shape)


class SeamBlended:
    """A Half whose Moon meets DE440's at the seam (see MOON_BLEND)."""

    def __init__(self, half, de440, seam):
        moon441 = lambda k: (k[3, 301].compute(seam) - k[3, 399].compute(seam),
                             k[3, 301].compute_and_differentiate(seam)[1] - k[3, 399].compute_and_differentiate(seam)[1])
        p440, v440 = moon441(de440)
        p441, v441 = moon441(half)
        self.dp, self.dv = p440 - p441, v440 - v441
        self.half = half
        self.start, self.end = half.start, half.end
        self.moon = _Blended(half[3, 301], seam, self.dp, self.dv)

    def __getitem__(self, key):
        return self.moon if key == (3, 301) else self.half[key]


def fit_chunked(K, name, deg, seg, t0, nseg):
    parts = [fit_series(K, name, deg, seg, t0 + lo * seg, min(CHUNK, nseg - lo)) for lo in range(0, nseg, CHUNK)]
    return np.concatenate(parts, axis=0)


def build_side(K, side, origin, count, target, report):
    """Fit and quantise every planetary series on one side of DE440."""
    fitted = {}
    for (sid, name, frame, seg, deg, dmin, knee) in PLANETARY:
        o, nseg = origin[sid], count[sid]
        t0 = time.time()
        c = fit_chunked(K, name, deg, seg, o, nseg)
        sh = shifts(deg, knee)
        step = target / RS * dmin * 1.5 / math.sqrt(deg + 1)
        for attempt in range(6):
            unit = step / (2.0 ** sh)
            q = np.round(c / unit).astype(np.int64)
            err_as = worst_error_km(K, name, deg, seg, o, nseg, q * unit) / dmin * RS
            if err_as <= target:
                break
            step *= 0.8
        else:
            sys.exit(f"{side} {name}: {err_as:.4f}″ misses {target}″ even at the finest step")
        vel = worst_velocity_rel(K, name, deg, seg, o, nseg, q * unit)
        fitted[sid] = (q, step, sh)
        report[name] = {"id": sid, "frame": frame, "segment_days": seg, "degree": deg, "knee": knee,
                        "step_km": step, "bytes": 0, "worst_arcsec": round(err_as, 5),
                        "worst_velocity_rel": float(f"{vel:.2e}")}
        print(f"{side:6s} {name:8s} {nseg:7d} seg  step {step*1000:9.3f} m  worst {err_as:.4f}″  "
              f"velocity {vel:.1e}  {time.time()-t0:6.1f} s", flush=True)
    return fitted


def write_files(files, fitted, origin, out, manifest, report):
    for (label, start, end) in files:
        blobs = []
        for (sid, name, frame, seg, deg, dmin, knee) in PLANETARY:
            q, step, sh = fitted[sid]
            o = origin[sid]
            i0 = max(0, int(math.floor((start - o) / seg + 1e-9)))
            i1 = min(len(q), int(math.ceil((end - o) / seg - 1e-9)))
            if i1 <= i0:
                continue
            payload = encode(q[i0:i1], sh)
            head = struct.pack("<HBBBddIdI", sid, frame, deg, knee, float(seg), float(o + i0 * seg), i1 - i0, step, len(payload))
            blobs.append(head + payload)
            report[name]["bytes"] += len(payload)
        body = b"NEPH" + struct.pack("<BBdd", 2, len(blobs), float(start), float(end)) + b"".join(blobs)
        fname = name_for(label)
        with open(os.path.join(out, fname), "wb") as f:
            f.write(body)
        manifest["files"].append({"name": fname, "century": label, "start_jd_tdb": start, "end_jd_tdb": end,
                                  "bytes": len(body), "sha256": hashlib.sha256(body).hexdigest()})
        print(f"{fname:16s} {len(body)/1024:7.1f} KB", flush=True)


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser()
    ap.add_argument("--de441", required=True)
    ap.add_argument("--de440", required=True, help="the DE440 file the near tables were built from: the Moon is joined to it at the seams")
    ap.add_argument("--near", default=os.path.join(here, "..", "data", "manifest.json"),
                    help="the DE440 tables' manifest, whose span the far tables meet")
    ap.add_argument("--out", default=os.path.join(here, "..", "data-far"))
    ap.add_argument("--target", type=float, default=0.01, help="arcseconds as seen from Earth")
    ap.add_argument("--side", choices=["past", "future", "both"], default="both")
    ap.add_argument("--trial-years", type=float, default=0, help="time a short span next to each seam; writes nothing")
    args = ap.parse_args()

    near = json.load(open(args.near, encoding="utf-8"))
    if near["target_arcsec"] != args.target:
        sys.exit(f"the DE440 tables were built to {near['target_arcsec']}″; build the far tables to the same")
    G, END = near["span"]["start_jd_tdb"], near["span"]["end_jd_tdb"]

    sha = hashlib.sha256()
    with open(args.de441, "rb") as f:
        for block in iter(lambda: f.read(1 << 24), b""):
            sha.update(block)
    spk = SPK.open(args.de441)
    de440 = SPK.open(args.de440)
    if hashlib.sha256(open(args.de440, "rb").read()).hexdigest() != near["source"]["sha256"]:
        sys.exit("--de440 is not the file the near tables were built from")
    # The half before 1969 and the half after, each with its Moon joined to DE440's at its seam.
    PAST, FUTURE = SeamBlended(Half(spk, G), de440, G), SeamBlended(Half(spk, END), de440, END)
    de_start, de_end = PAST.start, FUTURE.end
    for side, h in (("1550", PAST), ("2649", FUTURE)):
        print(f"Moon at the {side} seam: DE440 − DE441 = {np.linalg.norm(h.dp)*1000:.1f} m, "
              f"{np.linalg.norm(h.dv)*1000:.1f} m/day; blended over {MOON_BLEND:.0f} days", flush=True)
    os.makedirs(args.out, exist_ok=True)

    manifest = {
        "engine": "Naelos Ephemeris",
        "format": 2,
        "source": {"name": "JPL DE441", "file": os.path.basename(args.de441), "sha256": sha.hexdigest(),
                   "licence": "public domain (US Government work, NASA/JPL)"},
        "built": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "generator": "packages/naelos-ephemeris/tools/build_far_tables.py",
        "target_arcsec": args.target,
        "frame": "ICRF, km, TDB; barycentric, except the Moon (geocentric)",
        "moon_seam_blend": {
            "days": MOON_BLEND, "velocity_days": MOON_VBLEND,
            "1550": {"position_m": round(float(np.linalg.norm(PAST.dp)) * 1000, 1), "velocity_m_per_day": round(float(np.linalg.norm(PAST.dv)) * 1000, 1)},
            "2649": {"position_m": round(float(np.linalg.norm(FUTURE.dp)) * 1000, 1), "velocity_m_per_day": round(float(np.linalg.norm(FUTURE.dv)) * 1000, 1)},
        },
        "meets": {"name": near["source"]["name"], "sha256": near["source"]["sha256"],
                  "start_jd_tdb": G, "end_jd_tdb": END},
        "files": [],
        "series": {"past": {}, "future": {}},
        "span": {},
    }
    print(f"DE441 {de_start} .. {de_end}; DE440 tables {G} .. {END}", flush=True)

    if args.trial_years:
        de_start = max(de_start, G - args.trial_years * 365.25)
        de_end = min(de_end, END + args.trial_years * 365.25)
        args.out = os.path.join(args.out, "_trial")
        os.makedirs(args.out, exist_ok=True)

    if args.side in ("past", "both"):
        # Each grid ENDS at G: whole segments back from the seam, as far as DE441 goes.
        count = {sid: int(math.floor((G - (de_start + MARGIN)) / seg)) for (sid, _n, _f, seg, *_) in PLANETARY}
        origin = {sid: G - count[sid] * seg for (sid, _n, _f, seg, *_) in PLANETARY}
        start = max(origin.values())  # where every series has begun
        fitted = build_side(PAST, "past", origin, count, args.target, manifest["series"]["past"])
        y0 = 100 * math.floor((2000 + (start - 2451544.5) / 365.2425) / 100)
        edges = [(y, gregorian_jd(y)) for y in range(y0 + 100, 1600, 100) if start < gregorian_jd(y) < G]
        labels = [y0] + [y for y, _ in edges]
        bounds = [start] + [jd for _, jd in edges] + [G]
        files = [(labels[i], bounds[i], bounds[i + 1]) for i in range(len(labels))]
        write_files(files, fitted, origin, args.out, manifest, manifest["series"]["past"])
        manifest["span"]["past"] = {"start_jd_tdb": start, "end_jd_tdb": G}
        del fitted

    if args.side in ("future", "both"):
        # Each grid BEGINS at END: whole segments forward from the seam.
        count = {sid: int(math.floor((de_end - MARGIN - END) / seg)) for (sid, _n, _f, seg, *_) in PLANETARY}
        origin = {sid: END for (sid, *_) in PLANETARY}
        stop = min(END + count[sid] * seg for (sid, _n, _f, seg, *_) in PLANETARY)
        fitted = build_side(FUTURE, "future", origin, count, args.target, manifest["series"]["future"])
        y_end = 2000 + (END - 2451544.5) / 365.2425
        edges = [(y, gregorian_jd(y)) for y in range(100 * math.ceil(y_end / 100), 17300, 100) if END < gregorian_jd(y) < stop]
        labels = [math.floor(y_end)] + [y for y, _ in edges]
        bounds = [END] + [jd for _, jd in edges] + [stop]
        files = [(labels[i], bounds[i], bounds[i + 1]) for i in range(len(labels))]
        write_files(files, fitted, origin, args.out, manifest, manifest["series"]["future"])
        manifest["span"]["future"] = {"start_jd_tdb": END, "end_jd_tdb": stop}

    manifest["files"].sort(key=lambda f: f["start_jd_tdb"])
    with open(os.path.join(args.out, "manifest.json"), "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)
    total = sum(f["bytes"] for f in manifest["files"])
    print(f"wrote {len(manifest['files'])} files, {total/1024/1024:.1f} MB")


if __name__ == "__main__":
    main()
