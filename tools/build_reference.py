"""build_reference.py — the answers Naelos Ephemeris is tested against.

Computes, with Skyfield and DE440, the apparent places the engine must
reproduce, at 300 instants spread over the whole span of the tables (TT,
1550–2649; the asteroids where they have tables, 1600–2497):

  geo    apparent geocentric longitude/latitude, true ecliptic of date, and the
         longitude speed (central difference over ±1 minute)
  equ    apparent right ascension/declination, true equator of date
  helio  heliocentric longitude/latitude, true ecliptic of date, with light-time
         from the body to the Sun (the convention Swiss Ephemeris uses)
  gast   Greenwich apparent sidereal time, from ERFA's gst06a (IAU 2006/2000A),
         at UT = TT − 69.184 s (any fixed ΔT serves: both sides use the same)
  lunar  the true node and true apogee from DE440's exact lunar state, and the
         mean node and mean apogee from the published ELP2000-82B elements

    python -I tools/build_reference.py --de440 /path/to/de440.bsp --horizons /path/to/spk-folder

Skyfield applies IAU 2000A nutation and the engine IAU 2000B's terms; the two
differ by up to 2.6 milliarcseconds, well inside the engine's 0.01″ table
budget. Everything here comes from JPL data and ERFA and is safe to publish
with the engine.
"""
import argparse, gzip, json, os
import numpy as np
import erfa
from skyfield.api import load
from skyfield.framelib import ecliptic_frame

ap = argparse.ArgumentParser()
ap.add_argument("--de440", required=True)
ap.add_argument("--horizons", required=True, help="folder with the asteroid SPKs (see build_tables.py)")
ap.add_argument("--out", default=os.path.join(os.path.dirname(__file__), "..", "test", "reference.json.gz"))
args = ap.parse_args()

eph = load(args.de440)
ts = load.timescale(builtin=True)
# The engine's own span, from the tables' manifest.
span = json.load(open(os.path.join(os.path.dirname(__file__), "..", "data", "manifest.json"), encoding="utf-8"))["span"]
LO, HI = span["start_jd_tdb"], span["end_jd_tdb"]
A_LO, A_HI = span["asteroids"]["start_jd_tdb"], span["asteroids"]["end_jd_tdb"]
rng = np.random.default_rng(20261017)
jd = np.sort(np.concatenate([
    rng.uniform(LO + 1, HI - 1, 294),
    [LO + 0.5, HI - 0.5, 2451545.0, 2448871.946, 2415020.5, 2488069.5],
]))
t, tm, tp = ts.tt_jd(jd), ts.tt_jd(jd - 1 / 1440), ts.tt_jd(jd + 1 / 1440)
earth, sun = eph["earth"], eph["sun"]
names = {0: "sun", 1: "moon", 2: "mercury", 3: "venus", 4: "mars barycenter", 5: "jupiter barycenter",
         6: "saturn barycenter", 7: "uranus barycenter", 8: "neptune barycenter", 9: "pluto barycenter"}
r9 = lambda a: [None if v is None else round(float(v), 9) for v in a]

out = {"source": "JPL DE440 via Skyfield; see tools/build_reference.py", "jd_tt": r9(jd),
       "geo": {}, "equ": {}, "helio": {}}
for bid, name in names.items():
    app = earth.at(t).observe(eph[name]).apparent()
    lat, lon, _ = app.frame_latlon(ecliptic_frame)
    _, lm, _ = earth.at(tm).observe(eph[name]).apparent().frame_latlon(ecliptic_frame)
    _, lp, _ = earth.at(tp).observe(eph[name]).apparent().frame_latlon(ecliptic_frame)
    speed = ((lp.degrees - lm.degrees + 180) % 360 - 180) * 720
    ra, dec, _ = app.radec(epoch="date")
    out["geo"][bid] = {"lon": r9(lon.degrees), "lat": r9(lat.degrees), "speed": r9(speed)}
    out["equ"][bid] = {"ra": r9(ra.hours * 15), "dec": r9(dec.degrees)}
    if bid != 0:
        hlat, hlon, _ = sun.at(t).observe(eph[name]).frame_latlon(ecliptic_frame)
        out["helio"][bid] = {"lon": r9(hlon.degrees), "lat": r9(hlat.degrees)}
hlat, hlon, _ = sun.at(t).observe(earth).frame_latlon(ecliptic_frame)
out["helio"][14] = {"lon": r9(hlon.degrees), "lat": r9(hlat.degrees)}

ut = jd - 69.184 / 86400
out["gast"] = {"jd_ut": [round(float(v), 12) for v in ut],
               "deg": r9([np.degrees(erfa.gst06a(u, 0.0, j, 0.0)) for u, j in zip(ut, jd)])}

# The Moon's nodes and apogee. True points from DE440's exact lunar state
# (jplephem's analytic velocity), osculating elements in the inertial frame,
# then IAU 2006/2000A to the true ecliptic of date; mean points from the
# ELP2000-82B mean elements (Meeus ch. 47) plus 2000A nutation.
from jplephem.spk import SPK
K = SPK.open(args.de440)
GM = (398600.435507 + 4902.800118) * 86400.0**2
lunar = {"true_node": [], "true_apogee_lon": [], "true_apogee_lat": [], "mean_node": [],
         "mean_apogee_lon": [], "mean_apogee_lat": []}
for j in jd:
    p1, v1 = K[3, 301].compute_and_differentiate(j)
    p2, v2 = K[3, 399].compute_and_differentiate(j)
    r, v = p1 - p2, v1 - v2
    h = np.cross(r, v)
    e = np.cross(v, h) / GM - r / np.linalg.norm(r)
    dpsi, deps = erfa.nut06a(j, 0.0)
    eps = erfa.obl06(j, 0.0) + deps
    c, s = np.cos(eps), np.sin(eps)
    M = np.array([[1, 0, 0], [0, c, s], [0, -s, c]]) @ erfa.pnm06a(j, 0.0)  # ICRF → ecliptic of date
    n = np.cross([0, 0, 1], M @ h)
    lunar["true_node"].append(np.degrees(np.arctan2(n[1], n[0])) % 360)
    a = M @ (-e)
    lunar["true_apogee_lon"].append(np.degrees(np.arctan2(a[1], a[0])) % 360)
    lunar["true_apogee_lat"].append(np.degrees(np.arcsin(a[2] / np.linalg.norm(a))))
    T = (j - 2451545.0) / 36525
    om = 125.0445479 - 1934.1362891 * T + 0.0020754 * T**2 + T**3 / 467441 - T**4 / 60616000
    per = 83.3532465 + 4069.0137287 * T - 0.0103200 * T**2 - T**3 / 80053 + T**4 / 18999000
    lunar["mean_node"].append((om + np.degrees(dpsi)) % 360)
    u = np.radians(per + 180 - om)
    inc = np.radians(5.1453964)
    lunar["mean_apogee_lon"].append((om + np.degrees(np.arctan2(np.sin(u) * np.cos(inc), np.cos(u)) + dpsi)) % 360)
    lunar["mean_apogee_lat"].append(np.degrees(np.arcsin(np.sin(u) * np.sin(inc))))
out["lunar"] = {k: r9(v) for k, v in lunar.items()}

# Chiron and the asteroids, through NASA's SPICE toolkit instead of Skyfield
# (which cannot read Horizons' SPK type 21): SPICE for light-time ('LT'), then
# ERFA for the rest of the apparent place — the Sun's light deflection (ld)
# and relativistic aberration (ab) — and the frame of date (2006/2000A).
# Heliocentric places take light-time from the Sun and nothing else.
import spiceypy as sp
sp.furnsh(args.de440)
AST = {15: ("chiron", 20002060), 17: ("ceres", 20000001), 18: ("pallas", 20000002), 19: ("juno", 20000003), 20: ("vesta", 20000004)}
for bid, (name, spkid) in AST.items():
    sp.furnsh(os.path.join(args.horizons, f"{name}.bsp"))

def ecl_of_date(v, j):
    dpsi, deps = erfa.nut06a(j, 0.0)
    eps = erfa.obl06(j, 0.0) + deps
    c, s = np.cos(eps), np.sin(eps)
    w = np.array([[1, 0, 0], [0, c, s], [0, -s, c]]) @ erfa.pnm06a(j, 0.0) @ np.asarray(v)
    return np.degrees(np.arctan2(w[1], w[0])) % 360, np.degrees(np.arctan2(w[2], np.hypot(w[0], w[1])))

AU_KM, C_KMS = 149597870.7, 299792.458
for bid, (name, spkid) in AST.items():
    g, h = ([], []), ([], [])
    for j in jd:
        if not (A_LO + 1 < j < A_HI - 1):
            # Outside what Horizons publishes for the asteroids: no reference.
            for arr in (*g, *h):
                arr.append(None)
            continue
        et = (j - 2451545.0) * 86400.0  # TDB seconds; TT and TDB differ by under 2 ms
        p = np.array(sp.spkezr(str(spkid), et, "J2000", "LT", "399")[0][:3])
        earth_state = np.array(sp.spkezr("399", et, "J2000", "NONE", "0")[0])
        sun_pos = np.array(sp.spkezr("10", et, "J2000", "NONE", "0")[0][:3])
        e = earth_state[:3] - sun_pos
        em = np.linalg.norm(e) / AU_KM
        # ERFA's ld with the body's own direction from the Sun (q). ldsun would
        # take q = p, which holds for stars, not for a body a few AU away.
        q = earth_state[:3] + p - sun_pos
        pdir = erfa.ld(1.0, p / np.linalg.norm(p), q / np.linalg.norm(q), e / np.linalg.norm(e), em, 1e-9)
        v = earth_state[3:] / C_KMS
        app = erfa.ab(pdir, v, em, np.sqrt(1 - v @ v))
        lon, lat = ecl_of_date(app, j)
        g[0].append(lon); g[1].append(lat)
        lon, lat = ecl_of_date(sp.spkezr(str(spkid), et, "J2000", "LT", "10")[0][:3], j)
        h[0].append(lon); h[1].append(lat)
    out["geo"][bid] = {"lon": r9(g[0]), "lat": r9(g[1])}
    out["helio"][bid] = {"lon": r9(h[0]), "lat": r9(h[1])}

os.makedirs(os.path.dirname(args.out), exist_ok=True)
with gzip.open(args.out, "wt", encoding="utf-8") as f:
    json.dump(out, f, separators=(",", ":"))
print(f"wrote {args.out}: {len(jd)} instants")
