"""build_ltp_reference.py — ERFA's long-term precession at a few epochs, for the tests.

Writes test/ltp-reference.json: the bias-precession matrix of Vondrák,
Capitaine & Wallace (2011) as ERFA's ltpb computes it, and the mean obliquity
of date (the angle between ERFA's ltpequ and ltpecl poles), at epochs across
the DE441 tables. test/far.test.ts holds the engine's own implementation to
these, so a typo in the generated coefficients cannot pass.

    python -I tools/build_ltp_reference.py
"""
import json, os
import numpy as np
import erfa

J2000 = 2451545.0
JDS = [-3099000.0, -2000000.0, -1000000.0, 0.0, 1000000.0, 1600000.0, 2100000.0, 2260000.0,
       2720000.0, 3000000.0, 5000000.0, 7999000.0]
rows = []
for jd in JDS:
    epj = 2000.0 + (jd - J2000) / 365.25
    pb = erfa.ltpb(epj)
    eps = float(np.arccos(np.dot(erfa.ltpequ(epj), erfa.ltpecl(epj))))
    rows.append({"jd_tt": jd, "pb": [float(v) for v in pb.ravel()], "epsA": eps})
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "test", "ltp-reference.json")
with open(out, "w", encoding="utf-8", newline="\n") as f:
    json.dump({"source": f"ERFA {erfa.__version__} ltpb, ltpequ, ltpecl (Vondrák, Capitaine & Wallace 2011)", "rows": rows}, f, indent=1)
print(f"wrote {len(rows)} epochs to {out}")
