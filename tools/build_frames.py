"""build_frames.py — the precession and nutation constants, generated, never typed.

Writes src/generated/frames.ts:
  * IAU 2006 precession as Fukushima-Williams angles (gamma-bar, phi-bar,
    psi-bar) and the IAU 2006 mean obliquity: degree-5 polynomials in Julian
    centuries of TT, as published, checked here to reproduce ERFA's pfw06 and
    obl06 to 1e-12 rad across 1550–2650.
  * Nutation: the 77 luni-solar terms and planetary offsets of IAU 2000B
    (McCarthy & Luzum 2003), evaluated with the full IERS 2003 fundamental
    arguments. 2000B itself keeps only the linear part of each argument, which
    is all it needs across 1995–2050, the span it was designed for; by 1550
    or 2650 the powers it drops have moved the arguments far enough to cost
    14 mas. With them restored, the worst error in longitude against the
    complete IAU 2000A series is 2.6 mas across 1550–2650 and 1.5 mas across
    1995–2050, where nothing changes: about half of it the truncation to 77
    terms and half the planetary terms, which 2000B represents by two
    constants. (150 terms would halve it and double the cost of every frame;
    the tables themselves are fitted to 10 mas.) Coefficients taken from the
    IERS 2003 tables as shipped in Skyfield; the term table is checked
    against ERFA's nut00b to the last bit and the result against nut00a.

    python -I tools/build_frames.py

ERFA (BSD-3) and Skyfield (MIT) are used only to produce and check these
numbers; nothing of either ships in the engine.
"""
import os
import numpy as np
import erfa
import skyfield.nutationlib as sk

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "src", "generated", "frames.ts")
J2000 = 2451545.0
AS2R = np.pi / 648000.0

# The IAU 2006 coefficients as published (Capitaine & Wallace 2006), arcseconds,
# ascending powers of T. A least-squares fit to ERFA recovered exactly these;
# stating them as published keeps the file readable, and the check below holds
# them to ERFA across 1550–2650.
gamb = [-0.052928, 10.556378, 0.4932044, -0.00031238, -0.000002788, 0.0000000260]
phib = [84381.412819, -46.811016, 0.0511268, 0.00053289, -0.000000440, -0.0000000176]
psib = [-0.041775, 5038.481484, 1.5584175, -0.00018522, -0.000026452, -0.0000000148]
epsa = [84381.406, -46.836769, -0.0001831, 0.00200340, -0.000000576, -0.0000000434]

T = np.linspace(-4.51, 6.51, 1103)
jd = J2000 + T * 36525.0
fw = np.array([erfa.pfw06(j, 0.0) for j in jd])  # gamb, phib, psib, epsa (rad)
obl = np.array([erfa.obl06(j, 0.0) for j in jd])
for coeffs, ref in ((gamb, fw[:, 0]), (phib, fw[:, 1]), (psib, fw[:, 2]), (epsa, obl)):
    back = np.polynomial.polynomial.polyval(T, coeffs) * AS2R
    assert np.abs(back - ref).max() < 1e-12, np.abs(back - ref).max()

# Nutation: IAU 2000B's 77 terms, offsets in 0.1 µas.
N = 77
nals = sk.nals_t[:N].astype(int)
lon = sk.lunisolar_longitude_coefficients[:N]
obq = sk.lunisolar_obliquity_coefficients[:N]
OFF_PSI, OFF_EPS = -0.000135e7, 0.000388e7
# The fundamental arguments l, l', F, D, Omega (Simon et al. 1994, as adopted
# by IERS 2003): one degree-4 polynomial each, arcseconds, ascending powers.
DELAUNAY = np.hstack([sk.fa0, sk.fa1, sk.fa2, sk.fa3, sk.fa4])

def nutation(jd_tt, powers=5):
    t = (jd_tt - J2000) / 36525.0
    a = sum(DELAUNAY[:, k] * t**k for k in range(powers)) % 1296000.0 * AS2R
    arg = nals @ a
    s, c = np.sin(arg), np.cos(arg)
    dpsi = (lon[:, 0] + lon[:, 1] * t) @ s + lon[:, 2] @ c + OFF_PSI
    deps = (obq[:, 0] + obq[:, 1] * t) @ c + obq[:, 2] @ s + OFF_EPS
    return dpsi * 1e-7 * AS2R, deps * 1e-7 * AS2R

# The term table, transcribed right: with only the linear part of each
# argument, which is how 2000B is defined, this is ERFA's nut00b to the bit.
worst = 0.0
for j in np.linspace(J2000 - 3 * 36525, J2000 + 5 * 36525, 2001):
    mine = nutation(j, powers=2)
    ref = erfa.nut00b(j, 0.0)
    worst = max(worst, abs(mine[0] - ref[0]), abs(mine[1] - ref[1]))
assert worst < 1e-15, worst

# And what the engine evaluates, against the complete IAU 2000A series (678
# luni-solar and 687 planetary terms) across the whole span of the tables.
# In longitude this is the error in every ecliptic place: the true ecliptic
# of date does not depend on the nutation in obliquity at all.
MAS = AS2R / 1000.0
worst_psi = worst_eps = 0.0
for j in np.linspace(J2000 - 4.51 * 36525, J2000 + 6.51 * 36525, 4001):
    mine = nutation(j)
    ref = erfa.nut00a(j, 0.0)
    worst_psi = max(worst_psi, abs(mine[0] - ref[0]) / MAS)
    worst_eps = max(worst_eps, abs(mine[1] - ref[1]) / MAS)
assert worst_psi < 3.0 and worst_eps < 2.0, (worst_psi, worst_eps)

# ── Long-term precession, for the DE441 tables (13,200 BC – AD 17,191) ─────
# IAU 2006 is a polynomial fitted to a few centuries of observation, not meant
# for more than a few thousand years: by 13,000 BC it is 47 arcMINUTES out.
# Vondrák, Capitaine & Wallace (2011, corrigendum 2012) give the ecliptic pole
# and the equator pole as cubics plus long-period terms, valid for ±200,000
# years and agreeing with IAU 2006 at J2000. Their Tables 1 and 2: periods in
# centuries, then the cos and sin amplitudes of each coordinate, arcseconds.
LTP_ECL_POLY = [[5851.607687, -0.1189000, -0.00028913, 0.000000101],
                [-1600.886300, 1.1689818, -0.00000020, -0.000000437]]
LTP_ECL_PER = [[708.15, -5486.751211, -684.661560, 667.666730, -5523.863691],
               [2309.00, -17.127623, 2446.283880, -2354.886252, -549.747450],
               [1620.00, -617.517403, 399.671049, -428.152441, -310.998056],
               [492.20, 413.442940, -356.652376, 376.202861, 421.535876],
               [1183.00, 78.614193, -186.387003, 184.778874, -36.776172],
               [622.00, -180.732815, -316.800070, 335.321713, -145.278396],
               [882.00, -87.676083, 198.296701, -185.138669, -34.744450],
               [547.00, 46.140315, 101.135679, -120.972830, 22.885731]]
LTP_EQU_POLY = [[5453.282155, 0.4252841, -0.00037173, -0.000000152],
                [-73750.930350, -0.7675452, -0.00018725, 0.000000231]]
LTP_EQU_PER = [[256.75, -819.940624, 75004.344875, 81491.287984, 1558.515853],
               [708.15, -8444.676815, 624.033993, 787.163481, 7774.939698],
               [274.20, 2600.009459, 1251.136893, 1251.296102, -2219.534038],
               [241.45, 2755.175630, -1102.212834, -1257.950837, -2523.969396],
               [2309.00, -167.659835, -2660.664980, -2966.799730, 247.850422],
               [492.20, 871.855056, 699.291817, 639.744522, -846.485643],
               [396.10, 44.769698, 153.167220, 131.600209, -1393.124055],
               [288.90, -512.313065, -950.865637, -445.040117, 368.526116],
               [231.10, -819.415595, 499.754645, 584.522874, 749.045012],
               [1610.00, -538.071099, -145.188210, -89.756563, 444.704518],
               [620.00, -189.793622, 558.116553, 524.429630, 235.934465],
               [157.87, -402.922932, -23.923029, -13.549067, 374.049623],
               [220.30, 179.516345, -165.405086, -210.157124, -171.330180],
               [1200.00, -9.814756, 9.344131, -44.919798, -22.899655]]
LTP_EPS0 = 84381.406                     # the J2000 obliquity the ecliptic pole is referred through
LTP_BIAS = [-0.016617, -0.0068192, -0.0146]  # frame bias dx, de, dr (arcseconds)


def ltp_pole(poly, per, t):
    a = b = 0.0
    w = 2 * np.pi * t
    for (P, ca, cb, sa, sb) in per:
        s, c = np.sin(w / P), np.cos(w / P)
        a += c * ca + s * sa
        b += c * cb + s * sb
    for k in range(4):
        a += poly[0][k] * t**k
        b += poly[1][k] * t**k
    return a * AS2R, b * AS2R


def ltpb(epj):
    """ICRS → mean equator and equinox of date, long-term (as ERFA's ltpb)."""
    t = (epj - 2000.0) / 100.0
    p, q = ltp_pole(LTP_ECL_POLY, LTP_ECL_PER, t)
    w = np.sqrt(max(0.0, 1 - p * p - q * q))
    s, c = np.sin(LTP_EPS0 * AS2R), np.cos(LTP_EPS0 * AS2R)
    pecl = np.array([p, -q * c - w * s, -q * s + w * c])
    x, y = ltp_pole(LTP_EQU_POLY, LTP_EQU_PER, t)
    peq = np.array([x, y, np.sqrt(max(0.0, 1 - x * x - y * y))])
    v = np.cross(peq, pecl)
    eqx = v / np.linalg.norm(v)
    rp = np.array([eqx, np.cross(peq, eqx), peq])
    dx, de, dr = (b * AS2R for b in LTP_BIAS)
    return np.array([[r[0] - r[1] * dr + r[2] * dx, r[0] * dr + r[1] + r[2] * de, -r[0] * dx - r[1] * de + r[2]]
                     for r in rp])


# Transcribed right: the same matrix as ERFA's ltpb, across ±200 millennia.
worst_ltp = max(np.abs(ltpb(e) - erfa.ltpb(e)).max() for e in np.linspace(-198000, 202000, 2001))
assert worst_ltp < 1e-14, worst_ltp


def small_angle_mas(d):
    """The angle of a rotation matrix close to the identity, in mas."""
    a = np.array([d[1, 2] - d[2, 1], d[2, 0] - d[0, 2], d[0, 1] - d[1, 0]]) / 2
    return np.degrees(np.arcsin(min(1.0, np.linalg.norm(a)))) * 3.6e6


# Where the engine changes model. Inside the DE440 tables it keeps IAU 2006,
# so no modern chart moves; beyond them it blends to the long-term model over
# BLEND_YEARS, because the two already differ by 10–15 mas at the seams —
# under anything a chart shows, but more than the tables are fitted to, and a
# hard switch would be a step in every position at one instant.
SPAN = (2287194.5, 2688602.5)  # the DE440 tables, data/manifest.json (TDB ≈ TT)
BLEND_YEARS = 50
gaps = [small_angle_mas(ltpb(2000.0 + (j - J2000) / 365.25) @ erfa.pmat06(j, 0.0).T) for j in SPAN]

# ── Sidereal time, long-term ───────────────────────────────────────────────
# IAU 2006 sidereal time is the Earth rotation angle plus a polynomial for how
# far the equinox has drifted along the equator from the CIO — the equation
# of the origins. That polynomial is a fit, and diverges like the precession
# one does. Over 30,000 years it must be built from the long-term model: the
# CIO locator s integrated from its definition (Capitaine et al. 2000),
#     ds/dt = −(X dY/dt − Y dX/dt) / (1 + Z),
# with X, Y, Z the celestial pole of date in the ICRS (the third row of the
# bias-precession matrix), then EO = s − atan2(q, p) as in ERFA's eors. This
# is the MEAN equation of the origins (precession only); the engine adds the
# equation of the equinoxes for nutation as it does today. Integrated on a
# one-year grid by Simpson's rule (every term in the model has a period of
# 157 centuries or more), sampled every EO_STEP years, and checked against
# ERFA's gmst06 where IAU 2006 holds.
EO_STEP_YEARS = 10
EO_FROM, EO_TO = -13260.0, 17260.0  # Julian epochs; beyond the DE441 tables either side
epochs = np.arange(EO_FROM, EO_TO + 0.5, 1.0)
PB = erfa.ltpb(epochs)  # (n, 3, 3), the same as ltpb() above
X, Y, Z = PB[:, 2, 0], PB[:, 2, 1], PB[:, 2, 2]
dX, dY = np.gradient(X, epochs), np.gradient(Y, epochs)
rate = -(X * dY - Y * dX) / (1 + Z)  # rad per year
i0 = int(np.where(epochs == 2000.0)[0][0])
s = np.zeros_like(epochs)
# cumulative trapezoid, then the J2000 constant of IAU 2006's s (94 µas).
s[1:] = np.cumsum((rate[1:] + rate[:-1]) / 2 * np.diff(epochs))
s = s - s[i0] + 94e-6 * AS2R
# Plus the secular term the integral above cannot see: nutation sweeping its
# small ellipse around the moving pole adds 3,808.65 µas a century to s (the
# "crossed nutation-precession" term of IAU 2006's s06). The engine applies
# nutation separately, so this mean EO must carry it.
s = s + 3808.65e-6 * AS2R * (epochs - 2000.0) / 100.0
ax = X / (1 + Z)
xs, ys, zs = 1 - ax * X, -ax * Y, -X
p = PB[:, 0, 0] * xs + PB[:, 0, 1] * ys + PB[:, 0, 2] * zs
q = PB[:, 1, 0] * xs + PB[:, 1, 1] * ys + PB[:, 1, 2] * zs
eo = s - np.arctan2(q, p)
eo = np.unwrap(eo)  # radians, continuous

# Checked: GMST = ERA − EO must be IAU 2006's GMST where that polynomial holds —
# exactly at J2000 (0.003 mas over a decade), parting as the two precession
# models do, and inside the 10 mas the tables are fitted to across 1550–2650.
check = (epochs >= 1550) & (epochs <= 2650)
jd_chk = J2000 + (epochs[check] - 2000.0) * 365.25
gmst06 = np.array([erfa.gmst06(j, 0.0, j, 0.0) for j in jd_chk])
era = np.array([erfa.era00(j, 0.0) for j in jd_chk])
d = (era - eo[check]) - gmst06
d = (d + np.pi) % (2 * np.pi) - np.pi
worst_gmst_mas = np.abs(d).max() / MAS
assert worst_gmst_mas < 10.0, worst_gmst_mas

# Sampled every EO_STEP_YEARS and interpolated by cubic Hermite, as the engine
# does; the interpolation is held to a hundredth of a mas.
step_idx = slice(None, None, EO_STEP_YEARS)
eo_samples = eo[step_idx] / AS2R  # arcseconds
eo_epochs = epochs[step_idx]


def eo_hermite(epj):
    x = (epj - EO_FROM) / EO_STEP_YEARS
    i = np.floor(x).astype(int)
    u = x - i
    p0, p1, p2, p3 = eo_samples[i - 1], eo_samples[i], eo_samples[i + 1], eo_samples[i + 2]
    m1, m2 = (p2 - p0) / 2, (p3 - p1) / 2
    return ((2 * u**3 - 3 * u**2 + 1) * p1 + (u**3 - 2 * u**2 + u) * m1 + (-2 * u**3 + 3 * u**2) * p2 + (u**3 - u**2) * m2)


inner = (epochs > EO_FROM + 2 * EO_STEP_YEARS) & (epochs < EO_TO - 2 * EO_STEP_YEARS)
worst_eo_interp = np.abs(eo_hermite(epochs[inner]) - eo[inner] / AS2R).max() * 1000
assert worst_eo_interp < 0.01, worst_eo_interp

rows = []
for i in range(N):
    m = ", ".join(str(int(v)) for v in nals[i])
    co = ", ".join(f"{float(v):.10g}" for v in (*lon[i], *obq[i]))
    rows.append(f"  {m}, {co},")

# repr: the shortest decimal that reads back as exactly the same double.
delaunay = "\n".join("  [" + ", ".join(repr(float(v)) for v in row) + "]," for row in DELAUNAY)

os.makedirs(os.path.dirname(OUT), exist_ok=True)
body = "\n".join(rows)
with open(OUT, "w", encoding="utf-8", newline="\n") as f:
    f.write(f"""// GENERATED by tools/build_frames.py — do not edit. Rerun the tool instead.
//
// IAU 2006 precession (Fukushima-Williams angles) and mean obliquity, and IAU
// 2000B's nutation terms with the full fundamental arguments. Polynomials are
// in arcseconds, ascending powers of Julian centuries of TT from J2000.
// Checked against ERFA (pfw06, obl06, nut00b, nut00a). Across 1550–2650 the
// nutation is within {worst_psi:.1f} mas of IAU 2000A in longitude, {worst_eps:.1f} mas in obliquity.

/** gamma-bar, phi-bar, psi-bar and epsilon-A: IAU 2006, frame bias included. */
export const FW_GAMB = {gamb} as const;
export const FW_PHIB = {phib} as const;
export const FW_PSIB = {psib} as const;
export const OBL06 = {epsa} as const;

/**
 * The fundamental arguments l, l', F, D and Omega (IERS 2003): one polynomial
 * each, arcseconds. IAU 2000B as published keeps only the first two powers.
 */
export const DELAUNAY = [
{delaunay}
] as const;

/** Offsets standing in for the planetary terms, in 0.1 microarcseconds. */
export const NUT_OFFSET_PSI = {OFF_PSI};
export const NUT_OFFSET_EPS = {OFF_EPS};

/**
 * The 77 terms, 11 numbers each: multipliers of l, l', F, D, Omega, then
 * longitude (sin, sin*t, cos) and obliquity (cos, cos*t, sin), in 0.1 µas.
 */
// prettier-ignore
export const NUT_TERMS = [
{body}
] as const;
""")
print(f"wrote {OUT}")

# The long-term model goes in its own file, so frames.ts — what every chart
# inside the DE440 tables uses — is untouched by it.
OUT_LONG = os.path.join(HERE, "..", "src", "generated", "frames-long.ts")
rows5 = lambda table: "\n".join("  [" + ", ".join(repr(float(v)) for v in r) + "]," for r in table)
eo_values = ",".join(f"{v:.6f}" for v in eo_samples)
with open(OUT_LONG, "w", encoding="utf-8", newline="\n") as f:
    f.write(f"""// GENERATED by tools/build_frames.py — do not edit. Rerun the tool instead.
//
// Long-term precession (Vondrák, Capitaine & Wallace 2011, corrigendum 2012),
// valid for ±200,000 years, for moments beyond the DE440 tables. Checked
// against ERFA's ltpb to {worst_ltp:.0e}. At the DE440 tables' ends it differs
// from IAU 2006 by {gaps[0]:.1f} mas (1550) and {gaps[1]:.1f} mas (2650), so the engine blends
// the two over {BLEND_YEARS} years outside them rather than switching.

/** The DE440 tables' TDB span: IAU 2006 inside it, the long-term model beyond. */
export const SPAN_START = {SPAN[0]};
export const SPAN_END = {SPAN[1]};
export const BLEND_YEARS = {BLEND_YEARS};

/** Ecliptic pole: cubic in centuries from J2000 for P and Q, arcseconds. */
export const LTP_ECL_POLY = [
{rows5(LTP_ECL_POLY)}
] as const;
/** Ecliptic pole, long-period terms: period (centuries), cos P, cos Q, sin P, sin Q (arcseconds). */
export const LTP_ECL_PER = [
{rows5(LTP_ECL_PER)}
] as const;
/** Equator pole: cubic in centuries from J2000 for X and Y, arcseconds. */
export const LTP_EQU_POLY = [
{rows5(LTP_EQU_POLY)}
] as const;
/** Equator pole, long-period terms: period (centuries), cos X, cos Y, sin X, sin Y (arcseconds). */
export const LTP_EQU_PER = [
{rows5(LTP_EQU_PER)}
] as const;
/** The J2000 obliquity the ecliptic pole is referred through, and the frame bias dx, de, dr (arcseconds). */
export const LTP_EPS0 = {LTP_EPS0};
export const LTP_BIAS = {LTP_BIAS} as const;

/**
 * The mean equation of the origins, arcseconds, every {EO_STEP_YEARS} Julian years from
 * epoch {EO_FROM}: the CIO locator integrated through the long-term model.
 * GMST = ERA − EO. Matches IAU 2006 GMST within {worst_gmst_mas:.2f} mas across 1550–2650;
 * the Hermite interpolation is within {worst_eo_interp * 1000:.2f} µas of the integral.
 */
export const EO_FROM_EPOCH = {EO_FROM};
export const EO_STEP_YEARS = {EO_STEP_YEARS};
// prettier-ignore
export const EO_MEAN = [{eo_values}] as const;
""")
print(f"wrote {OUT_LONG}")
print(f"long-term precession vs ERFA ltpb: {worst_ltp:.1e}; vs IAU 2006 at the seams: {gaps[0]:.2f} mas, {gaps[1]:.2f} mas")
print(f"long-term GMST vs IAU 2006 across 1550–2650: {worst_gmst_mas:.3f} mas; EO interpolation {worst_eo_interp*1000:.3f} µas")
print("FW gamb", gamb)
print("FW phib", phib)
print("FW psib", psib)
print("OBL06", epsa)
print(f"nutation vs IAU 2000A, 1550–2650: {worst_psi:.3f} mas in longitude, {worst_eps:.3f} mas in obliquity")
