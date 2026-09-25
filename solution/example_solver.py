"""Example solver: probe the black box, then fit a hypothesised model.

Strategy (what an expert would do after exploratory probing):
  1. The output is >= 0, collapses to ~0 when value_a is large and rises steeply as
     value_a falls, value_b rises or value_c rises; it falls as value_d rises.
     d(output)/d(log value_a) is large and negative in the slow regime (Reerink-Overbeek
     behaviour) -> log10 of a colloidal stability ratio: value_a = salt (screening),
     value_b = surface potential, value_c = particle size, value_d = Hamaker constant.
  2. Barrier heights inferred from the output scale as value_c * value_b^2 -> constant
     potential, Derjaguin-scaled double-layer repulsion (HHF form).
  3. Textbook DLVO + Fuchs undershoots badly for large spheres with large value_d, and
     the gap grows with the separation of the barrier -> retarded attraction
     (Gregory form, one characteristic length fitted here).
  4. Even with a negligible barrier (large value_a, small value_c) the output is not
     what the plain Fuchs integral gives, and small spheres with weak attraction sit
     ~1 unit higher than predicted -> hydrodynamic resistance in both Fuchs integrals
     (Honig-Roebersen-Wiersema beta(u), McGown-Parfitt ratio definition).
  Fitted parameters: water permittivity eps_r and Gregory length lambda/14.

This solver is illustrative, not guaranteed optimal. Standard library only.
"""
import math

PROBE_BUDGET = 200

E, KB, NA, EPS0 = 1.602176634e-19, 1.380649e-23, 6.02214076e23, 8.8541878188e-12
TEMP = 298.15

_params = None


def _model(c, psi_mv, a_nm, a_zj, eps_r, lgreg, n=6000):
    kt = KB * TEMP
    a, psi, ham = a_nm * 1e-9, psi_mv * 1e-3, a_zj * 1e-21
    kappa = math.sqrt(2 * NA * 1000 * c * E * E / (eps_r * EPS0 * kt))
    s0, s1 = math.log(1e-9), math.log(1e7)
    ds = (s1 - s0) / n
    num, den = [], []
    for i in range(n + 1):
        u = math.exp(s0 + i * ds)
        h = u * a
        w = 0.5 if i in (0, n) else 1.0                       # trapezoid on ln u
        beta = (6 * u * u + 13 * u + 2) / (6 * u * u + 4 * u)
        base = math.log(w * beta * u) - 2 * math.log(2 + u)
        va = -ham * a / (12 * h) / (1 + h / lgreg) / kt
        vr = 2 * math.pi * eps_r * EPS0 * a * psi * psi * math.log1p(math.exp(-kappa * h)) / kt
        den.append(base + va)
        num.append(base + va + vr)
    tail = -math.log(2 + math.exp(s1)) - math.log(ds)

    def lse(v):
        m = max(v)
        return m + math.log(sum(math.exp(x - m) for x in v))

    return (lse(num + [tail]) - lse(den + [tail])) / math.log(10)


def _invert(target, f, lo, hi, increasing=True, iters=40):
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        if (f(mid) < target) == increasing:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def fit(query):
    global _params

    def q(a, b, c, d):
        return query({"value_a": a, "value_b": b, "value_c": c, "value_d": d})["output"]

    # Exploratory probes (one-at-a-time sweeps used to identify the family; see docstring).
    for c in (0.001, 0.003, 0.01, 0.03, 0.1):
        q(c, 30.0, 200.0, 10.0)
    for psi in (15.0, 25.0, 40.0):
        q(0.01, psi, 200.0, 10.0)
    for rad in (50.0, 150.0, 500.0):
        q(0.01, 30.0, rad, 10.0)
    for ham in (3.0, 10.0, 20.0):
        q(0.01, 30.0, 200.0, ham)

    # Calibration: y1 is insensitive to retardation (small spheres, weak attraction, thin
    # double layer), y2 is dominated by it.  Alternate 1-D inversions.
    y1 = q(0.1, 40.0, 50.0, 3.0)
    y2 = q(0.01, 30.0, 500.0, 20.0)
    eps_r, lgreg = 78.5, 7e-9
    for _ in range(4):
        eps_r = _invert(y1, lambda e: _model(0.1, 40.0, 50.0, 3.0, e, lgreg), 70.0, 90.0)
        # weaker attraction (shorter retardation length) -> higher W, so W decreases with lgreg
        lgreg = _invert(y2, lambda L: _model(0.01, 30.0, 500.0, 20.0, eps_r, L), 1e-9, 1e-7,
                        increasing=False)
    _params = (eps_r, lgreg)


def solve(inputs):
    if _params is None:
        raise RuntimeError("call fit(query) first")
    a, b, c, d = (float(inputs[k]) for k in ("value_a", "value_b", "value_c", "value_d"))
    return {"output": round(_model(a, b, c, d, *_params), 4)}
