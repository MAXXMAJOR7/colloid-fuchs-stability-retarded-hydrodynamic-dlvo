"""Shared utilities for task scripts and grader."""
import importlib.util
import json
import math
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GOLDEN_PATH = os.path.join(ROOT, "golden", "test_data.json")
ORACLE_PATH = os.path.join(ROOT, "oracle", "implement.py")

INPUT_KEYS = ("value_a", "value_b", "value_c", "value_d")
INPUT_RANGES = {
    "value_a": (0.001, 0.1),
    "value_b": (15.0, 40.0),
    "value_c": (50.0, 500.0),
    "value_d": (3.0, 20.0),
}
# A case passes when |got - expected| <= ABS_TOL + REL_TOL * |expected|.
# log10 W spans ~0 .. 100; the relative part keeps the bar at ~2e-4 kT per kT of barrier.
ABS_TOL = 0.001
REL_TOL = 1e-4


def tolerance(expected):
    return ABS_TOL + REL_TOL * abs(expected)


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def load_oracle():
    return load_module(ORACLE_PATH, "oracle_impl")


def load_golden(path=GOLDEN_PATH):
    with open(path, "r", encoding="utf-8") as fh:
        return json.load(fh)


def variant(c, psi_mv, a_nm, a_zj, retard=True, hydro=True, rounded=True):
    """Oracle with individual twists frozen, for twist-integrity checks and baselines.

    retard=False -> textbook unretarded Hamaker sphere-sphere attraction  -A a / 12 h   (T1 off)
    hydro=False  -> no hydrodynamic resistance, beta = 1 (original Fuchs integrals)  (T2 off)
    """
    o = load_oracle()
    kt = o.K_B * o.TEMP
    a, psi, ham = a_nm * o.NM, psi_mv * o.MV, a_zj * o.ZJ
    kappa = o.debye_kappa(c)
    s0, s1 = math.log(o.U_MIN), math.log(o.U_MAX)
    ds = (s1 - s0) / o.N_SIMPSON
    ln_num, ln_den = [], []
    for i in range(o.N_SIMPSON + 1):
        u = math.exp(s0 + i * ds)
        h = u * a
        weight = 1 if i in (0, o.N_SIMPSON) else (4 if i % 2 else 2)
        b = o.beta_hrw(u) if hydro else 1.0
        base = math.log(weight * b * u) - 2 * math.log(2 + u)
        va = (o.v_attr(h, a, ham) if retard else -ham * a / (12 * h)) / kt
        ln_den.append(base + va)
        ln_num.append(base + va + o.v_rep(h, a, psi, kappa) / kt)
    tail = -math.log(2 + o.U_MAX) - math.log(ds / 3)
    log_w = (o._log_sum_exp(ln_num + [tail]) - o._log_sum_exp(ln_den + [tail])) / math.log(10)
    return round(log_w, 4) if rounded else log_w
