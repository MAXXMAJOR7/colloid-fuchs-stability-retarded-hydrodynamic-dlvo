"""Automated GATE 2 checks. Exits non-zero on any failure."""
import ast
import json
import math
import random
import sys

from helpers import (GOLDEN_PATH, INPUT_KEYS, INPUT_RANGES, ORACLE_PATH, load_golden,
                     load_oracle, tolerance, variant)

# Every numeric literal allowed in oracle/implement.py, with its justification.
ALLOWED_CONSTANTS = {
    1.602176634e-19: "elementary charge (exact, SI 2019)",
    1.380649e-23: "Boltzmann constant (exact, SI 2019)",
    6.02214076e23: "Avogadro constant (exact, SI 2019)",
    8.8541878188e-12: "vacuum permittivity (CODATA 2022)",
    1000.0: "L per m^3 (exact)",
    298.15: "IUPAC standard temperature 25 degC (exact)",
    78.30: "Malmberg & Maryott 1956: eps_r of water at 25 degC",
    14.0: "Gregory 1981: 1 + 14 h / lambda",
    100e-9: "Gregory 1981: characteristic wavelength lambda = 100 nm",
    6.0: "Honig et al. 1971: 6u^2 (numerator and denominator)",
    13.0: "Honig et al. 1971: 13u",
    2.0: "Honig et al. 1971: + 2",
    4.0: "Honig et al. 1971: 4u",
    1e-3: "V per mV (exact); also lower bound of value_a domain (mol/L)",
    1e-9: "m per nm (exact)",
    1e-21: "J per zJ (exact)",
    1e-10: "quadrature lower limit on u = h/a (integrand < exp(-1e9) there)",
    1e6: "quadrature upper limit on u; analytic tail beyond",
    8192: "fixed Simpson panel count (determinism; converged to < 1e-6)",
    12: "Derjaguin equal spheres: A a / 12 h (= sphere-plate A a_p / 6 h with a_p = a/2)",
    2: "structural: 2 N_A (1:1 salt), 2 pi (HHF equal spheres), squares, (2 + u)^2, Simpson weight",
    4: "Simpson weight; 4-decimal output rounding (task spec)",
    3: "Simpson ds/3; lower bound of value_d domain (zJ)",
    1: "structural: 1 + ..., Simpson end weight",
    0: "structural: first grid index",
    10: "log10 base",
    0.1: "value_a domain upper bound (mol/L)",
    15: "value_b domain lower bound (mV)", 40: "value_b domain upper bound (mV)",
    50: "value_c domain lower bound (nm)", 500: "value_c domain upper bound (nm)",
    20: "value_d domain upper bound (zJ)",
}

failures = []


def check(name, cond, detail=""):
    print(f"[{'PASS' if cond else 'FAIL'}] {name}" + (f" -- {detail}" if detail else ""))
    if not cond:
        failures.append(name)


def numeric_literals(path):
    tree = ast.parse(open(path, encoding="utf-8").read())
    vals = []
    for node in ast.walk(tree):
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub) and isinstance(node.operand, ast.Constant):
            if isinstance(node.operand.value, (int, float)) and not isinstance(node.operand.value, bool):
                vals.append(-node.operand.value)
                node.operand.value = None
        elif isinstance(node, ast.Constant) and isinstance(node.value, (int, float)) and not isinstance(node.value, bool):
            vals.append(node.value)
    return vals


def sample(rng):
    lo, hi = (math.log10(x) for x in INPUT_RANGES["value_a"])
    return (10 ** rng.uniform(lo, hi), rng.uniform(*INPUT_RANGES["value_b"]),
            rng.uniform(*INPUT_RANGES["value_c"]), rng.uniform(*INPUT_RANGES["value_d"]))


def main():
    o = load_oracle()
    rng = random.Random(7)
    pts = [sample(rng) for _ in range(120)]
    ys = [o.compute(*p) for p in pts]

    # Determinism and schema
    check("deterministic", all(o.compute(*p) == y for p, y in zip(pts, ys)))
    out = o.oracle({"value_a": 0.01, "value_b": 30.0, "value_c": 200.0, "value_d": 10.0})
    check("output key is 'output' and scalar float", list(out) == ["output"] and isinstance(out["output"], float))
    check("4-decimal precision", all(round(y, 4) == y for y in ys))
    check("output is log10 of a ratio >= 1 (non-negative)", all(y >= 0 for y in ys))

    # No invented constants
    unknown = [v for v in numeric_literals(ORACLE_PATH) if v not in ALLOWED_CONSTANTS]
    check("every oracle numeric literal is justified", not unknown, f"unjustified: {unknown}" if unknown else "")

    # Literature anchors
    debye_nm = 1e9 / o.debye_kappa(0.001)
    check("Debye length at 1 mM 1:1 salt, 25 degC ~ 9.6 nm", abs(debye_nm - 9.61) < 0.02, f"{debye_nm:.3f} nm")
    check("HRW beta(u) -> 1 at large separation", abs(o.beta_hrw(1e6) - 1) < 2e-6)
    check("HRW beta(u) ~ 1/(2u) at contact", abs(o.beta_hrw(1e-8) * 2e-8 - 1) < 1e-6)
    ratio = o.v_attr(10e-9, 1e-7, 1e-20) / (-1e-20 * 1e-7 / (12 * 10e-9))
    check("Gregory factor at h = 10 nm is 1/2.4", abs(ratio - 1 / 2.4) < 1e-12, f"{ratio:.6f}")
    v0 = o.v_rep(0.0, 1e-7, 0.025, 1e8) / (2 * math.pi * o.EPS_R * o.EPS_0 * 1e-7 * 0.025 ** 2)
    check("HHF at contact reduces to 2 pi eps a psi^2 ln 2", abs(v0 - math.log(2)) < 1e-12)

    # Quadrature convergence (oracle grid vs a 4x finer grid)
    for p in [(0.001, 40.0, 500.0, 3.0), (0.02, 25.0, 60.0, 6.0), (0.1, 15.0, 50.0, 20.0)]:
        old = o.N_SIMPSON
        o.N_SIMPSON = 4 * old
        fine = variant(*p, rounded=False)
        o.N_SIMPSON = old
        check(f"quadrature converged at {p}", abs(fine - variant(*p, rounded=False)) < 1e-6)

    # Helper re-implementation agrees
    check("helper variant(all twists) == oracle", all(variant(*p) == y for p, y in zip(pts[:40], ys)))

    # Twist integrity
    d1 = [y - variant(*p, retard=False) for p, y in zip(pts, ys)]
    d2 = [y - variant(*p, hydro=False) for p, y in zip(pts, ys)]
    check("T1 retardation changes output (> 100x tolerance somewhere)", max(d1) > 0.5, f"max shift {max(d1):.3f}")
    check("T2 hydrodynamics changes output (> 100x tolerance somewhere)", max(d2) > 0.5, f"max shift {max(d2):.3f}")
    check("both twists raise W (weaker attraction / slower approach)", min(d1) >= 0 and min(d2) >= 0)
    both = sum(1 for a, b, y in zip(d1, d2, ys) if a > tolerance(y) and b > tolerance(y))
    check("both twists fire on most of the domain", both / len(pts) > 0.8, f"{both}/{len(pts)}")
    fast = [p for p, y in zip(pts, ys) if y <= 0.001] + [(0.1, 15.0, 50.0, 20.0), (0.1, 15.0, 500.0, 20.0)]
    check("both twists neutral (within tolerance) in the fast regime W ~ 1 (control)",
          all(abs(o.compute(*p) - variant(*p, retard=False, hydro=False)) <= tolerance(0.0) for p in fast),
          f"{len(fast)} fast points")
    # I/O isolation: T1 grows with particle size at fixed barrier conditions; T2 dominates small spheres
    s_small = o.compute(0.02, 25.0, 60.0, 6.0)
    check("T2 > T1 for small weakly attracting spheres",
          s_small - variant(0.02, 25.0, 60.0, 6.0, hydro=False) > s_small - variant(0.02, 25.0, 60.0, 6.0, retard=False))
    s_big = o.compute(0.01, 30.0, 500.0, 20.0)
    check("T1 > T2 for large strongly attracting spheres",
          s_big - variant(0.01, 30.0, 500.0, 20.0, retard=False) > s_big - variant(0.01, 30.0, 500.0, 20.0, hydro=False))

    # Golden data
    g = load_golden()
    cats = [c["category"] for c in g["cases"]]
    check(">=2 discriminating edge cases", cats.count("discriminating") >= 2)
    check(">=2 control edge cases", cats.count("control") >= 2)
    check(">=1 boundary edge case", cats.count("boundary") >= 1)
    check("golden inputs neutral and in range", all(
        sorted(c["input"]) == sorted(INPUT_KEYS)
        and all(INPUT_RANGES[k][0] <= c["input"][k] <= INPUT_RANGES[k][1] for k in INPUT_KEYS)
        for c in g["cases"]))
    check("golden expected outputs match oracle",
          all(o.oracle(c["input"]) == c["expected"] for c in g["cases"]), GOLDEN_PATH)

    print(json.dumps({"failures": failures}))
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
