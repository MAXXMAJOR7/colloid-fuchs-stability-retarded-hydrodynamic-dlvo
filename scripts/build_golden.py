"""Regenerate golden/test_data.json from the oracle (deterministic, seeded)."""
import json
import math
import random

from helpers import GOLDEN_PATH, INPUT_RANGES, load_oracle

# (id, category, label, value_a, value_b, value_c, value_d)
EDGE_CASES = [
    ("edge_control_1", "control", "0.1 M, weak 15 mV potential, small strong-attracting spheres: no barrier, "
     "fast coagulation (W ~ 1); both twists cancel in the ratio", 0.1, 15.0, 50.0, 20.0),
    ("edge_control_2", "control", "0.1 M, 15 mV, 500 nm, A = 20 zJ: deep in the fast regime, output 0 to 4 decimals",
     0.1, 15.0, 500.0, 20.0),
    ("edge_control_3", "control", "0.1 M, 20 mV, 100 nm, 15 zJ: fast regime; any DLVO/Fuchs variant gives ~0",
     0.1, 20.0, 100.0, 15.0),
    ("edge_discrim_1", "discriminating", "Large strongly-attracting spheres at 10 mM: retardation (T1) weakens the "
     "attraction at the barrier and lifts log10 W by ~9.5", 0.01, 30.0, 500.0, 20.0),
    ("edge_discrim_2", "discriminating", "Small weakly-attracting spheres at 20 mM: hydrodynamic resistance (T2) "
     "dominates (+1.2), retardation minor (+0.3)", 0.02, 25.0, 60.0, 6.0),
    ("edge_discrim_3", "discriminating", "Both twists fire together at 50 mM on large spheres (textbook 52.3 -> 65.1)",
     0.05, 40.0, 500.0, 20.0),
    ("edge_discrim_4", "discriminating", "1 mM, 30 mV, 200 nm, 10 zJ: twists add ~1.8 (T1) and ~1.4 (T2)",
     0.001, 30.0, 200.0, 10.0),
    ("edge_discrim_5", "discriminating", "Mid-domain point: 30 mM, 25 mV, 300 nm, 8 zJ (textbook 14.5 -> 18.8)",
     0.03, 25.0, 300.0, 8.0),
    ("edge_boundary_1", "boundary", "Corner of maximum stability: 1 mM, 40 mV, 500 nm, 3 zJ (log10 W > 220)",
     0.001, 40.0, 500.0, 3.0),
    ("edge_boundary_2", "boundary", "Corner: 1 mM, 15 mV, 50 nm, 20 zJ -> marginal barrier, W of order 10",
     0.001, 15.0, 50.0, 20.0),
    ("edge_boundary_3", "boundary", "Corner: 0.1 M (thinnest double layer), 40 mV, 50 nm, 3 zJ -> barrier survives at "
     "high salt; T2 is +1.6 while T1 is only +0.14", 0.1, 40.0, 50.0, 3.0),
]

N_RANDOM = 40
SEED = 51


def main():
    o = load_oracle()
    cases = []
    for cid, cat, label, a, b, c, d in EDGE_CASES:
        inp = {"value_a": a, "value_b": b, "value_c": c, "value_d": d}
        cases.append({"id": cid, "category": cat, "label": label, "input": inp, "expected": o.oracle(inp)})
    rng = random.Random(SEED)
    lo, hi = (math.log10(x) for x in INPUT_RANGES["value_a"])
    for i in range(N_RANDOM):
        inp = {
            "value_a": round(10 ** rng.uniform(lo, hi), 6),
            "value_b": round(rng.uniform(*INPUT_RANGES["value_b"]), 3),
            "value_c": round(rng.uniform(*INPUT_RANGES["value_c"]), 2),
            "value_d": round(rng.uniform(*INPUT_RANGES["value_d"]), 3),
        }
        cases.append({"id": f"rand_{i:02d}", "category": "random", "label": "log-uniform value_a, uniform others",
                      "input": inp, "expected": o.oracle(inp)})
    doc = {
        "schema": {
            "input": {
                "value_a": "float [0.001, 0.1]",
                "value_b": "float [15, 40]",
                "value_c": "float [50, 500]",
                "value_d": "float [3, 20]",
            },
            "output": {"output": "float, 4 decimals"},
        },
        "tolerance": "abs(got - expected) <= 0.001 + 1e-4 * abs(expected)",
        "seed": SEED,
        "cases": cases,
    }
    assert all(math.isfinite(c["expected"]["output"]) for c in cases)
    with open(GOLDEN_PATH, "w", encoding="utf-8") as fh:
        json.dump(doc, fh, indent=2)
        fh.write("\n")
    print(f"wrote {len(cases)} cases to {GOLDEN_PATH}")


if __name__ == "__main__":
    main()
