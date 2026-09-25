# Blackbox Task 51 - Fuchs Stability Ratio of Charged Spheres Under Retarded Dispersion and Hydrodynamic Resistance

**Domain:** Chemistry - colloid & interface chemistry (DLVO theory, coagulation kinetics, electrical double layers)
**Tier:** Med (floor: 3h, 2 genuine twists, 3 substantive steps)
**Honest counts:** 4 substantive steps (plus a trivial rounding step), 2 genuine twists

**Repo name:** `colloid-fuchs-stability-retarded-hydrodynamic-dlvo`
**GitHub description:** Reverse-engineer the log stability ratio of a charged colloid from salt, surface potential, particle size and Hamaker constant. Textbook DLVO plus Fuchs misses by up to 13 decades, because the dispersion attraction is retarded and the particles slow down hydrodynamically as they approach. Requires double-layer theory, Lifshitz-retardation approximations, and diffusion-limited coagulation kinetics.

---

## Task description (for reviewers)

The hidden function maps four real-valued inputs to the base-10 logarithm of a dimensionless kinetic ratio, the factor by which a colloidal dispersion coagulates more slowly than it would without an energy barrier. The computation first builds an interparticle pair potential. It then evaluates two diffusion-with-drift integrals over the whole separation range and takes their ratio. Two published refinements of the textbook treatment are built in. The first weakens the attraction at nanometre separations. The second slows the relative diffusion of the two particles near contact. A solver must identify the coagulation family, recover the pair-potential forms, find both refinements, and match the output to about 1e-4 relative. Expected solver knowledge: graduate colloid chemistry (Russel-Saville-Schowalter / Elimelech-Gregory level).

---

## 1. Domain & Algorithm

Colloid chemistry. The output is **log₁₀ W**, where W is the Fuchs stability ratio for doublet formation between equal spheres of radius a. The spheres carry a constant surface potential ψ in an aqueous 1:1 electrolyte of molar concentration c at 25 °C, and they attract with Hamaker constant A. W is the ratio of the barrier-free (fast) coagulation rate to the actual (slow) rate. Both rates include van der Waals attraction and hydrodynamic interaction (McGown-Parfitt definition). The integrals are evaluated to converged precision.

## 2. Core Method

- **Debye screening.** κ² = 2 N_A e² (1000 c) / (ε_r ε₀ k T) for a symmetric 1:1 electrolyte (Debye-Hückel). ε_r = 78.30 at 25 °C (Malmberg & Maryott 1956). T = 298.15 K.
- **Double-layer repulsion.** Hogg-Healy-Fuerstenau (1966) constant-potential formula for spheres, V = πεε₀ a_i a_j/(a_i+a_j) [2ψ_iψ_j ln((1+e^{−κh})/(1−e^{−κh})) + (ψ_i²+ψ_j²) ln(1−e^{−2κh})]. For equal spheres at equal potential it reduces exactly to V_R = 2π ε_r ε₀ a ψ² ln(1 + e^{−κh}).
- **Stability ratio.** Fuchs (1934), in the complete form of McGown & Parfitt (1967):
  W = ∫_{2a}^{∞} [β(r)/r²] e^{V_T/kT} dr / ∫_{2a}^{∞} [β(r)/r²] e^{V_A/kT} dr, with V_T = V_A + V_R. In u = h/a this becomes W = ∫₀^∞ β(u) e^{V_T/kT} (2+u)^−2 du / ∫₀^∞ β(u) e^{V_A/kT} (2+u)^−2 du.

## 3. Twists (2 genuine)

**T1 - Retarded dispersion attraction (Gregory 1981).** The textbook DLVO attraction between equal spheres is the unretarded Hamaker/Derjaguin form −A a/(12h). The oracle instead uses Gregory's retarded approximation, V_A = −(A a/12h)/(1 + 14h/λ) with λ = 100 nm. This is Gregory's sphere-plate expression −(A a_p/6h)/(1+14h/λ), carried to two equal spheres by the Derjaguin substitution a_p → a/2 (the same a_i a_j/(a_i+a_j) scaling that relates HHF eqs 7 and 8). At the barrier position (h ≈ 1-10 nm) retardation removes 10-60% of the attraction. That raises the barrier, and since W depends exponentially on it, log₁₀ W rises by up to 16 decades over the domain. The effect grows with the Hamaker constant, the radius, and the barrier separation (low salt).

**T2 - Hydrodynamic resistance in both Fuchs integrals (Honig, Roebersen & Wiersema 1971).** Classic Fuchs theory takes free Stokes-Einstein relative diffusion. The oracle multiplies the integrand by β(u) = (6u² + 13u + 2)/(6u² + 4u), which diverges as 1/(2u) at contact. β enters both integrals, so the fast reference rate also changes. The net effect is non-trivial. It is +0.7 to +2.2 decades where a barrier exists, is largest relative to the output for small, weakly attracting spheres, and vanishes when W → 1.

Both twists are I/O-isolable:
- **Fast regime (value_a → 0.1, value_b small).** Both twists are neutral within tolerance (controls).
- **Large value_c and value_d.** T1 dominates. At (0.01, 30, 500, 20), T1 = +9.5 and T2 = +1.8.
- **Small value_c, small value_d.** T2 dominates. At (0.02, 25, 60, 6), T2 = +1.19 and T1 = +0.32. At the high-salt corner (0.1, 40, 50, 3), T2 = +1.62 and T1 = +0.14.

## 4. Multi-Step Pipeline

| Step | Computation | Inputs involved |
|---|---|---|
| S1 | Debye screening κ(c) from ε_r(25 °C) and exact SI constants | value_a |
| S2 | Pair potential V_T(h) = V_R (HHF) + V_A (Gregory-retarded, T1) | value_a..d |
| S3 | Hydrodynamic resistance β(u) (T2) and assembly of both Fuchs integrands on a log-u grid | value_c (through u = h/a) |
| S4 | Both integrals (composite Simpson on ln u, log-sum-exp, analytic tail), ratio → log₁₀ W | all |
| (S5) | round to 4 decimals | trivial, not counted |

The steps cannot be collapsed. The output is a ratio of two Laplace-type integrals of an exponentiated non-separable potential. No closed form exists.

## 5. Input Schema

```json
{
  "value_a": "float [0.001, 0.1]",   // 1:1 electrolyte concentration, mol L^-1
  "value_b": "float [15, 40]",       // surface potential magnitude, mV
  "value_c": "float [50, 500]",      // particle radius, nm
  "value_d": "float [3, 20]"         // Hamaker constant, zJ (1e-21 J)
}
```

## 6. Output Schema

```json
{ "output": "float, 4 decimals" }     // log10 W; range 0 (fast coagulation) to ~225
```

Tolerance: a case passes when |got − expected| ≤ 0.001 + 1e-4·|expected|.

## 7. Hardness & Twists

These baselines are measured on the 51 golden cases (`.scratch/baselines.py`):

| Model | Pass | Median abs err | Max abs err |
|---|---|---|---|
| Textbook DLVO + Fuchs (no T1, no T2) | 4/51 | 4.24 | 12.76 |
| T2 only (unretarded attraction) | 4/51 | 2.48 | 10.57 |
| T1 only (no hydrodynamics) | 4/51 | 1.57 | 2.14 |
| Reerink-Overbeek barrier estimate, W ≈ e^{V_max/kT}/(2κa) | 4/51 | 3.67 | 11.95 |
| Unretarded model with best global effective Hamaker constant (A × 0.3-0.8, with or without β) | ≤ 4/51 | ≥ 0.80 | ≥ 6.0 |
| Full model with textbook ε_r = 78.54 | 4/51 | 0.16 | 0.76 |
| Example solver (16 probes; fits ε_r and λ/14) | 51/51 | 4e-5 | - |

The four cases that every model passes are the fast-regime controls and one fast random point, where W ≈ 1.

Why this is hard:
- **The family is recognisable, but the details are not.** The steep log W vs log c slope suggests DLVO + Fuchs to a colloid chemist. Every textbook version of that pipeline misses by several decades.
- **Two corrections act in the same direction.** Both raise W. A solver who adds only one of them, or who absorbs the difference into a fitted Hamaker constant, cannot match both the large-sphere and small-sphere regimes. Retardation scales with h/λ and hydrodynamics with h/a, so their input dependences differ.
- **The fast reference is non-trivial.** β also enters the denominator, and there the fast rate with hydrodynamics depends on A and a through a Laplace-type integral near contact. A solver who normalises by the Smoluchowski rate or by a hydrodynamics-free fast rate is off by a size- and Hamaker-dependent offset.
- **Exponential sensitivity.** Barriers reach ~500 kT, so the output needs the potential correct to ~1e-4 relative. Even the textbook ε_r = 78.54 instead of the measured 78.30 fails 47/51.
- **Linear and additive surrogates fail.** log W is ~0 over a region of input space and then rises almost linearly in the barrier height, with a smooth but strongly curved knee. The inputs interact multiplicatively (a·ψ² for V_R, A·a for V_A, κ with ε_r).
- **Cold read:** four neutral floats and a non-negative output that is often exactly 0. The schema names no field, unit, or physical quantity.

## 8. Constants & Citations

Every numeric literal in `oracle/implement.py` is whitelisted with its justification in `scripts/verify_task.py` (AST scan).

| Constant | Value | Status | Source |
|---|---|---|---|
| e | 1.602176634×10⁻¹⁹ C | EXACT | SI 2019 defining constant |
| k_B | 1.380649×10⁻²³ J K⁻¹ | EXACT | SI 2019 defining constant |
| N_A | 6.02214076×10²³ mol⁻¹ | EXACT | SI 2019 defining constant |
| ε₀ | 8.8541878188×10⁻¹² F m⁻¹ | CITED | CODATA 2022 (NIST) |
| T | 298.15 K | EXACT | IUPAC standard temperature (25 °C) |
| ε_r | 78.30 | CITED | Malmberg & Maryott 1956 (value at 25 °C) |
| 14, λ = 100 nm | - | CITED | Gregory 1981; λ = 100 nm standard choice |
| 12 | - | EXACT | Derjaguin: equal spheres = sphere-plate with a_p = a/2 (A a_p/6h → A a/12h) |
| 2π, ln(1+e^{−κh}) | - | EXACT | HHF equal-sphere reduction |
| 6, 13, 2, 6, 4 | - | CITED | Honig, Roebersen & Wiersema 1971 β(u) |
| 2 in κ² | - | EXACT | 1:1 electrolyte, Σ c_i z_i² = 2c |
| 1000, 1e-3, 1e-9, 1e-21 | - | EXACT | unit conversions (L m⁻³, V/mV, m/nm, J/zJ) |
| 1e-10, 1e6, 8192 | - | STRUCTURAL | quadrature range on u and fixed Simpson panel count (converged to < 1e-6 against a 4× finer grid; checked in verify_task.py) |
| domain bounds | - | STRUCTURAL | input-domain check |

### Citations fetched & verified (quote < 15 words, location)

| Item | Verbatim quote | Where fetched |
|---|---|---|
| Gregory retarded form | "described by Gregory (1981) as U = −(A a_P/6h)(1/(1 + 14h/λ))" | Jin et al., *Sci. Rep.* 5, 17747 (2015), text above eq (7) (nature.com PDF, p. 6) |
| λ = 100 nm | "the characteristic wavelength of the interactions, often assumed to be 100 nm" | Langmuir 2023 (PMC10173465), text after eq (20), citing Gregory 1981 |
| HHF | "Based on the classic Hogg-Healy-Fuerstenau expression" (eq 7: πεε₀a_ia_j/(a_i+a_j)[…]) | Jin et al. 2015, eq (7), p. 6 |
| β(u) (HRW) | "f_corr = (6u² + 4u)/(6u² + 13.0u + 2) … proposed by Honig et al." | Urbina-Villalba et al., *Int. J. Mol. Sci.* 2009 (PMC2672001), eq (81); β = 1/f_corr |
| W definition | "The complete formula of W was deduced later by McGown and Parfitt" | same, eq (8): W = ∫[f(r)/r²]e^{V_T/kT}dr / ∫[f(r)/r²]e^{V_A/kT}dr |
| ε_r(25 °C) | "At 25° C the dielectric constant was found to have the value 78.30" | Malmberg & Maryott 1956, abstract (archive.org jresv56n1p1 OCR) |
| ε₀ | "8.854 187 8188(14) x 10⁻¹² F m⁻¹" | physics.nist.gov CODATA value page (ep0) |
| Gregory validity | "for the sphere-plate system, a simpler result … up to … 20% of the sphere radius" | Gregory 1981 abstract (ScienceDirect) |

**Citation-code match:** `GREGORY_B = 14.0`, `LAMBDA = 100e-9`; `HRW_A2, HRW_A1, HRW_A0 = 6, 13, 2`, `HRW_B2, HRW_B1 = 6, 4`; `EPS_R = 78.30`; `EPS_0 = 8.8541878188e-12`; `v_rep = 2π ε_r ε₀ a ψ² ln1p(e^{−κh})` (HHF eq 7 with a_i = a_j, ψ_i = ψ_j). `verify_task.py` also checks the following:
- the Debye length at 1 mM is 9.61 nm
- the Gregory factor is exactly 1/2.4 at h = 10 nm
- HHF at contact equals 2πεaψ² ln 2
- β → 1 far apart and β → 1/(2u) at contact

**Honest notes (CONCERN-level, documented):**
- **Gregory validity range.** Jin et al. quote the sphere-plate form "if h < λ/2π" (≈ 16 nm). The oracle applies it at every separation. Beyond ~16 nm both V_A and V_R are ≪ kT everywhere in the domain, so this affects neither the barrier nor W beyond tolerance. It is still a modelling choice.
- **HHF linearisation.** HHF is derived from the linearised Poisson-Boltzmann equation and is quantitatively best for |ψ| ≲ 25 mV and κa ≫ 1. The domain reaches 40 mV, and κa ≥ 5.2 at its corner (1 mM, 50 nm). The oracle is a deterministic model, not a claim of experimental accuracy at 40 mV.
- **Integration down to h → 0.** The Fuchs integrals run from contact (u = 1e-10), as in McGown-Parfitt. Below ~0.1 nm the integrand is suppressed by the unretarded attraction (> 50 kT for the smallest A·a in the domain), so the result does not depend on a physical cut-off within tolerance.
- **Secondary sources.** The primary Gregory (1981) and Honig et al. (1971) papers are paywalled. Their formulas were verified in open-access papers that quote them with the constants shown. The Gregory abstract was read on ScienceDirect.

## 9. Edge Cases

| id | value_a | value_b | value_c | value_d | output | textbook | Rationale |
|---|---|---|---|---|---|---|---|
| edge_control_1 | 0.1 | 15 | 50 | 20 | 0.0006 | 0.0000 | Fast regime: W ≈ 1, twists cancel in the ratio |
| edge_control_2 | 0.1 | 15 | 500 | 20 | 0.0000 | 0.0000 | Deep fast regime |
| edge_control_3 | 0.1 | 20 | 100 | 15 | 0.0005 | 0.0000 | Fast regime, any DLVO/Fuchs variant ≈ 0 |
| edge_discrim_1 | 0.01 | 30 | 500 | 20 | 50.4927 | 39.0996 | T1 dominates (+9.5), T2 +1.8 |
| edge_discrim_2 | 0.02 | 25 | 60 | 6 | 4.6332 | 3.1046 | T2 dominates (+1.19), T1 +0.32 |
| edge_discrim_3 | 0.05 | 40 | 500 | 20 | 65.0564 | 52.2936 | Both twists together, +12.8 |
| edge_discrim_4 | 0.001 | 30 | 200 | 10 | 41.0781 | 37.8754 | T1 +1.76, T2 +1.39 |
| edge_discrim_5 | 0.03 | 25 | 300 | 8 | 18.8365 | 14.5350 | Mid-domain, +4.3 |
| edge_boundary_1 | 0.001 | 40 | 500 | 3 | 223.6015 | 219.8731 | Maximum-stability corner |
| edge_boundary_2 | 0.001 | 15 | 50 | 20 | 0.9014 | 0.2483 | Marginal barrier: W ≈ 8 vs textbook 1.8 |
| edge_boundary_3 | 0.1 | 40 | 50 | 3 | 14.8438 | 13.0666 | Thinnest double layer, barrier survives; T2 +1.62, T1 +0.14 |

The golden file adds 40 seeded random cases (seed 51; value_a log-uniform), for 51 cases in total. Random outputs span 0 to 116.

---

## GATE 1: Screening

| Test | Result | Notes |
|---|---|---|
| Algebraic collapse | PASS | The output is the log of a ratio of two integrals of exp(non-separable potential). In the slow regime it is ≈ V_max/kT·log₁₀e + corrections, but V_max has no closed form in (c, ψ, a, A), and both twists change the prefactor and the barrier differently. |
| Domain recall | PASS (CONCERN) | A colloid chemist may say "DLVO stability ratio" after probing. The specific combination (HHF + Gregory-retarded Derjaguin attraction + HRW resistance in both McGown-Parfitt integrals) is not named by the schema, and no single published method is this exact composition. |
| Twist survives | PASS | T1 is a separation-dependent factor inside an exponent. T2 is a separation-dependent weight inside both integrals. Neither can be absorbed into a constant: fitting an effective Hamaker constant fails, as the T1-only and T2-only baselines show. |
| Genuine twist | PASS (2) | Both are published refinements that the standard textbook DLVO/Fuchs calculation omits. HHF, Debye screening, and the McGown-Parfitt ratio are core method and are not counted. |
| Two-trap | PASS | value_c and value_d isolate T1 (large a·A) from T2 (small a, small A). value_a and value_b switch both off (fast regime). |
| Tier fit | PASS | 4 steps ≥ 3; 2 twists ≥ 2 (Med). |
| PhD authenticity | PASS | Retardation and hydrodynamic corrections are the known reasons classical DLVO/Fuchs predictions of W disagree with experiment by orders of magnitude, a standard topic in colloid-kinetics research. |

**Verdict: PASS** (one documented CONCERN on family-level domain recall)

## GATE 2: Verification

| Check | Result |
|---|---|
| Twist integrity (freeze-one tests) | PASS: over 120 random points, T1 max shift 7.8 and T2 max shift 2.2 decades; both fire beyond tolerance on 116/120 points (`verify_task.py`) |
| I/O isolation | PASS: both neutral in the fast regime; T1 > T2 for large, strongly attracting spheres; T2 > T1 for small, weakly attracting spheres |
| No invented constant | PASS: AST scan finds only whitelisted literals |
| Recall-reconstructability | PASS: the example solver recovers the function with 16 probes (≪ 50) once the family is hypothesised, fitting ε_r = 78.300 and λ/14 = 7.143 nm |
| Citations fetched & verified | PASS: all sources fetched (table above). The two paywalled primaries were verified through open papers quoting them. |
| Citation-code match | PASS |
| Output/schema neutralisation | PASS: key "output", float, value_a to value_d |
| Tier floor | PASS: 4 steps / 2 twists vs Med 3 / 2 |
| Edge cases ≥ 5 | PASS: 11 total (3 control, 5 discriminating, 3 boundary) |
| Determinism / numerics | PASS: fixed grid, bit-identical reruns, converged to < 1e-6 against a 4× finer grid |
| PhD-level QA | PASS: needs double-layer theory, Lifshitz-retardation approximations, and diffusion-limited coagulation with hydrodynamic interaction together |

**Verdict: PASS**

---

## Files

```
task-51/
├── oracle/implement.py          reference oracle (hidden)
├── solution/example_solver.py   probe + fit solver (16 probes, 51/51)
├── grader/grade_submission.py   tolerance 0.001 + 1e-4·|expected|, probe budget 200
├── golden/test_data.json        51 cases (11 edge + 40 random)
├── scripts/helpers.py           shared utilities, twist-freezing variant()
├── scripts/build_golden.py      regenerates golden data
├── scripts/verify_task.py       automated GATE 2 checks
├── scripts/run_tests.py         verify + oracle self-grade + solver grade
└── PROPOSAL.md
```

## References

- Hogg, R., Healy, T. W. & Fuerstenau, D. W. (1966). Mutual coagulation of colloidal dispersions. *Trans. Faraday Soc.* 62, 1638-1651. doi:10.1039/TF9666201638
- Gregory, J. (1981). Approximate expressions for retarded van der Waals interaction. *J. Colloid Interface Sci.* 83, 138-145. doi:10.1016/0021-9797(81)90018-7
- Honig, E. P., Roebersen, G. J. & Wiersema, P. H. (1971). Effect of hydrodynamic interaction on the coagulation rate of hydrophobic colloids. *J. Colloid Interface Sci.* 36, 97-109.
- Fuchs, N. (1934). Über die Stabilität und Aufladung der Aerosole. *Z. Phys.* 89, 736-743.
- McGown, D. N. L. & Parfitt, G. D. (1967). Improved theoretical calculation of the stability ratio for colloidal systems. *J. Phys. Chem.* 71, 449-450.
- Derjaguin, B. (1934). Untersuchungen über die Reibung und Adhäsion, IV. *Kolloid-Z.* 69, 155-164.
- Malmberg, C. G. & Maryott, A. A. (1956). Dielectric constant of water from 0° to 100 °C. *J. Res. Natl. Bur. Stand.* 56, 1-8. doi:10.6028/jres.056.001
- Elimelech, M., Gregory, J., Jia, X. & Williams, R. A. (1995). *Particle Deposition and Aggregation*. Butterworth-Heinemann.
- CODATA 2022 recommended values, NIST. https://physics.nist.gov/cuu/Constants/
- Open-access verification sources: Jin, C. et al. *Sci. Rep.* 5, 17747 (2015); PMC10173465 (*Langmuir* 2023); Urbina-Villalba, G. et al. *Int. J. Mol. Sci.* 10, 761 (2009), PMC2672001; archive.org jresv56n1p1.
