"""Reference oracle (hidden from solver).

log10 of the Fuchs stability ratio W of a monodisperse suspension of equal
charged spheres in an aqueous 1:1 electrolyte at 25 degC.

Pipeline (see PROPOSAL.md sections 2-4):
  S1  electrolyte screening   kappa^2 = 2 N_A e^2 (1000 c) / (eps_r eps_0 kT)        (Debye-Hueckel)
  S2  pair potential          V_R = 2 pi eps_r eps_0 a psi^2 ln(1 + exp(-kappa h))      (HHF 1966, equal spheres)
                              V_A = -(A a / 12 h) / (1 + 14 h / lambda), lambda = 100 nm
                                    (Gregory 1981 retarded form, Derjaguin a -> a/2)    (T1)
  S3  diffusional resistance  beta(u) = (6u^2 + 13u + 2) / (6u^2 + 4u), u = h / a
                                    (Honig, Roebersen & Wiersema 1971)                   (T2)
  S4  stability ratio         W = int beta e^{V_T/kT} (2+u)^-2 du / int beta e^{V_A/kT} (2+u)^-2 du
                                    (Fuchs 1934 as completed by McGown & Parfitt 1967)
                              output = log10 W

Inputs (neutral names):
    value_a : 1:1 electrolyte concentration, mol L^-1   [0.001, 0.1]
    value_b : surface potential magnitude, mV           [15, 40]
    value_c : particle radius, nm                       [50, 500]
    value_d : Hamaker constant, zJ (1e-21 J)            [3, 20]
Output:
    output  : log10 W, rounded to 4 decimals
"""
import math

# --- SI exact defining constants and CODATA 2022 --------------------------------
E_CHARGE = 1.602176634e-19      # C, elementary charge (exact, SI 2019)
K_B = 1.380649e-23              # J K^-1, Boltzmann constant (exact, SI 2019)
N_A = 6.02214076e23             # mol^-1, Avogadro constant (exact, SI 2019)
EPS_0 = 8.8541878188e-12        # F m^-1, vacuum permittivity (CODATA 2022)
L_PER_M3 = 1000.0               # L m^-3 (exact): mol L^-1 -> mol m^-3

# --- Solvent: water at 25 degC ------------------------------------------------------
TEMP = 298.15                   # K, IUPAC standard temperature (25 degC)
EPS_R = 78.30                   # Malmberg & Maryott (1956), J. Res. NBS 56, 1: eps_r at 25 degC

# --- Gregory (1981), J. Colloid Interface Sci. 83, 138: retarded sphere-plate form --
GREGORY_B = 14.0                # V = -(A a / 6h) / (1 + 14 h / lambda)
LAMBDA = 100e-9                 # m, characteristic wavelength of the dispersion interaction

# --- Honig, Roebersen & Wiersema (1971), J. Colloid Interface Sci. 36, 97 -----------
HRW_A2, HRW_A1, HRW_A0 = 6.0, 13.0, 2.0     # numerator   6u^2 + 13u + 2
HRW_B2, HRW_B1 = 6.0, 4.0                   # denominator 6u^2 + 4u

# --- Unit conversions (exact) --------------------------------------------------------
MV = 1e-3                       # V per mV
NM = 1e-9                       # m per nm
ZJ = 1e-21                      # J per zJ

# --- Quadrature (structural; fixed -> bit-for-bit determinism) ----------------------
U_MIN, U_MAX = 1e-10, 1e6       # u = h/a range on a logarithmic grid; integrand ~ 0 below,
N_SIMPSON = 8192                #   analytic tail 1/(2 + U_MAX) above (beta -> 1, V -> 0)


def debye_kappa(c_molar):
    """Inverse Debye length (m^-1) of a symmetric 1:1 electrolyte."""
    return math.sqrt(2 * N_A * L_PER_M3 * c_molar * E_CHARGE ** 2 / (EPS_R * EPS_0 * K_B * TEMP))


def v_attr(h, a, hamaker):
    """Retarded van der Waals energy (J) between equal spheres (T1)."""
    return -hamaker * a / (12 * h) / (1 + GREGORY_B * h / LAMBDA)


def v_rep(h, a, psi, kappa):
    """Constant-potential double-layer energy (J) between equal spheres (HHF)."""
    return 2 * math.pi * EPS_R * EPS_0 * a * psi ** 2 * math.log1p(math.exp(-kappa * h))


def beta_hrw(u):
    """Relative diffusional resistance of two approaching equal spheres (T2)."""
    return (HRW_A2 * u * u + HRW_A1 * u + HRW_A0) / (HRW_B2 * u * u + HRW_B1 * u)


def _log_sum_exp(terms):
    m = max(terms)
    return m + math.log(math.fsum(math.exp(t - m) for t in terms))


def compute(value_a, value_b, value_c, value_d):
    c = float(value_a)
    psi = float(value_b) * MV
    a = float(value_c) * NM
    hamaker = float(value_d) * ZJ
    if not (0.001 <= c <= 0.1 and 15 <= value_b <= 40 and 50 <= value_c <= 500 and 3 <= value_d <= 20):
        raise ValueError("inputs out of domain")

    # S1: screening
    kappa = debye_kappa(c)
    kt = K_B * TEMP

    # S2-S4: both Fuchs integrals on s = ln u with composite Simpson, in log space
    s0, s1 = math.log(U_MIN), math.log(U_MAX)
    ds = (s1 - s0) / N_SIMPSON
    ln_num, ln_den = [], []
    for i in range(N_SIMPSON + 1):
        u = math.exp(s0 + i * ds)
        h = u * a
        weight = 1 if i in (0, N_SIMPSON) else (4 if i % 2 else 2)
        # du = u ds ;  dr / r^2 = du / (a (2 + u)^2) ; the common factor a cancels in W
        base = math.log(weight * beta_hrw(u) * u) - 2 * math.log(2 + u)
        va = v_attr(h, a, hamaker) / kt
        ln_den.append(base + va)
        ln_num.append(base + va + v_rep(h, a, psi, kappa) / kt)
    tail = -math.log(2 + U_MAX) - math.log(ds / 3)     # int_{U_MAX}^inf du/(2+u)^2, Simpson-scaled
    log_w = _log_sum_exp(ln_num + [tail]) - _log_sum_exp(ln_den + [tail])
    return round(log_w / math.log(10), 4)


def oracle(inputs):
    return {"output": compute(inputs["value_a"], inputs["value_b"], inputs["value_c"], inputs["value_d"])}


if __name__ == "__main__":
    import json
    import sys
    print(json.dumps(oracle(json.loads(sys.stdin.read()))))
