import sys, math, statistics
sys.path.insert(0, 'scripts')
from helpers import load_golden, load_oracle, variant, tolerance
o = load_oracle(); cases = load_golden()['cases']

def ro(c, psi, a_nm, A):  # Reerink-Overbeek: W ~ exp(Vmax/kT) / (2 kappa a), textbook unretarded
    kt = o.K_B*o.TEMP; a=a_nm*o.NM; k=o.debye_kappa(c); vmax=0.0
    for i in range(1, 20000):
        h = i*0.005e-9
        v = (o.v_rep(h,a,psi*o.MV,k) - A*o.ZJ*a/(12*h))/kt
        vmax = max(vmax, v)
    return max(0.0, (vmax - math.log(2*k*a))/math.log(10))

def eps_alt(c, psi, a, A):
    old = o.EPS_R; o.EPS_R = 78.54
    try: return o.compute(c, psi, a, A)
    finally: o.EPS_R = old

models = {
  'Textbook DLVO + Fuchs (no T1, no T2)': lambda *p: variant(*p, retard=False, hydro=False),
  'T1 only (no hydrodynamics)': lambda *p: variant(*p, hydro=False),
  'T2 only (unretarded attraction)': lambda *p: variant(*p, retard=False),
  'Reerink-Overbeek barrier estimate': ro,
  'Full model, eps_r = 78.54 instead of 78.30': eps_alt,
}
for name, f in models.items():
    errs=[]; ok=0
    for c in cases:
        p = tuple(c['input'][k] for k in ('value_a','value_b','value_c','value_d'))
        y = c['expected']['output']; g = f(*p); e = abs(g-y); errs.append(e); ok += e <= tolerance(y)
    print(f"| {name} | {ok}/51 | {statistics.median(errs):.3f} | {max(errs):.2f} |")
