"""
Sensitivity and supplemental models reported in the manuscript (2023 RUCC).
Run after clean_panel.py, from the main repository folder:  python scripts/supplemental_analyses.py
Writes results/supplemental_results.txt

All models: county clinician counts, log(population) offset, calendar-year and Census-region
fixed effects, standard errors clustered by county, unless noted.
  A. Primary care physicians as the comparator (instead of all physicians)
  B. Negative binomial version of the main model                      (Supplemental Table S1)
  C. RUCC entered as a categorical variable                            (Supplemental Table S1)
  D. Year x RUCC x profession interaction (did gradients change?)
  E. Nonmetro counties only: urban population size and metro adjacency (Supplemental Table S2)
"""
import io
import os
import warnings

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf

warnings.filterwarnings("ignore")
buf = io.StringIO()


def out(*a):
    print(*a)
    print(*a, file=buf)


d = pd.read_csv("data/AHRF_allphys_NP_PA_panel_2018_2023_clean.csv", dtype={"fips": str})


def to_long(x, profs):
    """One row per county-year-profession."""
    return pd.concat([x.assign(prof=p, count=x[c]) for p, c in profs], ignore_index=True)


def poisson(formula, x):
    return smf.glm(formula, data=x, family=sm.families.Poisson(), offset=np.log(x.population)
                   ).fit(cov_type="cluster", cov_kwds={"groups": pd.factorize(x.fips)[0]})


def combo(m, terms):
    """exp(sum of coefficients) with 95% CI, using the model's (clustered) covariance."""
    est = sum(m.params[t] for t in terms)
    V = m.cov_params().loc[terms, terms].values
    se = np.sqrt(V.sum())
    return np.exp(est), np.exp(est - 1.96 * se), np.exp(est + 1.96 * se)


def fmt(rr, lo, hi):
    return f"{rr:.3f} ({lo:.3f}-{hi:.3f})"


ALL = [("phys", "phys_all"), ("np", "np"), ("pa", "pa")]
P = 'C(prof, Treatment("phys"))'
long = to_long(d, ALL)
out(f"Panel: {d.fips.nunique():,} counties x {d.year.nunique()} years. Long file: {len(long):,} rows.\n")

# A. Primary care physician comparator -------------------------------------------------
lp = to_long(d, [("pcp", "pcphys"), ("np", "np"), ("pa", "pa")])
m = poisson('count ~ C(prof, Treatment("pcp")) * rucc23 + C(year) + C(region)', lp)
out("A. Primary care physicians as comparator (Poisson, per one-category RUCC increase)")
out(f"   PCPs: rate ratio {fmt(*combo(m, ['rucc23']))}, {(np.exp(m.params['rucc23']) - 1) * 100:+.1f}%")
for p in ["np", "pa"]:
    k = f'C(prof, Treatment("pcp"))[T.{p}]:rucc23'
    rr = combo(m, ["rucc23", k])
    lo, hi = np.exp(m.conf_int().loc[k])
    out(f"   {p.upper()}s: rate ratio {fmt(*rr)}, {(rr[0] - 1) * 100:+.1f}%; "
        f"ratio vs PCPs {np.exp(m.params[k]):.3f} ({lo:.3f}-{hi:.3f}), p={m.pvalues[k]:.3f}")
y = d[d.year == 2023]
for lab, s in [("nonmetro", y[y.nonmetro]), ("metro", y[y.metro])]:
    out(f"   2023 NPs per primary care physician, {lab}: {s.np.sum() / s.pcphys.sum():.2f}")
out()

# B. Negative binomial --------------------------------------------------------------------
nb = smf.negativebinomial(f"count ~ {P} * rucc23 + C(year) + C(region)", data=long,
                          offset=np.log(long.population)
                          ).fit(cov_type="cluster", cov_kwds={"groups": pd.factorize(long.fips)[0]},
                                maxiter=200, disp=0)
out("B. Negative binomial (NB2), per one-category RUCC increase")
out(f"   Physicians: rate ratio {fmt(*combo(nb, ['rucc23']))}")
for p in ["np", "pa"]:
    k = f"{P}[T.{p}]:rucc23"
    lo, hi = np.exp(nb.conf_int().loc[k])
    out(f"   {p.upper()}s: rate ratio {fmt(*combo(nb, ['rucc23', k]))}; "
        f"ratio vs physicians {np.exp(nb.params[k]):.3f} ({lo:.3f}-{hi:.3f}), p={nb.pvalues[k]:.1e}")
out(f"   Dispersion (alpha): {nb.params['alpha']:.3f}\n")

# C. RUCC as categorical ----------------------------------------------------------------
m = poisson(f"count ~ {P} * C(rucc23) + C(year) + C(region)", long)
rows = []
for r in range(2, 10):
    t = f"C(rucc23)[T.{r}]"
    row = {"RUCC": r, "Physicians": fmt(*combo(m, [t]))}
    for p, lab in [("np", "NPs"), ("pa", "PAs")]:
        row[lab] = fmt(*combo(m, [t, f"{P}[T.{p}]:{t}"]))
    rows.append(row)
out("C. RUCC categorical (Poisson): rate ratio vs RUCC 1, by profession")
out(pd.DataFrame(rows).set_index("RUCC").to_string(), "\n")

# D. Did the gradients change over time? -------------------------------------------------
long["year_c"] = long.year - 2018
m = poisson(f"count ~ {P} * rucc23 * year_c + C(year) + C(region)", long)
out("D. Year x RUCC x profession (year centered at 2018)")
for p in ["np", "pa"]:
    k = f"{P}[T.{p}]:rucc23:year_c"
    out(f"   {p.upper()} vs physician gradient, change per year: ratio {np.exp(m.params[k]):.4f}, p={m.pvalues[k]:.3g}")
out()

# E. Nonmetro only: urban population size and adjacency -------------------------------------
n = d[d.nonmetro].copy()
n["size"] = n.rucc23.map({4: "20k+", 5: "20k+", 6: "5k-20k", 7: "5k-20k", 8: "<5k", 9: "<5k"})
n["nonadjacent"] = n.rucc23.isin([5, 7, 9]).astype(int)
ln = to_long(n, ALL)
S = 'C(size, Treatment("20k+"))'
m = poisson(f"count ~ {P} * {S} + {P} * nonadjacent + C(year) + C(region)", ln)
out("E. Nonmetro counties only: rate ratios for urban population size (ref 20,000+) and not adjacent to metro")
for t, lab in [(f"{S}[T.5k-20k]", "Urban pop 5,000-20,000"), (f"{S}[T.<5k]", "Urban pop <5,000"),
               ("nonadjacent", "Not adjacent to metro")]:
    parts = [f"physicians {fmt(*combo(m, [t]))}"]
    for p in ["np", "pa"]:
        k = f"{P}[T.{p}]:{t}"
        parts.append(f"{p.upper()}s {fmt(*combo(m, [t, k]))} [vs phys p={m.pvalues[k]:.3f}]")
    out(f"   {lab}: " + "; ".join(parts))

os.makedirs("results", exist_ok=True)
open("results/supplemental_results.txt", "w").write(buf.getvalue())
