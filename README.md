# PanelLocalProjections.jl

[![CI](https://github.com/andrerecio/PanelLocalProjections.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/andrerecio/PanelLocalProjections.jl/actions/workflows/CI.yml)
[![docs: stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://andrearecine.me/PanelLocalProjections.jl/stable/)
[![docs: dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://andrearecine.me/PanelLocalProjections.jl/dev/)

Panel local projections in Julia. For each horizon `h`, the package estimates a panel
regression of `y_{i,t+h}` on a shock, and gathers the estimators and inference
procedures of the recent panel-LP literature in one place:

| Method | Reference |
|---|---|
| Fixed-effects (within/LSDV) panel LP, with unit-, time- and two-way-clustered, heteroskedasticity-robust, Driscoll–Kraay and equal-weighted cosine (EWC) standard errors | Jordà (2005); Driscoll & Kraay (1998); Cameron, Gelbach & Miller (2011); Lazarus, Lewis, Stock & Watson (2018) |
| Split-panel jackknife (SPJ) correction of the Nickell bias of the FE estimator | Mei, Sheng & Shi (2026, *JIE*) |
| Cumulative-difference OLS (CD-OLS), which avoids the within transformation | Ugarte-Ruiz (2026, BBVA Research) |
| t-LAHR inference for aggregate (macro) shocks: time-clustered, lag-augmented, heteroskedasticity-robust, with an HC2 + Imbens–Kolesár refinement | Almuzara & Sancibrián (2026) |
| Panel LP-IV (2SLS, with lag augmentation and t-LAHR), with Kleibergen–Paap, Montiel Olea–Pflueger and Cragg–Donald first-stage diagnostics | Almuzara & Sancibrián (2026, §3.4); Kleibergen & Paap (2006); Montiel Olea & Pflueger (2013) |
| LP-DiD: difference-in-differences by local projections on clean-control samples, for absorbing, non-absorbing and one-off treatments, with reweighting, regression adjustment, pre-mean differencing, pooled estimates and a pre-trend test | Dube, Girardi, Jordà & Taylor (2025, *JAE*) |
| Pre-trends (negative horizons) and the response holding the later path of the treatment fixed, for a continuous, instrumented treatment | Ilzetzki (2024, *AER*) |

Version 0.2.2, for Julia 1.10 (LTS) and later. What comes after it (joint bands across
horizons, weak-IV-robust inference) is listed in [`TODO.md`](TODO.md);
the changes of each version are in [`CHANGELOG.md`](CHANGELOG.md). The
[documentation](#documentation) covers the methods, the conventions of each reference
implementation and how `panellp` reproduces them, the cross-checks, and audits of the
reference code.

## Installation

The package is not in the General registry. Install it from GitHub:

```julia
using Pkg
Pkg.add(url = "https://github.com/andrerecio/PanelLocalProjections.jl")                   # main
Pkg.add(url = "https://github.com/andrerecio/PanelLocalProjections.jl", rev = "v0.2.2")   # or a tagged release
```

It depends on StatsModels, Tables, StatsAPI, StatsFuns, RecipesBase and Julia's standard
library. Makie is optional (a package extension, for `irfplot`).

## Quick start

`panellp(formula, data; id, time, horizon)` returns one fit per horizon. `data` is any
[Tables.jl](https://github.com/JuliaData/Tables.jl) table in long format (a `DataFrame`, a
`NamedTuple` of vectors), one row per unit and period, in any order; `id` and `time` name
the unit and period columns. Periods are integers (consecutive integers are consecutive
periods), or dates with `frequency = :year`, `:halfyear`, `:quarter` or `:month`. Gaps and
missing values are allowed: an undefined lead or lag drops the row from that horizon's
sample.

This block runs as pasted:

```julia
using PanelLocalProjections, Random

# A simulated panel: 50 units, 40 periods, y_it = μ_i + 0.5 y_{i,t-1} + x_it + e_it
rng = Xoshiro(1)
N, T = 50, 40
id = repeat(1:N; inner = T)
t = repeat(1:T; outer = N)
x = randn(rng, N * T)
y = zeros(N * T)
for i in 1:N
    μ = randn(rng)
    for p in 1:T
        k = (i - 1) * T + p
        y[k] = μ + (p == 1 ? 0.0 : 0.5 * y[k - 1]) + x[k] + randn(rng)
    end
end
data = (; id, t, y, x)

# y_{i,t+h} on x_it for h = 0, …, 6, with unit fixed effects
r = panellp(@formula(leads(y) ~ x + lags(x, 1) + lags(y, 1) + fe(id)), data;
            id = :id, time = :t, horizon = 6)
s = summarize(r, Clustered(:unit))
```

```
Impulse response of leads(y) to x
Variance: Clustered(:unit; ssc = :stata); 95% bands
h        coef         se        lower      upper  df  nobs
0     0.98566  0.0219284     0.941593    1.02973  49  1950
1    0.511577  0.0310447     0.449191   0.573964  49  1900
2    0.232064  0.0432445     0.145161   0.318967  49  1850
3    0.104635  0.0403817    0.0234849   0.185785  49  1800
4   0.0660779  0.0376194  -0.00952113   0.141677  49  1750
5  -0.0173723  0.0394707   -0.0966917   0.061947  49  1700
6  -0.0250721  0.0444163     -0.11433  0.0641857  49  1650
```

The true response is `0.5^h`. The output is from Julia 1.13; the random stream, and so
the digits, change with the Julia version.

The formula helpers:

| Side | Terms |
|---|---|
| response | `leads(y)` (`y_{i,t+h}`), `ldiff(y)` (`y_{i,t+h} − y_{i,t−1}`), `anchor(y, z)` (`y_{i,t+h} − z_{i,t}`), `cumul(y)` (`Σ_{k≤h} y_{i,t+k}`), or a pre-built column with `horizon = 0` |
| regressors | columns and functions of columns, `lags(x, 1:p)`, `leads(x, 1:p)`, `firstdiff(x)`, interactions `s & X`, fixed effects `fe(a)` and `fe(a) & fe(b)` (absorbed), an IV block `(endo ~ inst)`, and `leads(x)` for the path `x_{i,t+1}, …, x_{i,t+h}` at horizon `h` |

`horizon = 6` means the horizons `0:6`; a range may start below zero (`horizon = -4:8`)
for pre-trends. `weights = :w` gives weighted least squares, and `pooled = true` one
regression of the response averaged over the horizons.

## Results

Variances are computed after estimation, from score sums stored at every horizon, so
several can be compared without re-estimating:

```julia
summarize(r, DriscollKraaySE(); level = 0.90)   # another variance, other bands
summarize(r, Clustered(:unit); term = "L1.y")   # another coefficient

coefnames(r)                           # "x", "L1.x", "L1.y"
coefpath(r)                            # the shock's coefficient at every horizon
coef(r, 2), nobs(r, 2)                 # all coefficients, and the observations, at h = 2
horizons(r), nobspath(r), esample(r, 2)

cv = vcov(Clustered(:unit, :time), r)  # variances at every horizon
vcov(cv, 2)                            # the matrix at h = 2
stderror(cv; term = "x")               # standard errors across horizons

s.horizon, s.coef, s.se, s.lower, s.upper, s.df, s.nobs
```

A summary is a Tables.jl table, and it plots with Makie or Plots.jl:

```julia
using DataFrames
DataFrame(s)

using CairoMakie      # loads the Makie extension
irfplot(s)            # a Figure; irfplot!(ax, s; color, bandalpha) draws on an axis

import Plots
Plots.plot(s)         # the same with Plots.jl (Makie and Plots both export `plot`)
```

## Estimators and variances

| `estimator =` | What it fits | What it needs |
|---|---|---|
| `WithinFE()` (default) | OLS after absorbing the `fe()` terms; pooled OLS with an intercept when there are none | — |
| `SplitPanelJackknife()` | the fixed-effects estimator with the split-panel jackknife correction | `fe(id)` in the formula; no IV block |
| `CDOLS()` | pooled OLS of the long difference on the shock and its first lag | the response `ldiff(y)` (or `cumul(g)` of a growth rate), `lags(shock, 1)` among the regressors, no unit fixed effects (time effects are allowed); no IV block. `CDOLS(validate = false)` skips the checks |
| `LagAugmented(rule)` | fixed effects, with `p(h)` lags of the outcome and of every shock added at horizon `h` | unit fixed effects. The default rule `ICLags()` needs the first shock to be an aggregate series, alone or times an exposure (`X`, `s & X`); `ConstantLags(p)` and `CappedLags(pmax)` work with any shock |
| an IV block `(endo ~ inst)` | two-stage least squares at every horizon; `firststage(r)` gives the diagnostics | `WithinFE()` or `LagAugmented(rule)` |
| `LPDiD(…)` | the long difference on the treatment switch, on the newly treated rows and the clean controls of each horizon | the formula `ldiff(y) ~ firstdiff(D) + fe(time)` with a 0/1 treatment `D`; no IV block. `nonabsorbing = L` for a treatment that can switch off |

| Variance | What it is | Works with |
|---|---|---|
| `Robust()` | heteroskedasticity-robust (HC1) | every estimator |
| `Clustered(:unit)`, `Clustered(:time)`, `Clustered(:unit, :time)`, `Clustered(:region)` | one- or two-way cluster-robust | every estimator. A column that is not the unit, the period or a fixed effect is declared at estimation: `panellp(…; clusters = (:region,))` |
| `DriscollKraaySE(bw)` | Driscoll–Kraay, Bartlett kernel | every estimator |
| `CosineHAR(B)` | equal-weighted cosine (EWC) | every estimator |
| `TLAHR(refine = false)` | clustered by period, no small-sample factor, normal critical values; the t-LAHR variance of an LP-IV | every estimator but `SplitPanelJackknife`; meant for `LagAugmented` |
| `TLAHR()` | t-LAHR with the HC2 + Imbens–Kolesár refinement | unweighted OLS fits (no IV block, not `SplitPanelJackknife`) whose shocks are an aggregate series, alone or times an exposure; meant for `LagAugmented` |
| `Homoskedastic()` | classical, `σ̂²(X̃'X̃)⁻¹`; a compatibility option (`regress`, `ivreg2` without `robust`) | least-squares and 2SLS fits (not `SplitPanelJackknife`, not regression-adjustment LP-DiD) |

Small-sample factors, critical values, the Driscoll–Kraay bandwidth and the repair of a
two-way variance are options of each specification (`ssc`, `psd`, `timelags`); see
[Inference](https://andrearecine.me/PanelLocalProjections.jl/stable/guide/inference/).

## Usage by method

`data` and `firms` stand for your panel.

```julia
# Mei–Sheng–Shi: fixed effects against the split-panel jackknife, unit-clustered standard errors
f     = @formula(leads(y) ~ x + lags(x, 1:4) + lags(y, 1:4) + fe(id))
r_fe  = panellp(f, data; id = :id, time = :year, horizon = 10)
r_spj = panellp(f, data; id = :id, time = :year, horizon = 10,
                estimator = SplitPanelJackknife())
summarize(r_spj, Clustered(:unit))

# Ugarte-Ruiz: cumulative-difference OLS
r_cd  = panellp(@formula(ldiff(y) ~ s + lags(s, 1)), data; id = :id, time = :year,
                horizon = 9, estimator = CDOLS())
summarize(r_cd, Clustered(:unit))

# Almuzara–Sancibrián: an aggregate shock mp interacted with an exposure lev, t-LAHR
# inference. `shock` is required when the formula has an interaction and other regressors
r_la  = panellp(@formula(ldiff(logk) ~ lev & mp + lev + fe(firm) + fe(sector) & fe(date)),
                firms; id = :firm, time = :date, frequency = :quarter, horizon = 12,
                shock = "lev & mp", estimator = LagAugmented(ICLags()))
summarize(r_la, TLAHR(); level = 0.90)
lagorder(r_la)                             # the lags p(h) added at each horizon

# Panel LP-IV: s × ffr instrumented by s × mp, lag-augmented
r_iv  = panellp(@formula(leads(y) ~ (s & ffr ~ s & mp) + fe(id) + fe(t)), data;
                id = :id, time = :t, horizon = 12, estimator = LagAugmented())
summarize(r_iv, TLAHR(refine = false))     # TLAHR() is not defined for IV
firststage(r_iv)                           # Kleibergen–Paap F at every horizon
```

In the third call the periods are dates, so `fe(date)` and clusters on `date` group by
quarter. The Montiel Olea–Pflueger effective F is reported for one endogenous regressor;
with `LagAugmented` the lags of the endogenous shock are endogenous too, and it is `NaN`.

### LP-DiD

A binary treatment `D` that units enter at different dates (Dube, Girardi, Jordà & Taylor
2025). The formula is the regression; the estimator picks, at each horizon, the rows
entering treatment and the clean controls (not yet treated at `t + h`):

```julia
f  = @formula(ldiff(y) ~ firstdiff(D) + fe(t))
es = panellp(f, data; id = :id, time = :t, horizon = -5:8, estimator = LPDiD())
summarize(es, Clustered(:unit))            # event study; h = −1 is zero by construction
pretrendtest(es)                           # joint test of the coefficients at h ≤ −2
averageeffect(es)                          # average over the post- and pre-treatment windows

# the equally weighted effect on the treated, with a covariate and lagged outcome growth
# (regression adjustment: its scores are the influence function, so no small-sample factor)
rw = panellp(@formula(ldiff(y) ~ firstdiff(D) + x + lags(firstdiff(y), 1) + fe(t)), data;
             id = :id, time = :t, horizon = -5:8, estimator = LPDiD(reweight = true))
summarize(rw, Clustered(:unit; ssc = :none))

# one pooled estimate over the post-treatment window; the same rows at every horizon
panellp(f, data; id = :id, time = :t, horizon = 0:8, pooled = true, estimator = LPDiD())
panellp(f, data; id = :id, time = :t, horizon = -5:8, sample = :common, estimator = LPDiD())

# a treatment that switches on and off, with effects assumed to stabilise within 5 periods
panellp(@formula(ldiff(y2) ~ firstdiff(D2) + fe(t)), data; id = :id, time = :t,
        horizon = -5:8, estimator = LPDiD(nonabsorbing = 5))
```

Other options: `control = :nevertreated` or `:notyet`, `firsttreat`, `oneoff` (a treatment
that lasts one period), `untreated_before`, `pmd = k` or `:max` (pre-mean differencing),
`weights`. It is written from the paper and reproduces the Stata command `lpdid` 1.1.0
option by option; the [LP-DiD page](https://andrearecine.me/PanelLocalProjections.jl/stable/methods/lpdid/) has the map between the
two, and the [audit](https://andrearecine.me/PanelLocalProjections.jl/stable/audit/lpdid/) the places where they differ (the
reweighted estimate on unbalanced panels, among them).

### Pre-trends and the path of the treatment

For a continuous, instrumented treatment, as in Ilzetzki (2024): horizons from −6, and at
`h > 0` the leads `x_{t+1}, …, x_{t+h}` as controls, so that the coefficient of `x_t` is
the response holding later treatment fixed.

```julia
r = panellp(@formula(ldiff(y) ~ (x ~ z) + lags(x, 1:6) + leads(x) + fe(id) + fe(t)), data;
            id = :id, time = :t, horizon = -6:18)
summarize(r, Clustered(:unit))
firststage(r; vcov = Clustered(:unit))     # Kleibergen–Paap, effective and Cragg–Donald F
```

The [replication page](https://andrearecine.me/PanelLocalProjections.jl/stable/replication/ilzetzki2024/) reproduces his Figure 9 on
the replication data (`singletons = :keep` and `Homoskedastic(ssc = :none)` for the
published numbers).

## Example: monetary policy and US states

How much more does personal income fall in states more exposed to manufacturing after a
rise in the 1-year Treasury rate? The rate is instrumented by monetary policy surprises,
a panel LP-IV in Almuzara & Sancibrián's setting. The data are public: BEA state personal
income, FRED's 1-year rate (`GS1`), and the Jarociński–Karadi surprises.
`examples/us_states/build.jl` downloads them and builds the panel
(`julia --project examples/us_states/build.jl`, from the repository root); the plot comes
from `examples/us_states/figure.jl`.

```julia
using PanelLocalProjections
include("examples/us_states/csv.jl")       # readtable
d = readtable("docs/data/us_states/processed/state_panel.csv")

# ldiff(lpi): 100 × (log income_{t+h} − log income_{t−1}); s: 1998 manufacturing share
# (standardised); gs1: 1-year rate; mp: surprise; state and quarter effects
iv = panellp(@formula(ldiff(lpi) ~ (s & firstdiff(gs1) ~ s & mp) + fe(fips) + fe(t)), d;
             id = :fips, time = :t, horizon = 12, estimator = LagAugmented())
summarize(iv, TLAHR(refine = false); level = 0.90)    # t-LAHR for LP-IV (A&S A.4)
summarize(iv, Clustered(:unit); level = 0.90)         # clustered by state, for comparison
firststage(iv).kpF                                    # Kleibergen–Paap F, clustered by period

using CairoMakie
irfplot(summarize(iv, TLAHR(refine = false); level = 0.90))
```

![LP-IV and reduced-form responses of US state personal income, with t-LAHR and state-clustered 90% bands](docs/src/assets/us_states_irf.png)

Clustered by state, the responses look significant at several horizons. The t-LAHR
standard errors are 1.5 to 2.7 times those clustered by state, and none of the t-LAHR
bands excludes zero. The first stage shows the same thing: its Kleibergen–Paap F is about
2 clustered by period, but 35 robust and 200,000 clustered by state. The shock is
aggregate, so the effective sample size is the 84 quarters, not the 4,284 state-quarters.
See the [example page](https://andrearecine.me/PanelLocalProjections.jl/stable/examples/us_states/) for the details.

Data: BEA Regional Economic Accounts; Board of Governors of the Federal Reserve System
(H.15), via FRED; Jarociński & Karadi (2020, *AEJ: Macroeconomics*), updated surprises
(CC BY 4.0). The data are not distributed with the package.

## Conventions

Leads, lags and differences are always computed on the (unit, calendar period) grid, so
gaps in a unit's time series are respected and nothing is shifted by row position. Every
convention on which reference implementations disagree is an explicit option, and each
is cross-checked against the reference that uses it: degrees-of-freedom factors,
singleton handling, two-way demeaning, the Driscoll–Kraay bandwidth and lag ordering, the
split-panel rule and critical values. The references are the Stata commands `xtlp`,
`reghdfe`, `xtreg`, `regress`, `xtscc`, `locproj`, `xtspj`, `ivreghdfe`, `ivreg2`,
`xtivreg2`, `weakivtest`, `harreg` and `lpdid`, R `fixest`, and the authors' own R code (Mei–Sheng–Shi's
`LP_panel`, Almuzara–Sancibrián's `panel_LP.R`). [`crosscheck/README.md`](crosscheck/README.md)
maps each reference command to the `panellp` options that reproduce it, and
[Deviations](https://andrearecine.me/PanelLocalProjections.jl/stable/deviations/) lists the defaults that differ from a reference.

## Documentation

The documentation is at <https://andrearecine.me/PanelLocalProjections.jl/stable/> for the latest release and <https://andrearecine.me/PanelLocalProjections.jl/dev/> for `main`:

- Guide: [getting started](https://andrearecine.me/PanelLocalProjections.jl/stable/guide/example/), [panel data](https://andrearecine.me/PanelLocalProjections.jl/stable/guide/data/),
  [formula](https://andrearecine.me/PanelLocalProjections.jl/stable/guide/formula/), [inference](https://andrearecine.me/PanelLocalProjections.jl/stable/guide/inference/).
- Methods: [fixed-effect absorption](https://andrearecine.me/PanelLocalProjections.jl/stable/methods/absorption/),
  [split-panel jackknife](https://andrearecine.me/PanelLocalProjections.jl/stable/methods/spj/), [CD-OLS](https://andrearecine.me/PanelLocalProjections.jl/stable/methods/cdols/),
  [t-LAHR](https://andrearecine.me/PanelLocalProjections.jl/stable/methods/tlahr/), [panel LP-IV](https://andrearecine.me/PanelLocalProjections.jl/stable/methods/lpiv/),
  [lead paths and negative horizons](https://andrearecine.me/PanelLocalProjections.jl/stable/methods/leads/),
  [LP-DiD](https://andrearecine.me/PanelLocalProjections.jl/stable/methods/lpdid/).
- Example and replications: [US states](https://andrearecine.me/PanelLocalProjections.jl/stable/examples/us_states/),
  [Mei, Sheng & Shi](https://andrearecine.me/PanelLocalProjections.jl/stable/replication/mss2026/),
  [Ugarte-Ruiz](https://andrearecine.me/PanelLocalProjections.jl/stable/replication/ugarte2026/),
  [Ilzetzki](https://andrearecine.me/PanelLocalProjections.jl/stable/replication/ilzetzki2024/) (on his replication data, not
  distributed),
  [firm-level studies](https://andrearecine.me/PanelLocalProjections.jl/stable/replication/firm_level_wip/) (call mappings, not run on
  the data).
- [Cross-checks](https://andrearecine.me/PanelLocalProjections.jl/stable/crosscheck/), [deviations](https://andrearecine.me/PanelLocalProjections.jl/stable/deviations/) and
  audits of the reference code, starting with [`xtlp`](https://andrearecine.me/PanelLocalProjections.jl/stable/audit/xtlp/) and [`lpdid`](https://andrearecine.me/PanelLocalProjections.jl/stable/audit/lpdid/).
- API: every exported name has a docstring (`?panellp`); the list is in
  the [API page](https://andrearecine.me/PanelLocalProjections.jl/stable/api/).

The sources are in [`docs/src`](docs/src). To build the HTML pages locally (Julia 1.11 or later), into `docs/build`:

```bash
julia --project=docs -e 'using Pkg; Pkg.instantiate()'
julia --project=docs docs/make.jl
```

## Acknowledgements and inspiration

The fixed-effect absorption, the OLS solver and every variance estimator are implemented
in this package. Their design takes inspiration from:

- [Regress.jl](https://github.com/gragusa/Regress.jl) and
  [CovarianceMatrices.jl](https://github.com/gragusa/CovarianceMatrices.jl) by Giuseppe
  Ragusa;
- [FixedEffectModels.jl](https://github.com/FixedEffects/FixedEffectModels.jl) by
  Matthieu Gomez;
- [LocalProjections.jl](https://github.com/gragusa/LocalProjections.jl) (and the
  [andrerecio fork](https://github.com/andrerecio/LocalProjections.jl)), whose formula
  interface (`leads`, `cumul`, `ldiff`, …) this package extends to panels;
- the Stata commands [`reghdfe`](https://github.com/sergiocorreia/reghdfe) (Sergio
  Correia) and [`xtlp`](https://github.com/shenshuuu/panel-local-projection-stata)
  (Mei, Sheng & Shi), and R [`fixest`](https://github.com/lrberge/fixest) (Laurent Bergé).

None of them is a dependency. Source files that follow one of these projects say so in a
comment.

## Citing

If you use the package, please cite it together with the paper whose method you apply;
see [`CITATION.bib`](CITATION.bib).

- Mei, Z., Sheng, L. & Shi, Z. (2026). Nickell bias in panel local projection: Financial
  crises are worse than you think. *Journal of International Economics*, 104210.
- Ugarte-Ruiz, A. (2026). Panel local projections without fixed-effects: A
  cumulative-difference OLS estimator. BBVA Research.
- Almuzara, M. & Sancibrián, V. (2026). Micro responses to macro shocks. Working paper,
  July 2026; earlier version: Federal Reserve Bank of New York Staff Report No. 1090.
- Dube, A., Girardi, D., Jordà, Ò. & Taylor, A. M. (2025). A local projections approach to
  difference-in-differences. *Journal of Applied Econometrics* 40(7), 741–758.
- Ilzetzki, E. (2024). Learning by necessity: Government demand, capacity constraints, and
  productivity growth. *American Economic Review*.

## License

MIT; see [`LICENSE`](LICENSE).
