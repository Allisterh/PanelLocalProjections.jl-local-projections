# crosscheck/

Scripts that compute the reference numbers the Julia tests compare against: Stata through
PyStata, and R. The scripts are committed. The local toolchain and the run outputs
(`.venv/`, `crosscheck/ado/`, `crosscheck/out/`) are not.

## Setup (once per machine)

```bash
uv venv .venv --python 3.13
uv pip install --python .venv/bin/python -r crosscheck/requirements.txt
.venv/bin/python crosscheck/run_do.py crosscheck/stata/setup_stata.do --log crosscheck/out/setup_stata.log
.venv/bin/python crosscheck/ado_manifest.py --check     # installed files match the committed checksums
.venv/bin/python crosscheck/run_do.py crosscheck/stata/smoke.do
Rscript crosscheck/R/smoke.R
```

`setup_stata.do` installs every package into the project-local tree `crosscheck/ado/`.
`common.do` points `sysdir PLUS/PERSONAL` there, so the user's own Stata setup is never
touched. Versions and pins are in `VERSIONS.md`.

## Running a cross-check

```bash
.venv/bin/python crosscheck/run_do.py crosscheck/stata/<script>.do --log crosscheck/out/<script>.log
Rscript crosscheck/R/<script>.R
```

Each do-file starts with `do "$PLP_ROOT/crosscheck/stata/common.do"`. `run_do.py` sets
`PLP_ROOT` to the repository root and runs the do-file from its own directory, in a fresh
Stata session each time (~1–2 s). The exit code is the Stata return code.

## Conventions every script follows

- `set type double` (in `common.do`) and `import delimited …, asdouble`. Stata's float
  storage otherwise costs about 1e-7 relative error.
- Inputs are the committed synthetic panels in `test/fixtures/synthetic/`, or local-only
  data under `docs/references/` whose outputs go to `test/fixtures/local/` (gitignored).
- Outputs are written at full precision: Stata `%21x` or `%21.0g`; R `format(x, digits = 17)`
  or `data.table::fwrite`. Committed outputs go to `test/fixtures/{stata,r}/` and are
  listed with their provenance in `test/fixtures/README.md`.
- A re-run rewrites the committed fixture in place, so `git diff --exit-code test/fixtures`
  shows any change. Every committed fixture has reproduced byte for byte so far.
- fixest ≥ 0.13: `ssc(K.adj = , G.adj = )`, explicit `fixef.rm`, explicit `vcov =`.
- MSS `LP_panel` needs a base `data.frame` (`as.data.frame()` first). It is reliable only in
  its `H = 0` mode with pre-built leads and lags.

## Command ↔ Julia call map

Filled in phase by phase. Each row gives the exact `panellp`/`vcov` options that reproduce
a reference command: `singletons`, `timelags`, `ssc`, `psd`, `absorb`, `dist`.

### Phase 1: fixed effects

Every row is a committed test (`test/test_crosscheck_p1.jl`), at horizons 0–6 on both
synthetic panels unless noted. The test requires each reference spec to have its exact
number of rows and a Julia mapping. It compares coefficients, standard errors,
observations, critical degrees of freedom, numbers of clusters or periods, and reghdfe's
`e(df_a)` and `e(df_a_initial)`. The formula is `leads(y) ~ x + lags(x, 1) + lags(y, 1) + w`
plus the fixed effects listed. "drop"/"keep" is the `singletons` keyword. Agreement is
the largest relative difference over all horizons, terms and panels.

| Reference command (Stata `xtset id t`) | `panellp` options | `vcov` spec | Agreement |
|---|---|---|---|
| `xtlp F#.y …, fe m(fe) vce(cluster id)` | `fe(id)`, drop | `Clustered(:unit)` | b 5e-13, se 5e-15 |
| `xtlp …, fe m(fe) vce(cluster id, ase)` | `fe(id)`, drop | `Clustered(:unit; ssc = :none)` (xtlp keeps t(G−1)) | se 6e-15 |
| `xtlp …, fe m(fe) vce(robust)` | `fe(id)`, drop | `Robust()` | se 5e-15 |
| `xtlp …, fe m(fe) vce(cluster t)` | `fe(id)`, drop | `Clustered(:time)` | se 5e-15 |
| `xtlp …, fe m(fe) vce(cluster id t [, ase])` | `fe(id)`, drop | `Clustered(:unit, :time [; ssc = :none])` | se 6e-15 |
| `xtlp …, fe m(fe) vce(dkraay [lag(3)])` | `fe(id)`, drop | `DriscollKraaySE(:nw94 \| 3; timelags = :position)` | se 6e-15 |
| `xtlp …, vce(dkraay lag(3), nodfadj \| ase)` | `fe(id)`, drop | `DriscollKraaySE(3; ssc = :xtscc \| :none, timelags = :position)` | se 5e-15 |
| `xtlp …, keepsingletons` | keep | as above | se 5e-15 |
| `xtlp …, tfe m(fe) vce(…)`, the same `vce()` options and `keepsingletons` | `fe(id) + fe(t)`, drop / keep | as above | b 7e-14, se 2e-11 (xtlp stops on a 1e-9 absolute step) |
| `xtlp y …, fe m(fe) hor(0 6) ytransf(level \| diff \| cmltdiff \| cmltsum)` | LHS `leads(y)`, `leads(firstdiff(y))`, `ldiff(y)`, `cumul(y)` | `Clustered(:unit)` | 1e-13 |
| `xtlp F#.y x w, fe \| tfe m(fe) vce(cluster id) [keepsingletons]` (no lags: units with one or two rows enter, and singleton drops cascade) | `leads(y) ~ x + w + fe(id) [+ fe(t)]` | `Clustered(:unit)` | se 3e-12 |
| `reghdfe …, absorb(id t) vce(cluster id t) tol(1e-14)` | `fe(id) + fe(t)`, drop | `Clustered(:unit, :time)` | se 5e-11 |
| `reghdfe …, absorb(id) vce(cluster reg)` | `fe(id)`, `clusters = (:reg,)` | `Clustered(:reg)` | se 9e-15 |
| `reghdfe …, absorb(id t grp#t) vce(cluster id)` | `fe(id) + fe(t) + fe(grp) & fe(t)` | `Clustered(:unit)` | se 2e-10 |
| `reghdfe …, absorb(id [t]) vce(dkraay 4)` | `fe(id) [+ fe(t)]` | `DriscollKraaySE(3; timelags = :position)` | se 5e-12 |
| `reghdfe …, absorb(id t) vce(robust)` | `fe(id) + fe(t)` | `Robust()` | se 2e-11 |
| `reghdfe y x w, absorb(id t) vce(cluster id)` (no lags) | `fe(id) + fe(t)`, drop | `Clustered(:unit)` | se 2e-11 |
| `reghdfe …, absorb(id) vce(cluster reg t)`: not PSD, the repair fires | `fe(id)` | `Clustered(:reg, :time)` | se 5e-15 |
| `reghdfe …, noabsorb vce(cluster reg t)`: pooled, `_cons` reported, the repair fires (slope block clipped on its own) | no fixed effects (intercept) | `Clustered(:reg, :time)` | se 1e-14 where reghdfe posts a V (*) |
| `xtreg …, fe vce(cluster id)`, also with `i.t` | `fe(id) [+ fe(t)]`, keep | `Clustered(:unit)` | se 2e-14 |
| `regress …, vce(cluster id)` | no fixed effects (intercept) | `Clustered(:unit)` | se 8e-15 |
| `xtscc …, fe lag(3)` (fx_balanced: xtscc refuses a period absent from the data for every unit) | `fe(id)`, keep | `DriscollKraaySE(3; ssc = :xtscc, timelags = :position)` | se 6e-15 |
| Derived panel `_t29` (t ≤ 29): `xtlp …, fe \| tfe m(fe) vce(dkraay)`, the default lags change with each horizon's T (3, then 2) | `fe(id) [+ fe(t)]`, drop | `DriscollKraaySE(:nw94; timelags = :position)` | se 2e-11 |
| Derived panel `fx_balanced_gap` (w missing at t = 20 for every unit, unit 30 at t = 10, 11 only): `xtscc …, fe lag(3)`; `xtlp …, fe m(fe) vce(dkraay lag(3))` | keep (xtscc), drop (xtlp) | `DriscollKraaySE(3; …, timelags = :position)` (calendar lags differ here) | se 6e-15 |
| Derived panel `fx_balanced_2blocks` (two connected components of unit and time effects): `reghdfe …, absorb(id t) vce(cluster id \| robust)` | `fe(id) + fe(t)`, drop | `Clustered(:unit)`, `Robust()` | se 4e-15; `e(df_a)` exact |
| R `feols(f(y,h) ~ … \| id [+ t], fixef.rm = "singleton" \| "none")` with `ssc(K.adj = FALSE, G.adj = FALSE)`: cluster, hetero, twoway | drop / keep | `Clustered(...; ssc = :none)`, `Robust(; ssc = :none)`, `Clustered(:unit, :time; ssc = :none, psd = :raw)` | 7e-14 |
| R `feols(...)`, default `ssc`, `cluster = ~id` | drop / keep | `Clustered(:unit)` | 2e-14 |
| R `feols(… \| id)` on units 1–4 (`fx_balanced_id4`), `cluster = ~id + t`: not PSD; default `vcov_fix = TRUE` / `vcov_fix = FALSE` | drop | `Clustered(:unit, :time; ssc = :none, psd = :raw \| false)` | 1e-10 (NaN where fixest's variance is negative) |
| R MSS `LP_panel(…, method = "FE", te, two_way_cluster, robust, eigens)` (H = 0, pre-built leads; h = 0–3, also on `fx_balanced_id4`, where `eigens = T` repairs) | keep; `te = T` ⇒ `fe(id) + fe(t)`, `absorb = OneSweep()` | `Clustered(:unit [, :time]; ssc = :none \| :mss, psd = :none \| :raw)` | 7e-15 (NaN where R's variance is negative) |

(*) At h = 0, 1 on fx_balanced Stata refuses reghdfe's repaired matrix when posting it
("variance matrix is nonsymmetric or highly singular") and reports zero standard errors.
The test then only requires that our repair fired. At h = 2 the repair fires and reghdfe
posts a valid matrix, which we match.

`xtlp …, tfe` counts the absorbed degrees of freedom as N + T − 1, which assumes the unit
and time effects form one connected component. `panellp` counts the components, as
reghdfe does. On a sample with disconnected blocks the two differ, so the `xtlp tfe`
mapping holds for connected samples only.

### Phase 2: split-panel jackknife

Committed tests in `test/test_crosscheck_p2.jl`, at horizons 0–3 on both synthetic panels
unless noted, with the same checks as Phase 1: row counts, every spec mapped,
coefficients, standard errors, observations, critical degrees of freedom (not under
`ase`) and numbers of clusters. The `estimator` is `SplitPanelJackknife()`; the formula
is that of Phase 1, or `leads(y) ~ x + w` ("no lags": units with two or three rows enter,
with single-row halves).

| Reference command | `panellp` options | `vcov` spec | Agreement |
|---|---|---|---|
| `xtlp …, fe m(spj) vce(…)`: `cluster id [ase]`, `cluster t`, `cluster reg`, `cluster id t [ase]`, `robust`, `dkraay` (default; `lag(2), ase`; `lag(2), nodfadj`), `keepsingletons`; no lags, drop / keep | `fe(id)`, drop / keep | as for `m(fe)` | b 1e-13, se 9e-15 |
| `xtlp …, tfe m(spj) vce(…)`: `cluster id [ase]`, `cluster t`, `cluster id t, ase`, `robust`, `dkraay lag(2), ase`, `keepsingletons`; no lags, drop / keep | `fe(id) + fe(t)` | as for `m(fe)` | b 1e-13, se 2e-11 |
| `xtspj …, model(regress) method(parm)` (fx_balanced, h = 1, 3: every unit with the same even number of periods) | `fe(id)` | — (point estimates) | 6e-14 |
| R MSS `LP_panel(…, method = "SPJ", te, two_way_cluster, robust, eigens)`, also on `fx_balanced_id4` | keep, as for FE | as for FE | 1e-14 (NaN where R's variance is negative) |
| R MSS notebook, four applications (local data and outputs: `crosscheck/R/mss_apps_local.R`; `test/test_mss_apps.jl`) | plain response, `horizon = 0`, keep (see `docs/src/replication/mss2026.md`) | per application | about 1e-12 |
| `xtlp …, tfe m(fe \| spj) vce(cluster id, ase) keepsingletons` on the Romer–Romer data (local: `crosscheck/stata/mss_apps_local.do`), the exact two-way columns of the replication page | `fe(id) + fe(t)`, keep, `AlternatingProjections()` | `Clustered(:unit; ssc = :none)` | b 1e-13, se 2e-10 |

`xtlp`'s `vce(unadjusted)` has no counterpart (`panellp` has no homoskedastic variance).
The variance of the split-panel jackknife is validated by these deterministic
cross-checks. The Monte Carlo test checks the bias pattern, and coverage within its
Monte Carlo error.

### Phase 3: cumulative-difference OLS

Committed tests in `test/test_crosscheck_p3.jl`, at horizons 0–9 on the three `xtarsim`
panels (`crosscheck/stata/p3_xtarsim.do`), `fx_cdols` and the unbalanced `fx_unbalanced`
(shock `x`; gaps, singleton units, a period absent for every unit, so no `xtscc` and no
`xtlp` rows there), with row counts and every spec mapped, coefficients and standard
errors, observations (not for `xtlp` rows), critical degrees of freedom (not for
`locproj` rows) and cluster counts. The formula is `ldiff(y) ~ s + lags(s, 1)` with
`estimator = CDOLS()` unless noted.

| Reference command (Stata `xtset id t`) | `panellp` options | `vcov` spec | Agreement |
|---|---|---|---|
| `regress F#.y−L.y s L.s, vce(cluster id \| robust)` | `CDOLS()` | `Clustered(:unit)`, `Robust()` | b 8e-14, se 3e-14 |
| `regress … LD.y` / `… L2.s` | `+ lags(firstdiff(y), 1)` / `lags(s, 1:2)` | `Clustered(:unit)` | b 1e-12, se 5e-14 |
| `regress … i.t, vce(cluster id \| robust)` | `+ fe(t)`, keep (`regress` keeps singleton periods) | `Clustered(:unit)`, `Robust()` | b 9e-14, se 3e-14 |
| `xtscc F#.y−L.y s L.s, lag(3)` (balanced panels) | `CDOLS()` | `DriscollKraaySE(3; ssc = :xtscc, timelags = :position)` | b 7e-14, se 4e-14 |
| `locproj y s, m(reg) tr(cmlt) sl(1) [yl(1)] hor(9) vce(cluster id \| robust)` (`e(irf)`) | `CDOLS()` [`+ lags(firstdiff(y), 1)`] | `Clustered(:unit)`, `Robust()` | b 8e-14, se 3e-14 |
| `xtreg F#.y−L.y s L.s, fe vce(cluster id)` (the paper's CMFE1) | `WithinFE()`, `fe(id)`, keep | `Clustered(:unit)` | b 3e-14, se 2e-14 |
| `xtlp y s, fe m(spj) hor(0 9) vce(cluster id)` (the jackknife of the paper's applications; its Monte Carlo used `locproj m(xtspj)`) | `leads(y) ~ s + fe(id)`, `SplitPanelJackknife()` | `Clustered(:unit)` | b 1e-14, se 1e-14 |
| R `feols(ld ~ s + l(s, 1))`, `cluster = ~id` (default `ssc`, and `ssc(K.adj = FALSE, G.adj = FALSE)`) / `"hetero"` | `CDOLS()` | `Clustered(:unit [; ssc = :none])` / `Robust()` | b 4e-14, se 1e-14 |

### Phase 4: t-LAHR (Almuzara & Sancibrián)

Committed tests in `test/test_crosscheck_p4.jl`, at horizons 0–6 on `fx_tlahr` (balanced,
N = 200, T = 40) and `fx_tlahr_unbal` (10% of the rows removed).

**The authors' `panel_LP.R`** (`crosscheck/R/p4_as_panelLP.R`; the code has no licence and
is only run, from the local copy in `docs/references/code`). Specs: `s_twfe` (s_i X_t,
`FE = cbind(id, t)`), `ones_unit` (s ≡ 1, `FE = cbind(id)`), `sv_twfe` (time-varying
exposure), `s_twfe_w` (`W = sv`), `s2_twfe` (exposures `s` and `sv`); `p_max ∈ {0, 2}`;
`cumulative` and `small_sample` on and off; 672 rows.

| `panel_LP(...)` | `panellp` options | `vcov` spec | Agreement |
|---|---|---|---|
| `y, X, s, FE = cbind(id, t), p_max = p` | `leads(y) ~ s & X + fe(id) + fe(t)`, `LagAugmented(CappedLags(p))`, `singletons = :keep` | `TLAHR(refine = false)` | b 2e-12, se 1e-13 (balanced); 2e-9 (unbalanced: fixest's tolerance) |
| `… small_sample = TRUE` | as above | `TLAHR(weights = :raw)` | se as above, df 3e-13 / 2e-9, except `s2_twfe`'s second df (their e₁ contrast) |
| `… cumulative = TRUE` | `cumul(y) ~ …` | as above | as above |
| `… small_sample = TRUE` (balanced panel: `s_twfe`, `s_twfe_w`, `ones_unit`, where ω_t is constant) | as above | `TLAHR()` (the paper's weights) | as above |

**Stata** (`crosscheck/stata/p4_tlahr.do`).

| Reference command (Stata `xtset id t`) | `panellp` options | `vcov` spec | Agreement |
|---|---|---|---|
| `reghdfe F#.y sX L(1/p).y L(1/p).sX, absorb(id t) vce(cluster t)`, p = min(h, 2), both panels | `leads(y) ~ s & X + fe(id) + fe(t)`, `LagAugmented(CappedLags(2))` | `Clustered(:time)` | b 2e-14 / 4e-12, se idem; N, df, G exact |
| synthetic regression, s ≡ 1: `collapse` means of `F#.y`, X, `L.X`, `L2.X`, each demeaned over the sample periods; `regress …, noconstant vce(hc2, dfadjust)` (se, `e(adj_df)`) | `leads(y) ~ X + lags(X, 1:2) + fe(id)` | `TLAHR()` | b, se, df 2e-14 |
| the same with ŝ_i = s_i − mean(s) weights | `leads(y) ~ s & X + lags(s & X, 1:2) + fe(id) + fe(t)` | `TLAHR()` | 3e-14 |
| `harreg …, noconstant estimator(ewc)` on both synthetic regressions (se, `e(df_fb)`) | as above | `CosineHAR()` | 3e-14 |

### Phase 6a: panel LP-IV (two-stage least squares)

Committed tests in `test/test_crosscheck_p6.jl`, at horizons 0–4 on `fx_lpiv` (balanced,
N = 60, T = 50) and `fx_lpiv_unbal` (staggered entry, period 25 absent, 8% of the rows
removed). The model: `leads(y) ~ (s & ffr ~ s & mp) + w + lags(y, 1) + fe(id) + fe(t)`
(Stata `Fy w Ly (sffr = smp)`), and with several endogenous regressors
`(g1 & ffr + g2 & ffr + g3 & ffr ~ g1 & mp + g2 & mp + g3 & mp) + w + fe(id)`.

| Reference command (Stata `xtset id t`) | `panellp` options | `vcov` spec | Agreement |
|---|---|---|---|
| `ivreghdfe …, absorb(id t) robust [small]` (`small` changes nothing: ivreghdfe always applies its correction) | `fe(id) + fe(t)`, drop | `Robust()` | b 4e-15, se 1e-11 |
| `ivreghdfe …, absorb(id t) cluster(id \| t) [small]`: a fixed effect is nested in the cluster | as above | `Clustered(:unit \| :time; ssc = :none)` × (N−1)/(N − K − df′)·G/(G−1), df′ = df_a minus the levels of the nested effect (two fewer than reghdfe's count, `ssc = :stata`) | se 2e-10 |
| `ivreghdfe …, absorb(id) cluster(t)`: none nested | `fe(id)` | `Clustered(:time)` | se 3e-15 |
| `ivreghdfe …, absorb(id t) cluster(id t) [small]` | `fe(id) + fe(t)` | `Clustered(:unit, :time; psd = false)` | se 6e-12 |
| `ivreghdfe …, absorb(id t) dkraay(3) [small]` (`dkraay(b)`: bandwidth b, b − 1 lags) | `fe(id) + fe(t)` | `DriscollKraaySE(2)` | se 8e-12 |
| `ivreghdfe Fy w (g1ffr g2ffr g3ffr = g1mp g2mp g3mp), absorb(id) cluster(t)` | `fe(id)` | `Clustered(:time)` | b 4e-15, se 3e-15 |
| `xtivreg2 …, fe cluster(id)` (xtivreg2 clusters only on groups nested in the panels) | `fe(id)` | `Clustered(:unit; ssc = :none)` | se 3e-15 |
| `xtivreg2 …, fe robust small` | `fe(id)` | `Robust()` | se 4e-15 |
| R `feols(Fy ~ w + Ly \| id + t \| sffr ~ smp, fixef.rm = "none")`, `ssc(K.adj = FALSE, G.adj = FALSE)`: hetero, cluster id / t, `~id + t` with `vcov_fix = FALSE`, `DK(2) ~ t` | keep | `Robust(; ssc = :none)`, `Clustered(…; ssc = :none)`, `Clustered(:unit, :time; ssc = :none, psd = false)`, `DriscollKraaySE(2; ssc = :none, timelags = :position)` (fixest lags over the periods present) | 1e-11 balanced, 5e-9 unbalanced (fixest's demeaning tolerance) |

The Kleibergen–Paap F of each ivreghdfe fit (`e(widstat)`) is recorded in the `F` column of
`p6_lpiv.csv`, for the first-stage diagnostics of Phase 6c.

### Phase 6b: lag-augmented panel LP-IV

Same panels and test file. At horizon h, p = min(h, 2) lags of y (exogenous) and of
`s & ffr` (endogenous, instrumented by the same lags of `s & mp`):
`LagAugmented(CappedLags(2))` with `leads(y) ~ (s & ffr ~ s & mp) + w + fe(id) [+ fe(t)]`.

| Reference command (Stata `xtset id t`) | `panellp` options | `vcov` spec | Agreement |
|---|---|---|---|
| `ivreghdfe Fy w L(1/p).y (sffr L(1/p).sffr = smp L(1/p).smp), absorb(id t) cluster(t)` | `fe(id) + fe(t)`, drop | `Clustered(:time; ssc = :none)` × ivreghdfe's nested-cluster factor (Phase 6a) | b within 1e-12, se 1e-9 |
| the same, `robust` / `dkraay(3)` | as above | `Robust()` / `DriscollKraaySE(2)` | as above |
| the same, `absorb(id) cluster(t)` | `fe(id)` | `Clustered(:time)` | as above |
| R `feols(Fy ~ w + Ly1 + … \| id + t \| sffr + Lsffr1 + … ~ smp + Lsmp1 + …, fixef.rm = "none")`, `cluster = ~t` and `"hetero"` with `ssc(K.adj = FALSE, G.adj = FALSE)` | keep | `TLAHR(refine = false)` (A&S supplement A.4), `Robust(; ssc = :none)` | 1e-11 balanced, 5e-9 unbalanced |

`ICLags` on IV fits equals the criterion of the OLS reduced form (the instruments as the
shocks), tested in `test/test_lpiv.jl`.

### Phase 6c: first-stage diagnostics

`firststage(r; vcov)` against ivreghdfe's `e(widstat)` (Kleibergen–Paap rk Wald F) and
weakivtest's `r(F_eff)` (Montiel Olea–Pflueger effective F). These are the just-identified
rows of `p6_lpiv.csv` (one and three endogenous regressors) and the over-identified ones of
`p6c_firststage.csv` (`s & mp + sv & mp`). Their small-sample factor counts the first
stage's regressors: the exogenous ones plus the excluded instruments, which is
`firststage`'s own count. With a fixed effect nested in the cluster, ivreghdfe's df′ is
max(df_a − its levels, 1).

| Reference command | `firststage` vcov | Agreement |
|---|---|---|
| `ivreghdfe …, absorb(id t) robust`, `e(widstat)` | `Robust()` | 2e-11 |
| `ivreghdfe …, absorb(id t) cluster(id \| t)` | `Clustered(:unit \| :time; ssc = :none)` ÷ ivreghdfe's nested factor | 7e-11 |
| `ivreghdfe …, absorb(id t) cluster(id t)` / `dkraay(3)` | `Clustered(:unit, :time; psd = false)` / `DriscollKraaySE(2)` | 3e-11 |
| `ivreghdfe …, absorb(id) cluster(t)`, three endogenous regressors | `Clustered(:time)` | 2e-15 |
| `xtivreg2 …, fe cluster(id)` | `Clustered(:unit; ssc = :none)` ÷ (N−1)/(N−K)·G/(G−1) | 2e-14 |
| `xtivreg2 …, fe robust small` | `Robust()` | 1e-15 |
| `ivreghdfe …, absorb(id) cluster(id) \| robust`, then `weakivtest`: `r(F_eff)` | `Clustered(:unit; ssc = :none)` ÷ nested factor / `Robust()`, `effF` | 4e-14 |
| two endogenous regressors, three instruments (`sffr svffr = smp svmp wmp`): robust, `cluster(t)`, `absorb(id) cluster(id)` | as above | 2e-11 |
| `ivreghdfe Fy (g1ffr g2ffr g3ffr = g1mp g2mp g3mp), absorb(id) robust` (diagonal Π̂, tied singular values) | `Robust()` | 4e-15 |

### Phase 7: LP-DiD, negative horizons, lead paths, weights

**LP-DiD** against `lpdid` 1.1.0 (Busch & Girardi; GPL-3, run as a black box) on the
`fx_lpdid` panels: `crosscheck/stata/p7_lpdid.do` → `test/fixtures/stata/p7_lpdid.csv`
(1424 rows: 47 option sets × 2 panels, with the event-study, pooled, aggregate and
pre-trend tables), compared in `test/test_crosscheck_p7.jl`. The formula is
`ldiff(y) ~ firstdiff(d) + fe(t)` with `estimator = LPDiD(…, compat = :lpdid)`,
`horizon = -5:8`. `compat = :lpdid` selects the command's rules where the package's
defaults differ (`docs/src/audit/lpdid.md`): it matters on `fx_lpdid_unbal` only; on
`fx_lpdid` the test also checks that the defaults give the same numbers.

| `lpdid y, unit(id) time(t) treat(d) pre_window(5) post_window(8) …` | `panellp` | `vcov` spec | Agreement |
|---|---|---|---|
| (baseline); `cluster(st)`; `[pweight = w]` | `LPDiD()`; `clusters = (:st,)`; `weights = :w` | `Clustered(:unit)`; `Clustered(:st)` | b 5e-15, se 5e-11 |
| `ylags(2) controls(x)`, `dylags(2)`, `controls(x xb)`, `absorb(reg)` | `lags(y, 1:2) + x`, `lags(firstdiff(y), 1:2)`, `x + xb`, `fe(reg)` | `Clustered(:unit)` | as above |
| `pmd(max)`, `pmd(3)`; `nevertreated`; `nocomp` | `LPDiD(pmd = :max \| 3)`; `control = :nevertreated`; `sample = :common` | as above | as above |
| `pooled`, `post_pooled(2 6) pre_pooled(3 5)` | `pooled = true` with `horizon = 0:8 \| -5:-2`, `2:6 \| -5:-3` | as above | 1e-15 |
| `rw` (also with weights, `nocomp`, `pmd(max) nevertreated`, `pooled`) | `LPDiD(reweight = true, compat = :lpdid)` | as above | 5e-15 |
| a hand-written `reghdfe [pweight]` with the factor from the regression rows (`d_rw_hand`) | `LPDiD(reweight = true)` | as above | 4e-15 |
| `rw` with `controls(x)`, `ylags(1)`, `dylags(2)`, `absorb(reg)`, `[pweight = w]`, `cluster(st)`, `pooled` (regression adjustment through `listreg`) | `LPDiD(reweight = true)` with the same terms | `Clustered(…; ssc = :none)`, normal bands | 2e-15 |
| `nonabsorbing(5 \| 3)`, `(5, notyet)`, `(5, firsttreat)`, `(, firsttreat notyet)`, `(5, firsttreat notyet) pmd(2)`, with `nevertreated`, `untreated_before`, `pmd(2 \| max)`, `rw`, `nocomp`, `rw dylags(1)`; treatment `d2` | `LPDiD(nonabsorbing = L \| true, control = :notyet \| :nevertreated, firsttreat, untreated_before, …)` | `Clustered(:unit)` | 1e-14 |
| `nonabsorbing(3, oneoff)` and its `notyet`, `firsttreat`, `nevertreated`, `untreated_before`, `pmd(2)`, `rw` variants; treatment `d3` | `LPDiD(nonabsorbing = 3, oneoff = true, …)` | as above | 1e-14 |
| `aggregate_average` (rows VWATT / ATT, Pre_avg) | `averageeffect(r; cluster)` | — | 1e-12 |
| `pretrend_test` (F, p, df) | `pretrendtest(r; cluster)` | — | 3e-11 |

Row `pre{k}` of `e(results)` is h = −k and `tau{h}` is h. `pre1` is h = −1, which `lpdid`
reports as (0, missing) and `panellp` as a zero coefficient with a zero standard error. A
pre-treatment difference spanned by the lagged outcomes is zero by construction in both.
The wild bootstrap is not compared (not implemented).

**General additions**: `crosscheck/stata/p7_general.do` → `p7_general.csv` (112 rows).

| Reference command (Stata `xtset id t`) | `panellp` | `vcov` spec | Agreement |
|---|---|---|---|
| `reghdfe Δ_h y x L(1/2).x [F(1/h).x], absorb(id t) vce(cluster id)`, h = −3..4, with Δ_h y = `F{h}.y − L.y` (`L{−h}.y − L.y` at h < 0); `fx_balanced`, `fx_unbalanced` | `ldiff(y) ~ x + lags(x, 1:2) + leads(x) + fe(id) + fe(t)`, `horizon = -3:4` | `Clustered(:unit)` | b 3e-14, se 4e-11 |
| the same, `vce(unadjusted)` | as above | `Homoskedastic()` | as above |
| `reghdfe Δ_h y x F(1/2).x, …` | `leads(x, 1:2)` | `Clustered(:unit)` | as above |
| `ivreg2 Δ_h y (sffr = smp) L.sffr [F(1/h).sffr] w i.id i.t`, no option / `small`, h = −2..4; `fx_lpiv`, `fx_lpiv_unbal` | `ldiff(y) ~ (s & ffr ~ s & mp) + lags(s & ffr, 1) + leads(s & ffr) + w + fe(id) + fe(t)`, `singletons = :keep` | `Homoskedastic(ssc = :none)` / `Homoskedastic()`; `e(cdf)` is `firststage(r).cdF` | b, se and F 6e-15 |
| `reghdfe Fy x L.y [pweight = w], absorb(id t) vce(robust \| cluster id \| cluster st)`; `absorb(t)`; `regress … [pweight = w], vce(cluster id)`; `fx_lpdid_unbal` | `weights = :w` | `Robust()`, `Clustered(:unit \| :st)` | b 5e-15, se 2e-11 |
| `reghdfe … [aweight = w], vce(unadjusted)` | `weights = :w` | `Homoskedastic()` | as above |
| `ivreghdfe Fy (x = L.x) [pweight = w], absorb(id t) robust` | `(x ~ lags(x, 1))`, `weights = :w` | `Robust()` | as above |
| `reghdfe ȳ x L(1/2).x, absorb(id t) vce(cluster id)`, ȳ the mean of Δ_h y over h = 0..3 (and over −3..−2) | `pooled = true`, `horizon = 0:3` (`-3:-2`) | `Clustered(:unit)` | within the test tolerance (b 1e-10, se 1e-9) |

At h = −1 the long difference is identically zero: `reghdfe` reports a zero coefficient;
its classical variance of an all-zero response is not compared.

**Local** (third-party data): `crosscheck/stata/ilzetzki_local.do` builds the data of
Ilzetzki (2024) with the replication package's own `WWII_Planes.do` and re-runs
`ivreg2 Δ_h y (x = z) L(1/6).x i.Time i.PlantByModelID controls [F(1/h).x] if …` for his
Figures 9(c), 9(d), 6(a–b) and A.4, h = −6..18. `test/test_ilzetzki.jl` (tag `:local_ei`)
compares coefficients and classical standard errors (1e-8), N, the degrees of freedom of
the bands and the Cragg–Donald F. See `docs/src/replication/ilzetzki2024.md`.
