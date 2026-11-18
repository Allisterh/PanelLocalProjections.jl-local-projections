# Deviations from the reference implementations

Where this package's defaults differ from a reference implementation, the table gives
the reason and, where a published number depends on it, the switch that reproduces the
reference. The audit pages document each reference in detail.

## Mei, Sheng & Shi (`LP_panel_all.r`, `xtlp`)

| topic | reference | `panellp` default | switch |
|---|---|---|---|
| two-way demeaning | one unit-then-time sweep (R code): exact only on balanced samples | alternating projections to convergence | `absorb = OneSweep()` |
| split-panel jackknife halves | first ``\lceil n_i/2 \rceil`` rows (R code, `xtlp`); the paper writes ``\lfloor n_i/2 \rfloor`` | `split = :ceil`, as the code | `SplitPanelJackknife(split = :floor)` |
| small-sample factors (R code, `robust = T`) | ``N/(N-1) \cdot (n-1)/(n-K)``; two-way ``G/(G-1) \cdot (n-1)/(n-N-K)``, ``G = \min(N, T)`` | `reghdfe`'s factors (`ssc = :stata`) | `Clustered(:unit or :unit, :time; ssc = :mss)` |
| repair of a two-way variance that is not positive semi-definite | R code: clip the eigenvalues of V (`eigens = T`); `xtlp`: clip in the metric of the regressors' standard deviations | `psd = :stata`, as `xtlp` and `reghdfe` | `psd = :raw` |
| Driscoll–Kraay lags (`xtlp`, `xtscc`) | by position among the sample's periods | by calendar period | `timelags = :position` |
| `xtlp`'s `e(irf)` bands | ±1.96 standard errors whatever the variance | the specification's critical values | `summarize(…; dist = :normal)` |
| transformed outcomes (`xtlp`) | stored as floats unless `set type double` | double | — |
| a horizon whose variance is undefined (`xtlp`) | coefficients posted as missing | coefficients kept, variance NaN | — |

## Ugarte-Ruiz (2026)

| topic | paper | `panellp` | switch |
|---|---|---|---|
| specification rules | stated in the text | checked (`CDOLS()`) | `CDOLS(validate = false)` |
| variance | "robust", formula not stated (heteroskedasticity-robust, presumably) | any; unit clusters are valid, robust ones are not | `Robust()` |

## Almuzara & Sancibrián (`panel_LP.m`/`.R`)

| topic | reference code | `panellp` default | switch |
|---|---|---|---|
| lags | ``\min(h, p_{\max})``: no lag at ``h = 0``; no information criterion | the paper's ``\max\{1, \min(h, \hat p)\}`` with AIC (`ICLags()`) | `LagAugmented(CappedLags(p_max))` |
| the HC2/Imbens–Kolesár refinement | off by default | on (the paper's recommendation) | `TLAHR(refine = false)` |
| synthetic regression | raw exposure ``s_{it}``; periods unweighted | the exposure residualised per period (the paper's ``\hat s_i``, and on the exposures of the other regressors of the same aggregate series); periods weighted by ``\omega_t`` (fn 24: CR2 by period) | `TLAHR(weights = :raw)` |
| degrees-of-freedom contrast | ``e_1`` for every coefficient (a bug) | each coefficient's own | — |
| singletons | kept | dropped, as `reghdfe` and `xtlp` | `singletons = :keep` |
| no fixed effects | no intercept | an intercept | `0 +` |
| `Date` time index (R) | broken (day counts) | calendar periods with `frequency` | — |

## Dube, Girardi, Jordà & Taylor (`lpdid` 1.1.0)

| topic | `lpdid` | `panellp` (`LPDiD`) | switch |
|---|---|---|---|
| reweighting factor of the weighted regression (`rw`) | the treated share of the period among clean rows with the outcome, including rows whose ``\Delta D`` is undefined and that the regression drops | the treated share among the rows of the regression (the equally weighted average of the cohorts' comparisons) | `reweight = :lpdid` |
| a treatment that switches off, without `nonabsorbing()` | kept, with a regressor of −1 | an error | `nonabsorbing = L` |
| `nocomp` with a one-sided window | also requires the pre-treatment horizon 2 (or horizon 0) | the requested horizons (`sample = :common`) | `horizon = -2:H` |
| pooled regressions under `nocomp` | the flags of both pooled windows | the horizons of the window | — |
| ``h = -1`` under `pmd()` | reported as zero | estimated | — |
| standard errors of regression adjustment | no small-sample factor, normal critical values | the specification's | `Clustered(:unit; ssc = :none)` |
| wild bootstrap | `bootstrap()` through `boottest` | not implemented | — |

## Ilzetzki (2024, `ivreg2` with dummies)

| topic | replication package | `panellp` default | switch |
|---|---|---|---|
| variance | classical, ``RSS/N`` (no option to `ivreg2`) | none: the variance is chosen at `summarize` | `Homoskedastic(ssc = :none)` |
| critical values | Student t with ``N - K``, on a standard error without correction | normal for `ssc = :none` | the degrees of freedom of `Homoskedastic()` |
| singletons of the line and month dummies | kept | dropped | `singletons = :keep` |
| first-stage F | Cragg–Donald | Kleibergen–Paap with the chosen variance | `firststage(r).cdF` |

## Variance conventions shared with other software

| topic | reference | `panellp` | switch |
|---|---|---|---|
| Driscoll–Kraay bandwidth | `xtscc`, `xtlp`: ``\lfloor 4(T/100)^{2/9} \rfloor`` lags; Almuzara & Sancibrián: kernel bandwidth ``\lceil 4((T-h)/100)^{2/9} \rceil`` | `DriscollKraaySE(:nw94)` | `:as2026` (``M - 1`` lags), `:mss2026`, an integer |
| Driscoll–Kraay factor | `xtscc` and `xtlp …, nodfadj`: ``(n-1)/(n-K-1) \cdot T/(T-1)`` | `xtlp`'s default (`ssc = :stata`) | `ssc = :xtscc`; `ssc = :none` (`ase`) |
| `fixest`'s two-way variance | eigenvalues floored at 1e-16 (`vcov_fix`) | `psd = :stata` | `psd = :raw` (clips at 0) |
| equal-weighted cosine variance | `harreg`: ``B = \lfloor 0.41\,T^{2/3} \rfloor``, t(``B``) | the same (`CosineHAR()`) | an integer ``B`` |
| absorption stopping rule | `xtlp`: absolute step below 1e-9; the Almuzara–Sancibrián code: 1e-8 (MATLAB: absolute norm of the change; R: `fixest::demean`) | relative to the column after its first sweep, 1e-12 | `AlternatingProjections(tol = …)` |
