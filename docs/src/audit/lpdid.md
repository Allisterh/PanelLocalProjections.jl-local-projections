# Audit: `lpdid`

`lpdid` (Busch & Girardi; SSC S459273), version 1.1.0 of 20 September 2026, is the Stata
implementation of LP-DiD by one of the paper's authors. It is GPL-3, so it was run as a
black box and nothing of it is ported: `LPDiD` is written from Dube, Girardi, Jordà &
Taylor (2025) and from the command's documented behaviour, and compared with its outputs
on the synthetic panels `fx_lpdid` and `fx_lpdid_unbal` (`crosscheck/stata/p7_lpdid.do`,
47 option sets on each panel; `test/test_crosscheck_p7.jl`).

## What it computes

Per horizon, `reghdfe` of the long difference on `D.treat`, the covariates and time
effects, on the clean-control sample, clustered by unit; with `rw` and covariates,
regression adjustment through `listreg`; optionally pooled regressions, an average of the
event-study coefficients, a joint pre-trend test, and a wild bootstrap through `boottest`.
The [LP-DiD page](../methods/lpdid.md) maps every option to `panellp`.

## Agreement

Event-study coefficients, standard errors, confidence bounds and numbers of observations
agree to round-off (coefficients 1e-14, standard errors 5e-11) for every option set:
covariates, lagged outcomes, extra fixed effects, weights, clusters, `pmd`,
`nevertreated`, `nocomp`, the non-absorbing and one-off rules with `notyet`, `firsttreat`
and `untreated_before`, `rw` by weighted regression and by regression adjustment; so do
the pooled estimates, `aggregate_average` and `pretrend_test`.

## Findings

1. **`rw` on unbalanced panels is not the equally weighted estimate.** The weighted
   regression's factor ``1/(1 - p_t)`` should use the treated share ``p_t`` of the rows of
   the regression. `lpdid` 1.1.0 computes it on the rows that are clean and have the
   outcome, which includes rows whose ``\Delta D`` is undefined (a unit's first row, a row
   after a gap) and which the regression then drops: the controls of those periods are
   overcounted. On `fx_lpdid_unbal` the coefficients differ from the equally weighted
   average of the cohorts' comparisons by up to 0.011 (0.2%); `LPDiD(reweight = true)`
   equals that average to round-off (a test), and equals a hand-written `reghdfe` with the
   weights from the regression rows (`d_rw_hand` in the do-file). On balanced panels, and
   with regression adjustment (which uses the regression rows), there is no difference.
   `LPDiD(reweight = :lpdid)` reproduces the command.
2. **An absorbing treatment is not checked.** Without `nonabsorbing()`, a treatment that
   switches off stays in the sample with a regressor of −1. `LPDiD()` refuses it and
   points to `nonabsorbing`.
3. **`nonabsorbing(0)` is `nonabsorbing(1)`.** The same here (the rule needs
   ``\Delta D_{t-1}``), stated in the table of clean controls.
4. **`nocomp` with a one-sided window** still requires the pre-treatment horizon 2 (and,
   with only a pre window, horizon 0), outcomes included. `sample = :common` uses the
   requested horizons: pass `horizon = -2:H` for the command's sample.
5. **The pooled regressions under `nocomp`** require the clean-control flags of both
   pooled windows. `pooled = true` with `sample = :common` uses the horizons of its own
   window.
6. **Mixed conventions in the output.** The event-study standard errors carry `reghdfe`'s
   small-sample factor and t(G − 1) critical values, the regression-adjustment ones no
   factor and normal critical values, and `aggregate_average` no factor with t(M − 1).
   `averageeffect` and `pretrendtest` follow the command; for the event study choose with
   `Clustered(…; ssc = …)`.
7. **The pre-treatment row `pre1` is always (0, missing)**, also under `pmd()`, where the
   difference between ``y_{t-1}`` and the pre-treatment mean is not zero by construction.
   `panellp` estimates ``h = -1`` like any other horizon then.
8. **The reported `t` and `p` are rounded** (2 and 4 decimals); compare coefficients,
   standard errors and bounds.
9. A pre-treatment difference that the lagged outcomes span (`ylags(k)`: horizons up to
   ``-k``; `dylags(k)`: up to ``-(k+1)``) is zero by construction: `reghdfe` returns
   round-off, regression adjustment (0, 0). Here the coefficient is round-off as well, and
   `pretrendtest` leaves those horizons out, as the command does.

Not compared: the wild bootstrap (`bootstrap()`), which is not implemented.
