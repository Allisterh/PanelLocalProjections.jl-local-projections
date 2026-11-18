# Lead paths and negative horizons

Four general options of `panellp` come from the event-study use of local projections.
They are what Ilzetzki (2024) uses on a continuous, instrumented treatment
([replication](../replication/ilzetzki2024.md)) and what [LP-DiD](lpdid.md) builds on.

## Negative horizons

`horizon = -6:12` estimates the regression at horizons before the shock as well: the
response is dated before the regressors, so a coefficient away from zero is a pre-trend.

- `leads(y)` is ``y_{i,t+h}`` and `ldiff(y)` is ``y_{i,t+h} - y_{i,t-1}`` for every integer
  ``h``. At ``h = -1`` the long difference is identically zero: the coefficient and its
  standard error are 0 (the horizon is the origin of the response, as in the figures of
  Ilzetzki 2024 and the `pre1` row of `lpdid`).
- The right-hand side at ``h < 0`` is that of ``h = 0``.
- `cumul(y)` has no negative horizons, and `LagAugmented` sets its lags ``p(h)`` for
  ``h \ge 0`` only: both are errors. The other estimators accept them.
- An integer `horizon = H` still means `0:H`.

## The path of a regressor: `leads(x)`

A local projection of ``y_{t+h}`` on ``x_t`` includes the effect of the later values of
``x`` that ``x_t`` predicts. Adding ``x_{t+1}, \dots, x_{t+h}`` as controls gives the
response to ``x_t`` *holding its later values fixed*:

```math
y_{i,t+h} - y_{i,t-1} = \beta_h x_{it} + \sum_{j=1}^{h} \phi_{h,j}\, x_{i,t+j} + \dots
```

This is the regression of Ilzetzki (2024, Figure 9), and what conditioning on future
treatment does in LP-DiD (DGJT Section 4.3).

```julia
panellp(@formula(ldiff(y) ~ x + lags(x, 1:6) + leads(x) + fe(id) + fe(t)), data;
        id = :id, time = :t, horizon = -6:18)
```

- `leads(x)` without a range is the path: the columns `F1.x`, …, `Fh.x` at horizon ``h``
  and none at ``h \le 0``. A row needs the leads of its horizon only. A lead beyond the
  horizon is not in the regression (a NaN coefficient, not an omitted regressor).
- `leads(x, a:b)` adds the same leads at every horizon.
- In an IV fit the path of the endogenous regressor is exogenous: it is included in both
  stages and is not instrumented (as in Ilzetzki's `ivreg2` call). A path term cannot be
  inside the `(endo ~ inst)` block nor nested in another term.

## Observation weights

`weights = :w` is weighted least squares, as Stata's `[pweight = w]`: the fixed effects
are absorbed with weighted means, the scores are ``w_{it}\,\tilde x_{it}\,\hat e_{it}``, and
the degrees-of-freedom factors count rows. Rows with a missing or zero weight leave the
sample; a negative weight is an error. With `WithinFE`, `CDOLS` and `LPDiD` (and with an IV
block). `TLAHR(refine = true)` is not defined for weighted fits.

## Pooled fits

`pooled = true` runs one regression of the response averaged over the horizons, on the
rows where every horizon's response is defined: the average response over the window.
The horizons must be on one side of the shock (`0:8`, or `-5:-2`). The result holds a
single fit, stored at the far end of the window (`r.horizons`), and `r.pooled` is the
window. With the same regressors at every horizon, on a common sample without singletons
to drop, it equals the mean of the horizons' coefficients.

## The classical variance and the Cragg–Donald F

[`Homoskedastic`](@ref) is ``\hat\sigma^2 (\tilde X'\tilde X)^{-1}``, with
``(\hat X'\hat X)^{-1}`` and the structural residuals for LP-IV:

| | ``\hat\sigma^2`` | critical values | as |
|---|---|---|---|
| `Homoskedastic()` | ``RSS/(n - K - df_a)`` | ``t(n - K - df_a)`` | `regress`, `reghdfe, vce(unadjusted)`, `ivreg2 …, small` |
| `Homoskedastic(ssc = :none)` | ``RSS/n`` | normal | `ivreg2` without options |

It is the variance of the figures in Ilzetzki (2024); it is not robust to
heteroskedasticity or to correlation within units, so it is a compatibility option rather
than a recommendation. `firststage` reports, next to the Kleibergen–Paap and effective F,
the **Cragg–Donald** Wald F (`cdF`), the homoskedastic weak-identification statistic that
`ivreg2` prints without `robust` or `cluster`.
