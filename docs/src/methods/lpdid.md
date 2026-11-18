# LP-DiD

LP-DiD is the local-projections approach to difference-in-differences of Dube, Girardi,
Jordà & Taylor (2025, *Journal of Applied Econometrics*; DGJT below). With a binary
treatment ``D_{it}`` that units enter at different dates, it regresses the long difference
of the outcome on the treatment switch and time effects,

```math
y_{i,t+h} - y_{i,t-1} = \beta_h\, \Delta D_{it} + \delta^h_t + e^h_{it},
\qquad h = -Q, \dots, H,
```

on the rows that are either **newly treated** (``\Delta D_{it} = 1``) or **clean controls**
(``D_{i,t+h} = 0``). Every comparison is then between units entering treatment at ``t`` and
units not yet treated at ``t + h``, which avoids the comparisons with already-treated units
that bias two-way fixed-effects estimators when effects are heterogeneous.

```julia
r = panellp(@formula(ldiff(y) ~ firstdiff(D) + fe(t)), data;
            id = :id, time = :t, horizon = -5:8, estimator = LPDiD())
summarize(r, Clustered(:unit))        # event study: h = −5, …, 8
```

The formula is the regression: the response is `ldiff(y)`, the shock is `firstdiff(D)` with
`D` the 0/1 treatment column, and the time effects are written out (`fe(t)`; an error
without them). Covariates, lags of the outcome and further fixed effects are ordinary
terms. Negative horizons are the pre-treatment coefficients; ``h = -1`` is zero by
construction (the long difference is identically zero), with a zero standard error.

The method is implemented from the paper. The Stata command `lpdid` (Busch & Girardi,
GPL-3) served only as a reference to compare outputs with: every option below reproduces
it to round-off on the cross-check panels (`crosscheck/README.md`, Phase 7), and the
[audit page](../audit/lpdid.md) lists where the two differ.

## What the coefficient estimates

With ``N^1_g`` units entering at ``g`` and ``N^0_{g,h}`` clean controls,
``\hat\beta_h = \sum_g \omega_{g,h} \hat\tau_{g,h}`` (DGJT eq. 9–10), where
``\hat\tau_{g,h}`` is the difference between the mean long differences of the two groups
and

```math
\omega_{g,h} \propto \frac{N^1_g N^0_{g,h}}{N^1_g + N^0_{g,h}} .
```

The weights are positive: a variance-weighted average effect on the treated. With
`reweight = true` every treated row has the same weight, and ``\hat\beta_h`` is the
equally weighted average effect on the treated, ``\sum_g N^1_g \hat\tau_{g,h} / \sum_g
N^1_g`` (Section 3.3; with not-yet-treated controls it equals the Callaway–Sant'Anna
estimator). Both identities are tests of the package.

## Clean controls

| `LPDiD(...)` | rows in the sample of horizon ``h`` | `lpdid` |
|---|---|---|
| default (absorbing treatment) | ``\Delta D_{it} = 1``, or ``D_{i,t+h} = 0`` (``D_{it} = 0`` at ``h < 0``) | default |
| `control = :nevertreated` | controls only from units never treated in the data | `nevertreated` |
| `nonabsorbing = L` | ``\Delta D_{it} \in \{0, 1\}``, no switch in ``t-L, \dots, t-1`` nor in ``t+1, \dots, t+h``; at ``h < 0`` the same at ``t, t-1, \dots, t-\lvert h\rvert+1``. `nonabsorbing = 0` is 1 (the rule needs ``\Delta D_{i,t-1}``) | `nonabsorbing(L)` |
| `nonabsorbing = L, oneoff = true` | the same with events (``D = 1``) in place of switches | `nonabsorbing(L, oneoff)` |
| `nonabsorbing = L, control = :notyet` | controls only from units untreated from their first row through ``t + h`` | `nonabsorbing(L, notyet)` |
| `nonabsorbing = L, firsttreat = true` | rows with at most one switch (event) through ``t + h`` | `nonabsorbing(L, firsttreat)` |
| `nonabsorbing = true, firsttreat = true, control = :notyet` (or `:nevertreated`) | the absorbing rule on the first treatment spell: no stabilisation window is needed | `nonabsorbing(, firsttreat notyet)` |
| `untreated_before = true` | a unit counts as untreated before its first row | `untreated_before` |

- A condition on a value that is not observed (a missing treatment, an absent period) is
  false: the row is left out.
- Switches are counted between consecutive *observed* treatment values. A unit observed
  untreated and later treated has switched, also when the switch fell in an absent period
  or behind a missing value: it is no longer a not-yet-treated control, and a later
  switch is not its first.
- A control is a row with ``\Delta D_{it} = 0``: a row where the treatment switches off is
  neither treated nor a control, under every rule.
- With the default rule a treatment that ever switches off, in a gap too, is an error
  pointing to `nonabsorbing`. (`lpdid` runs, with a regressor of −1 on those rows.)
- With a non-absorbing treatment (DGJT Section 4.2) the estimand is the effect of a
  treatment event under the assumption that effects stabilise within ``L`` periods; with
  `firsttreat` and not-yet- or never-treated controls it is the effect of entering
  treatment for the first time.
- The rules are evaluated on all the data, before `subset`: a row can use treatment values
  of rows that `subset` excludes, as with every lead and lag in `panellp`.

`sample = :common` keeps the same rows at every horizon (no composition change; `lpdid`'s
`nocomp`): the intersection of the samples of the requested horizons, pre-treatment ones
included.

## Reweighting and covariates

- `reweight = true` without covariates and with only time effects is a weighted
  regression: each period's rows are weighted by ``1/(1 - p_t)``, ``p_t`` the share of
  treated rows of the period in the regression (times the observation weights, if any).
  Periods without a treated row or without a control row leave the sample.
- With covariates or further fixed effects, direct inclusion (`reweight = false`) keeps the
  variance weights and assumes that effects do not vary with the covariates (Section
  4.1.2). `reweight = true` is then **regression adjustment** (Section 4.1.1): the outcome
  model is fitted on the clean controls and

  ```math
  \hat\beta_h = \frac{1}{N^1} \sum_{\text{treated}} \bigl(\Delta_h y_{it} - Z_{it}'\hat\gamma_h\bigr),
  ```

  with ``Z`` the covariates and the fixed effects. The fit reports the treatment
  coefficient only, and stores its influence function as the score, so the variance of
  the estimator is `Clustered(:unit; ssc = :none)`, with normal critical values (as
  `lpdid`). A treated row whose period (or other fixed-effect level, or value of a
  two-valued covariate) has no control row cannot be predicted and is left out. With
  several fixed-effect dimensions the prediction also needs the row's combination of
  levels to be connected to the controls', as any two-way prediction does.
- The time effects may be group × time effects, `fe(region) & fe(t)` (comparisons within a
  region and period); `reweight = true` is then regression adjustment as well.
- `compat = :lpdid` applies the command's own rules where they differ from the ones on
  this page; see [Reproducing `lpdid`](@ref) below.

## Base period, pooled estimates, averages

- `pmd = k` (pre-mean differencing, Section 3.4) replaces the base ``y_{i,t-1}`` by the
  mean of ``y`` over ``t-k, \dots, t-1`` (all required); `pmd = :max` by the mean over
  every earlier row of the unit. With `nonabsorbing = L` the base must be clean as well:
  no switch (event) in the ``k + L - 1`` periods before ``t``, so that the outcomes in the
  base window are themselves ``L`` periods away from a switch; `:max` then means
  ``k = L``. With `pmd`, ``h = -1`` is estimated like
  any other horizon (`lpdid` reports it as zero).
- `pooled = true` gives the pooled estimate (Section 3.5): one regression of the outcome
  averaged over the horizons, less the base, on the rows where the whole window is
  observed, with the clean-control rule of the far end of the window.

  ```julia
  post = panellp(f, data; id, time, horizon = 0:8, pooled = true, estimator = LPDiD())
  pre = panellp(f, data; id, time, horizon = -5:-2, pooled = true, estimator = LPDiD())
  ```
- [`averageeffect`](@ref) averages the event-study coefficients over the post-treatment
  horizons (and over ``h \le -2``), weighting each by its number of treated rows: the
  average effect over all treated observations of the window, with the standard error of
  its influence function (which includes the estimated weights). Unlike the pooled
  estimate it does not need complete windows.
- [`pretrendtest`](@ref) is the joint test that the pre-treatment coefficients are zero,
  from their covariance across horizons (the score sums by cluster are stored with global
  codes, so no regression is stacked).

## Inference

The paper and `lpdid` cluster by unit: `Clustered(:unit)`, with `reghdfe`'s small-sample
factor and t(G − 1) critical values, is `lpdid`'s output. Any other specification applies
to the weighted-regression and direct-inclusion fits. The wild bootstrap of `lpdid` is not
implemented.

## Reproducing `lpdid`

On a balanced panel without missing values whose units all start untreated, the rules
above are the command's. Otherwise they differ in three places, and `compat = :lpdid`
switches to the command's:

| | `LPDiD(...)` | `compat = :lpdid` (`lpdid` 1.1.0) |
|---|---|---|
| controls of a period in the weighted regression of `reweight` | the rows of the regression | also clean rows whose ``\Delta D`` is undefined, which the regression drops |
| switches, for `control = :notyet`, `firsttreat` and the check of an absorbing treatment | changes between consecutive observed values | the observed first differences: a switch in a gap is not seen |
| first treatment with not-yet- or never-treated controls | rows where the treatment switches off are left out | kept, with a regressor of −1, for a unit treated from its first row |

`reweight = :lpdid` (from 0.2.1) is short for `reweight = true, compat = :lpdid`. The
[audit page](../audit/lpdid.md) has the details and the size of the differences.

## Map from `lpdid`

| `lpdid` | `panellp` |
|---|---|
| `lpdid y, unit(id) time(t) treat(D) pre_window(K) post_window(H)` | `ldiff(y) ~ firstdiff(D) + fe(t)`, `horizon = -K:H`, `estimator = LPDiD()` |
| `ylags(p)`, `dylags(p)`, `controls(x)`, `absorb(v)` | `lags(y, 1:p)`, `lags(firstdiff(y), 1:p)`, `x`, `fe(v)` in the formula |
| `[pweight = w]`, `weights(w)` | `weights = :w` |
| `cluster(c)` | `Clustered(c)` (declare `clusters = (c,)` for a column that is not a fixed effect) |
| `rw` | `LPDiD(reweight = true)` |
| the command's numbers on unbalanced panels | `LPDiD(…, compat = :lpdid)` |
| `pmd(k)`, `pmd(max)` | `LPDiD(pmd = k)`, `LPDiD(pmd = :max)` |
| `nocomp` | `sample = :common` |
| `pooled`, `post_pooled(a b)`, `pre_pooled(a b)` | `pooled = true` with `horizon = a:b`, `horizon = -b:-a` |
| `aggregate_average`, `pretrend_test` | `averageeffect(r)`, `pretrendtest(r)` |
| `level(90)` | `summarize(r, spec; level = 0.90)` |
| `bootstrap()` | not implemented |

## References

- Dube, A., Girardi, D., Jordà, Ò. & Taylor, A. M. (2025). A local projections approach to
  difference-in-differences. *Journal of Applied Econometrics* 40(7), 741–758.
- Busch, A. & Girardi, D. `lpdid`: Stata module implementing local projections
  difference-in-differences. Statistical Software Components S459273.
