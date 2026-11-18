# Ilzetzki (2024): Learning by Necessity

Ilzetzki (2024, *American Economic Review*) estimates the response of aircraft plants'
productivity to demand in World War II by panel LP-IV: 143 production lines, monthly,
1940–1945, with the output of a line instrumented by the output of the other lines of its
broad aircraft type. Two features of his regressions are general options of `panellp`
since v0.2.1 ([Lead paths and negative horizons](../methods/leads.md)):

- **pre-trends**: the loop runs from a negative horizon (``h = -6`` in his Figure A.4,
  ``-1`` elsewhere), with the outcome ``y_{i,t+h} - y_{i,t-1}`` at every ``h``;
- **the path of future output as controls** (his "diff-in-diff" files `LPIV_DiD.do` and
  `LPIV3D_DiD.do`, Figure 9c–d): at ``h > 0`` the leads ``x_{t+1}, \dots, x_{t+h}`` of the
  treatment enter as exogenous regressors, so the coefficient of ``x_t`` is the response
  to demand at ``t`` holding later demand fixed.

There is no binary treatment and no clean-control sample here: it is not the
[LP-DiD](../methods/lpdid.md) estimator, only its outcome and the idea of conditioning on
later treatment.

The data and code are the author's replication package (code BSD-3, data CC BY 4.0); they
are not distributed with this package. `crosscheck/stata/ilzetzki_local.do` builds the data
with his `WWII_Planes.do`, unmodified, and re-runs the regressions with `ivreg2` exactly as
his files issue them; `test/test_ilzetzki.jl` (tag `:local_ei`, skipped when the data are
absent) compares.

## The regression

His call is `ivreg2 Δy (x = z) L(1/6).x i.Time i.PlantByModelID controls F(1/h).x if
Model != "Total"`, with no variance option. In `panellp`:

```julia
r = panellp(@formula(ldiff(TFPFloorCU) ~ (logOutputAMPR ~ logLOOInstrument) +
                     lags(logOutputAMPR, 1:6) + leads(logOutputAMPR) +
                     mainProducer50 + mainProducer25 +
                     fe(PlantByModelID) + fe(Month)),
            planes; id = :PlantByModelID, time = :Month, horizon = -1:18,
            subset = :lines, singletons = :keep)
s = summarize(r, Homoskedastic(ssc = :none))     # ivreg2 without options
firststage(r).cdF                                # its first-stage F (Cragg–Donald)
```

| His code | `panellp` |
|---|---|
| `f{h}.y - l.y`, `l{h}.y - l.y` (`$CumGrowth`) | `ldiff(y)`, `horizon = -1:18` |
| `(x = z)` | `(x ~ z)` |
| `l(1/6).x` | `lags(x, 1:6)` |
| `f(1/h).x` at ``h > 0`` (`LPIV_DiD.do`) | `leads(x)` |
| `i.Time i.PlantBy` (dummies; singletons stay) | `fe(Month) + fe(PlantByModelID)`, `singletons = :keep` |
| `if Model != "Total"` | `subset` |
| no variance option: ``\hat\sigma^2 = RSS/N``, no correction | `Homoskedastic(ssc = :none)` |
| bands with `invt(e(N) - e(df_m) - 1, ·)` | Student t with the degrees of freedom of `Homoskedastic()`, ``N - K`` counting the dummies |
| `e(widstat)` | `firststage(r).cdF` |
| 3D files: `(x D·x = z D·z)`, plotting the interaction | `(x + D & x ~ z + D & z)`, `term = "D & x"` |

## Results

Against `ivreg2` on his data, for Figure 9(c) and 9(d), Figures 6(a–b)/A.4 (``h = -6,
\dots, 18``) and an OLS variant: coefficients and standard errors within 1e-8 relative
(4e-12 absolute), the number of observations, the degrees of freedom of the bands and the
first-stage F exact (the F to 5e-15). Figure 9(c), the response of TFP (adjusted for
capacity utilisation) to output, holding later output fixed:

| ``h`` | ``\hat\beta_h`` | se (classical) | 95% band, t(N − K) | se (clustered by line) | Cragg–Donald F | KP F (by line) | N |
|---|---|---|---|---|---|---|---|
| 0 | 0.1026 | 0.0441 | [0.016, 0.189] | 0.0612 | 36.6 | 18.3 | 1860 |
| 5 | 0.1409 | 0.1011 | [−0.057, 0.339] | 0.2006 | 22.1 | 18.4 | 1331 |
| 9 | 0.4604 | 0.1831 | [0.101, 0.820] | 0.3341 | 12.3 | 9.0 | 997 |
| 12 | 0.2125 | 0.1481 | [−0.078, 0.503] | 0.2932 | 9.4 | 7.7 | 779 |
| 18 | 0.2131 | 0.2617 | [−0.301, 0.728] | 0.5254 | 3.6 | 4.3 | 453 |

The first four columns and the Cragg–Donald F are the paper's numbers (its Figure 9c
prints the 5-month F, 22). The other two are ours.

## Remarks

- The published bands use the classical variance, which assumes errors that are
  homoskedastic and uncorrelated within a production line. Clustered by line
  (`Clustered(:unit)`), the standard errors of Figure 9(c) are 1.4 to 2.2 times as large:
  the classical 95% bands exclude zero at ``h = 0`` and ``7``–``10``, the clustered ones at
  no horizon.
- The bands mix conventions: a standard error without a small-sample factor and a
  Student-t quantile with ``N - K`` degrees of freedom.
- At ``h = -1`` the outcome is identically zero, so the coefficient and its standard error
  are zero: the plotted origin. `ivreg2` fails there with `robust` or `cluster()`
  (r(504)); `panellp` returns the zero with any variance.
- The package stores variables as floats (Stata's default), so its own outputs differ from
  a double-precision run by about 5e-8.
- Not reproduced here (read, not run): the files with horizon-specific controls built from
  long differences of other variables (growth controls, spillovers), the cumulative-output
  treatment (`LP_MultiIV_Reg.do`), and the tables.

## Reference

Ilzetzki, E. (2024). Learning by necessity: Government demand, capacity constraints, and
productivity growth. *American Economic Review*. Replication package: AEA Data and Code
Repository (AEAREP-5450).
