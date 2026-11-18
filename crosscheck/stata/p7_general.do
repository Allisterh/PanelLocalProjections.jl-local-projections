* Phase 7 cross-check of the general additions to panellp: negative horizons, the path of
* a regressor between the shock and the response (leads 1..h as controls, as Ilzetzki
* 2024), the classical variance, the Cragg-Donald F, observation weights and pooled fits.
*
* A. fx_balanced, fx_unbalanced: reghdfe of y_{t+h} - y_{t-1} on x_t, two lags of x and, at
*    h > 0, the leads x_{t+1..t+h}; unit and time effects; h = -3..4; clustered by unit and
*    classical (vce(unadjusted)).
* B. fx_lpiv, fx_lpiv_unbal: ivreg2 with unit and time dummies (singletons kept) of the
*    long difference on (s ffr = s mp), its lag, its leads 1..h and w; h = -2..4; without
*    options (classical variance without correction, Cragg-Donald F) and with `small`.
* C. fx_lpdid_unbal: weighted regressions, reghdfe [pweight = w] (robust, clustered) and
*    [aweight = w] with vce(unadjusted); ivreghdfe [pweight = w] with x instrumented by its
*    lag; h = 0..3.
* D. fx_unbalanced: the pooled regression, the mean over h = 0..3 of y_{t+h} - y_{t-1}.
* Writes test/fixtures/stata/p7_general.csv; spec labels map to panellp/vcov calls in
* test/test_crosscheck_p7.jl.
do "$PLP_ROOT/crosscheck/stata/common.do"

capture mkdir "$PLP_ROOT/test/fixtures/stata"
tempname fh
file open `fh' using "$PLP_ROOT/test/fixtures/stata/p7_general.csv", write replace
file write `fh' "dataset,spec,h,term,b,se,df,N,F" _n

capture program drop _row
program define _row
    args fh ds spec h term
    local df = cond(missing(e(df_r)), ., e(df_r))
    local F = cond(missing(e(cdf)), ., e(cdf))
    file write `fh' "`ds',`spec',`h',`term'," %25.17e (_b[`term']) "," ///
        %25.17e (_se[`term']) "," %12.0f (`df') "," %12.0f (e(N)) "," %25.17e (`F') _n
end

* y_{t+h} - y_{t-1} for any integer h
capture program drop _ldiff
program define _ldiff
    args new y h
    if `h' >= 0 {
        quietly generate double `new' = F`h'.`y' - L.`y'
    }
    else {
        local k = -`h'
        quietly generate double `new' = L`k'.`y' - L.`y'
    }
end

* --- A. negative horizons and lead paths, OLS
foreach ds in fx_balanced fx_unbalanced {
    import delimited using "$PLP_ROOT/test/fixtures/synthetic/`ds'.csv", clear asdouble case(preserve)
    xtset id t
    forvalues h = -3/4 {
        tempvar dy
        _ldiff `dy' y `h'
        local leads ""
        if `h' > 0 local leads "F(1/`h').x"
        quietly reghdfe `dy' x L(1/2).x `leads', absorb(id t) vce(cluster id)
        _row `fh' `ds' ols_cl_id `h' x
        quietly reghdfe `dy' x L(1/2).x `leads', absorb(id t) vce(unadjusted)
        _row `fh' `ds' ols_classical `h' x
        if `h' >= 2 {
            quietly reghdfe `dy' x L(1/2).x `leads', absorb(id t) vce(cluster id)
            _row `fh' `ds' ols_cl_id `h' F2.x
        }
        * leads(x, 1:2) at every horizon
        quietly reghdfe `dy' x F(1/2).x, absorb(id t) vce(cluster id)
        _row `fh' `ds' fixedleads_cl_id `h' x
        drop `dy'
    }
}

* --- B. IV with lead paths, classical variance and the Cragg-Donald F
foreach ds in fx_lpiv fx_lpiv_unbal {
    import delimited using "$PLP_ROOT/test/fixtures/synthetic/`ds'.csv", clear asdouble case(preserve)
    xtset id t
    generate double sffr = s * ffr
    generate double smp = s * mp
    forvalues h = -2/4 {
        tempvar dy
        _ldiff `dy' y `h'
        local leads ""
        if `h' > 0 local leads "F(1/`h').sffr"
        quietly ivreg2 `dy' (sffr = smp) L.sffr `leads' w i.id i.t
        _row `fh' `ds' iv_classical `h' sffr
        quietly ivreg2 `dy' (sffr = smp) L.sffr `leads' w i.id i.t, small
        _row `fh' `ds' iv_classical_small `h' sffr
        drop `dy'
    }
}

* --- C. observation weights
import delimited using "$PLP_ROOT/test/fixtures/synthetic/fx_lpdid_unbal.csv", clear asdouble case(preserve)
xtset id t
forvalues h = 0/3 {
    tempvar fy
    quietly generate double `fy' = F`h'.y
    quietly reghdfe `fy' x L.y [pweight = w], absorb(id t) vce(robust)
    _row `fh' fx_lpdid_unbal pw_robust `h' x
    quietly reghdfe `fy' x L.y [pweight = w], absorb(id t) vce(cluster id)
    _row `fh' fx_lpdid_unbal pw_cl_id `h' x
    quietly reghdfe `fy' x L.y [pweight = w], absorb(id t) vce(cluster st)
    _row `fh' fx_lpdid_unbal pw_cl_st `h' x
    quietly reghdfe `fy' x L.y [aweight = w], absorb(id t) vce(unadjusted)
    _row `fh' fx_lpdid_unbal aw_classical `h' x
    quietly reghdfe `fy' x L.y [pweight = w], absorb(t) vce(cluster id)
    _row `fh' fx_lpdid_unbal pw_tfe_cl_id `h' x
    quietly regress `fy' x L.y [pweight = w], vce(cluster id)
    _row `fh' fx_lpdid_unbal pw_nofe_cl_id `h' x
    quietly ivreghdfe `fy' (x = L.x) [pweight = w], absorb(id t) robust
    _row `fh' fx_lpdid_unbal pw_iv_robust `h' x
    drop `fy'
}

* --- D. pooled regression
import delimited using "$PLP_ROOT/test/fixtures/synthetic/fx_unbalanced.csv", clear asdouble case(preserve)
xtset id t
generate double ybar = ((y - L.y) + (F1.y - L.y) + (F2.y - L.y) + (F3.y - L.y)) / 4
quietly reghdfe ybar x L(1/2).x, absorb(id t) vce(cluster id)
_row `fh' fx_unbalanced pooled_cl_id 3 x
generate double ypre = ((L3.y - L.y) + (L2.y - L.y)) / 2
quietly reghdfe ypre x L(1/2).x, absorb(id t) vce(cluster id)
_row `fh' fx_unbalanced pooled_cl_id -3 x

file close `fh'
display as result "p7_general: done"
