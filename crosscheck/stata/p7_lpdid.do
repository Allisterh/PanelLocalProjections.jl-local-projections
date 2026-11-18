* Phase 7 cross-check: LP-DiD (Dube, Girardi, Jordà & Taylor 2025) on the fx_lpdid panels,
* against the Stata command lpdid 1.1.0 (Busch & Girardi; SSC), run as a black box.
*
* A. Absorbing treatment `d`, outcome `y`: baseline, cluster, weights, covariates, extra
*    fixed effects, pre-mean differencing, never-treated controls, no composition change,
*    pooled estimates, reweighting (weighted regression and regression adjustment), the
*    aggregate average and the pre-trend test.
* B. Non-absorbing persistent treatment `d2`, outcome `y2`: stabilisation window, not-yet-
*    treated and never-treated controls, first treatment, untreated_before.
* C. One-off treatment `d3`, outcome `y3`.
* D. The reweighted regression with the weights computed on the regression rows (the
*    definition panellp uses; lpdid 1.1.0 also counts clean rows whose treatment switch is
*    undefined), by hand with reghdfe.
* Writes test/fixtures/stata/p7_lpdid.csv; the case labels map to panellp calls in
* test/test_crosscheck_p7.jl.
do "$PLP_ROOT/crosscheck/stata/common.do"

capture mkdir "$PLP_ROOT/test/fixtures/stata"
tempname fh
file open `fh' using "$PLP_ROOT/test/fixtures/stata/p7_lpdid.csv", write replace
file write `fh' "dataset,case,table,row,b,se,lo,hi,obs" _n

capture program drop _num
program define _num, rclass
    * the whole argument is one expression (it may contain spaces)
    if missing(`0') {
        return local s ""
    }
    else {
        return local s = strtrim(string(`0', "%25.17e"))
    }
end

* one row of a results matrix: coefficient, se, ci_low, ci_high, obs (columns 1, 2, 5, 6, 7)
capture program drop _matrow
program define _matrow
    args fh ds case table row M i
    local line "`ds',`case',`table',`row'"
    foreach c in 1 2 5 6 7 {
        _num `M'[`i', `c']
        local line "`line',`r(s)'"
    }
    file write `fh' "`line'" _n
end

* everything lpdid stored: the event study, the pooled estimates, the aggregate average
* and the pre-trend test (F, p, numerator and denominator df in the b, se, lo, hi columns)
capture program drop _out
program define _out
    args fh ds case
    tempname R
    foreach tb in results pooled_results aggregate {
        capture confirm matrix e(`tb')
        if _rc continue
        matrix `R' = e(`tb')
        local table = cond("`tb'" == "results", "event", cond("`tb'" == "pooled_results", "pooled", "aggregate"))
        local rn : rownames `R'
        local i 0
        foreach r of local rn {
            local ++i
            _matrow `fh' `ds' `case' `table' `r' `R' `i'
        }
    }
    if !missing(e(pretrend_df)) {
        local line "`ds',`case',pretrend,F"
        foreach s in pretrend_F pretrend_p pretrend_df pretrend_df_r {
            _num e(`s')
            local line "`line',`r(s)'"
        }
        file write `fh' "`line'," _n
    }
end

capture program drop _lp
program define _lp
    gettoken fh 0 : 0
    gettoken ds 0 : 0
    gettoken case 0 : 0
    capture noisily lpdid `0' nograph
    if _rc {
        display as error "p7_lpdid: case `case' on `ds' failed with r(" _rc ")"
        exit _rc
    }
    _out `fh' `ds' `case'
end

foreach ds in fx_lpdid fx_lpdid_unbal {
    import delimited using "$PLP_ROOT/test/fixtures/synthetic/`ds'.csv", clear asdouble case(preserve)
    xtset id t
    local U "unit(id) time(t)"
    local W "pre_window(5) post_window(8)"

    * --- A. absorbing treatment
    local A "y, `U' treat(d) `W'"
    _lp `fh' `ds' a_base       `A'
    _lp `fh' `ds' a_cl_st      `A' cluster(st)
    _lp `fh' `ds' a_pw         y [pweight = w], `U' treat(d) `W'
    _lp `fh' `ds' a_ylags      `A' ylags(2) controls(x)
    _lp `fh' `ds' a_dylags     `A' dylags(2)
    _lp `fh' `ds' a_controls   `A' controls(x xb)
    _lp `fh' `ds' a_absorb     `A' absorb(reg)
    _lp `fh' `ds' a_pmdmax     `A' pmd(max)
    _lp `fh' `ds' a_pmd3       `A' pmd(3)
    _lp `fh' `ds' a_never      `A' nevertreated
    _lp `fh' `ds' a_nocomp     `A' nocomp
    _lp `fh' `ds' a_pooled     `A' pooled
    _lp `fh' `ds' a_pooled2    `A' post_pooled(2 6) pre_pooled(3 5) ylags(1)
    _lp `fh' `ds' a_agg        `A' aggregate_average pretrend_test
    _lp `fh' `ds' a_agg_st     `A' aggregate_average pretrend_test cluster(st) dylags(1)
    _lp `fh' `ds' a_agg_pw     y [pweight = w], `U' treat(d) `W' aggregate_average pretrend_test absorb(reg)
    * reweighted, weighted regression
    _lp `fh' `ds' a_rw         `A' rw pooled aggregate_average pretrend_test
    _lp `fh' `ds' a_rw_pw      y [pweight = w], `U' treat(d) `W' rw
    _lp `fh' `ds' a_rw_nocomp  `A' rw nocomp
    _lp `fh' `ds' a_rw_pmd     `A' rw pmd(max) nevertreated
    * reweighted, regression adjustment
    _lp `fh' `ds' a_ra_x       `A' rw controls(x) pooled aggregate_average pretrend_test
    _lp `fh' `ds' a_ra_ylags   `A' rw ylags(1)
    _lp `fh' `ds' a_ra_dylags  `A' rw dylags(2)
    _lp `fh' `ds' a_ra_absorb  `A' rw absorb(reg)
    _lp `fh' `ds' a_ra_xb_pw   y [pweight = w], `U' treat(d) `W' rw controls(x xb)
    _lp `fh' `ds' a_ra_st      `A' rw controls(x) absorb(reg) cluster(st)

    * --- B. non-absorbing, persistent treatment
    local B "y2, `U' treat(d2) `W'"
    _lp `fh' `ds' b_L5          `B' nonabsorbing(5) aggregate_average pretrend_test
    _lp `fh' `ds' b_L3          `B' nonabsorbing(3)
    _lp `fh' `ds' b_L5_notyet   `B' nonabsorbing(5, notyet)
    _lp `fh' `ds' b_L5_first    `B' nonabsorbing(5, firsttreat)
    _lp `fh' `ds' b_first_ny    `B' nonabsorbing(, firsttreat notyet)
    _lp `fh' `ds' b_L5_first_ny `B' nonabsorbing(5, firsttreat notyet) pmd(2)
    _lp `fh' `ds' b_L5_never    `B' nonabsorbing(5) nevertreated
    _lp `fh' `ds' b_first_never `B' nonabsorbing(, firsttreat) nevertreated
    _lp `fh' `ds' b_L5_ub       `B' nonabsorbing(5) untreated_before
    _lp `fh' `ds' b_L5_pmd2     `B' nonabsorbing(5) pmd(2)
    _lp `fh' `ds' b_L3_pmdmax   `B' nonabsorbing(3) pmd(max) untreated_before
    _lp `fh' `ds' b_L5_rw       `B' nonabsorbing(5) rw pooled
    _lp `fh' `ds' b_L5_nocomp   `B' nonabsorbing(5) nocomp
    _lp `fh' `ds' b_L5_ra       `B' nonabsorbing(5) rw dylags(1)

    * --- C. one-off treatment
    local C "y3, `U' treat(d3) `W'"
    _lp `fh' `ds' c_L3          `C' nonabsorbing(3, oneoff) aggregate_average pretrend_test
    _lp `fh' `ds' c_L3_notyet   `C' nonabsorbing(3, oneoff notyet)
    _lp `fh' `ds' c_L3_first    `C' nonabsorbing(3, oneoff firsttreat)
    _lp `fh' `ds' c_first_ny    `C' nonabsorbing(, firsttreat notyet oneoff)
    _lp `fh' `ds' c_first_never `C' nonabsorbing(, firsttreat oneoff) nevertreated
    _lp `fh' `ds' c_L3_ub       `C' nonabsorbing(3, oneoff) untreated_before
    _lp `fh' `ds' c_L3_pmd2     `C' nonabsorbing(3, oneoff) pmd(2)
    _lp `fh' `ds' c_L3_rw       `C' nonabsorbing(3, oneoff) rw

    * --- D. reweighted regression by hand, weights from the regression rows:
    *     w_it = (N1_t + N0_t) / N0_t over the rows of period t in the regression
    forvalues h = 0/8 {
        tempvar dy dd smp n1 n0 wt
        quietly generate double `dy' = F`h'.y - L.y
        quietly generate double `dd' = D.d
        quietly generate byte `smp' = (`dd' == 1 | F`h'.d == 0) & !missing(`dy') & !missing(`dd')
        quietly bysort t: egen double `n1' = total(`smp' * (`dd' == 1))
        quietly bysort t: egen double `n0' = total(`smp' * (`dd' == 0))
        quietly generate double `wt' = (`n1' + `n0') / `n0' if `smp' & `n1' > 0 & `n0' > 0
        sort id t
        quietly reghdfe `dy' `dd' [pweight = `wt'] if `smp', absorb(t) vce(cluster id)
        local line "`ds',d_rw_hand,event,tau`h'"
        _num _b[`dd']
        local line "`line',`r(s)'"
        _num _se[`dd']
        local line "`line',`r(s)'"
        _num _b[`dd'] - invttail(e(df_r), 0.025) * _se[`dd']
        local line "`line',`r(s)'"
        _num _b[`dd'] + invttail(e(df_r), 0.025) * _se[`dd']
        local line "`line',`r(s)'"
        _num e(N)
        file write `fh' "`line',`r(s)'" _n
        drop `dy' `dd' `smp' `n1' `n0' `wt'
    }
}
file close `fh'
display as result "p7_lpdid: done"
