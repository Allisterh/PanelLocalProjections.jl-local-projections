* Ilzetzki (2024, AER), "Learning by Necessity": the panel LP-IV figures with pre-trends
* (negative horizons) and with the path of future output as controls (Figure 9c-d,
* LPIV_DiD.do and LPIV3D_DiD.do), re-run with ivreg2 exactly as the replication package
* issues them (unit and month dummies, no variance option).
*
* Inputs are the local, gitignored replication package under docs/references/code/aer_ei
* (code BSD-3, data CC BY 4.0): its WWII_Planes.do builds the data, unmodified. The
* outputs are local too:
*   docs/data/ilzetzki/ilzetzki_panel.csv      the estimation panel (numeric columns)
*   test/fixtures/local/ilzetzki_lp.csv        coefficient, classical se, N, rank, t df,
*                                              Cragg-Donald F per horizon
* The loop and the bookkeeping below are ours; test/test_ilzetzki.jl compares.
do "$PLP_ROOT/crosscheck/stata/common.do"

local pkg "$PLP_ROOT/docs/references/code/aer_ei"
capture confirm file "`pkg'/data/AircraftData.xlsx"
if _rc {
    display as text "Ilzetzki replication package absent (local-only): skipped"
    exit
}
* the build needs ereplace (SSC); the package stores floats, as Stata's default
capture which ereplace
if _rc ssc install ereplace, replace
set type float
global data "`pkg'/data"
global code "`pkg'/code/Code"
global figures "$PLP_ROOT/crosscheck/out"
quietly do "$code/WWII_Planes.do"
set type double

* ---------------------------------------------------------------- estimation panel
generate byte cond = Model != "Total"
local yv log_y TFPFloor TFPFloorCU
local xv logOutputAMPR logLOOInstrument mainProducer50 mainProducer25 initialCapacityUHigh
keep PlantByModelID Month Time cond `yv' `xv'
* rows that can never enter a sample, as themselves or as a lead or lag
egen byte anyobs = rownonmiss(`yv' logOutputAMPR)
summarize Month if anyobs > 0, meanonly
keep if inrange(Month, r(min) - 18, r(max) + 18)
drop anyobs
recast double `yv' `xv'
sort PlantByModelID Month
xtset PlantByModelID Month
capture mkdir "$PLP_ROOT/docs/data"
capture mkdir "$PLP_ROOT/docs/data/ilzetzki"
format `yv' `xv' %25.17e
* the month as its integer code (months since 1960m1), not as a %tm label
format Month %12.0g
export delimited using "$PLP_ROOT/docs/data/ilzetzki/ilzetzki_panel.csv", replace datafmt

* ---------------------------------------------------------------- regressions
generate double interaction = initialCapacityUHigh * logOutputAMPR
generate double InstInteract = initialCapacityUHigh * logLOOInstrument
capture mkdir "$PLP_ROOT/test/fixtures/local"
tempname fh
file open `fh' using "$PLP_ROOT/test/fixtures/local/ilzetzki_lp.csv", write replace
file write `fh' "spec,h,term,b,se,N,K,tdf,cdF" _n

capture program drop runspec
program define runspec
    syntax , fh(string) spec(string) y(string) pre(integer) [did(integer 0) iv(integer 1) threed(integer 0)]
    local x logOutputAMPR
    local z logLOOInstrument
    local rx "`x'"
    local coef `x'
    if `threed' {
        local rx "(`x' interaction = `z' InstInteract)"
        local coef interaction
    }
    else if `iv' local rx "(`x' = `z')"
    forvalues h = `pre'/18 {
        capture drop lhs
        if `h' < 0 {
            local a = -`h'
            quietly generate double lhs = L`a'.`y' - L.`y'
        }
        else quietly generate double lhs = F`h'.`y' - L.`y'
        local fut ""
        if `did' & `h' > 0 local fut "F(1/`h').`x'"
        quietly ivreg2 lhs `rx' L(1/6).`x' i.Time i.PlantByModelID mainProducer50 mainProducer25 `fut' if cond
        local F = cond(missing(e(cdf)), ., e(cdf))
        file write `fh' "`spec',`h',`coef'," %25.17e (_b[`coef']) "," %25.17e (_se[`coef']) "," ///
            %12.0f (e(N)) "," %12.0f (e(rankxx)) "," %12.0f (e(N) - e(df_m) - 1) "," %25.17e (`F') _n
    }
end

* Figure 9(c), 9(d): the response holding later output fixed
runspec, fh(`fh') spec(did_TFPCU)   y(TFPFloorCU) pre(-1) did(1)
runspec, fh(`fh') spec(did3d_TFPCU) y(TFPFloorCU) pre(-1) did(1) threed(1)
* Figure 6(b) and A.4(b): the baseline LP-IV with six pre-treatment horizons
runspec, fh(`fh') spec(reg_TFPCU)   y(TFPFloorCU) pre(-6)
runspec, fh(`fh') spec(reg_logy)    y(log_y)      pre(-6)
* not in the paper: the OLS version of Figure 9(c) with pre-treatment horizons
runspec, fh(`fh') spec(did_TFPCU_ols) y(TFPFloorCU) pre(-6) did(1) iv(0)
file close `fh'
display as result "ilzetzki_local: done"
