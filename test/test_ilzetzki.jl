# Ilzetzki (2024, AER), "Learning by Necessity": the panel LP-IV with pre-treatment horizons
# and with the path of future output held fixed (Figure 9c–d of the paper), against ivreg2
# run as the replication package runs it (crosscheck/stata/ilzetzki_local.do). The data
# are third-party and local-only: the item is skipped when they are absent
# (test/runtests.jl).
#
# Map: `ivreg2 lhs (x = z) L(1/6).x i.Time i.PlantByModelID controls F(1/h).x if cond` is
#     ldiff(y) ~ (x ~ z) + lags(x, 1:6) + leads(x) + controls + fe(unit) + fe(month)
# with horizon = -1:18, subset = cond, singletons = :keep (dummies keep them),
# `Homoskedastic(ssc = :none)` (ivreg2 without options), t(N − K) bands from
# `Homoskedastic()`'s degrees of freedom, and firststage's `cdF` (ivreg2's `widstat`).
# JuliaFormatter (sciml) mangles `@testitem ... setup = [...]`; this file is formatted by hand.
#! format: off

@testitem "Ilzetzki (2024): LP-IV with pre-trends and the path of future output" tags = [:local_ei, :lpdid] begin
    using StatsModels
    include(joinpath(@__DIR__, "support", "PLPDGPs.jl"))
    using .PLPDGPs: read_csv
    root = dirname(@__DIR__)
    d = read_csv(joinpath(root, "docs", "data", "ilzetzki", "ilzetzki_panel.csv"))
    data = merge(d, (; sub = d.cond .== 1))
    num(s) = (s = strip(s); s == "." || isempty(s) ? NaN : parse(Float64, s))
    lines = readlines(joinpath(@__DIR__, "fixtures", "local", "ilzetzki_lp.csv"))
    hdr = Tuple(Symbol.(split(lines[1], ",")))
    ref = [NamedTuple{hdr}(Tuple(strip.(split(l, ",")))) for l in lines[2:end]]
    @test length(ref) == 115
    kw = (; id = :PlantByModelID, time = :Month, subset = :sub, singletons = :keep)
    specs = Dict(
        "did_TFPCU" => (@formula(ldiff(TFPFloorCU) ~ (logOutputAMPR ~ logLOOInstrument) +
                lags(logOutputAMPR, 1:6) + leads(logOutputAMPR) + mainProducer50 +
                mainProducer25 + fe(PlantByModelID) + fe(Month)), -1:18, "logOutputAMPR"),
        "did3d_TFPCU" => (@formula(ldiff(TFPFloorCU) ~
                (logOutputAMPR + initialCapacityUHigh & logOutputAMPR ~
                 logLOOInstrument + initialCapacityUHigh & logLOOInstrument) +
                lags(logOutputAMPR, 1:6) + leads(logOutputAMPR) + mainProducer50 +
                mainProducer25 + fe(PlantByModelID) + fe(Month)), -1:18,
            "initialCapacityUHigh & logOutputAMPR"),
        "reg_TFPCU" => (@formula(ldiff(TFPFloorCU) ~ (logOutputAMPR ~ logLOOInstrument) +
                lags(logOutputAMPR, 1:6) + mainProducer50 + mainProducer25 +
                fe(PlantByModelID) + fe(Month)), -6:18, "logOutputAMPR"),
        "reg_logy" => (@formula(ldiff(log_y) ~ (logOutputAMPR ~ logLOOInstrument) +
                lags(logOutputAMPR, 1:6) + mainProducer50 + mainProducer25 +
                fe(PlantByModelID) + fe(Month)), -6:18, "logOutputAMPR"),
        "did_TFPCU_ols" => (@formula(ldiff(TFPFloorCU) ~ logOutputAMPR +
                lags(logOutputAMPR, 1:6) + leads(logOutputAMPR) + mainProducer50 +
                mainProducer25 + fe(PlantByModelID) + fe(Month)), -6:18, "logOutputAMPR"))
    for (spec, (f, hs, term)) in specs
        r = panellp(f, data; horizon = hs, kw...)
        s0 = summarize(r, Homoskedastic(ssc = :none); term)
        s1 = summarize(r, Homoskedastic(); term)
        fs = isempty(r.endogenous) ? nothing : firststage(r; vcov = Robust())
        for x in filter(x -> x.spec == spec, ref)
            k = parse(Int, x.h) - first(hs) + 1
            @test s0.nobs[k] == num(x.N)
            @test s0.coef[k] ≈ num(x.b) rtol = 1e-8 atol = 1e-12
            @test s0.se[k] ≈ num(x.se) rtol = 1e-8 atol = 1e-12
            @test s1.df[k] == num(x.tdf)                 # N − rank of the dummy regression
            fs === nothing || @test fs.cdF[k] ≈ num(x.cdF) rtol = 1e-8
        end
        # h = −1: the long difference is identically zero
        @test coef(r, -1)[1] == 0 && s0.se[findfirst(==(-1), hs)] == 0
    end
end
