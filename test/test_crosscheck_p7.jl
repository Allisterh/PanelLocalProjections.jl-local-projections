# Cross-checks of Phase 7 against Stata.
#
# LP-DiD (crosscheck/stata/p7_lpdid.do → test/fixtures/stata/p7_lpdid.csv): the command
# lpdid 1.1.0 (Busch & Girardi), run as a black box on the fx_lpdid panels. Map
# (crosscheck/README.md), with the formula ldiff(y) ~ firstdiff(d) + fe(t):
# - pre_window(K) post_window(H): horizon = -K:H (lpdid's row pre{k} is h = −k, tau{h} is h;
#   pre1, zero by construction, is h = −1);
# - ylags(k): lags(y, 1:k); dylags(k): lags(firstdiff(y), 1:k); controls(): regressors;
#   absorb(v): fe(v); [pweight = w]: weights = :w; cluster(c): Clustered(c);
# - nocomp: sample = :common; pooled (post_pooled(a b), pre_pooled(a b)): pooled = true
#   with horizon = a:b, -b:-a;
# - pmd(k | max): LPDiD(pmd = k | :max); nevertreated: control = :nevertreated;
#   nonabsorbing(L, notyet firsttreat oneoff): nonabsorbing = L, control = :notyet,
#   firsttreat = true, oneoff = true; untreated_before: untreated_before = true;
# - rw: reweight = :lpdid (lpdid's count of the controls) without covariates, reweight =
#   true with them (regression adjustment: Clustered(c; ssc = :none), normal bands);
# - aggregate_average: averageeffect(r; cluster); pretrend_test: pretrendtest(r; cluster).
#
# The general additions (crosscheck/stata/p7_general.do → p7_general.csv): reghdfe, ivreg2
# and ivreghdfe on negative horizons, lead paths, the classical variance, the Cragg–Donald
# F, observation weights and pooled fits.
# JuliaFormatter (sciml) mangles `@testitem ... setup = [...]`; this file is formatted by hand.
#! format: off

@testitem "Stata: lpdid 1.1.0" tags = [:crosscheck, :stata, :lpdid] begin
    using StatsModels
    include(joinpath(@__DIR__, "support", "PLPDGPs.jl"))
    using .PLPDGPs: read_csv
    const FIX = joinpath(@__DIR__, "fixtures")
    num(s) = isempty(s) ? NaN : parse(Float64, s)
    lines = readlines(joinpath(FIX, "stata", "p7_lpdid.csv"))
    ref = [(f = split(l, ","; keepempty = true);
            (ds = f[1], case = f[2], table = f[3], row = f[4], b = num(f[5]), se = num(f[6]),
                lo = num(f[7]), hi = num(f[8]), obs = num(f[9]))) for l in lines[2:end]]
    @test length(ref) == 1424

    A = @formula(ldiff(y) ~ firstdiff(d) + fe(t))
    B = @formula(ldiff(y2) ~ firstdiff(d2) + fe(t))
    C = @formula(ldiff(y3) ~ firstdiff(d3) + fe(t))
    st = (; clusters = (:st,))
    # case => (formula, estimator, panellp keywords, cluster)
    cases = Dict(
        "a_base" => (A, LPDiD(), (;), :unit),
        "a_cl_st" => (A, LPDiD(), st, :st),
        "a_pw" => (A, LPDiD(), (; weights = :w), :unit),
        "a_ylags" => (@formula(ldiff(y) ~ firstdiff(d) + lags(y, 1:2) + x + fe(t)), LPDiD(), (;), :unit),
        "a_dylags" => (@formula(ldiff(y) ~ firstdiff(d) + lags(firstdiff(y), 1:2) + fe(t)), LPDiD(), (;), :unit),
        "a_controls" => (@formula(ldiff(y) ~ firstdiff(d) + x + xb + fe(t)), LPDiD(), (;), :unit),
        "a_absorb" => (@formula(ldiff(y) ~ firstdiff(d) + fe(t) + fe(reg)), LPDiD(), (;), :unit),
        "a_pmdmax" => (A, LPDiD(pmd = :max), (;), :unit),
        "a_pmd3" => (A, LPDiD(pmd = 3), (;), :unit),
        "a_never" => (A, LPDiD(control = :nevertreated), (;), :unit),
        "a_nocomp" => (A, LPDiD(), (; sample = :common), :unit),
        "a_pooled" => (A, LPDiD(), (;), :unit),
        "a_pooled2" => (@formula(ldiff(y) ~ firstdiff(d) + lags(y, 1) + fe(t)), LPDiD(), (;), :unit),
        "a_agg" => (A, LPDiD(), (;), :unit),
        "a_agg_st" => (@formula(ldiff(y) ~ firstdiff(d) + lags(firstdiff(y), 1) + fe(t)), LPDiD(), st, :st),
        "a_agg_pw" => (@formula(ldiff(y) ~ firstdiff(d) + fe(t) + fe(reg)), LPDiD(), (; weights = :w), :unit),
        "a_rw" => (A, LPDiD(reweight = :lpdid), (;), :unit),
        "a_rw_pw" => (A, LPDiD(reweight = :lpdid), (; weights = :w), :unit),
        "a_rw_nocomp" => (A, LPDiD(reweight = :lpdid), (; sample = :common), :unit),
        "a_rw_pmd" => (A, LPDiD(reweight = :lpdid, pmd = :max, control = :nevertreated), (;), :unit),
        "d_rw_hand" => (A, LPDiD(reweight = true), (;), :unit),
        "a_ra_x" => (@formula(ldiff(y) ~ firstdiff(d) + x + fe(t)), LPDiD(reweight = true), (;), :unit),
        "a_ra_ylags" => (@formula(ldiff(y) ~ firstdiff(d) + lags(y, 1) + fe(t)), LPDiD(reweight = true), (;), :unit),
        "a_ra_dylags" => (@formula(ldiff(y) ~ firstdiff(d) + lags(firstdiff(y), 1:2) + fe(t)), LPDiD(reweight = true), (;), :unit),
        "a_ra_absorb" => (@formula(ldiff(y) ~ firstdiff(d) + fe(t) + fe(reg)), LPDiD(reweight = true), (;), :unit),
        "a_ra_xb_pw" => (@formula(ldiff(y) ~ firstdiff(d) + x + xb + fe(t)), LPDiD(reweight = true), (; weights = :w), :unit),
        "a_ra_st" => (@formula(ldiff(y) ~ firstdiff(d) + x + fe(t) + fe(reg)), LPDiD(reweight = true), st, :st),
        "b_L5" => (B, LPDiD(nonabsorbing = 5), (;), :unit),
        "b_L3" => (B, LPDiD(nonabsorbing = 3), (;), :unit),
        "b_L5_notyet" => (B, LPDiD(nonabsorbing = 5, control = :notyet), (;), :unit),
        "b_L5_first" => (B, LPDiD(nonabsorbing = 5, firsttreat = true), (;), :unit),
        "b_first_ny" => (B, LPDiD(nonabsorbing = true, firsttreat = true, control = :notyet), (;), :unit),
        "b_L5_first_ny" => (B, LPDiD(nonabsorbing = 5, firsttreat = true, control = :notyet, pmd = 2), (;), :unit),
        "b_L5_never" => (B, LPDiD(nonabsorbing = 5, control = :nevertreated), (;), :unit),
        "b_first_never" => (B, LPDiD(nonabsorbing = true, firsttreat = true, control = :nevertreated), (;), :unit),
        "b_L5_ub" => (B, LPDiD(nonabsorbing = 5, untreated_before = true), (;), :unit),
        "b_L5_pmd2" => (B, LPDiD(nonabsorbing = 5, pmd = 2), (;), :unit),
        "b_L3_pmdmax" => (B, LPDiD(nonabsorbing = 3, pmd = :max, untreated_before = true), (;), :unit),
        "b_L5_rw" => (B, LPDiD(nonabsorbing = 5, reweight = :lpdid), (;), :unit),
        "b_L5_nocomp" => (B, LPDiD(nonabsorbing = 5), (; sample = :common), :unit),
        "b_L5_ra" => (@formula(ldiff(y2) ~ firstdiff(d2) + lags(firstdiff(y2), 1) + fe(t)), LPDiD(nonabsorbing = 5, reweight = true), (;), :unit),
        "c_L3" => (C, LPDiD(nonabsorbing = 3, oneoff = true), (;), :unit),
        "c_L3_notyet" => (C, LPDiD(nonabsorbing = 3, oneoff = true, control = :notyet), (;), :unit),
        "c_L3_first" => (C, LPDiD(nonabsorbing = 3, oneoff = true, firsttreat = true), (;), :unit),
        "c_first_ny" => (C, LPDiD(nonabsorbing = true, oneoff = true, firsttreat = true, control = :notyet), (;), :unit),
        "c_first_never" => (C, LPDiD(nonabsorbing = true, oneoff = true, firsttreat = true, control = :nevertreated), (;), :unit),
        "c_L3_ub" => (C, LPDiD(nonabsorbing = 3, oneoff = true, untreated_before = true), (;), :unit),
        "c_L3_pmd2" => (C, LPDiD(nonabsorbing = 3, oneoff = true, pmd = 2), (;), :unit),
        "c_L3_rw" => (C, LPDiD(nonabsorbing = 3, oneoff = true, reweight = :lpdid), (;), :unit))
    # pooled windows (post, pre) of the cases with pooled rows
    pooledwin = Dict("a_pooled" => (0:8, -5:-2), "a_pooled2" => (2:6, -5:-3),
        "a_rw" => (0:8, -5:-2), "a_ra_x" => (0:8, -5:-2), "b_L5_rw" => (0:8, -5:-2))
    @test Set(x.case for x in ref) == Set(keys(cases))

    isra(r) = any(f -> f.did !== nothing && f.did.ra, r.fits)
    spec(r, cl) = isra(r) ? Clustered(cl; ssc = :none) : Clustered(cl)
    for ds in ("fx_lpdid", "fx_lpdid_unbal"), (case, (f, est, kw, cl)) in cases
        data = read_csv(joinpath(FIX, "synthetic", "$ds.csv"))
        mine = [x for x in ref if x.ds == ds && x.case == case]
        # the hand-written regression has the post-treatment horizons only
        hs = case == "d_rw_hand" ? (0:8) : (-5:8)
        r = panellp(f, data; id = :id, time = :t, horizon = hs, estimator = est, kw...)
        s = summarize(r, spec(r, cl))
        for x in mine
            x.table == "event" || continue
            x.row == "pre1" && continue
            h = startswith(x.row, "pre") ? -parse(Int, x.row[4:end]) : parse(Int, x.row[4:end])
            k = h - first(hs) + 1
            @test s.nobs[k] == x.obs
            if isnan(x.lo)
                # zero by construction (the difference is spanned by the lagged outcomes):
                # lpdid reports (0, 0) under regression adjustment
                @test abs(s.coef[k]) < 1e-12 && abs(x.b) < 1e-12
                continue
            end
            @test s.coef[k] ≈ x.b rtol = 1e-10 atol = 1e-12
            @test s.se[k] ≈ x.se rtol = 1e-9 atol = 1e-12
            @test s.lower[k] ≈ x.lo rtol = 1e-9 atol = 1e-10
            @test s.upper[k] ≈ x.hi rtol = 1e-9 atol = 1e-10
        end
        # h = −1 is zero by construction
        case == "d_rw_hand" || est.pmd != 0 || @test coef(r, -1)[1] == 0
        if haskey(pooledwin, case)
            for (row, win) in zip(("Post", "Pre"), pooledwin[case])
                x = only(x for x in mine if x.table == "pooled" && x.row == row)
                rp = panellp(f, data; id = :id, time = :t, horizon = win, estimator = est,
                    pooled = true, kw...)
                sp = summarize(rp, spec(rp, cl))
                @test rp.pooled == win && length(rp.fits) == 1
                @test sp.nobs[1] == x.obs
                @test sp.coef[1] ≈ x.b rtol = 1e-10
                @test sp.se[1] ≈ x.se rtol = 1e-9
                @test sp.lower[1] ≈ x.lo rtol = 1e-9 atol = 1e-10
            end
        end
        agg = [x for x in mine if x.table == "aggregate"]
        if !isempty(agg)
            a = averageeffect(r; cluster = cl)
            for x in agg
                k = findfirst(==(x.row == "Pre_avg" ? :pre : :post), a.window)
                @test a.ntreated[k] == x.obs
                @test a.coef[k] ≈ x.b rtol = 1e-10
                @test a.se[k] ≈ x.se rtol = 1e-9
                @test a.lower[k] ≈ x.lo rtol = 1e-9 atol = 1e-10
                @test a.upper[k] ≈ x.hi rtol = 1e-9 atol = 1e-10
            end
            x = only(x for x in mine if x.table == "pretrend")
            pt = pretrendtest(r; cluster = cl)
            @test pt.F ≈ x.b rtol = 1e-8
            @test pt.p ≈ x.se rtol = 1e-8
            @test pt.df == x.lo && pt.df_r == x.hi
        end
    end
end

@testitem "Stata: negative horizons, lead paths, classical variance, weights, pooled" tags = [:crosscheck, :stata, :lpdid] begin
    using StatsModels
    include(joinpath(@__DIR__, "support", "PLPDGPs.jl"))
    using .PLPDGPs: read_csv
    const FIX = joinpath(@__DIR__, "fixtures")
    num(s) = (s = strip(s); s == "." || isempty(s) ? NaN : parse(Float64, s))
    lines = readlines(joinpath(FIX, "stata", "p7_general.csv"))
    hdr = Tuple(Symbol.(split(lines[1], ",")))
    ref = [NamedTuple{hdr}(Tuple(strip.(split(l, ",")))) for l in lines[2:end]]
    @test length(ref) == 112
    data(ds) = read_csv(joinpath(FIX, "synthetic", "$ds.csv"))
    function check(r, cv, x; term = r.shocks[1], dfcheck = true)
        h = parse(Int, x.h)
        k = h - first(r.horizons) + 1
        j = findfirst(==(term), r.coefnames)
        @test nobs(r, h) == num(x.N)
        @test coef(r, h)[j] ≈ num(x.b) rtol = 1e-10 atol = 1e-13
        @test sqrt(cv.V[k][j, j]) ≈ num(x.se) rtol = 1e-9 atol = 1e-13
        dfcheck && isfinite(cv.df[k][j]) && @test cv.df[k][j] == num(x.df)
    end
    sel(ds, spec) = [x for x in ref if x.dataset == ds && x.spec == spec]

    # A. negative horizons and the path of x up to the horizon
    for ds in ("fx_balanced", "fx_unbalanced")
        d = data(ds)
        r = panellp(@formula(ldiff(y) ~ x + lags(x, 1:2) + leads(x) + fe(id) + fe(t)), d;
            id = :id, time = :t, horizon = -3:4)
        @test r.coefnames == ["x", "L1.x", "L2.x", "F1.x", "F2.x", "F3.x", "F4.x"]
        # the path has h leads at horizon h and none before
        @test [count(f.active) for f in r.fits] == [3, 3, 3, 3, 4, 5, 6, 7]
        @test coef(r, -1)[1] == 0 && all(isnan, coef(r, 1)[5:7])
        cl, cs = vcov(Clustered(:unit), r), vcov(Homoskedastic(), r)
        for x in sel(ds, "ols_cl_id")
            check(r, cl, x; term = x.term)
        end
        for x in sel(ds, "ols_classical")
            x.h == "-1" && continue        # reghdfe's variance of an all-zero response
            check(r, cs, x)
        end
        rf = panellp(@formula(ldiff(y) ~ x + leads(x, 1:2) + fe(id) + fe(t)), d;
            id = :id, time = :t, horizon = -3:4)
        clf = vcov(Clustered(:unit), rf)
        foreach(x -> check(rf, clf, x), sel(ds, "fixedleads_cl_id"))
    end

    # B. IV with the leads of the endogenous regressor as exogenous controls; ivreg2 with
    # dummies keeps singletons, and its Cragg–Donald F is firststage's cdF
    for ds in ("fx_lpiv", "fx_lpiv_unbal")
        r = panellp(@formula(ldiff(y) ~ (s & ffr ~ s & mp) + lags(s & ffr, 1) + leads(s & ffr) +
                w + fe(id) + fe(t)), data(ds); id = :id, time = :t, horizon = -2:4,
            singletons = :keep)
        c0, c1 = vcov(Homoskedastic(ssc = :none), r), vcov(Homoskedastic(), r)
        fs = firststage(r; vcov = Robust())
        for x in sel(ds, "iv_classical")
            check(r, c0, x)
            k = parse(Int, x.h) + 3
            @test fs.cdF[k] ≈ num(x.F) rtol = 1e-9
        end
        foreach(x -> check(r, c1, x), sel(ds, "iv_classical_small"))
        # normal critical values without the correction (NaN for the leads beyond h)
        @test all(isinf, filter(!isnan, c0.df[1]))
    end

    # C. observation weights
    d = data("fx_lpdid_unbal")
    kw = (; id = :id, time = :t, horizon = 0:3, weights = :w)
    r2 = panellp(@formula(leads(y) ~ x + lags(y, 1) + fe(id) + fe(t)), d; clusters = (:st,), kw...)
    r1 = panellp(@formula(leads(y) ~ x + lags(y, 1) + fe(t)), d; kw...)
    r0 = panellp(@formula(leads(y) ~ x + lags(y, 1)), d; kw...)
    riv = panellp(@formula(leads(y) ~ (x ~ lags(x, 1)) + fe(id) + fe(t)), d; kw...)
    @test r2.weights === :w
    for (spec, r, v) in (("pw_robust", r2, Robust()), ("pw_cl_id", r2, Clustered(:unit)),
        ("pw_cl_st", r2, Clustered(:st)), ("aw_classical", r2, Homoskedastic()),
        ("pw_tfe_cl_id", r1, Clustered(:unit)), ("pw_nofe_cl_id", r0, Clustered(:unit)),
        ("pw_iv_robust", riv, Robust()))
        cv = vcov(v, r)
        rows = sel("fx_lpdid_unbal", spec)
        @test length(rows) == 4
        foreach(x -> check(r, cv, x; term = "x"), rows)
    end

    # D. pooled fits
    d = data("fx_unbalanced")
    for x in sel("fx_unbalanced", "pooled_cl_id")
        win = x.h == "3" ? (0:3) : (-3:-2)
        r = panellp(@formula(ldiff(y) ~ x + lags(x, 1:2) + fe(id) + fe(t)), d;
            id = :id, time = :t, horizon = win, pooled = true)
        @test r.pooled == win && r.horizons == parse(Int, x.h):parse(Int, x.h)
        check(r, vcov(Clustered(:unit), r), x)
    end
end
