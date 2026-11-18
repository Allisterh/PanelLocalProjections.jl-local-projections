# LP-DiD (Dube, Girardi, Jordà & Taylor 2025) and the general additions of Phase 7:
# negative horizons, lead paths, observation weights, pooled fits, the classical variance.
# Oracles written with plain LinearAlgebra; the Stata cross-checks are in
# test_crosscheck_p7.jl.
# JuliaFormatter (sciml) mangles `@testitem ... setup = [...]`; this file is formatted by hand.
#! format: off

@testitem "LP-DiD: specification checks" tags = [:lpdid, :api] begin
    using StatsModels
    include(joinpath(@__DIR__, "support", "PLPDGPs.jl"))
    using .PLPDGPs: read_csv
    d = read_csv(joinpath(@__DIR__, "fixtures", "synthetic", "fx_lpdid.csv"))
    kw = (; id = :id, time = :t, horizon = -3:4)
    f = @formula(ldiff(y) ~ firstdiff(d) + fe(t))
    @test panellp(f, d; estimator = LPDiD(), kw...) isa PanelLocalProjection{LPDiD}
    err(f, est; k...) = @test_throws ArgumentError panellp(f, d; estimator = est, kw..., k...)
    # the regression is the long difference on the treatment switch, with time effects
    err(@formula(leads(y) ~ firstdiff(d) + fe(t)), LPDiD())
    err(@formula(ldiff(y) ~ d + fe(t)), LPDiD())
    err(@formula(ldiff(y) ~ firstdiff(d) + fe(id)), LPDiD())
    err(@formula(ldiff(y) ~ firstdiff(d)), LPDiD())
    err(@formula(ldiff(y) ~ firstdiff(x) + fe(t)), LPDiD())            # not binary
    err(@formula(ldiff(y) ~ (firstdiff(d) ~ x) + fe(t)), LPDiD())      # no IV version
    # d2 switches off: not absorbing
    e = try
        panellp(@formula(ldiff(y2) ~ firstdiff(d2) + fe(t)), d; estimator = LPDiD(), kw...)
    catch e
        e
    end
    @test e isa ArgumentError && occursin("not absorbing", e.msg) && occursin("nonabsorbing = L", e.msg)
    # rules of a non-absorbing treatment need the stabilisation window
    err(f, LPDiD(firsttreat = true))
    err(f, LPDiD(oneoff = true))
    err(f, LPDiD(untreated_before = true))
    err(f, LPDiD(control = :notyet))
    err(@formula(ldiff(y2) ~ firstdiff(d2) + fe(t)), LPDiD(nonabsorbing = true))
    err(f, LPDiD(control = :other))
    @test_throws ArgumentError LPDiD(pmd = 0)
    @test_throws ArgumentError LPDiD(reweight = "yes")
    @test_throws ArgumentError LPDiD(nonabsorbing = 2.5)
    # show names the estimator and its options
    r = panellp(@formula(ldiff(y2) ~ firstdiff(d2) + fe(t)), d;
        estimator = LPDiD(nonabsorbing = 3, control = :notyet, pmd = 2), weights = :w, kw...)
    s = sprint(show, MIME"text/plain"(), r)
    @test occursin("LPDiD(control = :notyet, nonabsorbing = 3, pmd = 2)", s)
    @test occursin("weights:  w", s) && occursin("horizons -3:4", s)
    # the average effect and the pre-trend test need what they average or test
    @test_throws ArgumentError averageeffect(panellp(
        @formula(ldiff(y) ~ x + fe(id) + fe(t)), d; kw...))
    rp = panellp(f, d; estimator = LPDiD(), id = :id, time = :t, horizon = 0:4, pooled = true)
    @test_throws ArgumentError averageeffect(rp)
    @test occursin("pooled over horizons 0:4", sprint(show, MIME"text/plain"(), rp))
    @test_throws ArgumentError pretrendtest(r; horizons = -9:-2)
    @test_throws ArgumentError pretrendtest(r; cluster = :st)
end

@testitem "LP-DiD: clean-control comparisons by hand" tags = [:lpdid] begin
    using StatsModels, LinearAlgebra, Statistics
    include(joinpath(@__DIR__, "support", "PLPDGPs.jl"))
    using .PLPDGPs: read_csv
    d = read_csv(joinpath(@__DIR__, "fixtures", "synthetic", "fx_lpdid.csv"))
    N, T = 60, 30
    Y = permutedims(reshape(d.y, T, N))        # unit × period (balanced, sorted)
    Dm = permutedims(reshape(d.d, T, N))
    f = @formula(ldiff(y) ~ firstdiff(d) + fe(t))
    kw = (; id = :id, time = :t, horizon = 0:6)
    rw = panellp(f, d; estimator = LPDiD(reweight = true), kw...)
    vw = panellp(f, d; estimator = LPDiD(), kw...)
    for h in 0:6
        # cohort by cohort: the mean long difference of the units entering at t less that
        # of the units still untreated at t + h (DGJT eq. 9), for the periods with both
        τ, n1, n0 = Float64[], Int[], Int[]
        for t in 2:(T - h)
            tr = [i for i in 1:N if Dm[i, t] == 1 && Dm[i, t - 1] == 0]
            ct = [i for i in 1:N if Dm[i, t + h] == 0]
            (isempty(tr) || isempty(ct)) && continue
            dy(i) = Y[i, t + h] - Y[i, t - 1]
            push!(τ, mean(dy, tr) - mean(dy, ct))
            push!(n1, length(tr))
            push!(n0, length(ct))
        end
        # the reweighted estimate weights every treated row equally (Section 3.3)
        @test coef(rw, h)[1] ≈ sum(n1 .* τ) / sum(n1) rtol = 1e-12
        @test rw.fits[h + 1].did.ntreated == sum(n1)
        # the default has the variance weights N₁N₀/(N₁ + N₀) (eq. 10)
        ω = n1 .* n0 ./ (n1 .+ n0)
        @test coef(vw, h)[1] ≈ sum(ω .* τ) / sum(ω) rtol = 1e-12
    end
    # with every row on the regression's rows, lpdid's count of the controls is the same
    rl = panellp(f, d; estimator = LPDiD(reweight = :lpdid), kw...)
    @test coefpath(rl) ≈ coefpath(rw) rtol = 1e-13
    # never-treated controls only
    rn = panellp(f, d; estimator = LPDiD(control = :nevertreated, reweight = true), kw...)
    never = [i for i in 1:N if all(==(0), Dm[i, :])]
    h = 3
    # the equally weighted average of the cohorts' comparisons, for a base period `base`
    # and a control group `controls(t)`
    function att(base, controls)
        num = den = 0.0
        for t in 2:(T - h)
            tr = [i for i in 1:N if Dm[i, t] == 1 && Dm[i, t - 1] == 0]
            isempty(tr) && continue
            dy(i) = Y[i, t + h] - base(i, t)
            num += length(tr) * (mean(dy, tr) - mean(dy, controls(t)))
            den += length(tr)
        end
        return num / den
    end
    @test coef(rn, h)[1] ≈ att((i, t) -> Y[i, t - 1], t -> never) rtol = 1e-12
    # pre-mean differencing with all earlier periods as the base
    rp = panellp(f, d; estimator = LPDiD(reweight = true, pmd = :max), kw...)
    @test coef(rp, h)[1] ≈ att((i, t) -> mean(Y[i, 1:(t - 1)]),
        t -> [i for i in 1:N if Dm[i, t + h] == 0]) rtol = 1e-12

    # group × time effects: comparisons within a region and period; reweighting is then
    # regression adjustment, the equally weighted average of the cells' comparisons
    reg = d.reg[1:T:end]
    fr = @formula(ldiff(y) ~ firstdiff(d) + fe(reg) & fe(t))
    rr = panellp(fr, d; estimator = LPDiD(reweight = true), kw...)
    @test all(f -> f.did.ra, rr.fits)
    num = den = 0.0
    for t in 2:(T - h), g in 1:4
        tr = [i for i in 1:N if reg[i] == g && Dm[i, t] == 1 && Dm[i, t - 1] == 0]
        ct = [i for i in 1:N if reg[i] == g && Dm[i, t + h] == 0]
        (isempty(tr) || isempty(ct)) && continue
        dy(i) = Y[i, t + h] - Y[i, t - 1]
        global num += length(tr) * (mean(dy, tr) - mean(dy, ct))
        global den += length(tr)
    end
    @test coef(rr, h)[1] ≈ num / den rtol = 1e-10
    @test panellp(fr, d; estimator = LPDiD(), kw...) isa PanelLocalProjection

    # On an unbalanced panel the weights come from the rows of the regression: the estimate
    # is still the equally weighted average of the cohorts' comparisons. lpdid 1.1.0 also
    # counts clean rows whose ΔD is undefined (a unit's first row) among the controls of a
    # period, and its estimate is not that average (`reweight = :lpdid`).
    du = read_csv(joinpath(@__DIR__, "fixtures", "synthetic", "fx_lpdid_unbal.csv"))
    row = Dict((i, t) => k for (k, (i, t)) in enumerate(zip(du.id, du.t)))
    val(v, i, t) = (k = get(row, (i, t), 0); k == 0 || ismissing(v[k]) ? NaN : Float64(v[k]))
    function attu(h)
        num = den = 0.0
        for t in 1:T
            tr, ct = Float64[], Float64[]
            for i in 1:N
                haskey(row, (i, t)) || continue
                dy = val(du.y, i, t + h) - val(du.y, i, t - 1)
                dd = val(du.d, i, t) - val(du.d, i, t - 1)
                (isfinite(dy) && isfinite(dd)) || continue
                dd == 1 ? push!(tr, dy) : (val(du.d, i, t + h) == 0 && push!(ct, dy))
            end
            (isempty(tr) || isempty(ct)) && continue
            num += length(tr) * (mean(tr) - mean(ct))
            den += length(tr)
        end
        return num / den
    end
    ru = panellp(f, du; estimator = LPDiD(reweight = true), kw...)
    rl = panellp(f, du; estimator = LPDiD(reweight = :lpdid), kw...)
    @test coefpath(ru) ≈ [attu(h) for h in 0:6] rtol = 1e-12
    @test maximum(abs, coefpath(rl) .- coefpath(ru)) > 1e-3
    @test nobspath(rl) == nobspath(ru)
end

@testitem "LP-DiD: regression adjustment against explicit dummies" tags = [:lpdid] begin
    using StatsModels, LinearAlgebra
    include(joinpath(@__DIR__, "support", "PLPDGPs.jl"))
    using .PLPDGPs: read_csv
    for ds in ("fx_lpdid", "fx_lpdid_unbal"), wts in (nothing, :w)
        d = read_csv(joinpath(@__DIR__, "fixtures", "synthetic", "$ds.csv"))
        r = panellp(@formula(ldiff(y) ~ firstdiff(d) + x + lags(firstdiff(y), 1) + fe(t) + fe(reg)),
            d; id = :id, time = :t, horizon = 0:3, estimator = LPDiD(reweight = true),
            weights = wts)
        @test r.coefnames == ["D.d", "x", "L1.D.y"]
        cv = vcov(Clustered(:unit; ssc = :none), r)
        # hand-built columns on the (unit, period) grid
        row = Dict((i, t) => k for (k, (i, t)) in enumerate(zip(d.id, d.t)))
        val(v, i, t) = (k = get(row, (i, t), 0); k == 0 || ismissing(v[k]) ? NaN : Float64(v[k]))
        for h in 0:3
            fit = r.fits[h + 1]
            rows = findall(fit.esample)
            n = length(rows)
            y = [val(d.y, d.id[k], d.t[k] + h) - val(d.y, d.id[k], d.t[k] - 1) for k in rows]
            T1 = [val(d.d, d.id[k], d.t[k]) - val(d.d, d.id[k], d.t[k] - 1) for k in rows]
            ly = [val(d.y, d.id[k], d.t[k] - 1) - val(d.y, d.id[k], d.t[k] - 2) for k in rows]
            @test all(isfinite, y) && all(x -> x == 0 || x == 1, T1)
            ts = sort!(unique(d.t[rows]))
            Z = hcat(Float64.(coalesce.(d.x[rows], NaN)), ly,
                [Float64(d.t[k] == t) for k in rows, t in ts],
                [Float64(d.reg[k] == g) for k in rows, g in 1:4])
            a = wts === nothing ? ones(n) : d.w[rows]
            c, tr = findall(==(0), T1), findall(==(1), T1)
            G = pinv(Z[c, :]' * (a[c] .* Z[c, :]))
            γ = G * (Z[c, :]' * (a[c] .* y[c]))
            res = y .- Z * γ
            A1 = sum(a[tr])
            β = sum(a[tr] .* res[tr]) / A1
            @test coef(r, h)[1] ≈ β rtol = 1e-10
            @test all(isnan, coef(r, h)[2:3])          # only the treatment coefficient
            @test fit.did.ra && fit.did.ntreated == length(tr) && fit.dof.nobs == n
            # the influence function of the estimator, summed by unit
            m = Z[tr, :]' * a[tr] ./ A1
            ψ = zeros(n)
            ψ[tr] = a[tr] .* (res[tr] .- β) ./ A1
            ψ[c] = -a[c] .* res[c] .* (Z[c, :] * (G * m))
            su = zeros(60)
            for (k, rr) in enumerate(rows)
                su[d.id[rr]] += ψ[k]
            end
            @test cv.V[h + 1][1, 1] ≈ sum(abs2, su) rtol = 1e-8
        end
        # without covariates and further fixed effects the same estimand is a weighted
        # regression: regression adjustment with only time effects gives its coefficient
        ra = panellp(@formula(ldiff(y) ~ firstdiff(d) + fe(t) + fe(reg)), d; id = :id,
            time = :t, horizon = 0:3, estimator = LPDiD(reweight = true), weights = wts)
        @test all(f -> f.did.ra, ra.fits)
        @test_throws ArgumentError vcov(Homoskedastic(), ra)
    end
end

@testitem "negative horizons and the path of a regressor" tags = [:lpdid, :formula] begin
    using StatsModels, LinearAlgebra
    include(joinpath(@__DIR__, "support", "PLPDGPs.jl"))
    using .PLPDGPs: read_csv
    d = read_csv(joinpath(@__DIR__, "fixtures", "synthetic", "fx_unbalanced.csv"))
    kw = (; id = :id, time = :t)
    # y_{t+h} at h = −2 is the regression of the second lag of y, on the same rows
    r = panellp(@formula(leads(y) ~ x + lags(x, 1) + fe(id) + fe(t)), d; horizon = -2:1, kw...)
    r0 = panellp(@formula(leads(lags(y, 2:2)) ~ x + lags(x, 1) + fe(id) + fe(t)), d; horizon = 0, kw...)
    @test coef(r, -2) ≈ coef(r0, 0) rtol = 1e-12
    @test nobs(r, -2) == nobs(r0, 0)
    @test horizons(r) == -2:1
    # the long difference at h = −1 is identically zero: a zero coefficient and variance
    rl = panellp(@formula(ldiff(y) ~ x + fe(id) + fe(t)), d; horizon = -3:0, kw...)
    @test coef(rl, -1) == [0.0]
    s = summarize(rl, Clustered(:unit))
    @test s.se[3] == 0 && s.horizon == [-3, -2, -1, 0]
    @test_throws ArgumentError panellp(@formula(ldiff(y) ~ x + fe(id)), d; horizon = -2, kw...)
    @test_throws ArgumentError panellp(@formula(cumul(y) ~ x + fe(id)), d; horizon = -2:2, kw...)
    @test_throws ArgumentError panellp(@formula(leads(y) ~ x + fe(id)), d; horizon = -2:2,
        estimator = LagAugmented(ConstantLags(1)), kw...)
    # the split-panel jackknife and CD-OLS have no rule that depends on the horizon
    rs = panellp(@formula(ldiff(y) ~ x + lags(x, 1) + fe(id)), d; horizon = -2:1,
        estimator = SplitPanelJackknife(), kw...)
    @test all(f -> f.status === :ok, rs.fits)

    # leads(x): x_{t+1}, …, x_{t+h} at horizon h; leads(x, a:b) at every horizon
    rp = panellp(@formula(leads(y) ~ x + leads(x) + fe(id) + fe(t)), d; horizon = -1:3, kw...)
    @test rp.coefnames == ["x", "F1.x", "F2.x", "F3.x"]
    for h in -1:3
        fit = rp.fits[h + 2]
        @test findall(fit.active) == collect(1:(1 + max(h, 0)))
        if h >= 1
            rf = panellp(@formula(leads(y) ~ x + leads(x, 1:3) + fe(id) + fe(t)), d;
                horizon = h:h, kw...)
            h == 3 && @test coef(rp, h) ≈ coef(rf, h) rtol = 1e-12
        else
            rf = panellp(@formula(leads(y) ~ x + fe(id) + fe(t)), d; horizon = h:h, kw...)
            @test coef(rp, h)[1] ≈ coef(rf, h)[1] rtol = 1e-12
            @test nobs(rp, h) == nobs(rf, h)
        end
    end
    # a lead beyond the horizon is not an omitted regressor
    @test isempty(PanelLocalProjections._fitnotes(rp))
    # programmatic terms
    t = StatsModels.term
    rq = panellp(leads(t(:y)) ~ t(:x) + leads(t(:x)) + fe(t(:id)) + fe(t(:t)), d;
        horizon = -1:3, kw...)
    @test coef(rq, 3) ≈ coef(rp, 3)
    @test leads(t(:x), 2) isa PanelLocalProjections.LeadsTerm
    # a path term stands alone
    @test_throws ArgumentError panellp(@formula(leads(y) ~ x + lags(leads(x), 1) + fe(id)), d;
        horizon = 2, kw...)
    @test_throws ArgumentError panellp(@formula(leads(y) ~ (x ~ leads(w)) + fe(id)), d;
        horizon = 2, kw...)
    @test_throws ArgumentError panellp(@formula(leads(y) ~ x + leads(x, 1, 2) + fe(id)), d;
        horizon = 2, kw...)
    # lag augmentation: the criterion's h = 0 regression has no path, and the fit of
    # horizon h has the lags of p(h) and the leads up to h
    db = read_csv(joinpath(@__DIR__, "fixtures", "synthetic", "fx_balanced.csv"))
    la = (; id = :id, time = :t, shock = "s & X")
    ra = panellp(@formula(leads(y) ~ s & X + w + leads(w) + fe(id) + fe(t)), db;
        horizon = 0:3, estimator = LagAugmented(ICLags()), la...)
    rb = panellp(@formula(leads(y) ~ s & X + w + fe(id) + fe(t)), db;
        horizon = 0:3, estimator = LagAugmented(ICLags()), la...)
    @test ra.lags.p == rb.lags.p && ra.lags.ic ≈ rb.lags.ic
    rc = panellp(@formula(leads(y) ~ s & X + w + leads(w, 1:2) + fe(id) + fe(t)), db;
        horizon = 2:2, estimator = LagAugmented(ConstantLags(lagorder(ra)[3])), la...)
    @test coef(ra, 2)[1] ≈ coef(rc, 2)[1] rtol = 1e-10
    @test nobs(ra, 2) == nobs(rc, 2)
    # with the jackknife, and at horizons without leads
    rj = panellp(@formula(leads(y) ~ x + leads(x) + fe(id)), d; horizon = 0:2,
        estimator = SplitPanelJackknife(), kw...)
    @test [count(f.keep) for f in rj.fits] == [1, 2, 3]
end

@testitem "observation weights, pooled fits and the classical variance" tags = [:lpdid, :vcov] begin
    using StatsModels, LinearAlgebra
    include(joinpath(@__DIR__, "support", "PLPDGPs.jl"))
    using .PLPDGPs: read_csv
    d0 = read_csv(joinpath(@__DIR__, "fixtures", "synthetic", "fx_unbalanced.csv"))
    n = length(d0.id)
    f = @formula(leads(y) ~ x + lags(y, 1) + fe(id) + fe(t))
    kw = (; id = :id, time = :t, horizon = 0:2)
    r = panellp(f, d0; kw...)
    # equal weights are no weights, whatever their level
    rc = panellp(f, merge(d0, (; wt = fill(3.5, n))); weights = :wt, kw...)
    @test coefpath(rc) ≈ coefpath(r) rtol = 1e-12
    for v in (Robust(), Clustered(:unit), Clustered(:unit, :time), DriscollKraaySE(2), Homoskedastic())
        @test stderror(vcov(v, rc); term = "x") ≈ stderror(vcov(v, r); term = "x") rtol = 1e-9
    end
    # a zero or missing weight takes the row out, as `subset` does
    drop = [mod(k, 7) == 0 for k in 1:n]
    wz = Union{Missing, Float64}[drop[k] ? (iseven(k) ? 0.0 : missing) : 1.0 for k in 1:n]
    rz = panellp(f, merge(d0, (; wt = wz)); weights = :wt, kw...)
    rs = panellp(f, d0; subset = .!drop, kw...)
    @test coefpath(rz) ≈ coefpath(rs) rtol = 1e-12
    @test nobspath(rz) == nobspath(rs)
    # weighted least squares against explicit dummies
    w = [0.5 + mod(7k, 11) / 4 for k in 1:n]
    dw = merge(d0, (; wt = w))
    rw = panellp(@formula(leads(y) ~ x + fe(id) + fe(t)), dw; weights = :wt, id = :id,
        time = :t, horizon = 0, singletons = :keep)
    rows = findall(rw.fits[1].esample)
    ids, ts = sort!(unique(d0.id[rows])), sort!(unique(d0.t[rows]))
    X = hcat(d0.x[rows], [Float64(d0.id[k] == i) for k in rows, i in ids],
        [Float64(d0.t[k] == t) for k in rows, t in ts])
    W = Diagonal(w[rows])
    b = pinv(X' * W * X) * (X' * W * d0.y[rows])
    @test coef(rw, 0)[1] ≈ b[1] rtol = 1e-9
    e = d0.y[rows] .- X * b
    @test rw.fits[1].rss ≈ sum(w[rows] .* e .^ 2) rtol = 1e-9
    @test_throws ArgumentError panellp(f, merge(d0, (; wt = -w)); weights = :wt, kw...)
    @test_throws ArgumentError panellp(f, dw; weights = :wt, estimator = SplitPanelJackknife(), kw...)
    @test_throws ArgumentError panellp(f, dw; weights = :wt,
        estimator = LagAugmented(ConstantLags(1)), kw...)
    @test_throws ArgumentError panellp(f, dw; weights = :nope, kw...)

    # the classical variance: σ̂² (X'X)⁻¹ of the dummy regression
    r1 = panellp(@formula(leads(y) ~ x + fe(id) + fe(t)), d0; id = :id, time = :t,
        horizon = 0, singletons = :keep)
    rows = findall(r1.fits[1].esample)
    ids, ts = sort!(unique(d0.id[rows])), sort!(unique(d0.t[rows]))
    X = hcat(d0.x[rows], [Float64(d0.id[k] == i) for k in rows, i in ids],
        [Float64(d0.t[k] == t) for k in rows, t in ts])
    b = pinv(X' * X) * (X' * d0.y[rows])
    rss = sum(abs2, d0.y[rows] .- X * b)
    V = pinv(X' * X)[1, 1]
    nn, rk = length(rows), rank(X)
    c1, c0 = vcov(Homoskedastic(), r1), vcov(Homoskedastic(ssc = :none), r1)
    @test c1.V[1][1, 1] ≈ rss / (nn - rk) * V rtol = 1e-9
    @test c0.V[1][1, 1] ≈ rss / nn * V rtol = 1e-9
    @test c1.df[1][1] == nn - rk && isinf(c0.df[1][1])
    @test_throws ArgumentError vcov(Homoskedastic(), panellp(
        @formula(leads(y) ~ x + fe(id)), d0; estimator = SplitPanelJackknife(), kw...))
    @test_throws ArgumentError Homoskedastic(ssc = :mss)

    # a pooled fit is the regression of the mean response on the rows with every horizon
    fp = @formula(ldiff(y) ~ x + lags(x, 1) + fe(id) + fe(t))
    rp = panellp(fp, d0; id = :id, time = :t, horizon = 0:3, pooled = true)
    rc4 = panellp(fp, d0; id = :id, time = :t, horizon = 0:3, sample = :common, singletons = :keep)
    rpk = panellp(fp, d0; id = :id, time = :t, horizon = 0:3, pooled = true, singletons = :keep)
    # on a common sample the pooled coefficient is the mean of the horizons' coefficients
    @test coef(rpk, 3) ≈ sum(coef(rc4, h) for h in 0:3) ./ 4 rtol = 1e-10
    @test rp.horizons == 3:3 && rp.pooled == 0:3 && length(rp.fits) == 1
    @test summarize(rp, Clustered(:unit)).horizon == [3]
    @test_throws ArgumentError panellp(fp, d0; id = :id, time = :t, horizon = -2:3, pooled = true)
    @test_throws ArgumentError panellp(fp, d0; id = :id, time = :t, horizon = 0:3,
        pooled = true, estimator = SplitPanelJackknife())
end

@testitem "bug fixes of the 2026-11 review" tags = [:core, :formula, :absorb] begin
    using StatsModels, LinearAlgebra
    include(joinpath(@__DIR__, "support", "PLPDGPs.jl"))
    using .PLPDGPs: read_csv
    const PLP = PanelLocalProjections
    FIX = joinpath(@__DIR__, "fixtures", "synthetic")
    du = read_csv(joinpath(FIX, "fx_unbalanced.csv"))
    db = read_csv(joinpath(FIX, "fx_balanced.csv"))
    kw = (; id = :id, time = :t)

    # 1. nested shifts compose: they never go through an intermediate row
    cols = Dict(:x => Float64.(du.x), :y => Float64.(du.y))
    ctx = PLP._EvalContext(PLP.PanelIndex(du.id, du.t), cols)
    x = PLP.VarTerm(:x)
    same(a, b) = all(isequal.(a, b))
    @test same(PLP._cols(PLP.LagsTerm(PLP.LagsTerm(x, 1:1), 1:1), ctx), PLP._cols(PLP.LagsTerm(x, 2:2), ctx))
    @test same(PLP._cols(PLP.FirstDiffTerm(PLP.LagsTerm(x, 2:2)), ctx),
        PLP._cols(PLP.LagsTerm(PLP.FirstDiffTerm(x), 2:2), ctx))
    @test same(PLP._cols(PLP.LeadsTerm(PLP.LagsTerm(x, 1:1), 1:1), ctx), PLP._cols(x, ctx))
    # ldiff(lags(y, 1)) ~ lags(x, 1:2) is ldiff(y) ~ x + lags(x, 1) dated one period later:
    # the same rows and coefficients (it used to need the row of t + h as well)
    a = panellp(@formula(ldiff(lags(y, 1)) ~ lags(x, 1:2)), db; horizon = 1:3, estimator = CDOLS(),
        shock = "L1.x", kw...)
    b = panellp(@formula(ldiff(y) ~ x + lags(x, 1)), db; horizon = 1:3, estimator = CDOLS(), kw...)
    @test nobspath(a) == nobspath(b) == [1740, 1710, 1680]
    @test all(h -> coef(a, h) ≈ coef(b, h), 1:3)
    # an absent period and a row of missing values give the same estimates (the formula
    # uses x at t, so the row of missing values is not an observation either)
    gap = [!(du.id[k] == 3 && du.t[k] == 20) for k in eachindex(du.id)]
    dg = map(v -> v[gap], du)
    miss = [du.id[k] == 3 && du.t[k] == 20 for k in eachindex(du.id)]
    blank(v) = Union{Missing, eltype(v)}[miss[k] ? missing : v[k] for k in eachindex(v)]
    dm = merge(du, (; y = blank(du.y), x = blank(du.x), w = blank(du.w)))
    fn = @formula(leads(y) ~ x + lags(lags(x, 1), 1) + firstdiff(lags(x, 2:2)) + fe(id))
    rg, rm = panellp(fn, dg; horizon = 2, kw...), panellp(fn, dm; horizon = 2, kw...)
    @test nobspath(rg) == nobspath(rm)
    @test all(h -> coef(rg, h) ≈ coef(rm, h), 0:2)

    # 2. a missing id is an error, also as NaN or nothing
    @test_throws ArgumentError PLP.PanelIndex([1.0, NaN, 2.0], [1, 2, 1])
    @test_throws ArgumentError PLP.PanelIndex(Any[1, nothing, 2], [1, 2, 1])

    # 3. a regressor constant in the sample is absorbed whatever its level
    n = length(db.id)
    for c in (0.1, 0.35, 0.25, 1 / 3, 1e6 + 0.1)
        dc = merge(db, (; tax = fill(c, n)))
        @test_throws ArgumentError panellp(@formula(leads(y) ~ tax + fe(id)), dc; horizon = 1, kw...)
        r = panellp(@formula(leads(y) ~ x + tax + fe(id)), dc; horizon = 1, kw...)
        @test all(f -> f.absorbed[2] && !f.keep[2], r.fits)
        r0 = panellp(@formula(leads(y) ~ x + tax), dc; horizon = 1, kw...)
        @test all(f -> f.absorbed[3], r0.fits)
    end
    # a shock that stops varying at later horizons is :singular there, at any scale
    for s in (1.0, 10.0)
        tariff = [db.t[k] >= 59 ? 0.25 * s : 0.1 * s for k in 1:n]
        r = panellp(@formula(leads(y) ~ tariff + fe(id)), merge(db, (; tariff)); horizon = 0:4, kw...)
        @test [f.status for f in r.fits] == [:ok, :ok, :singular, :singular, :singular]
    end

    # 4. never more than maxiter sweeps
    for m in 2:7
        r = panellp(@formula(leads(y) ~ x + fe(id) + fe(t)), du; horizon = 0,
            absorb = AlternatingProjections(maxiter = m), kw...)
        @test r.fits[1].iterations <= m
    end

    # 5. t-LAHR does not depend on the units of another exposure of the same aggregate
    dt = read_csv(joinpath(FIX, "fx_tlahr.csv"))
    s2 = [1.0 + mod(7 * dt.id[k], 13) / 5 for k in eachindex(dt.id)]
    ft = @formula(leads(y) ~ s & X + s2 & X + fe(id) + fe(t))
    function tl(scale)
        r = panellp(ft, merge(dt, (; s2 = scale .* s2)); id = :id, time = :t, horizon = 0:3,
            shock = "s & X")
        cv = vcov(TLAHR(), r)
        return stderror(cv; term = "s & X"), [d[1] for d in cv.df]
    end
    se1, df1 = tl(1.0)
    for scale in (1e3, 1e4, 1e6)
        se, df = tl(scale)
        @test se ≈ se1 rtol = 1e-8
        @test df ≈ df1 rtol = 1e-6
    end

    # 6. nearly collinear instruments: undefined diagnostics are NaN, never an exception
    dv = read_csv(joinpath(FIX, "fx_lpiv.csv"))
    noise = [sin(13.7 * k) for k in eachindex(dv.id)]
    for ε in (1e-2, 1e-6, 3e-9, 1e-8)
        z2 = dv.s .* dv.mp .+ ε .* noise
        r = panellp(@formula(leads(y) ~ (s & ffr ~ s & mp + z2) + fe(id) + fe(t)),
            merge(dv, (; z2)); id = :id, time = :t, horizon = 0:1)
        fs = firststage(r; vcov = Robust())
        @test all(k -> isnan(fs.kpF[k]) == isnan(fs.effF[k]), 1:2)
        ε == 1e-2 && @test all(isfinite, fs.kpF) && all(isfinite, fs.cdF)
    end
end
