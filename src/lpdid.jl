# LP-DiD (Dube, Girardi, Jordà & Taylor 2025, "A local projections approach to
# difference-in-differences", Journal of Applied Econometrics 40(7)): the long difference
# y_{i,t+h} − y_{i,t−1} regressed on the treatment switch ΔD_it and time effects, on the rows
# that enter treatment at t or are clean controls at horizon h.
#
# Written from the paper (Sections 3–4 and the online appendix) and from the documented
# behaviour of the Stata command `lpdid` (Busch & Girardi; GPL-3), which only served as a
# reference to compare outputs with: no code of it is translated here.

function LPDiD(;
        control::Symbol = :clean, nonabsorbing = nothing, firsttreat::Bool = false,
        oneoff::Bool = false, untreated_before::Bool = false, reweight = false,
        pmd = nothing, compat::Symbol = :none)
    L = nonabsorbing === nothing || nonabsorbing === false ? nothing :
        nonabsorbing === true ? -1 :
        nonabsorbing isa Integer ? Int(nonabsorbing) :
        throw(ArgumentError(
        "nonabsorbing is the number of periods L after which effects stabilise (or `true` " *
        "with firsttreat and control = :notyet or :nevertreated); got $(repr(nonabsorbing))"))
    # `reweight = :lpdid` (0.2.1) is `reweight = true, compat = :lpdid`
    reweight === :lpdid && (compat = :lpdid)
    rw = reweight === false ? :none :
         reweight === true || reweight === :lpdid ? :equal :
         throw(ArgumentError("reweight must be true or false; got $(repr(reweight))"))
    k = pmd === nothing || pmd === false ? 0 :
        pmd === :max ? -1 :
        pmd isa Integer && pmd >= 1 ? Int(pmd) :
        throw(ArgumentError("pmd must be a number of periods k ≥ 1 or :max; got $(repr(pmd))"))
    return LPDiD(control, L, firsttreat, oneoff, untreated_before, rw, k, compat)
end

function _label(e::LPDiD)
    opts = String[]
    e.control === :clean || push!(opts, "control = :$(e.control)")
    e.nonabsorbing === nothing ||
        push!(opts, "nonabsorbing = $(e.nonabsorbing < 0 ? true : e.nonabsorbing)")
    e.firsttreat && push!(opts, "firsttreat = true")
    e.oneoff && push!(opts, "oneoff = true")
    e.untreated_before && push!(opts, "untreated_before = true")
    e.reweight === :none || push!(opts, "reweight = true")
    e.pmd == 0 || push!(opts, "pmd = $(e.pmd < 0 ? ":max" : e.pmd)")
    e.compat === :none || push!(opts, "compat = :$(e.compat)")
    return "LPDiD($(join(opts, ", ")))"
end

# Treatment is handled as absorbing by the clean-control rule D_{i,t+h} = 0: the default,
# and the first-treatment estimator with not-yet-treated or never-treated controls (DGJT
# Section 4.2.2), where everything after a unit's first entry is left out.
function _pseudoabsorbing(e::LPDiD)
    e.nonabsorbing !== nothing && e.firsttreat &&
        e.control !== :clean
end
_absorbingrule(e::LPDiD) = e.nonabsorbing === nothing || _pseudoabsorbing(e)

# ---------------------------------------------------------------------------------------
# Treatment histories
# ---------------------------------------------------------------------------------------

# For each data row, the first period of its unit in the data.
function _firstperiod(p::PanelIndex)
    tmin = fill(typemax(Int), nunits(p))
    for r in eachindex(p.unit)
        tmin[p.unit[r]] = min(tmin[p.unit[r]], p.period[r])
    end
    return tmin
end

# Running sum of x over each unit's rows in time order (a NaN adds nothing).
function _runsum(p::PanelIndex, x::Vector{Float64})
    out = Vector{Float64}(undef, length(x))
    for i in axes(p.rowof, 1)
        acc = 0.0
        for j in axes(p.rowof, 2)
            r = p.rowof[i, j]
            r == 0 && continue
            isnan(x[r]) || (acc += x[r])
            out[r] = acc
        end
    end
    return out
end

# Number of switches of each unit up to each of its rows, counted between consecutive
# observed values: a change hidden in an absent period or behind a missing value counts.
function _switches(p::PanelIndex, D::Vector{Float64})
    out = Vector{Float64}(undef, length(D))
    for i in axes(p.rowof, 1)
        last, n = NaN, 0
        for j in axes(p.rowof, 2)
            r = p.rowof[i, j]
            r == 0 && continue
            if !isnan(D[r])
                (!isnan(last) && D[r] != last) && (n += 1)
                last = D[r]
            end
            out[r] = n
        end
    end
    return out
end

# Number of units whose observed treatment values fall from 1 to 0 at some point.
function _exits(p::PanelIndex, D::Vector{Float64})
    n = 0
    for i in axes(p.rowof, 1)
        last = NaN
        for j in axes(p.rowof, 2)
            r = p.rowof[i, j]
            (r == 0 || isnan(D[r])) && continue
            if last == 1 && D[r] == 0
                n += 1
                break
            end
            last = D[r]
        end
    end
    return n
end

# x_{i,t−1} subtracted from x_{i,t} (NaN when either is undefined).
_delta(p::PanelIndex, x::Vector{Float64}) = x .- _gather(x, _shift(p, -1))

# Mean of the outcome over the k periods before t (all required), or over every earlier
# row of the unit when k < 0 (`pmd = :max`; undefined when the row of t − 1 is absent or no
# earlier outcome is observed). `term` is the outcome's term, evaluated at each lag.
function _premean(p::PanelIndex, term, ctx::_EvalContext, k::Int)
    if k > 0
        acc = zeros(nrows(p))
        for l in 1:k
            acc .+= _col1(term, ctx, -l)
        end
        return acc ./ k
    end
    y = _col1(term, ctx)
    out = fill(NaN, length(y))
    for i in axes(p.rowof, 1)
        s, n = 0.0, 0
        for j in axes(p.rowof, 2)
            r = p.rowof[i, j]
            if r != 0
                # the mean of the rows before this one, when the previous period is a row
                j > 1 && p.rowof[i, j - 1] != 0 && n > 0 && (out[r] = s / n)
                isnan(y[r]) || (s += y[r]; n += 1)
            end
        end
    end
    return out
end

"""
    _did_resolve(est::LPDiD, d, ctx, horizons) -> _DiD

The treatment histories of an LP-DiD call: the treatment D (cumulated for a one-off
treatment under the first-treatment estimator), its switch ΔD, the counters of the
not-yet-treated and first-treatment rules, the base flag C₀ of the non-absorbing rules, and
the base period of the long difference (y_{t−1}, or a pre-treatment mean with `pmd`).
"""
function _did_resolve(est::LPDiD, d::PanelDesign, ctx::_EvalContext,
        horizons::UnitRange{Int}, tvar::Symbol)
    p = d.pidx
    n = nrows(p)
    D = copy(d.columns[tvar])
    if _pseudoabsorbing(est) && est.oneoff
        ev = _runsum(p, D)
        D = [ev[r] > 0 ? 1.0 : D[r] for r in 1:n]
    end
    dD = _delta(p, D)
    tmin = _firstperiod(p)
    first_t(r) = tmin[p.unit[r]]
    # never treated: D = 0 at every row of the unit where it is observed (and at some)
    mx = fill(-Inf, nunits(p))
    for r in 1:n
        isnan(D[r]) || (mx[p.unit[r]] = max(mx[p.unit[r]], D[r]))
    end
    never = BitVector(mx[p.unit[r]] == 0 for r in 1:n)
    # not yet treated at t: no switch up to t, and untreated at the unit's first row.
    # Switches are changes between consecutive observed values; `lpdid` sums the observed
    # first differences, which misses a switch hidden in a gap.
    P = est.compat === :lpdid ? _runsum(p, abs.(dD)) : _switches(p, D)
    entry = fill(NaN, nunits(p))
    for r in 1:n
        p.period[r] == first_t(r) && (entry[p.unit[r]] = D[r])
    end
    E0 = BitVector(entry[p.unit[r]] == 0 for r in 1:n)
    # first treatment: the number of events (one-off) or of switches up to t
    Q = est.oneoff && !_pseudoabsorbing(est) ? _runsum(p, D) : P
    L = est.nonabsorbing === nothing ? 0 : max(est.nonabsorbing, 0)
    L1 = max(L, 1)
    # the unit was untreated k periods ago, as far as the rule of the treatment type can
    # tell: no switch then (persistent) or no event then (one-off); with `untreated_before`
    # also when that period precedes the unit's first row
    function quiet(k::Int)
        idx = _shift(p, -k)
        out = falses(n)
        for r in 1:n
            j = idx[r]
            s = p.period[r] - k
            if est.oneoff
                out[r] = (j != 0 && D[j] == 0) ||
                         (est.untreated_before && s < first_t(r))
            else
                out[r] = (j != 0 && dD[j] == 0) ||
                         (est.untreated_before && s - 1 < first_t(r) &&
                          (s < first_t(r) || (j != 0 && D[j] == 0)))
            end
        end
        return out
    end
    C0 = falses(n)
    if !_absorbingrule(est)
        for r in 1:n
            C0[r] = est.oneoff ? ((D[r] == 1 && dD[r] == 1) || (D[r] == 0 && dD[r] == 0)) :
                    (dD[r] == 0 || dD[r] == 1)
        end
        for k in 1:L1
            C0 .&= quiet(k)
        end
    end
    term = d.response.term
    y = _col1(term, ctx)
    k = est.pmd
    # with a non-absorbing treatment and no first-treatment rule, the pre-treatment mean
    # is over the stabilisation window
    k < 0 && !_absorbingrule(est) && (k = L)
    k == 0 && est.pmd != 0 &&
        throw(ArgumentError(
            "LPDiD(pmd = :max) with nonabsorbing = 0 has no pre-treatment window; give pmd a " *
            "number of periods"))
    base = k == 0 ? _col1(term, ctx, -1) : _premean(p, term, ctx, k)
    if k > 0 && est.nonabsorbing !== nothing
        for kk in 1:(k + L - 1)
            q = quiet(kk)
            for r in 1:n
                q[r] || (base[r] = NaN)
            end
        end
    end
    respok = trues(n)
    for h in horizons
        respok .&= isfinite.(_col1(term, ctx, h) .- base)
    end
    return _DiD(est, horizons, D, dD, never, P, E0, Q, C0, y, base, respok, false)
end

"""
    _cleanflag(s::_DiD, pidx, h) -> BitVector

The rows of the LP-DiD sample of horizon h: those entering treatment at t and the clean
controls of horizon h (DGJT Sections 3.1 and 4.2), before any missing-value rule. A
condition on an unobserved value is false.

- Absorbing rule (the default; first treatment with not-yet-treated or never-treated
  controls): ΔD_t = 1, or ΔD_t = 0 and D_{t+h} = 0 (D_t = 0 at pre-treatment horizons).
  `compat = :lpdid` does not require ΔD_t = 0 of a control, which lets in the rows where a
  unit treated from its first row switches off.
- `nonabsorbing = L`: ΔD_t ∈ {0, 1} and no switch in the L periods before t, nor in
  t + 1, …, t + h; at h < 0, the same at t, t − 1, …, t − |h| + 1. With `oneoff`, events
  (D = 1) in place of switches.
- `control = :nevertreated` keeps as controls only the rows of never-treated units;
  `:notyet` those of units untreated from their first row through t + h;
  `firsttreat` keeps the rows with at most one switch (event) through t + h.
"""
function _cleanflag(s::_DiD, p::PanelIndex, h::Int)
    e = s.spec
    n = nrows(p)
    hp = max(h, 0)
    lead = _shift(p, hp)
    flag = falses(n)
    if _absorbingrule(e)
        lp = e.compat === :lpdid
        for r in 1:n
            j = lead[r]
            flag[r] = s.dD[r] == 1 || ((lp || s.dD[r] == 0) && j != 0 && s.D[j] == 0)
        end
    elseif h >= 0
        flag .= s.C0
        for k in 1:h
            idx = _shift(p, k)
            x = e.oneoff ? s.D : s.dD
            for r in 1:n
                j = idx[r]
                flag[r] &= j != 0 && x[j] == 0
            end
        end
    else
        flag .= s.C0
        for k in 1:(-h - 1)
            idx = _shift(p, -k)
            for r in 1:n
                j = idx[r]
                flag[r] &= j != 0 && s.C0[j]
            end
        end
    end
    for r in 1:n
        flag[r] || continue
        j = lead[r]
        control = s.dD[r] == 0
        if e.control === :nevertreated && control && !s.never[r]
            flag[r] = false
        elseif e.control === :notyet && control &&
               !(j != 0 && s.P[j] == 0 && s.E0[r])
            flag[r] = false
        end
        e.firsttreat && !(j != 0 && s.Q[j] <= 1) && (flag[r] = false)
    end
    return flag
end

# ---------------------------------------------------------------------------------------
# Specification checks
# ---------------------------------------------------------------------------------------

# The treatment column of an LP-DiD formula: the shock is its first difference.
function _did_treatment(d::PanelDesign)
    length(d.shockcols) == 1 || throw(ArgumentError(
        "LPDiD has one treatment: the shock must be the single column firstdiff(D)"))
    j = d.shockcols[1]
    ct = d.colterm[j]
    if ct !== nothing
        t = ct[1]
        t isa FirstDiffTerm && t.term isa VarTerm && return t.term.sym
    end
    throw(ArgumentError(
        "LPDiD regresses the long difference on the treatment switch: write the shock as " *
        "firstdiff(D), with D the 0/1 treatment column (e.g. ldiff(y) ~ firstdiff(D) + " *
        "fe(time)); got $(repr(d.coefnames[j]))"))
end

function _validate(est::LPDiD, d::PanelDesign, f)
    _check_option(est.control, (:clean, :nevertreated, :notyet), "control")
    _check_option(est.compat, (:none, :lpdid), "compat")
    d.response isa LDiffResponse || throw(ArgumentError(
        "LPDiD estimates the long difference y_{t+h} − y_{t−1}: write the response as " *
        "ldiff(y)"))
    tvar = _did_treatment(d)
    D = d.columns[tvar]
    all(x -> isnan(x) || x == 0 || x == 1, D) || throw(ArgumentError(
        "LPDiD needs a binary treatment: column :$tvar must hold 0, 1 or missing"))
    isempty(d.timefe) && throw(ArgumentError(
        "LPDiD needs time effects (it compares units within a period): add " *
        "fe(<time column>), or group × time effects such as fe(region) & fe(<time " *
        "column>), to the formula"))
    if est.nonabsorbing === nothing
        (est.firsttreat || est.oneoff || est.untreated_before) && throw(ArgumentError(
            "firsttreat, oneoff and untreated_before are rules for a non-absorbing " *
            "treatment: set nonabsorbing = L (the periods after which effects stabilise)"))
        est.control === :notyet && throw(ArgumentError(
            "control = :notyet is a rule for a non-absorbing treatment (nonabsorbing = L); " *
            "with an absorbing treatment the default clean controls are the not-yet-treated " *
            "and the never-treated"))
        # also a switch-off hidden in a gap or behind a missing value (which `lpdid`
        # does not see)
        nexit = est.compat === :lpdid ? count(==(-1), _delta(d.pidx, D)) :
                _exits(d.pidx, D)
        nexit == 0 || throw(ArgumentError(
            "the treatment :$tvar is not absorbing: it switches off for $nexit unit(s) " *
            "(or row(s)). Set LPDiD(nonabsorbing = L), with L the periods after which " *
            "effects are assumed to stabilise (DGJT Section 4.2)"))
    else
        est.nonabsorbing < 0 && !_pseudoabsorbing(est) &&
            throw(ArgumentError(
                "nonabsorbing = true (no stabilisation window) is for the first-treatment " *
                "estimator: firsttreat = true with control = :notyet or :nevertreated; " *
                "otherwise give the number of periods L"))
    end
    return nothing
end

function _augment(est::LPDiD, d::PanelDesign, ctx::_EvalContext, base::BitVector,
        horizons::UnitRange{Int}, singletons::Symbol, absorb::AbsorptionMethod)
    return d, _did_resolve(est, d, ctx, horizons, _did_treatment(d)), nothing
end

# Whether a fixed-effect dimension is the period itself (not a finer group × period).
function _timeeffect(codes::Vector{Int32}, p::PanelIndex)
    seen = Dict{Int32, Int}()
    back = Dict{Int, Int32}()
    for r in eachindex(codes)
        c = codes[r]
        c == 0 && continue
        t = p.period[r]
        get!(seen, c, t) == t || return false
        get!(back, t, c) == c || return false
    end
    return !isempty(seen)
end

# ---------------------------------------------------------------------------------------
# Fits
# ---------------------------------------------------------------------------------------

# Counts of the treated rows and of all rows of a fit by group, per score dimension.
function _didinfo(d::PanelDesign, s::_DiD, rows::Vector{Int}, ra::Bool)
    treated = [zeros(Int32, G) for G in d.ngroups]
    nr = [zeros(Int32, G) for G in d.ngroups]
    n1 = 0
    for r in rows
        tr = s.dD[r] == 1
        n1 += tr
        for g in eachindex(d.groups)
            c = d.groups[g][r]
            nr[g][c] += Int32(1)
            tr && (treated[g][c] += Int32(1))
        end
    end
    return DiDInfo(n1, treated, nr, ra)
end

# Whether the equally weighted fit is a regression adjustment (covariates, or fixed effects
# besides the time effects) rather than a weighted regression.
function _usera(s::_DiD, d::PanelDesign)
    s.spec.reweight === :none && return false
    # the weighted regression reweights by period: it is the equally weighted estimate
    # only when the time effects are the single fixed effect
    return size(d.X, 2) > 1 || length(d.fe) > 1 || !_timeeffect(d.fe[1], d.pidx)
end

"""
    _rwweights(d, s, mask, extra) -> (weights, mask)

Weights of the equally weighted LP-DiD regression (DGJT Section 3.3): within each period,
every row is weighted by 1/(1 − p_t), p_t the (weighted) share of treated rows of the
period in the sample, which turns the variance weights N_t p_t (1 − p_t) of the cohorts
into their numbers of treated rows. Periods without a treated row or without a control
row leave the sample. `extra` marks rows outside the regression that `compat = :lpdid`
counts among the controls of their period (as `lpdid` 1.1.0 does: clean rows whose ΔD is
undefined).
"""
function _rwweights(d::PanelDesign, s::_DiD, mask::BitVector, extra::BitVector)
    p = d.pidx
    T = nspan(p)
    A1, A0 = zeros(T), zeros(T)
    a(r) = isempty(d.weights) ? 1.0 : d.weights[r]
    for r in findall(mask)
        k = p.period[r] - p.pmin + 1
        s.dD[r] == 1 ? (A1[k] += a(r)) : (A0[k] += a(r))
    end
    if s.spec.compat === :lpdid
        for r in findall(extra)
            A0[p.period[r] - p.pmin + 1] += a(r)
        end
    end
    w = zeros(nrows(p))
    out = copy(mask)
    for r in findall(mask)
        k = p.period[r] - p.pmin + 1
        if A1[k] > 0 && A0[k] > 0
            w[r] = a(r) * (A1[k] + A0[k]) / A0[k]
        else
            out[r] = false
        end
    end
    return w, out
end

function _fit_horizon(s::_DiD, d::PanelDesign, h::Int, y::Vector{Float64},
        mask::BitVector, singletons::Int, absorb::AbsorptionMethod)
    cols = _horizoncols(d, h)
    if _usera(s, d)
        return _fit_ra(s, d, h, y, mask, absorb, cols)
    end
    w = d.weights
    if s.spec.reweight !== :none
        extra = s.spec.compat === :lpdid ? _rwextra(s, d, h, y) : falses(0)
        w, mask = _rwweights(d, s, mask, extra)
    end
    fit = _fit_within(d, h, y, mask, singletons, absorb, cols, 0; weights = w)
    fit.status === :ok || return fit
    return _withdid(fit, _didinfo(d, s, findall(fit.esample), false))
end

# Rows that `lpdid` 1.1.0 counts among the controls of the reweighting although they are
# not in the regression: clean at horizon h, with the outcome defined, and ΔD undefined
# (with a common sample: clean and with the outcome defined at every horizon).
function _rwextra(s::_DiD, d::PanelDesign, h::Int, y::Vector{Float64})
    flag = _cleanflag(s, d.pidx, h)
    if s.common
        for hh in s.horizons
            flag .&= _cleanflag(s, d.pidx, hh)
        end
        flag .&= s.respok
    end
    out = falses(length(flag))
    for r in eachindex(flag)
        out[r] = flag[r] && d.subset[r] && isfinite(y[r]) && isnan(s.dD[r]) &&
                 (isempty(d.weights) || isfinite(d.weights[r])) &&
                 all(c -> c[r] != 0, d.fe) && all(c -> c[r] != 0, d.groups)
    end
    return out
end

# Solve the normal equations of a regression on the dummies of the fixed effects `fs`
# (weights `a`) for the right-hand sides `rhs` (one vector per dimension, by level):
# coefficients c with Σ_{i ∈ level} a_i Σ_k c_k[level_k(i)] = rhs. Block Gauss–Seidel;
# exact in one pass with a single dimension.
function _fe_solve(fs::_FESample, a::Vector{Float64}, rhs::Vector{Vector{Float64}};
        tol::Float64 = 1e-12, maxiter::Int = 10_000)
    D = length(fs.codes)
    n = length(a)
    c = [zeros(length(r)) for r in rhs]
    wsum = [zeros(length(r)) for r in rhs]
    for k in 1:D, i in 1:n

        wsum[k][fs.codes[k][i]] += a[i]
    end
    fitted = zeros(n)
    acc = zeros(maximum(length, rhs; init = 0))
    for it in 1:maxiter
        delta = 0.0
        scale = 0.0
        for k in 1:D
            ck, codes = c[k], fs.codes[k]
            fill!(view(acc, 1:length(ck)), 0.0)
            @inbounds for i in 1:n
                acc[codes[i]] += a[i] * (fitted[i] - ck[codes[i]])
            end
            @inbounds for l in eachindex(ck)
                new = wsum[k][l] > 0 ? (rhs[k][l] - acc[l]) / wsum[k][l] : 0.0
                delta = max(delta, abs(new - ck[l]))
                scale = max(scale, abs(new))
                acc[l] = new - ck[l]
                ck[l] = new
            end
            @inbounds for i in 1:n
                fitted[i] += acc[codes[i]]
            end
        end
        (D == 1 || delta <= tol * scale) && return c, fitted, true
    end
    return c, fitted, false
end

"""
    _fit_ra(s, d, h, y, mask, absorb, cols) -> HorizonFit

Regression-adjustment LP-DiD (DGJT Section 4.1.1), the equally weighted estimator with
covariates or further fixed effects: the outcome model, the long difference on the
covariates and the fixed effects Z, is fitted on the clean controls (weighted least
squares with the observation weights a), and the estimate is the weighted mean over the
treated rows of the long difference less its prediction,

    β̂ = Σ_{treated} a_i (y_i − Z_i'γ̂) / A₁,   A₁ = Σ_{treated} a_i.

Its influence function is a_i (r_i − β̂)/A₁ for a treated row and −a_i r_i p_i for a
control, with r the residuals and p_i = Z_i'(Z₀'AZ₀)⁻ m, m = Σ_{treated} a_i Z_i / A₁, the
weight of the control in the predictions for the treated. The fit stores it as its score
(with bread 1), so `Clustered(…; ssc = :none)` is the variance of the estimator; only the
treatment coefficient is reported.

A treated row whose period or other fixed-effect level, or whose value of a two-valued
covariate, has no control row cannot be predicted and leaves the sample; so do the control
rows of such cells (a single pass, as `lpdid`).
"""
function _fit_ra(s::_DiD, d::PanelDesign, h::Int, y::Vector{Float64}, mask::BitVector,
        absorb::AbsorptionMethod, cols::Vector{Int})
    K = length(d.coefnames)
    j = d.shockcols[1]
    active = falses(K)
    active[j] = true
    mask = _ra_guard(s, d, mask, cols)
    rows = findall(mask)
    tr = Int[i for (i, r) in enumerate(rows) if s.dD[r] == 1]
    ct = Int[i for (i, r) in enumerate(rows) if s.dD[r] == 0]
    (isempty(tr) || isempty(ct)) && return _empty_fit(d, h, mask, 0; active)
    a = isempty(d.weights) ? ones(length(rows)) : d.weights[rows]
    xc = Int[c for c in cols if c != j]
    nx = length(xc)
    # the outcome model on the controls
    crow = rows[ct]
    ac = a[ct]
    A = Matrix{Float64}(undef, length(ct), 1 + nx)
    A[:, 1] = view(y, crow)
    for (c, col) in enumerate(xc)
        A[:, c + 1] = view(d.X, crow, col)
    end
    raw = _colstats(view(A, :, 2:(1 + nx)), ac)
    X0 = A[:, 2:end]
    y0 = A[:, 1]
    fs = _fe_sample(d.fe, crow)
    it, conv = _absorb_matrix!(A, fs, absorb, ac)
    sq = sqrt.(ac)
    A .*= sq
    γ = zeros(nx)
    B = zeros(nx, nx)
    xkeep = falses(nx)
    resid = A[:, 1]
    if nx > 0
        absd = _absorbed_cols(view(A, :, 2:(1 + nx)), raw, true, false)
        if !all(absd)
            o = _ols(view(A, :, 2:(1 + nx)), A[:, 1], raw, true, false, Int[];
                absorbed = absd)
            xkeep = o.keep
            γ[xkeep] = o.β
            B[xkeep, xkeep] = o.B
            resid = o.resid
        end
    end
    resid ./= sq                                  # r_i of the controls
    # fixed effects of the outcome model: y₀ − X₀γ̂ less its residual is their fitted sum
    u0 = y0 .- X0 * γ
    rhs = [zeros(length(cnt)) for cnt in fs.counts]
    for k in eachindex(rhs), i in eachindex(ct)

        rhs[k][fs.codes[k][i]] += ac[i] * u0[i]
    end
    cfe, _, conv1 = _fe_solve(fs, ac, rhs)
    # levels of the treated rows in the controls' coding
    trow = rows[tr]
    at = a[tr]
    A1 = sum(at)
    tcodes = _codes_like(d.fe, crow, fs, trow)
    rt = Vector{Float64}(undef, length(tr))
    for (i, r) in enumerate(trow)
        pred = 0.0
        for c in 1:nx
            pred += d.X[r, xc[c]] * γ[c]
        end
        for k in eachindex(cfe)
            pred += cfe[k][tcodes[k][i]]
        end
        rt[i] = y[r] - pred
    end
    β = sum(at .* rt) / A1
    # p: the element of the span of Z₀ with Z₀'A p = m
    md = [zeros(length(cnt)) for cnt in fs.counts]
    for k in eachindex(md), i in eachindex(tr)

        md[k][tcodes[k][i]] += at[i] / A1
    end
    _, u, conv2 = _fe_solve(fs, ac, md)
    p = u
    if any(xkeep)
        kx = findall(xkeep)
        mx = zeros(length(kx))
        for (c, col) in enumerate(kx), (i, r) in enumerate(trow)

            mx[c] += at[i] * d.X[r, xc[col]] / A1
        end
        g = mx .- X0[:, kx]' * (ac .* u)
        q = B[kx, kx] * g
        # X̃₀ q, with the weighted within-transformed covariates (unscaled)
        p = u .+ (view(A, :, 1 .+ kx) * q) ./ sq
    end
    ψ = zeros(length(rows), 1)
    for (i, k) in enumerate(tr)
        ψ[k, 1] = at[i] * (rt[i] - β) / A1
    end
    for (i, k) in enumerate(ct)
        ψ[k, 1] = -ac[i] * resid[i] * p[i]
    end
    keep = falses(K)
    keep[j] = true
    groups = [view(g, rows) for g in d.groups]
    stats, present = _score_stats(ψ, keep, groups, d.ngroups)
    D = length(d.fe)
    dof = DofInfo(length(rows), 1, 0, zeros(Int, D), zeros(Int, D, D), present,
        Matrix{Bool}(falses(D, length(d.groups))), 0)
    coef = fill(NaN, K)
    coef[j] = β
    bread = fill(NaN, K, K)
    bread[j, j] = 1.0
    xsd = fill(NaN, K)
    xsd[j] = 1.0
    rss = sum(at .* (rt .- β) .^ 2) + sum(ac .* resid .^ 2)
    tss = sum(a .* view(y, rows) .^ 2)
    return HorizonFit(h, :ok, coef, keep, falses(K), bread, xsd, mask,
        _present(groups[2], d.ngroups[2]), dof, stats, it, conv && conv1 && conv2, active,
        0, SyntheticSeries[], falses(size(d.Z, 2)), nothing, rss, tss,
        _didinfo(d, s, rows, true))
end

# The codes of the rows `trow` in the coding that `fs` gives the rows `crow` (every level
# of `trow` is among those of `crow`).
function _codes_like(fe::Vector{Vector{Int32}}, crow::Vector{Int}, fs::_FESample,
        trow::Vector{Int})
    out = Vector{Vector{Int32}}(undef, length(fe))
    for (k, c) in enumerate(fe)
        remap = zeros(Int32, maximum(c; init = Int32(0)))
        for (i, r) in enumerate(crow)
            remap[c[r]] = fs.codes[k][i]
        end
        out[k] = Int32[remap[c[r]] for r in trow]
    end
    return out
end

# The rows of `mask` whose cells all hold a control row: every fixed-effect level, and
# every value of a covariate that takes exactly two values on the `subset` rows.
function _ra_guard(s::_DiD, d::PanelDesign, mask::BitVector, cols::Vector{Int})
    rows = findall(mask)
    ctl = [r for r in rows if s.dD[r] == 0]
    out = copy(mask)
    for c in d.fe
        has = falses(maximum(c; init = Int32(0)))
        for r in ctl
            has[c[r]] = true
        end
        for r in rows
            has[c[r]] || (out[r] = false)
        end
    end
    j = d.shockcols[1]
    for col in cols
        col == j && continue
        x = view(d.X, :, col)
        vals = Float64[]
        for r in eachindex(x)
            (d.subset[r] && isfinite(x[r])) || continue
            x[r] in vals && continue
            push!(vals, x[r])
            length(vals) > 2 && break
        end
        length(vals) == 2 || continue
        has = [any(r -> x[r] == v, ctl) for v in vals]
        for r in rows
            has[findfirst(==(x[r]), vals)] || (out[r] = false)
        end
    end
    return out
end

# ---------------------------------------------------------------------------------------
# Joint inference across horizons: the average effect and the pre-trend test
# ---------------------------------------------------------------------------------------

# Cluster sums of the influence function of coefficient j at one horizon, by group of the
# score dimension g: S_g B e_j.
function _influence(fit::HorizonFit, g::Int, j::Int)
    kidx = findall(fit.keep)
    return fit.scores.S[g][:, kidx] * fit.bread[kidx, j]
end

# A horizon whose response is a linear combination of the regressors (the long difference
# at h = −1; a pre-treatment difference spanned by lags of the outcome): its coefficient
# is zero by construction.
_degenerate(fit::HorizonFit) = fit.rss <= 1e-20 * fit.tss || fit.tss == 0

function _clusterdim(r::PanelLocalProjection, cluster::Symbol)
    g = _dimindex(r, cluster)
    return g
end

"""
    pretrendtest(r; cluster = :unit, horizons = h ≤ −2, term = first shock)

Joint test that the coefficients of `term` are zero at the pre-treatment horizons: the
Wald statistic b'V⁻¹b/P, with V the covariance of the P coefficients across horizons
clustered by `cluster` (from the stored score sums; no small-sample factor), referred to
F(P, M − 1), M the clusters in the samples of the tested horizons. As `lpdid`'s
`pretrend_test`. Horizons that were not estimated, and those whose coefficient is zero by
construction (h = −1 with `ldiff`; a difference spanned by the lagged outcomes), are left
out. Returns `(F, p, df, df_r, horizons)`; `F` is NaN when V is singular or M ≤ P.
"""
function pretrendtest(r::PanelLocalProjection; cluster::Symbol = :unit,
        horizons = nothing, term = r.shocks[1])
    g = _clusterdim(r, cluster)
    j = _termindex(r.coefnames, term)
    hs = horizons === nothing ? [h for h in r.horizons if h <= -2] :
         collect(Int, horizons)
    all(h -> h in r.horizons, hs) ||
        throw(ArgumentError("horizons $hs are not all among the estimated $(r.horizons)"))
    fits = [f
            for f in (_fit(r, h) for h in hs)
            if f.status === :ok && f.keep[j] && !_degenerate(f)]
    used = [f.h for f in fits]
    P = length(fits)
    P == 0 && return (; F = NaN, p = NaN, df = 0, df_r = 0, horizons = used)
    Ψ = reduce(hcat, [_influence(f, g, j) for f in fits])
    M = _clusters(r, fits, g)
    b = [f.coef[j] for f in fits]
    V = LinearAlgebra.Symmetric(Ψ' * Ψ)
    ev = LinearAlgebra.eigvals(V)
    if M <= P || first(ev) <= 1e-12 * last(ev)
        return (; F = NaN, p = NaN, df = P, df_r = M - 1, horizons = used)
    end
    F = (b' * (V \ b)) / P
    return (; F, p = StatsFuns.fdistccdf(P, M - 1, F), df = P, df_r = M - 1,
        horizons = used)
end

# Number of groups of score dimension g with a row in the sample of any of the fits.
function _clusters(r::PanelLocalProjection, fits::Vector{HorizonFit}, g::Int)
    seen = falses(r.ngroups[g])
    for f in fits
        if f.did !== nothing
            seen .|= f.did.rows[g] .> 0
        elseif g == 2
            seen .|= f.timepresent
        else
            # without row counts, a group is in the sample when its score sum is not zero
            S = f.scores.S[g]
            for c in axes(S, 1)
                seen[c] |= any(!iszero, view(S, c, :))
            end
        end
    end
    return count(seen)
end

"""
    averageeffect(r; cluster = :unit, level = 0.95)

The average of the LP-DiD event-study coefficients over the post-treatment horizons
(h ≥ 0) and over the pre-treatment ones (h ≤ −2), each horizon weighted by its number of
treated rows N₁ₕ: θ̂ = Σₕ N₁ₕ β̂ₕ / Σₕ N₁ₕ, the average effect over all treated observations
of the window (the equally weighted ATT with `reweight`, its variance-weighted counterpart
without). As `lpdid`'s `aggregate_average`: the standard error is that of the influence
function of θ̂, which includes the estimated weights,

    A_i = [Σₕ N₁ₕ ψ_ih + Σₕ (β̂ₕ − θ̂) 1{i is a treated row of horizon h}] / Σₕ N₁ₕ,

summed by `cluster`, with no small-sample factor and t(M − 1) critical values, M the
clusters in the samples of the window. A window with a horizon that was not estimated is
NaN. Returns a column table with rows `:post` and `:pre`.
"""
function averageeffect(r::PanelLocalProjection; cluster::Symbol = :unit,
        level::Real = 0.95)
    r.estimator isa LPDiD ||
        throw(ArgumentError("averageeffect needs an LP-DiD result (estimator = LPDiD(…))"))
    isempty(r.pooled) || throw(ArgumentError(
        "averageeffect averages event-study coefficients; this is a pooled estimate"))
    0 < level < 1 || throw(ArgumentError("level must be in (0, 1)"))
    g = _clusterdim(r, cluster)
    j = _termindex(r.coefnames, r.shocks[1])
    windows = (:post => [h for h in r.horizons if h >= 0],
        :pre => [h for h in r.horizons if h <= -2])
    out = (window = Symbol[], coef = Float64[], se = Float64[], lower = Float64[],
        upper = Float64[], df = Float64[], ntreated = Int[], horizons = UnitRange{Int}[])
    for (name, hs) in windows
        isempty(hs) && continue
        fits = [_fit(r, h) for h in hs]
        push!(out.window, name)
        push!(out.horizons, first(hs):last(hs))
        if !all(f -> f.status === :ok && f.keep[j] && f.did !== nothing, fits)
            foreach(v -> push!(v, NaN), (out.coef, out.se, out.lower, out.upper, out.df))
            push!(out.ntreated, 0)
            continue
        end
        n1 = [f.did.ntreated for f in fits]
        S = sum(n1)
        θ = sum(n1 .* [f.coef[j] for f in fits]) / S
        A = zeros(r.ngroups[g])
        for (f, n) in zip(fits, n1)
            A .+= n .* _influence(f, g, j) .+ (f.coef[j] - θ) .* f.did.treated[g]
        end
        se = sqrt(sum(abs2, A)) / S
        M = _clusters(r, fits, g)
        q = M > 1 ? StatsFuns.tdistinvcdf(M - 1, 1 - (1 - level) / 2) : NaN
        push!(out.coef, θ)
        push!(out.se, se)
        push!(out.lower, θ - q * se)
        push!(out.upper, θ + q * se)
        push!(out.df, M - 1)
        push!(out.ntreated, S)
    end
    return out
end
