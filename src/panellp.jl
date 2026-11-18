"""
    panellp(f, data; id, time, horizon, estimator = WithinFE(), shock = nothing,
            frequency = nothing, subset = nothing, sample = :horizon, singletons = :drop,
            absorb = AlternatingProjections(), clusters = (), weights = nothing,
            pooled = false)

Estimate a panel local projection: for each horizon `h`, regress the response defined by
the left-hand side of `f` on its right-hand side, absorbing the `fe()` terms.

# Formula
- Left-hand side: `leads(y)` (y_{i,t+h}), `ldiff(y)` (y_{i,t+h} − y_{i,t−1}),
  `anchor(y, z)` (y_{i,t+h} − z_{i,t}), `cumul(y)` (Σ_{k≤h} y_{i,t+k}), or a plain column
  (a pre-built outcome) with `horizon = 0`.
- Right-hand side: columns and functions of columns, interactions `s & X`, `lags(x, p)`,
  `lags(x, a:b)`, `leads(x, a:b)`, `firstdiff(x)`, and fixed effects `fe(v)`,
  `fe(a) & fe(b)`. An intercept is included when there are no fixed effects (unless
  `0 +`). `leads(x)` without a range adds the path `x_{i,t+1}, …, x_{i,t+h}` at horizon
  `h` (nothing at `h ≤ 0`): the response to `x_{i,t}` holding its later values fixed.
- Panel LP-IV: one block `(endo ~ inst)`, e.g. `(s & ffr ~ s & mp)`, makes the fit two-stage
  least squares (Almuzara & Sancibrián 2026, §3.4) on each horizon's fixed-effects
  design: the endogenous regressors on the left (listed first among the coefficients,
  and the default `shock`), the excluded instruments on the right; the other regressors
  are exogenous and instrument themselves. Every variance applies to the IV scores
  X̂ ξ̂ (ξ̂ the structural residual); `TLAHR(refine = false)` is the paper's LP-IV variance.
  An instrument absorbed by the fixed effects or collinear with the exogenous regressors
  (or with earlier instruments) is left out of the first stage, with a note; a horizon with
  fewer instruments than endogenous regressors left is `:underidentified`. With
  `WithinFE` or [`LagAugmented`](@ref) (whose lags of an endogenous shock are instrumented
  by lags of the instruments).

# Keywords
- `id`, `time`: the unit and period columns. Integer periods are consecutive integers;
  `Date`/`DateTime` periods need `frequency` (`:year`, `:halfyear`, `:quarter`, `:month`).
  `fe(time)` and clusters on `time` group by period, not by date.
- `horizon`: `H` (horizons `0:H`) or a range of consecutive horizons such as `1:6`.
  Negative horizons put the response before the regressors (pre-trends): with
  `horizon = -6:12`, `ldiff(y)` is `y_{i,t−|h|} − y_{i,t−1}` at `h < 0`, identically zero
  at `h = −1` (coefficient and standard error 0). Not with `cumul` or
  [`LagAugmented`](@ref).
- `estimator`: [`WithinFE`](@ref) (fixed effects; pooled OLS without them),
  [`SplitPanelJackknife`](@ref), [`CDOLS`](@ref), [`LagAugmented`](@ref) or
  [`LPDiD`](@ref) (difference-in-differences with clean controls).
- `shock`: coefficient name(s) of the shock; by default the columns of the first
  regressor (the endogenous ones with an IV block). Required when the formula has
  interaction terms and other regressors.
- `subset`: a Bool column or vector restricting the estimation sample *after* leads and
  lags are computed (they may read outside it).
- `sample`: `:horizon` (each horizon its own complete cases) or `:common` (the
  intersection over horizons).
- `singletons`: `:drop` (iteratively, as `reghdfe` and `xtlp`) or `:keep`.
- `absorb`: [`AlternatingProjections`](@ref) or [`OneSweep`](@ref).
- `clusters`: extra columns (or tuples of columns, for intersections) whose score sums are
  stored for clustered variances, e.g. `clusters = (:region, (:sector, :quarter))`; unit,
  time and the fixed-effect dimensions always are.
- `weights`: a column of observation weights (weighted least squares, as Stata's
  `[pweight]`/`[aweight]`: weighted fixed-effect means, scores `w x̃ ê`, and the number of
  rows in every degrees-of-freedom factor). Rows with a missing or zero weight leave the
  sample. With `WithinFE`, `CDOLS` and `LPDiD`.
- `pooled = true`: one regression of the response averaged over the horizons (every one
  required), instead of one per horizon: the average response over the window, on the rows
  where the whole window is observed (the pooled LP-DiD estimate of Dube, Girardi, Jordà &
  Taylor 2025, Section 3.5, with `LPDiD`). The horizons must be all ≥ 0 or all < 0; the
  result has one fit, at the far end of the window, and `r.pooled` holds the window. With
  `WithinFE` and `LPDiD`.

Leads, lags and differences are taken on the (unit, calendar period) grid: an absent
period makes them undefined, and rows where the response or any regressor is undefined
(missing, `NaN`, or a missing input of a function term) leave that horizon's sample.

Horizons without observations, horizons where the shock is absorbed by the fixed effects
or collinear in the horizon's sample (status `:singular`), omitted regressors and
non-convergence are reported in one warning each per call, and their coefficients are
NaN; it is an error when no horizon identifies the shock.
"""
function panellp(f::StatsModels.FormulaTerm, data; id::Symbol, time::Symbol, horizon,
        estimator::PanelLPEstimator = WithinFE(), shock = nothing, frequency = nothing,
        subset = nothing, sample::Symbol = :horizon, singletons::Symbol = :drop,
        absorb::AbsorptionMethod = AlternatingProjections(), clusters = (),
        weights = nothing, pooled::Bool = false)
    _check_option(sample, (:horizon, :common), "sample")
    _check_option(singletons, (:drop, :keep), "singletons")
    (id === :time || time === :unit || id === time) && throw(ArgumentError(
        "`:unit` and `:time` name the panel dimensions (e.g. in Clustered(:unit, :time)): " *
        "the id column cannot be called `time`, nor the time column `unit`"))
    horizons = _horizons(horizon)
    weights === nothing || _check_weights(estimator)
    pooled && _check_pooled(estimator, horizons)
    d = _build_design(f, data; id, time, frequency, subset, clusters, shock, weights,
        hmax = max(last(horizons), 0))
    isempty(d.endo) || _check_iv(estimator)
    _check_horizons(estimator, d, horizons)
    _validate(estimator, d, f)
    d.response isa PlainResponse && horizons != 0:0 &&
        throw(ArgumentError(
            "a plain response is a pre-built outcome and is estimated at horizon 0 only; use " *
            "leads(y), ldiff(y), anchor(y, z) or cumul(y) to estimate several horizons"))
    ctx = _EvalContext(d.pidx, d.columns)
    base = d.rhsok .& d.subset
    for codes in Iterators.flatten((d.fe, d.groups))
        base .&= codes .!= 0
    end
    # LagAugmented: p(h) and the lag columns; LPDiD: the treatment histories (the other
    # estimators pass through)
    d, est, lags = _augment(estimator, d, ctx, base, horizons, singletons, absorb)
    est = _withsample(est, sample)
    # the responses are built again when each horizon is fitted, not kept for all
    masks = [_sample(est, d, ctx, base, h) for h in horizons]
    if pooled
        fits = [_pooled_fit(est, d, ctx, horizons, masks, sample, singletons, absorb)]
        window = horizons
        horizons = fits[1].h:fits[1].h
    else
        if sample === :common
            common = _common_sample(masks, horizons)
            masks = [copy(common) for _ in masks]
        end
        fits = _estimate(est, d, ctx, horizons, masks, singletons, absorb)
        window = 1:0
    end
    estimator isa LPDiD && !any(f -> f.status === :ok, fits) &&
        !all(f -> f.status === :empty, fits) &&
        throw(ArgumentError(
            "LP-DiD: no horizon has both rows entering treatment and clean controls in " *
            "the same period (the treatment switch is constant in every sample). Check " *
            "the treatment column and the clean-control options: $(_label(estimator))"))
    _check_identified(fits, d, time)
    p = d.pidx
    panel = PanelSummary(
        id, time, p.frequency, nunits(p), nspan(p), length(unique(p.period)),
        nrows(p), p.pmin, p.pmax)
    r = PanelLocalProjection(estimator, absorb, f, d.responsename, d.coefnames,
        d.intercept, d.coefnames[d.shockcols], horizons, fits, d.groupnames, d.ngroups,
        d.fenames, panel, sample, singletons, lags, d.owner, d.lagorder,
        d.coefnames[d.endo], d.instnames, d.zlagorder, window,
        weights === nothing ? nothing : Symbol(weights))
    foreach(note -> @warn(note), _fitnotes(r))
    return r
end

function _horizons(horizon)
    if horizon isa Integer
        horizon >= 0 || throw(ArgumentError(
            "horizon = H means the horizons 0:H; for negative horizons give a range, e.g. " *
            "horizon = $horizon:0"))
        return 0:Int(horizon)
    elseif horizon isa AbstractVector{<:Integer} && !isempty(horizon)
        hs = collect(Int, horizon)
        a, b = first(hs), last(hs)
        hs == a:b && return a:b
    end
    throw(ArgumentError(
        "horizon must be an integer H ≥ 0 (horizons 0:H) or a range of consecutive " *
        "horizons a:b; got $(repr(horizon))"))
end

# Negative horizons (the response before the regressors) need a response defined there
# and an estimator whose rules do not depend on the horizon.
function _check_horizons(est::PanelLPEstimator, d::PanelDesign, horizons::UnitRange{Int})
    first(horizons) < 0 || return nothing
    d.response isa CumulResponse && throw(ArgumentError(
        "cumul(y) sums the response from t to t + h: it has no negative horizons; use " *
        "ldiff of the level, or horizons from 0"))
    est isa LagAugmented && throw(ArgumentError(
        "LagAugmented sets the lags p(h) for horizons h ≥ 0 (Almuzara & Sancibrián 2026): " *
        "estimate negative horizons with WithinFE"))
    return nothing
end

_check_weights(::Union{WithinFE, CDOLS, LPDiD}) = nothing
function _check_weights(est::PanelLPEstimator)
    throw(ArgumentError(
        "observation weights are implemented for WithinFE, CDOLS and LPDiD, not for " *
        "$(nameof(typeof(est)))"))
end

function _check_pooled(est::PanelLPEstimator, horizons::UnitRange{Int})
    est isa Union{WithinFE, LPDiD} || throw(ArgumentError(
        "pooled = true is implemented for WithinFE and LPDiD, not for " *
        "$(nameof(typeof(est)))"))
    (first(horizons) >= 0 || last(horizons) < 0) || throw(ArgumentError(
        "pooled = true averages the response over a window on one side of the shock: the " *
        "horizons must be all ≥ 0 or all < 0; got $horizons"))
    return nothing
end

# The response of horizon h for every data row.
function _response(::PanelLPEstimator, d::PanelDesign, ctx::_EvalContext, h::Int)
    _lhs_at(d.response, ctx, h)
end
function _response(s::_DiD, d::PanelDesign, ctx::_EvalContext, h::Int)
    return _col1(d.response.term, ctx, h) .- s.base
end

# Estimator-specific restrictions of the sample of horizon h.
function _estmask!(mask::BitVector, est::PanelLPEstimator, d::PanelDesign, h::Int)
    _lagmask!(mask, est, d, h)
end
function _estmask!(mask::BitVector, s::_DiD, d::PanelDesign, h::Int)
    mask .&= _cleanflag(s, d.pidx, h)
end

# The rows where the response and the regressors of horizon h are defined.
function _sample(est::PanelLPEstimator, d::PanelDesign, ctx::_EvalContext,
        base::BitVector, h::Int)
    mask = base .& isfinite.(_response(est, d, ctx, h))
    _leadmask!(mask, d, h)
    return _estmask!(mask, est, d, h)
end

# The estimator told whether the horizons share one sample (LP-DiD's reweighting counts
# rows of that sample).
_withsample(est::PanelLPEstimator, sample::Symbol) = est
function _withsample(s::_DiD, sample::Symbol)
    _DiD((getfield(s, k) for k in 1:(fieldcount(_DiD) - 1))..., sample === :common)
end

# Whether singleton groups of the fixed effects are dropped before the fit.
_dropsingletons(::PanelLPEstimator, d::PanelDesign) = true
_dropsingletons(s::_DiD, d::PanelDesign) = !_usera(s, d)

# The pooled fit: the response averaged over the horizons, on the rows where every one is
# defined, with the regressors (and the clean-control rule) of the far end of the window;
# `sample = :common` also requires the samples of every horizon of the window.
function _pooled_fit(est::PanelLPEstimator, d::PanelDesign, ctx::_EvalContext,
        horizons::UnitRange{Int}, masks::Vector{BitVector}, sample::Symbol,
        singletons::Symbol, absorb::AbsorptionMethod)
    far = first(horizons) >= 0 ? last(horizons) : first(horizons)
    y = zeros(nrows(d.pidx))
    for h in horizons
        y .+= _response(est, d, ctx, h)
    end
    y ./= length(horizons)
    mask = copy(masks[far - first(horizons) + 1])
    if sample === :common
        foreach(m -> mask .&= m, masks)
    end
    mask .&= isfinite.(y)
    dropped = singletons === :drop && _dropsingletons(est, d) ?
              _drop_singletons!(mask, d.fe) : 0
    return _fit_horizon(est, d, far, y, mask, dropped, absorb)
end

# Intersection of the per-horizon samples; an error naming the horizon that empties it.
function _common_sample(masks::Vector{BitVector}, horizons::UnitRange{Int})
    common = copy(masks[1])
    for (k, m) in enumerate(masks)
        common .&= m
        any(common) && continue
        what = any(m) ?
               "the intersection of the samples of horizons " *
               "$(first(horizons))–$(horizons[k]) is empty" :
               "horizon $(horizons[k]) has no observations"
        throw(ArgumentError(
            "sample = :common: $what, so no row has the response and every regressor " *
            "defined at all horizons. Use fewer horizons or sample = :horizon."))
    end
    return common
end

# Function barrier: concrete types from here on.
function _estimate(est::PanelLPEstimator, d::PanelDesign, ctx::_EvalContext,
        horizons::UnitRange{Int}, masks::Vector{BitVector}, singletons::Symbol,
        absorb::AbsorptionMethod)
    fits = Vector{HorizonFit}(undef, length(horizons))
    for (k, h) in enumerate(horizons)
        mask = copy(masks[k])
        dropped = singletons === :drop && _dropsingletons(est, d) ?
                  _drop_singletons!(mask, d.fe) : 0
        fits[k] = _fit_horizon(est, d, h, _response(est, d, ctx, h), mask, dropped,
            absorb)
    end
    return fits
end

# An error when no horizon has observations, or when the shock is identified at none.
function _check_identified(fits::Vector{HorizonFit}, d::PanelDesign, time::Symbol)
    all(f -> f.status === :empty, fits) && throw(ArgumentError(
        "no horizon has observations: at every row the response or a regressor is " *
        "undefined, or only singletons are left (check `horizon`, `subset`, the lags and " *
        "missing values)"))
    any(f -> f.status === :ok, fits) && return nothing
    any(f -> f.status === :spj_half_singular, fits) && throw(ArgumentError(
        "no horizon can be estimated: at every horizon the shock is not identified or a " *
        "split-panel half is rank-deficient ($(_halfcause(fits, d.coefnames)))"))
    any(f -> f.status === :singular, fits) || throw(ArgumentError(
        "the model is underidentified at every horizon: once the fixed effects and the " *
        "exogenous regressors are partialled out, fewer excluded instruments than " *
        "endogenous regressors are left, or the instruments do not predict every " *
        "endogenous regressor (instruments $(join(repr.(d.instnames), ", ")))"))
    f = fits[findfirst(f -> f.status === :singular, fits)]
    absorbed = [d.coefnames[j] for j in d.shockcols if f.absorbed[j]]
    if isempty(absorbed)
        throw(ArgumentError(
            "the shock $(join(repr.(d.coefnames[d.shockcols]), ", ")) is collinear with the " *
            "other regressors at every horizon"))
    end
    timefe = any(nm -> string(time) in split(nm, '&'), d.fenames)
    throw(ArgumentError(
        "the shock $(join(repr.(absorbed), ", ")) is absorbed by the " *
        "$(isempty(d.fe) ? "constant" : "fixed effects") at every horizon" *
        (timefe ?
         " (an aggregate shock with time effects? Interact it with a unit-level " *
         "exposure, or drop the time effects)" : "")))
end

# Why split-panel halves are rank-deficient: the regressors the fixed effects absorb inside
# a half (marked in the fits' `absorbed`), or else too few rows.
function _halfcause(fits::Vector{HorizonFit}, coefnames::Vector{String})
    hf = [f for f in fits if f.status === :spj_half_singular]
    names = unique(coefnames[j] for f in hf for j in findall(f.absorbed .& f.keep))
    isempty(names) &&
        return "fewer rows than regressors, or collinear regressors, inside a half: too " *
               "few periods per unit?"
    # every regressor absorbed: units with a single row in a half
    allabs = any(f -> any(f.keep) && all(f.absorbed[f.keep]), hf)
    return "$(join(names, ", ")) absorbed by the fixed effects inside a half" *
           (allabs ? ": too few periods per unit?" : "")
end

# "0–3, 5, 7–9"
function _hlist(hs::Vector{Int})
    parts = String[]
    k = 1
    while k <= length(hs)
        j = k
        while j < length(hs) && hs[j + 1] == hs[j] + 1
            j += 1
        end
        push!(parts,
            j == k ? string(hs[k]) :
            j == k + 1 ? "$(hs[k]), $(hs[j])" :
            hs[k] < 0 ? "$(hs[k]) to $(hs[j])" : "$(hs[k])–$(hs[j])")
        k = j + 1
    end
    return join(parts, ", ")
end

# What went wrong at some horizons: reported by `panellp` (one warning each) and `show`.
function _fitnotes(r::PanelLocalProjection)
    notes = String[]
    hs(pred) = [f.h for f in r.fits if pred(f)]
    estimated = hs(f -> f.status === :ok)
    h = hs(f -> f.status === :empty)
    isempty(h) || push!(notes,
        "no observations at h = $(_hlist(h)) (the response or a regressor is undefined " *
        "at every row, or only singletons are left); coefficients are NaN")
    h = hs(f -> f.status === :singular)
    isempty(h) || push!(notes,
        "the shock is absorbed by the fixed effects or collinear in the sample of " *
        "h = $(_hlist(h)); coefficients are NaN")
    h = hs(f -> f.status === :underidentified)
    isempty(h) || push!(notes,
        "underidentified at h = $(_hlist(h)) (once the fixed effects and the exogenous " *
        "regressors are partialled out, fewer excluded instruments than endogenous " *
        "regressors are left, or the instruments do not predict every endogenous " *
        "regressor); coefficients are NaN")
    h = hs(f -> f.status === :spj_half_singular)
    isempty(h) || push!(notes,
        "a split-panel half is rank-deficient at h = $(_hlist(h)) " *
        "($(_halfcause(r.fits, r.coefnames))); coefficients are NaN")
    absorbed = (f, j) -> f.absorbed[j]
    # a lag LagAugmented added that the other regressors span is redundant, not omitted
    collinear = (f, j) -> f.active[j] && !f.keep[j] && !f.absorbed[j] && r.addedlag[j] == 0
    for (label, pred) in (("absorbed by the fixed effects or the constant", absorbed),
        ("collinear with the other regressors", collinear))
        items = String[]
        for (j, nm) in enumerate(r.coefnames)
            h = hs(f -> f.status === :ok && pred(f, j))
            isempty(h) && continue
            push!(items, h == estimated ? nm : "$nm (h = $(_hlist(h)))")
        end
        isempty(items) || push!(notes, "omitted, $label: $(join(items, "; "))")
    end
    items = String[]
    for (l, nm) in enumerate(r.instruments)
        # an added lag beyond p(h) is not in the fit, not left out
        h = hs(f -> f.status === :ok && r.instlag[l] <= f.nlags && !f.instkeep[l])
        isempty(h) && continue
        push!(items, h == estimated ? nm : "$nm (h = $(_hlist(h)))")
    end
    isempty(items) || push!(notes,
        "instruments left out of the first stage (absorbed by the fixed effects, or " *
        "collinear with the exogenous regressors or the other instruments): " *
        join(items, "; "))
    h = hs(f -> !f.converged)
    isempty(h) || push!(notes,
        "fixed-effect absorption did not converge at h = $(_hlist(h)): the estimates there " *
        "are not at the requested accuracy, and a regressor that the fixed effects absorb " *
        "may go undetected; raise `maxiter` in AlternatingProjections")
    return notes
end
