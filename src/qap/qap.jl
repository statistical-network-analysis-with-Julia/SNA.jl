"""
QAP inference: permutation tests and network regression.

Provides graph correlation (`gcor`, `gcov`), the quadratic assignment
procedure (QAP) test for arbitrary graph-level statistics, and OLS/logistic
network regression with QAP null-hypothesis testing (Dekker's
double-semi-partialing by default), following R `sna::gcor`, `sna::qaptest`,
`sna::netlm`, and `sna::netlogit`.
"""

using Distributions: TDist, Normal, cdf, ccdf, quantile
using LinearAlgebra
using Logging: NullLogger, with_logger
using Random
using Statistics

# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

# Apply the missing-dyad policy to every dyadic argument of a QAP routine.
# `Network` arguments are guarded with `require_observed`; raw matrices carry
# no mask, so nothing can be masked in them — but the policy itself is still
# validated so that a typo (`missing=:faec`) never passes silently.
function _require_observed_dyads(policy::Symbol, context::AbstractString,
                                 args...)
    policy in MISSING_POLICIES ||
        throw(ArgumentError("invalid missing-dyad policy $(repr(policy)); " *
                            "expected one of " *
                            "$(join(map(repr, MISSING_POLICIES), ", "))"))
    for a in args
        a isa AbstractNetwork && require_observed(a, policy; context=context)
    end
    return nothing
end

# What the missing-dyad policy actually AMOUNTED TO for these arguments, for the
# shared result-metadata protocol. Call only after `_require_observed_dyads` has
# passed, so `:error` + a mask has already thrown: what remains is either no mask
# at all (`:none`) or an explicit opt-in to face values (`:condition_on_face`).
function _missing_method_of(policy::Symbol, args...)
    masked = any(a -> a isa AbstractNetwork && n_missing_dyads(a) > 0, args)
    masked || return :none
    return :condition_on_face      # policy === :face, guard already let it through
end

# Apply a random vertex permutation to rows and columns (sna::rmperm). With a
# two-mode partition p, actors 1:p and p+1:n are permuted separately, so a
# two-mode network keeps its modes.
function _rmperm(rng::Random.AbstractRNG, A::Matrix{Float64}, partition::Union{Nothing,Int}=nothing)
    n = size(A, 1)
    p = partition === nothing ? randperm(rng, n) :
        vcat(randperm(rng, partition), partition .+ randperm(rng, n - partition))
    return A[p, p]
end

_is_rectangular(x) = x isa AbstractMatrix && size(x, 1) != size(x, 2)

# The two-mode partition (number of mode-1 actors) of the arguments, or
# `nothing` for one-mode data. A rectangular matrix is a two-mode incidence
# matrix whose rows are mode 1.
function _qap_partition(args...)
    partitions = Int[]
    for a in args
        if a isa BipartiteNetwork
            push!(partitions, a.n_mode1)
        elseif a isa Network && a.bipartite !== nothing
            push!(partitions, a.bipartite)
        elseif _is_rectangular(a)
            push!(partitions, size(a, 1))
        end
    end
    isempty(partitions) && return nothing
    all(==(first(partitions)), partitions) ||
        throw(ArgumentError("two-mode QAP arguments must have the same actor partition"))
    return first(partitions)
end

# Square sociomatrix of any QAP argument: a network's (two-mode networks are
# expanded to all actors), a square matrix as given, and a p×q incidence
# matrix embedded symmetrically in the (p+q)×(p+q) sociomatrix.
function _qap_sociomatrix(x)
    _is_rectangular(x) || return _sociomatrix(x)
    p, q = size(x)
    B = Matrix{Float64}(x)
    all(isfinite, B) || throw(ArgumentError("an incidence matrix must contain finite values"))
    A = zeros(p + q, p + q)
    A[1:p, (p+1):end] = B
    A[(p+1):end, 1:p] = B'
    return A
end

# Whether the dyadic data are directed: from `mode`, else from `y` when it is
# a network; a rectangular incidence matrix is undirected two-mode data, and
# square raw matrices are undirected only when `y` and every predictor are
# symmetric (an asymmetric predictor of a symmetric response is directed data).
function _qap_directed(y, xs, mode::Symbol)
    mode in (:auto, :digraph, :graph) ||
        throw(ArgumentError("mode must be :auto, :digraph, or :graph"))
    mode == :digraph && return true
    mode == :graph && return false
    y isa AbstractNetwork && return is_directed(y)
    _is_rectangular(y) && return false
    _symmetric_dyads(x) = x isa AbstractNetwork ? !is_directed(x) :
                          _is_rectangular(x) || issymmetric(x)
    return !(issymmetric(y) && all(_symmetric_dyads, xs))
end

# The dyads a QAP analysis observes, as (row, column) index pairs into the
# square sociomatrix, in sna's column-major order: every off-diagonal ordered
# pair for directed data, the lower triangle (i > j) for undirected data
# (sna::gvectorize). With a two-mode partition p only cross-mode pairs are
# dyads; the within-mode pairs are structurally impossible, and counting them
# as observed zeros would manufacture association. For directed two-mode data
# the mode-2 → mode-1 pairs are dyads only when some argument has a tie there
# (`reverse=true`); otherwise that block is structurally empty as well.
function _dyad_indices(n::Int, directed::Bool, partition::Union{Nothing,Int}=nothing;
                       reverse::Bool=true)
    idx = Tuple{Int,Int}[]
    for j in 1:n
        i_range = directed ? (1:n) : ((j+1):n)
        for i in i_range
            i == j && continue
            if partition !== nothing
                forward = i <= partition < j
                backward = j <= partition < i
                if directed
                    (forward || (reverse && backward)) || continue
                else
                    backward || continue           # lower triangle: i in mode 2, j in mode 1
                end
            end
            push!(idx, (i, j))
        end
    end
    return idx
end

# Whether any sociomatrix has a tie from a mode-2 to a mode-1 actor.
_has_reverse_ties(partition::Int, Ms) =
    any(M -> any(!iszero, @view M[(partition+1):end, 1:partition]), Ms)

# Whether some sociomatrix's mode-2 → mode-1 block says something its
# mode-1 → mode-2 block does not (it is neither empty nor the transpose).
_reverse_informative(partition::Int, Ms) = any(Ms) do M
    R = @view M[(partition+1):end, 1:partition]
    any(!iszero, R) && R != transpose(@view M[1:partition, (partition+1):end])
end

# Vectorize the selected dyads of a sociomatrix
_gvectorize(A::Matrix{Float64}, idx::Vector{Tuple{Int,Int}}) =
    Float64[A[i, j] for (i, j) in idx]

# Each replicate owns its RNG and workspace; scheduling never changes the draws.
function _qap_replicates(f, rng::AbstractRNG, n_sim::Int, threaded::Bool)
    seeds = rand(rng, UInt64, n_sim)
    errors = Vector{Any}(nothing, n_sim)
    function replicate(r)
        try
            f(r, Xoshiro(seeds[r]))
        catch e   # a fresh local: never the outer variable, so never a shared Core.Box
            e isa InterruptException && rethrow()
            errors[r] = e
        end
    end
    if threaded && Threads.nthreads() > 1
        Threads.@threads for r in 1:n_sim
            replicate(r)
        end
    else
        for r in 1:n_sim
            replicate(r)
            errors[r] === nothing || break
        end
    end
    failed = findfirst(!isnothing, errors)
    if failed !== nothing
        failure = errors[failed]
        if failure isa ArgumentError
            throw(ArgumentError("QAP replicate $failed failed: $(failure.msg). " *
                                "No permutation p-value was computed; every requested replicate is required."))
        end
        throw(failure)
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Graph correlation and covariance
# ---------------------------------------------------------------------------

# The paired dyad vectors of two graphs, as sna::gcor/gcov select them.
function _paired_dyads(g1, g2, mode::Symbol, diag::Bool, fname::String, policy::Symbol)
    _require_observed_dyads(policy, fname, g1, g2)
    partition = _qap_partition(g1, g2)
    A = _qap_sociomatrix(g1)
    B = _qap_sociomatrix(g2)
    size(A) == size(B) ||
        throw(ArgumentError("$fname: the two graphs must have the same number of vertices"))
    directed = _qap_directed(g1, (g2,), mode)
    n = size(A, 1)
    idx = partition === nothing ? _dyad_indices(n, directed) :
          _dyad_indices(n, directed, partition; reverse=_has_reverse_ties(partition, (A, B)))
    if diag && partition === nothing
        idx = vcat(idx, [(i, i) for i in 1:n])
    end
    return _gvectorize(A, idx), _gvectorize(B, idx)
end

"""
    gcor(g1, g2; mode=:auto, diag=false, missing=:error) -> Float64

Graph correlation, R `sna::gcor`: the Pearson correlation of the two
networks' tie values over the dyads they share. As in sna the diagonal is
left out unless `diag=true`; an undirected pair counts once
(`mode=:graph`), an ordered pair of a directed network once in each
direction (`mode=:digraph`). `mode=:auto` takes undirected networks and
symmetric matrices as undirected data.

The graphs are `Network`s or adjacency matrices. For two-mode data (two-mode
networks, or `p × q` incidence matrices) only cross-mode pairs are compared.

`missing=:error` refuses masked dyads; `missing=:face` uses their stored
face values.

# Example
```julia
using SNA
flo = load_dataset(:florentine_marriage)
biz = load_dataset(:florentine_business)
gcor(flo, biz)                       # 0.3719, as R sna::gcor
```
"""
function gcor(g1, g2; mode::Symbol=:auto, diag::Bool=false, missing::Symbol=:error)
    x, y = _paired_dyads(g1, g2, mode, diag, "gcor", missing)
    return cor(x, y)
end

"""
    gcov(g1, g2; mode=:auto, diag=false, missing=:error) -> Float64

Graph covariance, R `sna::gcov`: the covariance of the two networks' tie
values over the dyads they share, selected as in [`gcor`](@ref).

# Example
```julia
using SNA
flo = load_dataset(:florentine_marriage)
biz = load_dataset(:florentine_business)
gcov(flo, biz)
```
"""
function gcov(g1, g2; mode::Symbol=:auto, diag::Bool=false, missing::Symbol=:error)
    x, y = _paired_dyads(g1, g2, mode, diag, "gcov", missing)
    return cov(x, y)
end

# ---------------------------------------------------------------------------
# qaptest
# ---------------------------------------------------------------------------

"""
    QAPTestResult

Result of a [`qaptest`](@ref) permutation test, with the fields of R's
`qaptest` object.

# Fields
- `testval::Float64`: Observed value of the test statistic
- `dist::Vector{Float64}`: Monte Carlo null distribution (one value per
  replicate)
- `pgreq::Float64`: Proportion of null draws `>=` the observed value
- `pleeq::Float64`: Proportion of null draws `<=` the observed value
- `reps::Int`: Number of replicates (`n_sim`)

# Example
```julia
using SNA, Random
flo = load_dataset(:florentine_marriage)
biz = load_dataset(:florentine_business)
qt = qaptest(gcor, flo, biz; n_sim=200, rng=Xoshiro(1))
qt isa QAPTestResult, qt.pgreq
```
"""
struct QAPTestResult
    testval::Float64
    dist::Vector{Float64}
    pgreq::Float64
    pleeq::Float64
    reps::Int
end

function Base.show(io::IO, result::QAPTestResult)
    println(io, "QAP Test Results")
    println(io, "================")
    println(io, "Replications: $(result.reps)")
    println(io)
    println(io, "Test value:            $(round(result.testval, digits=6))")
    println(io, "Estimated p-values:")
    println(io, "  p(f(perm) >= f(d)):  $(NetworkCore.format_pvalue(result.pgreq; floor=1/result.reps))")
    print(io, "  p(f(perm) <= f(d)):  $(NetworkCore.format_pvalue(result.pleeq; floor=1/result.reps))")
end

"""
    qaptest(f, g1, g2; n_sim=1000, threaded=true, missing=:error,
            rng=Random.default_rng()) -> QAPTestResult

Quadratic assignment procedure (QAP) test of the graph-level statistic `f`
for the network pair `(g1, g2)`, R `sna::qaptest`.

`f(A1, A2)` is computed on the adjacency matrices of the two networks and
must return a scalar, for instance [`gcor`](@ref). The observed value is
compared with its distribution under `n_sim` independent uniform vertex
permutations of `g1` (rows and columns together) with `g2` held fixed: the
QAP null hypothesis of no association between the two structures,
conditional on both. The p-values are the proportions of permuted values at
least (`pgreq`) or at most (`pleeq`) the observed one, as in sna, so they can
be exactly 0.

`g1` and `g2` are `Network`s or adjacency matrices of equal size.

**Two-mode data.** For two-mode networks, or `p × q` incidence matrices, `f`
is called on the `p × q` incidence matrices (rows: mode 1, columns: mode 2),
and each permutation reorders the mode-1 actors and the mode-2 actors
separately. A directed two-mode network with ties from mode 2 to mode 1 is
refused; pass its square sociomatrix as a raw matrix to test a statistic of
the full matrix.

`missing=:error` refuses a network argument with masked (unobserved) dyads,
because a masked dyad read at face value would enter both the statistic and
the null distribution; `missing=:face` computes `f` from the stored face
values. Raw matrices carry no mask.

Replicates use independently seeded random number generators drawn from
`rng`, so results are reproducible and do not depend on the number of
threads; `threaded=false` runs them serially. `f` must be thread-safe when
threading is on.

# Example
```julia
using SNA, Random
flo = load_dataset(:florentine_marriage)
biz = load_dataset(:florentine_business)
qt = qaptest(gcor, flo, biz; n_sim=1000, rng=Xoshiro(1))
qt.testval, qt.pgreq                 # 0.3719, about 0.001
```
"""
function qaptest(f, g1, g2; n_sim::Int=1000, threaded::Bool=true, missing::Symbol=:error,
                 rng::Random.AbstractRNG=Random.default_rng())
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    _require_observed_dyads(policy, "qaptest", g1, g2)
    n_sim > 0 || throw(ArgumentError("n_sim must be positive"))
    partition = _qap_partition(g1, g2)
    A = _qap_sociomatrix(g1)
    B = _qap_sociomatrix(g2)
    size(A) == size(B) ||
        throw(ArgumentError("g1 and g2 must have the same number of vertices"))
    view_of = if partition === nothing
        identity
    else
        _reverse_informative(partition, (A, B)) &&
            throw(ArgumentError("qaptest: a directed two-mode network with ties from mode 2 " *
                                "to mode 1 has no single incidence matrix; pass the square " *
                                "sociomatrices as raw matrices to test a statistic of both blocks"))
        M -> M[1:partition, (partition+1):end]
    end

    testval = Float64(f(view_of(A), view_of(B)))
    Bv = view_of(B)
    dist = Vector{Float64}(undef, n_sim)
    _qap_replicates(rng, n_sim, threaded) do r, rrng
        dist[r] = Float64(f(view_of(_rmperm(rrng, A, partition)), Bv))
    end

    pgreq = count(>=(testval), dist) / n_sim
    pleeq = count(<=(testval), dist) / n_sim

    return QAPTestResult(testval, dist, pgreq, pleeq, n_sim)
end

# ---------------------------------------------------------------------------
# netlm
# ---------------------------------------------------------------------------

"""
    NetLMResult

Result of a [`netlm`](@ref) network regression.

# Fields
- `coefficients::Vector{Float64}`: OLS coefficient estimates
- `names::Vector{String}`: Coefficient names
- `tstat::Vector{Float64}`: t-statistics
- `covariance::Matrix{Float64}`: iid-dyad OLS covariance
- `loglik::Float64`: maximized iid Gaussian log likelihood
- `pleeq::Vector{Float64}`: `p(stat_perm <= stat_obs)` per coefficient
- `pgreq::Vector{Float64}`: `p(stat_perm >= stat_obs)` per coefficient
- `pgreqabs::Vector{Float64}`: Two-sided `p(|stat_perm| >= |stat_obs|)`
  (the p-value reported by `show`)
- `dist::Union{Nothing,Matrix{Float64}}`: Null distribution of the test
  statistics (`n_sim × k`; `nothing` for `nullhyp = :classical`)
- `r_squared::Float64`: Coefficient of determination
- `adj_r_squared::Float64`: Adjusted R²
- `n::Int`: Number of dyadic observations
- `df_residual::Int`: Residual degrees of freedom
- `nullhyp::Symbol`: Null hypothesis actually used
- `reps::Int`: Number of Monte Carlo replicates (0 for `:classical`)
- `intercept::Bool`: Whether an intercept column was included
- `directed::Bool`: Whether the dyads were treated as directed

# Example
```julia
using SNA, Random
flo = load_dataset(:florentine_marriage)
biz = load_dataset(:florentine_business)
fit = netlm(flo, biz; n_sim=200, rng=Xoshiro(1))
fit isa NetLMResult, coef(fit)
```
"""
struct NetLMResult
    coefficients::Vector{Float64}
    names::Vector{String}
    tstat::Vector{Float64}
    covariance::Matrix{Float64}
    loglik::Float64
    pleeq::Vector{Float64}
    pgreq::Vector{Float64}
    pgreqabs::Vector{Float64}
    dist::Union{Nothing,Matrix{Float64}}
    r_squared::Float64
    adj_r_squared::Float64
    n::Int
    df_residual::Int
    nullhyp::Symbol
    reps::Int
    intercept::Bool
    directed::Bool
    # The missing-dyad policy the fit ran under, so the result can report it
    # through the shared metadata protocol (`missing_method`) instead of
    # answering ":unspecified". `:none` when no argument carried a mask,
    # `:condition_on_face` when the caller opted into face values.
    missing_method::Symbol
end

function Base.show(io::IO, result::NetLMResult)
    println(io, "Network Regression (OLS + QAP)")
    println(io, "==============================")
    println(io, "Null hypothesis: $(result.nullhyp)" *
                (result.dist === nothing ? "" : " ($(result.reps) replications)"))
    println(io, "Dyadic observations: $(result.n) " *
                "($(result.directed ? "directed" : "undirected") dyads)")
    println(io)
    show(io, coeftable(result))
    println(io)
    print(io, "Multiple R-squared: $(round(result.r_squared, digits=4)), " *
              "Adjusted R-squared: $(round(result.adj_r_squared, digits=4))")
end

# t-values of the OLS regression of y on X
function _ols_tvals(X::Matrix{Float64}, y::Vector{Float64})
    coef = X \ y
    resid = y - X * coef
    rdf = size(X, 1) - size(X, 2)
    resvar = sum(abs2, resid) / rdf
    se = sqrt.(max.(diag(inv(Symmetric(X' * X))), 0.0) .* resvar)
    return coef ./ se
end

# Resolve the requested null hypothesis (mirrors sna: :qap is an alias for
# :qapspp, and semi-partialing degenerates to y-permutation with a single
# regressor)
function _resolve_nullhyp(nullhyp::Symbol, nx::Int, fname::String)
    nullhyp == :qap && (nullhyp = :qapspp)
    nullhyp in (:qapspp, :qapy, :qapx, :classical) ||
        throw(ArgumentError("Unknown nullhyp for $fname: $nullhyp (use " *
                            ":qap/:qapspp, :qapy, :qapx, or :classical)"))
    (nullhyp == :qapspp && nx == 1) && (nullhyp = :qapy)
    return nullhyp
end

# Build the dyadic design: returns (yv, X, Y, Gx, idx, names, directed, partition).
# `policy` is the missing-dyad policy applied to `y` and every predictor.
function _qap_design(y, xs, intercept::Bool, mode::Symbol, fname::String,
                     policy::Symbol)
    _require_observed_dyads(policy, fname, y, xs...)
    isempty(xs) && !intercept &&
        throw(ArgumentError("$fname needs at least one predictor"))
    partition = _qap_partition(y, xs...)
    directed = _qap_directed(y, xs, mode)
    Y = _qap_sociomatrix(y)
    n = size(Y, 1)

    Gx = Matrix{Float64}[]
    names = String[]
    if intercept
        push!(Gx, ones(n, n))
        push!(names, "(intercept)")
    end
    for (k, x) in enumerate(xs)
        M = _qap_sociomatrix(x)
        size(M) == (n, n) ||
            throw(ArgumentError("predictor x$k must have the same number " *
                                "of vertices as y"))
        push!(Gx, M)
        push!(names, "x$k")
    end

    idx = partition === nothing ? _dyad_indices(n, directed) :
          _dyad_indices(n, directed, partition;
                        reverse=_has_reverse_ties(partition, (Y, Gx[(intercept ? 2 : 1):end]...)))
    yv = _gvectorize(Y, idx)
    X = Matrix{Float64}(undef, length(idx), length(Gx))
    for (j, M) in enumerate(Gx)
        X[:, j] = _gvectorize(M, idx)
    end
    return yv, X, Y, Gx, idx, names, directed, partition
end

# Dekker double-semi-partialing residual matrix for predictor column i:
# the residual of x_i on the remaining predictors, written back into dyad
# positions (and symmetrized for undirected data)
function _dsp_residual_matrix(X::Matrix{Float64}, Gx::Vector{Matrix{Float64}},
                              idx::Vector{Tuple{Int,Int}}, i::Int,
                              others::Vector{Int}, directed::Bool)
    Xo = X[:, others]
    ei = X[:, i] - Xo * (Xo \ X[:, i])
    E = copy(Gx[i])
    for (k, (r, c)) in enumerate(idx)
        E[r, c] = ei[k]
    end
    if !directed
        for (r, c) in idx
            E[c, r] = E[r, c]
        end
    end
    return E
end

# Permutation p-values from an n_sim × k null distribution matrix
function _perm_pvalues(dist::Matrix{Float64}, tstat::Vector{Float64})
    n_sim = size(dist, 1)
    pleeq = [count(<=(tstat[i]), @view dist[:, i]) / n_sim for i in eachindex(tstat)]
    pgreq = [count(>=(tstat[i]), @view dist[:, i]) / n_sim for i in eachindex(tstat)]
    pgreqabs = [count(x -> abs(x) >= abs(tstat[i]), @view dist[:, i]) / n_sim
                for i in eachindex(tstat)]
    return pleeq, pgreq, pgreqabs
end

"""
    netlm(y, xs; intercept=true, nullhyp=:qapspp, n_sim=1000, mode=:auto,
          threaded=true, missing=:error, rng=Random.default_rng()) -> NetLMResult

Linear (OLS) regression of the network `y` on one or more predictor
networks `xs`, with QAP null-hypothesis testing, R `sna::netlm`.

The networks are vectorized over their dyads (all ordered pairs for
directed data, one entry per unordered pair for undirected data; the
diagonal is always excluded) and `y` is regressed on the predictors by OLS.
Coefficient significance is assessed by comparing the observed
t-statistics with a permutation null distribution.

# Arguments
- `y`: Dependent network (`Network` or adjacency matrix)
- `xs`: Predictor network(s): a single network/matrix, or a vector/tuple
  of them, all of the same order as `y`
- `intercept::Bool=true`: Include an intercept column
- `nullhyp::Symbol=:qapspp`: Null hypothesis for the coefficient tests
    - `:qap`/`:qapspp`: Dekker's double-semi-partialing QAP (default, as in
      modern R `sna`): each predictor is residualized on the remaining
      predictors, the residual matrix is repeatedly row/column-permuted,
      and the t-statistic of the permuted residual (refit alongside the
      other predictors) forms the null distribution. Robust to
      multicollinearity among the predictors. With a single regressor this
      degenerates to `:qapy`.
    - `:qapy`: Classical y-permutation QAP: permute the rows/columns of `y`
      and refit
    - `:qapx`: Permute each predictor matrix separately and refit
    - `:classical`: Parametric t-tests (no permutation), for reference
      only: dyadic dependence typically makes them badly anti-conservative
- `n_sim::Int=1000`: Number of permutation replicates (R's `reps`)
- `mode::Symbol=:auto`: Dyad set: `:digraph` (every ordered pair), `:graph`
  (every unordered pair), or `:auto`, which takes the directedness of `y` when
  it is a network and, for raw matrices, treats the data as undirected when
  `y` and every predictor are symmetric. (sna's default is `"digraph"` for
  any matrix, which counts each pair of a symmetric matrix twice and makes
  the classical standard errors √2 too small.)
- `missing::Symbol=:error`: Missing-dyad policy for `y` and the predictors
  (`NetworkCore.require_observed`). A masked (unobserved) dyad read at face
  value becomes a fabricated row of the design matrix, biasing both the
  coefficients and the QAP null, so the default `:error` rejects a masked
  network argument. `:face` is the explicit opt-in to regressing on the
  stored face values; `netlm` does **not** implement listwise deletion of
  unobserved dyads or any missing-data estimator
  (`NetworkCore.supports_missing(netlm) == false`)
- `rng`: Random number generator

**Two-mode data.** For two-mode networks (or `p × q` incidence matrices)
the observations are the cross-mode dyads only: within-mode pairs cannot
carry a tie, and treating them as observed zeros would add a spurious
association (and can flip the sign of a slope). Permutations reorder each
mode separately. For a directed two-mode network the mode-2 → mode-1 pairs
are included only when some argument has a tie there.

The two-sided permutation p-value `pgreqabs` is reported by `show`;
one-sided `pleeq`/`pgreq` are also stored. As in sna, a p-value is the
proportion of permuted statistics at least as extreme as the observed one
(no `+1` correction), so it can be exactly 0; `show` prints it as `< 1/n_sim`.

Replicates use independently seeded random number generators drawn from
`rng` and are reproducible across thread counts; `threaded=false` disables
threading.

# Example
```julia
using SNA, Random
flo = load_dataset(:florentine_marriage)
biz = load_dataset(:florentine_business)
fit = netlm(flo, biz; n_sim=500, rng=Xoshiro(1))
```
"""
function netlm(y, xs::Union{Tuple,AbstractVector}; intercept::Bool=true,
               nullhyp::Symbol=:qapspp, n_sim::Int=1000,
               mode::Symbol=:auto,
               threaded::Bool=true,
               missing::Symbol=:error,
               rng::Random.AbstractRNG=Random.default_rng())
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    yv, X, Y, Gx, idx, names, directed, partition =
        _qap_design(y, xs, intercept, mode, "netlm", policy)
    miss = _missing_method_of(policy, y, xs...)
    N, nx = size(X)
    N > nx || throw(ArgumentError("more predictors than dyadic observations"))
    rank(X) == nx || throw(ArgumentError("predictors are linearly dependent; remove redundant columns"))

    coef = X \ yv
    resid = yv - X * coef
    df_residual = N - nx
    resvar = sum(abs2, resid) / df_residual
    covariance = Matrix(inv(Symmetric(X' * X))) * resvar
    se = sqrt.(diag(covariance))
    loglik = -N / 2 * (log(2π) + 1 + log(sum(abs2, resid) / N))
    tstat = coef ./ se

    sse = sum(abs2, resid)
    sst = intercept ? sum(abs2, yv .- mean(yv)) : sum(abs2, yv)
    r_squared = sst > 0 ? 1 - sse / sst : 0.0
    adj_r_squared = 1 - (1 - r_squared) * (N - (intercept ? 1 : 0)) / df_residual

    nullhyp = _resolve_nullhyp(nullhyp, nx, "netlm")

    if nullhyp == :classical
        tdist = TDist(df_residual)
        pleeq = cdf.(tdist, tstat)
        pgreq = ccdf.(tdist, tstat)
        pgreqabs = 2 .* ccdf.(tdist, abs.(tstat))
        return NetLMResult(coef, names, tstat, covariance, loglik, pleeq, pgreq, pgreqabs,
                           nothing, r_squared, adj_r_squared, N, df_residual,
                           nullhyp, 0, intercept, directed, miss)
    end

    n_sim > 0 || throw(ArgumentError("n_sim must be positive"))
    dist = Matrix{Float64}(undef, n_sim, nx)
    if nullhyp == :qapy
        _qap_replicates(rng, n_sim, threaded) do r, rrng
            dist[r, :] = _ols_tvals(X, _gvectorize(_rmperm(rrng, Y, partition), idx))
        end
    elseif nullhyp == :qapx
        for i in 1:nx
            _qap_replicates(rng, n_sim, threaded) do r, rrng
                Xp = copy(X)
                Xp[:, i] = _gvectorize(_rmperm(rrng, Gx[i], partition), idx)
                dist[r, i] = _ols_tvals(Xp, yv)[i]
            end
        end
    else  # :qapspp — Dekker double semi-partialing
        for i in 1:nx
            others = [j for j in 1:nx if j != i]
            E = _dsp_residual_matrix(X, Gx, idx, i, others, directed)
            template = hcat(X[:, others], zeros(N))
            _qap_replicates(rng, n_sim, threaded) do r, rrng
                Xp = copy(template)
                Xp[:, end] = _gvectorize(_rmperm(rrng, E, partition), idx)
                dist[r, i] = _ols_tvals(Xp, yv)[end]
            end
        end
    end

    pleeq, pgreq, pgreqabs = _perm_pvalues(dist, tstat)
    return NetLMResult(coef, names, tstat, covariance, loglik, pleeq, pgreq, pgreqabs, dist,
                       r_squared, adj_r_squared, N, df_residual, nullhyp,
                       n_sim, intercept, directed, miss)
end

netlm(y, x::Union{AbstractNetwork,AbstractMatrix}; kwargs...) =
    netlm(y, (x,); kwargs...)


# ---------------------------------------------------------------------------
# netlogit
# ---------------------------------------------------------------------------

"""
    NetLogitResult

Result of a [`netlogit`](@ref) network logistic regression.

# Fields
- `coefficients::Vector{Float64}`: Logit coefficient estimates
- `names::Vector{String}`: Coefficient names
- `se::Vector{Float64}`: Standard errors (inverse Fisher information)
- `covariance::Matrix{Float64}`: iid-dyad inverse-information covariance
- `converged::Bool`: Whether the shared Newton fit converged
- `tstat::Vector{Float64}`: the test statistic of each coefficient: the
  signed root likelihood-ratio statistic for `statistic = :lr`, the Wald z
  (`coef ./ se`) for `statistic = :wald` and `nullhyp = :classical`
- `pleeq::Vector{Float64}`: `p(stat_perm <= stat_obs)` per coefficient
- `pgreq::Vector{Float64}`: `p(stat_perm >= stat_obs)` per coefficient
- `pgreqabs::Vector{Float64}`: Two-sided `p(|stat_perm| >= |stat_obs|)`
- `dist::Union{Nothing,Matrix{Float64}}`: Null distribution of the test
  statistics (`nothing` for `nullhyp = :classical`)
- `deviance::Float64`: Residual deviance
- `null_deviance::Float64`: Deviance of the model with every probability
  1/2 (`2N log 2`), as sna reports it; not the intercept-only deviance
- `aic::Float64`, `bic::Float64`: Information criteria
- `n::Int`: Number of dyadic observations
- `df_residual::Int`: Residual degrees of freedom
- `nullhyp::Symbol`: Null hypothesis actually used
- `reps::Int`: Number of Monte Carlo replicates (0 for `:classical`)
- `intercept::Bool`: Whether an intercept column was included
- `directed::Bool`: Whether the dyads were treated as directed
- `statistic::Symbol`: The QAP test statistic, `:lr` or `:wald`
- `separation::NetworkCore.SeparationVerdict`: the separation verdict on the
  observed design (`NetworkCore.logistic_separation`)
- `separated::Vector{String}`: the coefficients whose maximum-likelihood
  value is infinite (the separated terms); empty when the observed design
  has a finite maximiser. When it is not empty, `converged` is `false` and
  the test statistics, p-values and confidence intervals are `NaN`

# Example
```julia
using SNA, Random
flo = load_dataset(:florentine_marriage)
biz = load_dataset(:florentine_business)
fit = netlogit(flo, biz; n_sim=200, rng=Xoshiro(1))
fit isa NetLogitResult, fit.statistic
```
"""
struct NetLogitResult
    coefficients::Vector{Float64}
    names::Vector{String}
    se::Vector{Float64}
    covariance::Matrix{Float64}
    converged::Bool
    tstat::Vector{Float64}
    pleeq::Vector{Float64}
    pgreq::Vector{Float64}
    pgreqabs::Vector{Float64}
    dist::Union{Nothing,Matrix{Float64}}
    deviance::Float64
    null_deviance::Float64
    aic::Float64
    bic::Float64
    n::Int
    df_residual::Int
    nullhyp::Symbol
    reps::Int
    intercept::Bool
    directed::Bool
    # See NetLMResult: the missing-dyad policy the fit ran under.
    missing_method::Symbol
    statistic::Symbol
    separation::NetworkCore.SeparationVerdict
    separated::Vector{String}
end

function Base.show(io::IO, result::NetLogitResult)
    println(io, "Network Logit Model (QAP)")
    println(io, "=========================")
    println(io, "Null hypothesis: $(result.nullhyp)" *
                (result.dist === nothing ? "" : " ($(result.reps) replications)"))
    result.dist === nothing ||
        println(io, "Test statistic: " *
                    (result.statistic === :lr ? "signed root likelihood ratio (LR z)" :
                                                "Wald z"))
    println(io, "Dyadic observations: $(result.n) " *
                "($(result.directed ? "directed" : "undirected") dyads)")
    println(io)
    show(io, coeftable(result))
    println(io)
    println(io, "Null deviance (all probabilities 1/2, as sna): " *
                "$(round(result.null_deviance, digits=2)), " *
                "Residual deviance: $(round(result.deviance, digits=2))")
    print(io, "AIC: $(round(result.aic, digits=2)), " *
              "BIC: $(round(result.bic, digits=2))")
    caveat = NetworkCore.separation_caveat(result.separation, result.names)
    caveat === nothing || print(io, "\nWarning: ", caveat)
end

# The logistic fit of `y` on `X` with the ecosystem's separation verdict
# (`NetworkCore.logistic_separation`, a property of the data): on a separated
# design no finite maximum exists, the coefficients are where Newton stopped
# on the asymptote, and `converged` is false. The optimizer's own warning
# about the singular information there is silenced, because the separation
# warning the caller emits says what is wrong.
function _logit_fit(X::Matrix{Float64}, y::Vector{Float64};
                    maxiter::Int=100, tol::Float64=1e-8)
    yb = Bool.(y)
    verdict = NetworkCore.logistic_separation(X, yb)
    derivatives = NetworkCore.logistic_derivatives(X, yb)
    fit = if verdict.separated
        with_logger(NullLogger()) do
            NetworkCore.newton_fit(derivatives, zeros(size(X, 2)); maxiter, tol)
        end
    else
        NetworkCore.newton_fit(derivatives, zeros(size(X, 2)); maxiter, tol)
    end
    return (coef=fit.θ, se=fit.se, covariance=fit.vcov, deviance=-2 * fit.loglik,
            converged=fit.converged && !verdict.separated, verdict=verdict)
end


# The minimized logistic deviance of y on X (the infimum when the data are
# separated: Newton then follows the separating direction until the
# deviance stops decreasing, which leaves it within the optimizer's
# tolerance of the infimum). No identification is required, so it is the
# building block of the likelihood-ratio QAP statistic, which exists on
# separated permutations where the Wald statistic does not. Returns the
# deviance and the coefficients (whose signs give the signed root).
function _logit_deviance(X::AbstractMatrix{Float64}, y::AbstractVector{Float64})
    size(X, 2) == 0 && return (2 * length(y) * log(2), Float64[])
    derivatives = NetworkCore.logistic_derivatives(X, Bool.(y))
    # A separated permutation has singular information; newton_fit says so,
    # but here that is expected and only the deviance is used.
    fit = with_logger(NullLogger()) do
        NetworkCore.newton_fit(derivatives, zeros(size(X, 2)); maxiter=200, tol=1e-10)
    end
    return (-2 * fit.loglik, fit.θ)
end

# Signed root of the likelihood-ratio statistic for one coefficient.
_signed_root(d_reduced, d_full, θk) = sign(θk) * sqrt(max(d_reduced - d_full, 0.0))

"""
    netlogit(y, xs; intercept=true, nullhyp=:qapspp, statistic=:lr, n_sim=1000,
             mode=:auto, threaded=true, missing=:error,
             rng=Random.default_rng()) -> NetLogitResult

Logistic regression of the binary network `y` on one or more predictor
networks `xs`, with QAP null-hypothesis testing, R `sna::netlogit`.

The networks are vectorized over dyads exactly as in [`netlm`](@ref)
(including its treatment of `mode`, two-mode data and missing dyads), and a
binomial GLM with logit link is fitted by the shared Newton optimizer. The
supported null hypotheses (`:qap`/`:qapspp`, Dekker double semi-partialing,
the default; `:qapy`; `:qapx`; `:classical`) have the same meaning as in
[`netlm`](@ref). As in R `sna`, the semi-partialing residuals are computed by
*linear* regression of each predictor on the others, then permuted and
refitted in the logistic model.

`y` must be dichotomous (all dyad values 0 or 1).

# The test statistic

A QAP test compares a statistic with its permutation distribution, and any
statistic gives a valid test. sna uses the Wald z, `coef / se`. On sparse
binary data that statistic often does not exist for some permutation: the
permuted design is separated, the maximum-likelihood coefficient is
infinite, and sna's `glm` returns a meaningless z near 0 for it. On
Florentine marriage ~ business ties about 6 % of the permutations are
separated.

- `statistic=:lr` (default): the signed root likelihood-ratio statistic
  `sign(β̂ₖ) √(D₋ₖ − D)`, where `D` is the model's deviance and `D₋ₖ` the
  deviance without coefficient `k`. It is finite on separated permutations
  too (the deviance has an infimum), so every replicate counts and none is
  dropped. In large samples it is close to the Wald z.
- `statistic=:wald`: sna's Wald z. A separated permutation has no Wald z,
  so the permutation p-values of every coefficient with a separated
  replicate are withheld (`NaN`, with a warning saying how many replicates
  were separated) rather than computed from a fabricated z or from the
  surviving replicates only.

`nullhyp=:classical` fits only the observed design and reports sna's Wald
test against a t distribution; its p-values assume independent dyads, which
dyadic dependence typically violates badly, and are not a substitute for
QAP inference.

If the *observed* design is completely or quasi-completely separated there
is no finite estimate. `netlogit` then follows the ecosystem's separation
policy (`NetworkCore.logistic_separation`): it warns, returns the fit with
`converged == false`, names the separated coefficients in `fit.separated`,
and withholds inference. The test statistics, p-values and confidence
intervals are `NaN`, and no permutation is run. R's `glm` warns "fitted
probabilities numerically 0 or 1 occurred" for the same design, and
`sna::netlogit` reports its meaningless estimates and p-values.

The reported null deviance is that of the model with every probability
1/2, `2N log 2`, as sna computes it (sna fits without a separate
intercept); a pseudo-R² computed from it is larger than one computed from
the intercept-only deviance. Missing dyads are handled exactly as in
[`netlm`](@ref): `netlogit` implements no missing-data estimator
(`NetworkCore.supports_missing(netlogit) == false`).

`n_sim` replicates (R's `reps`) use independently
seeded random number generators drawn from `rng` and are reproducible
across thread counts; `threaded=false` disables threading.

# Example
```julia
using SNA, Random
flo = load_dataset(:florentine_marriage)
biz = load_dataset(:florentine_business)
fit = netlogit(flo, biz; n_sim=1000, rng=Xoshiro(1))   # DSP QAP, LR statistic
fit.pgreqabs                         # business ties: p < 0.001
```
"""
function netlogit(y, xs::Union{Tuple,AbstractVector}; intercept::Bool=true,
                  nullhyp::Symbol=:qapspp, statistic::Symbol=:lr, n_sim::Int=1000,
                  mode::Symbol=:auto,
                  threaded::Bool=true,
                  missing::Symbol=:error,
                  rng::Random.AbstractRNG=Random.default_rng())
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    statistic in (:lr, :wald) ||
        throw(ArgumentError("netlogit: statistic must be :lr or :wald; got :$statistic"))
    yv, X, Y, Gx, idx, names, directed, partition =
        _qap_design(y, xs, intercept, mode, "netlogit", policy)
    miss = _missing_method_of(policy, y, xs...)
    N, nx = size(X)
    N > nx || throw(ArgumentError("more predictors than dyadic observations"))
    rank(X) == nx || throw(ArgumentError("predictors are linearly dependent; remove redundant columns"))
    all(v -> v == 0.0 || v == 1.0, yv) ||
        throw(ArgumentError("netlogit requires a dichotomous (0/1) " *
                            "dependent network"))

    base = _logit_fit(X, yv)
    separation = base.verdict
    separated = names[separation.terms]
    NetworkCore.warn_separation("netlogit", separation, names;
        note="R's glm warns \"fitted probabilities numerically 0 or 1 " *
             "occurred\" for the same design, and sna::netlogit reports the " *
             "estimates and p-values regardless.")
    base.converged || separation.separated ||
        @warn "netlogit did not converge; inspect fit.converged"
    coef, se = base.coef, base.se
    deviance = base.deviance
    # sna calls glm.fit(intercept = FALSE) (the intercept is an explicit
    # ones column), so the null model is mu = 1/2
    null_deviance = 2 * N * log(2)
    aic = deviance + 2 * nx
    bic = deviance + nx * log(N)
    df_residual = N - nx

    nullhyp = _resolve_nullhyp(nullhyp, nx, "netlogit")

    if separation.separated
        # No finite maximum: inference is withheld (the ecosystem's policy),
        # and the permutation test, whose statistic is read off the same
        # non-existent maximum, is not run.
        nan = fill(NaN, nx)
        return NetLogitResult(coef, names, se, base.covariance, false, nan, copy(nan),
                              copy(nan), copy(nan), nothing, deviance, null_deviance,
                              aic, bic, N, df_residual, nullhyp, 0, intercept, directed,
                              miss, nullhyp == :classical ? :wald : statistic,
                              separation, separated)
    end

    if nullhyp == :classical
        # As in sna::netlogit, classical p-values use the t distribution
        # with the residual degrees of freedom
        tstat = coef ./ se
        tdist = TDist(df_residual)
        pleeq = cdf.(tdist, tstat)
        pgreq = ccdf.(tdist, tstat)
        pgreqabs = 2 .* ccdf.(tdist, abs.(tstat))
        return NetLogitResult(coef, names, se, base.covariance, base.converged, tstat, pleeq, pgreq, pgreqabs,
                              nothing, deviance, null_deviance, aic, bic, N,
                              df_residual, nullhyp, 0, intercept, directed, miss, :wald,
                              separation, separated)
    end

    n_sim > 0 || throw(ArgumentError("n_sim must be positive"))
    # sna's Wald z of a replicate; on a separated replicate it does not
    # exist (NaN), and the coefficients it touches get no permutation p-value
    function zvals(Xp, yp)
        f = _logit_fit(Xp, yp)
        f.verdict.separated && return fill(NaN, size(Xp, 2))
        f.converged && all(isfinite, f.se) ||
            throw(ArgumentError("a QAP logistic replicate failed to converge or identify its coefficients"))
        return f.coef ./ f.se
    end
    without(M, k) = M[:, [j for j in 1:size(M, 2) if j != k]]
    # Deviances of the observed design without each coefficient: the reduced
    # models of the observed LR statistics, and (for :qapx and :qapspp, which
    # leave the other columns and y fixed) of every replicate too.
    reduced = statistic === :lr ? [_logit_deviance(without(X, k), yv)[1] for k in 1:nx] :
              Float64[]
    tstat = if statistic === :lr
        d_obs = _logit_deviance(X, yv)[1]
        [_signed_root(reduced[k], d_obs, coef[k]) for k in 1:nx]
    else
        coef ./ se
    end

    dist = Matrix{Float64}(undef, n_sim, nx)
    if nullhyp == :qapy
        _qap_replicates(rng, n_sim, threaded) do r, rrng
            yp = _gvectorize(_rmperm(rrng, Y, partition), idx)
            if statistic === :lr
                d, θ = _logit_deviance(X, yp)
                for k in 1:nx
                    dist[r, k] = _signed_root(_logit_deviance(without(X, k), yp)[1], d, θ[k])
                end
            else
                dist[r, :] = zvals(X, yp)
            end
        end
    elseif nullhyp == :qapx
        for i in 1:nx
            _qap_replicates(rng, n_sim, threaded) do r, rrng
                Xp = copy(X)
                Xp[:, i] = _gvectorize(_rmperm(rrng, Gx[i], partition), idx)
                if statistic === :lr
                    d, θ = _logit_deviance(Xp, yv)
                    dist[r, i] = _signed_root(reduced[i], d, θ[i])
                else
                    dist[r, i] = zvals(Xp, yv)[i]
                end
            end
        end
    else  # :qapspp — Dekker double semi-partialing
        for i in 1:nx
            others = [j for j in 1:nx if j != i]
            E = _dsp_residual_matrix(X, Gx, idx, i, others, directed)
            template = hcat(X[:, others], zeros(N))
            _qap_replicates(rng, n_sim, threaded) do r, rrng
                Xp = copy(template)
                Xp[:, end] = _gvectorize(_rmperm(rrng, E, partition), idx)
                if statistic === :lr
                    # the reduced model y ~ X[:, others] is the observed one
                    d, θ = _logit_deviance(Xp, yv)
                    dist[r, i] = _signed_root(reduced[i], d, θ[end])
                else
                    dist[r, i] = zvals(Xp, yv)[end]
                end
            end
        end
    end

    pleeq, pgreq, pgreqabs = _perm_pvalues(dist, tstat)
    # A Wald replicate with no z (separated) leaves its coefficient without a
    # permutation distribution: withhold that p-value rather than compute it
    # from the surviving replicates.
    withheld = [k for k in 1:nx if any(isnan, view(dist, :, k))]
    if !isempty(withheld)
        nsep = count(r -> any(isnan, view(dist, r, :)), 1:n_sim)
        @warn "netlogit: $nsep of $n_sim permutation replicates are separated, so " *
              "sna's Wald statistic does not exist for them; the permutation " *
              "p-values of $(join(("`" * names[k] * "`" for k in withheld), ", ")) " *
              "are withheld (NaN). The default statistic=:lr exists on separated " *
              "replicates."
        for k in withheld
            pleeq[k] = pgreq[k] = pgreqabs[k] = NaN
        end
    end
    return NetLogitResult(coef, names, se, base.covariance, base.converged, tstat, pleeq, pgreq, pgreqabs,
                          dist, deviance, null_deviance, aic, bic, N,
                          df_residual, nullhyp, n_sim, intercept, directed, miss, statistic,
                          separation, separated)
end

netlogit(y, x::Union{AbstractNetwork,AbstractMatrix}; kwargs...) =
    netlogit(y, (x,); kwargs...)

# ---------------------------------------------------------------------------
# The shared result-metadata protocol (NetworkCore.jl `src/results.jl`)
# ---------------------------------------------------------------------------
#
# `fit_metadata(fit)` collects these accessors. The point they have to make for
# the QAP regressions is that the ESTIMATOR and the INFERENCE come from
# different places: the coefficients are a plain dyad-independent OLS/logit fit,
# and everything that makes the result a *network* method lives in the
# permutation null the p-values are read off.

estimand(::NetLMResult) = :network_regression
estimand(::NetLogitResult) = :network_logit_regression

"""
    objective(::NetLMResult) -> Symbol

`:least_squares` — the coefficients are ordinary least squares of the vectorized
dyads of `y` on the vectorized dyads of the predictors. Nothing about the QAP
null enters the point estimates.
"""
objective(::NetLMResult) = :least_squares

"""
    objective(::NetLogitResult) -> Symbol

`:likelihood` — the coefficients maximize the binomial (logit-link) likelihood of
the vectorized dyads, fitted by the shared Newton optimizer.
"""
objective(::NetLogitResult) = :likelihood

"""
    is_exact(::NetLMResult) -> Bool

`true`: OLS solves its objective in closed form, and that objective is the exact
Gaussian likelihood of the model actually fitted — a regression treating the
dyads as independent observations. It is **not** a statement that the dyads are
independent; that assumption is exactly what the QAP permutation null in
`NetworkCore.approximations` is there to work around.
"""
is_exact(::NetLMResult) = true

"""
    is_exact(r::NetLogitResult) -> Bool

`true` when the Newton fit converged to a finite maximum: it then maximizes
the exact binomial likelihood of the dyad-independent logit model that is
being fitted. As with [`netlm`](@ref), this says the objective is not
approximated — not that the independence assumption holds. `false` for an
unconverged fit and for a separated design, which has no maximum.
"""
is_exact(r::NetLogitResult) = r.converged

"""
    se_method(::NetLMResult) -> Symbol

`:ols` — homoskedastic iid-dyad OLS covariance. These standard errors
studentize the QAP statistic; they do not account for dyadic dependence.
"""
se_method(::NetLMResult) = :ols

"""
    se_method(::NetLogitResult) -> Symbol

`:fisher` — the inverse Fisher information of the binomial GLM (`inv(X'WX)` at
convergence), which treats every dyad as an independent observation and is
therefore anticonservative under dyadic dependence. The reported p-values do NOT
come from these standard errors (unless `nullhyp = :classical`): they come from
the QAP permutation null.
"""
se_method(::NetLogitResult) = :fisher

# Record the actual missing-dyad treatment of the fitted response/predictors.
missing_method(r::NetLMResult) = r.missing_method
missing_method(r::NetLogitResult) = r.missing_method

# Shared caveats for both QAP regressions.
function _qap_approximations(nullhyp::Symbol, reps::Int, se_note::String)
    out = String[
        "the dyads are treated as independent observations by the estimator; " *
        "the dyadic dependence is addressed only by the null distribution the " *
        "p-values are read off",
        se_note,
    ]
    if nullhyp === :classical
        push!(out, "nullhyp = :classical: the p-values are parametric t/z tests " *
                   "that assume independent dyads — for reference only, since " *
                   "dyadic dependence typically invalidates them. No permutation " *
                   "was performed")
    else
        push!(out, "p-values are Monte-Carlo QAP permutation p-values " *
                   "(nullhyp = :$nullhyp, $reps replications): they carry " *
                   "simulation error and, as in sna, are proportions without " *
                   "a +1 correction, so they can be exactly 0 or 1")
    end

    return out
end

approximations(r::NetLMResult) =
    _qap_approximations(r.nullhyp, r.reps,
        "the standard errors are homoskedastic " *
        "iid-dyad OLS standard errors and serve only as the test statistic " *
        "compared against the null distribution")

function approximations(r::NetLogitResult)
    out = _qap_approximations(r.nullhyp, r.reps,
        "the standard errors are the inverse Fisher information of a binomial " *
        "GLM that treats the dyads as independent: they are expected " *
        "anticonservative under dyadic dependence, and the reported p-values " *
        "are not derived from them")
    if r.nullhyp !== :classical && r.statistic === :lr
        push!(out, "the QAP statistic is the signed root likelihood ratio, not sna's " *
                   "Wald z: it exists on separated permutations, so no replicate is dropped")
    end
    caveat = NetworkCore.separation_caveat(r.separation, r.names)
    caveat === nothing || push!(out, caveat)
    if r.dist !== nothing && any(isnan, r.pgreqabs)
        push!(out, "some permutation replicates are separated, so their Wald z does not " *
                   "exist; the permutation p-values of " *
                   join(("`" * r.names[k] * "`" for k in findall(isnan, r.pgreqabs)), ", ") *
                   " are withheld")
    end
    r.converged || caveat !== nothing || push!(out, "logistic fit did not converge")
    return out
end


# These covariance matrices and intervals describe the dyad-independent point
# model; QAP inference resides in pgreqabs, never in a permutation-derived SE.
const _NetworkRegressionResult = Union{NetLMResult,NetLogitResult}
coef(r::_NetworkRegressionResult) = r.coefficients
coefnames(r::_NetworkRegressionResult) = copy(r.names)
vcov(r::_NetworkRegressionResult) = r.covariance
stderror(r::_NetworkRegressionResult) = sqrt.(diag(vcov(r)))
nobs(r::_NetworkRegressionResult) = r.n
dof(r::NetLMResult) = length(coef(r)) + 1 # Gaussian residual variance
dof(r::NetLogitResult) = length(coef(r))
loglikelihood(r::NetLMResult) = r.loglik
loglikelihood(r::NetLogitResult) = -r.deviance / 2
aic(r::_NetworkRegressionResult) = -2loglikelihood(r) + 2dof(r)
bic(r::_NetworkRegressionResult) = -2loglikelihood(r) + log(nobs(r)) * dof(r)
function confint(r::_NetworkRegressionResult; level::Real=0.95)
    0 < level < 1 || throw(ArgumentError("level must be between 0 and 1"))
    # a separated logit has no finite maximum: no interval (the policy)
    r isa NetLogitResult && r.separation.separated &&
        return fill(NaN, length(coef(r)), 2)
    q = quantile(r isa NetLMResult ? TDist(r.df_residual) : Normal(), (1 + level) / 2)
    delta = q * stderror(r)
    return hcat(coef(r) - delta, coef(r) + delta)
end
function coeftable(r::_NetworkRegressionResult)
    stat = r isa NetLMResult ? "t" : r.statistic === :lr ? "LR z" : "z"
    return NetworkCore.CoefficientTable(r.names, coef(r), stderror(r);
        z_values=r.tstat, p_values=r.pgreqabs,
        header=("Estimate", "Std. Error", "$stat value", "Pr(>=|$stat|)"),
        p_floor=r.reps > 0 ? 1 / r.reps : 1e-16)
end
