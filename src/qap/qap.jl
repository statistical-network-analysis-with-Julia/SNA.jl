"""
QAP inference: permutation tests and network regression.

Provides the quadratic assignment procedure (QAP) test for arbitrary
graph-level statistics, and OLS/logistic network regression with QAP
null-hypothesis testing (Dekker's double-semi-partialing by default),
following R `sna::qaptest`, `sna::netlm`, and `sna::netlogit`.
"""

using Distributions
using LinearAlgebra
using Random
using Statistics

# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

# Coerce a network or matrix argument to a dense Float64 sociomatrix
_sociomatrix(x::AbstractNetwork) = Float64.(as_matrix(x))
_sociomatrix(x::AbstractMatrix) = Matrix{Float64}(x)

# Apply a random vertex permutation to rows and columns (sna::rmperm)
function _rmperm(rng::Random.AbstractRNG, A::Matrix{Float64})
    p = randperm(rng, size(A, 1))
    return A[p, p]
end

# Off-diagonal dyad index list: all ordered pairs for directed data, the
# lower triangle (i > j) for undirected data (as in sna::gvectorize)
function _dyad_indices(n::Int, directed::Bool)
    idx = Tuple{Int,Int}[]
    for j in 1:n
        i_range = directed ? (1:n) : ((j+1):n)
        for i in i_range
            i == j && continue
            push!(idx, (i, j))
        end
    end
    return idx
end

# Vectorize the selected dyads of a sociomatrix
_gvectorize(A::Matrix{Float64}, idx::Vector{Tuple{Int,Int}}) =
    Float64[A[i, j] for (i, j) in idx]

# Whether the dyadic data should be treated as directed: taken from the
# network itself when possible, from the `mode` keyword for raw matrices
function _qap_directed(y, mode::Symbol)
    mode in (:auto, :digraph, :graph) ||
        throw(ArgumentError("mode must be :auto, :digraph, or :graph"))
    mode == :digraph && return true
    mode == :graph && return false
    return y isa AbstractNetwork ? is_directed(y) : true
end

# Significance stars for permutation/classical p-values (R conventions)
function _signif_code(p::Float64)
    return p < 0.001 ? "***" :
           p < 0.01 ? "**" :
           p < 0.05 ? "*" :
           p < 0.1 ? "." : ""
end

function _print_coef_table(io::IO, names::Vector{String},
                           coefficients::Vector{Float64},
                           tstat::Vector{Float64}, stat_label::String,
                           pgreqabs::Vector{Float64})
    println(io, rpad("", 14), lpad("Estimate", 12), lpad(stat_label, 12),
            lpad("Pr(>=|stat|)", 14))
    println(io, "-"^56)
    for i in eachindex(names)
        println(io, rpad(names[i], 14),
                lpad(round(coefficients[i], digits=6), 12),
                lpad(round(tstat[i], digits=4), 12),
                lpad(round(pgreqabs[i], digits=4), 14), " ",
                _signif_code(pgreqabs[i]))
    end
    println(io, "-"^56)
    println(io, "Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1")
end

# ---------------------------------------------------------------------------
# qaptest
# ---------------------------------------------------------------------------

"""
    QAPTestResult

Result of a [`qaptest`](@ref) permutation test.

# Fields
- `testval::Float64`: Observed value of the test statistic
- `dist::Vector{Float64}`: Monte Carlo null distribution (one value per
  replication)
- `pgreq::Float64`: Proportion of null draws `>=` the observed value
- `pleeq::Float64`: Proportion of null draws `<=` the observed value
- `reps::Int`: Number of replications
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
    println(io, "  p(f(perm) >= f(d)):  $(round(result.pgreq, digits=6))")
    print(io, "  p(f(perm) <= f(d)):  $(round(result.pleeq, digits=6))")
end

"""
    qaptest(f, g1, g2; reps=1000, rng=Random.default_rng()) -> QAPTestResult

Perform a quadratic assignment procedure (QAP) test for the graph-level
statistic `f` on the network pair `(g1, g2)`, following R `sna::qaptest`.

`f` is called as `f(A1, A2)` on the adjacency matrices of the two networks
and must return a scalar (e.g. graph correlation). The observed value
`f(A1, A2)` is compared against the null distribution obtained by applying
`reps` independent uniform vertex permutations to `g1` (simultaneous row
and column permutation) while holding `g2` fixed — the QAP null hypothesis
of no association between the two structures conditional on both.

`g1` and `g2` may be `Network` objects or adjacency matrices of equal size.

# Example
```julia
flo = load_dataset(:florentine_marriage)
biz = load_dataset(:florentine_business)
gcor(a, b) = cor(vec(a), vec(b))
qt = qaptest(gcor, flo, biz; reps=1000)
```
"""
function qaptest(f, g1, g2; reps::Int=1000,
                 rng::Random.AbstractRNG=Random.default_rng())
    reps > 0 || throw(ArgumentError("reps must be positive"))
    A = _sociomatrix(g1)
    B = _sociomatrix(g2)
    size(A) == size(B) ||
        throw(ArgumentError("g1 and g2 must have the same number of vertices"))

    testval = Float64(f(A, B))
    dist = Vector{Float64}(undef, reps)
    for r in 1:reps
        dist[r] = Float64(f(_rmperm(rng, A), B))
    end

    pgreq = count(>=(testval), dist) / reps
    pleeq = count(<=(testval), dist) / reps

    return QAPTestResult(testval, dist, pgreq, pleeq, reps)
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
- `pleeq::Vector{Float64}`: `p(stat_perm <= stat_obs)` per coefficient
- `pgreq::Vector{Float64}`: `p(stat_perm >= stat_obs)` per coefficient
- `pgreqabs::Vector{Float64}`: Two-sided `p(|stat_perm| >= |stat_obs|)`
  (the p-value reported by `show`)
- `dist::Union{Nothing,Matrix{Float64}}`: Null distribution of the test
  statistics (`reps × k`; `nothing` for `nullhyp = :classical`)
- `r_squared::Float64`: Coefficient of determination
- `adj_r_squared::Float64`: Adjusted R²
- `n::Int`: Number of dyadic observations
- `df_residual::Int`: Residual degrees of freedom
- `nullhyp::Symbol`: Null hypothesis actually used
- `reps::Int`: Number of Monte Carlo replications (0 for `:classical`)
- `intercept::Bool`: Whether an intercept column was included
- `directed::Bool`: Whether the dyads were treated as directed
"""
struct NetLMResult
    coefficients::Vector{Float64}
    names::Vector{String}
    tstat::Vector{Float64}
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
end

function Base.show(io::IO, result::NetLMResult)
    println(io, "Network Regression (OLS + QAP)")
    println(io, "==============================")
    println(io, "Null hypothesis: $(result.nullhyp)" *
                (result.dist === nothing ? "" : " ($(result.reps) replications)"))
    println(io, "Dyadic observations: $(result.n) " *
                "($(result.directed ? "directed" : "undirected") dyads)")
    println(io)
    _print_coef_table(io, result.names, result.coefficients, result.tstat,
                      "t-value", result.pgreqabs)
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

# Build the dyadic design: returns (yv, X, Gx, idx, names, directed)
function _qap_design(y, xs, intercept::Bool, mode::Symbol, fname::String)
    isempty(xs) && !intercept &&
        throw(ArgumentError("$fname needs at least one predictor"))
    directed = _qap_directed(y, mode)
    Y = _sociomatrix(y)
    n = size(Y, 1)

    Gx = Matrix{Float64}[]
    names = String[]
    if intercept
        push!(Gx, ones(n, n))
        push!(names, "(intercept)")
    end
    for (k, x) in enumerate(xs)
        M = _sociomatrix(x)
        size(M) == (n, n) ||
            throw(ArgumentError("predictor x$k must have the same number " *
                                "of vertices as y"))
        push!(Gx, M)
        push!(names, "x$k")
    end

    idx = _dyad_indices(n, directed)
    yv = _gvectorize(Y, idx)
    X = Matrix{Float64}(undef, length(idx), length(Gx))
    for (j, M) in enumerate(Gx)
        X[:, j] = _gvectorize(M, idx)
    end
    return yv, X, Gx, idx, names, directed
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

# Permutation p-values from a reps × k null distribution matrix
function _perm_pvalues(dist::Matrix{Float64}, tstat::Vector{Float64})
    reps = size(dist, 1)
    pleeq = [count(<=(tstat[i]), @view dist[:, i]) / reps for i in eachindex(tstat)]
    pgreq = [count(>=(tstat[i]), @view dist[:, i]) / reps for i in eachindex(tstat)]
    pgreqabs = [count(x -> abs(x) >= abs(tstat[i]), @view dist[:, i]) / reps
                for i in eachindex(tstat)]
    return pleeq, pgreq, pgreqabs
end

"""
    netlm(y, xs; intercept=true, nullhyp=:qapspp, reps=1000, mode=:auto,
          rng=Random.default_rng()) -> NetLMResult

Linear (OLS) regression of the network `y` on one or more predictor
networks `xs`, with QAP null-hypothesis testing, following R `sna::netlm`.

The networks are vectorized over off-diagonal dyads (all ordered pairs for
directed data, one entry per unordered pair for undirected data; the
diagonal is always excluded) and `y` is regressed on the predictors by OLS.
Coefficient significance is assessed by comparing the observed
t-statistics against a permutation null distribution.

# Arguments
- `y`: Dependent network (`Network` or adjacency matrix)
- `xs`: Predictor network(s) — a single network/matrix, or a vector/tuple
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
    - `:classical`: Parametric t-tests (no permutation) — for reference
      only; dyadic dependence typically invalidates these
- `reps::Int=1000`: Number of permutation replications
- `mode::Symbol=:auto`: Dyad set — `:auto` (from `is_directed(y)`, or
  directed for raw matrices), `:digraph`, or `:graph`
- `rng`: Random number generator

The two-sided permutation p-value `pgreqabs` is reported by `show`;
one-sided `pleeq`/`pgreq` are also stored.

# Example
```julia
flo = load_dataset(:florentine_marriage)
biz = load_dataset(:florentine_business)
fit = netlm(flo, biz)
```
"""
function netlm(y, xs::Union{Tuple,AbstractVector}; intercept::Bool=true,
               nullhyp::Symbol=:qapspp, reps::Int=1000, mode::Symbol=:auto,
               rng::Random.AbstractRNG=Random.default_rng())
    yv, X, Gx, idx, names, directed = _qap_design(y, xs, intercept, mode,
                                                  "netlm")
    N, nx = size(X)
    N > nx || throw(ArgumentError("more predictors than dyadic observations"))

    coef = X \ yv
    resid = yv - X * coef
    df_residual = N - nx
    resvar = sum(abs2, resid) / df_residual
    se = sqrt.(max.(diag(inv(Symmetric(X' * X))), 0.0) .* resvar)
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
        return NetLMResult(coef, names, tstat, pleeq, pgreq, pgreqabs,
                           nothing, r_squared, adj_r_squared, N, df_residual,
                           nullhyp, 0, intercept, directed)
    end

    reps > 0 || throw(ArgumentError("reps must be positive"))
    dist = Matrix{Float64}(undef, reps, nx)
    if nullhyp == :qapy
        Y = _sociomatrix(y)
        for r in 1:reps
            dist[r, :] = _ols_tvals(X, _gvectorize(_rmperm(rng, Y), idx))
        end
    elseif nullhyp == :qapx
        for i in 1:nx
            Xp = copy(X)
            for r in 1:reps
                Xp[:, i] = _gvectorize(_rmperm(rng, Gx[i]), idx)
                dist[r, i] = _ols_tvals(Xp, yv)[i]
            end
        end
    else  # :qapspp — Dekker double semi-partialing
        for i in 1:nx
            others = [j for j in 1:nx if j != i]
            E = _dsp_residual_matrix(X, Gx, idx, i, others, directed)
            Xp = hcat(X[:, others], zeros(N))
            for r in 1:reps
                Xp[:, end] = _gvectorize(_rmperm(rng, E), idx)
                dist[r, i] = _ols_tvals(Xp, yv)[end]
            end
        end
    end

    pleeq, pgreq, pgreqabs = _perm_pvalues(dist, tstat)
    return NetLMResult(coef, names, tstat, pleeq, pgreq, pgreqabs, dist,
                       r_squared, adj_r_squared, N, df_residual, nullhyp,
                       reps, intercept, directed)
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
- `tstat::Vector{Float64}`: z-statistics (`coef ./ se`)
- `pleeq::Vector{Float64}`: `p(stat_perm <= stat_obs)` per coefficient
- `pgreq::Vector{Float64}`: `p(stat_perm >= stat_obs)` per coefficient
- `pgreqabs::Vector{Float64}`: Two-sided `p(|stat_perm| >= |stat_obs|)`
- `dist::Union{Nothing,Matrix{Float64}}`: Null distribution of the test
  statistics (`nothing` for `nullhyp = :classical`)
- `deviance::Float64`: Residual deviance
- `null_deviance::Float64`: Null deviance
- `aic::Float64`, `bic::Float64`: Information criteria
- `n::Int`: Number of dyadic observations
- `df_residual::Int`: Residual degrees of freedom
- `nullhyp::Symbol`: Null hypothesis actually used
- `reps::Int`: Number of Monte Carlo replications (0 for `:classical`)
- `intercept::Bool`: Whether an intercept column was included
- `directed::Bool`: Whether the dyads were treated as directed
"""
struct NetLogitResult
    coefficients::Vector{Float64}
    names::Vector{String}
    se::Vector{Float64}
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
end

function Base.show(io::IO, result::NetLogitResult)
    println(io, "Network Logit Model (QAP)")
    println(io, "=========================")
    println(io, "Null hypothesis: $(result.nullhyp)" *
                (result.dist === nothing ? "" : " ($(result.reps) replications)"))
    println(io, "Dyadic observations: $(result.n) " *
                "($(result.directed ? "directed" : "undirected") dyads)")
    println(io)
    _print_coef_table(io, result.names, result.coefficients, result.tstat,
                      "z-value", result.pgreqabs)
    println(io)
    println(io, "Null deviance: $(round(result.null_deviance, digits=2)), " *
                "Residual deviance: $(round(result.deviance, digits=2))")
    print(io, "AIC: $(round(result.aic, digits=2)), " *
              "BIC: $(round(result.bic, digits=2))")
end

_binomial_deviance(y::Vector{Float64}, mu::Vector{Float64}) =
    -2 * sum(y .* log.(clamp.(mu, 1e-16, 1.0)) .+
             (1 .- y) .* log.(clamp.(1 .- mu, 1e-16, 1.0)))

# IRLS logistic regression (binomial GLM with logit link, no implicit
# intercept). Returns coefficients, standard errors from the inverse Fisher
# information, deviance, and fitted probabilities.
function _logit_fit(X::Matrix{Float64}, y::Vector{Float64};
                    maxiter::Int=100, tol::Float64=1e-8)
    N, k = size(X)
    # glm.fit-style start: shrink observed proportions toward 1/2
    mu = (y .+ 0.5) ./ 2
    eta = log.(mu ./ (1 .- mu))
    coef = zeros(k)
    dev = Inf

    for _ in 1:maxiter
        mu = 1 ./ (1 .+ exp.(-eta))
        w = clamp.(mu .* (1 .- mu), 1e-10, Inf)
        z = eta .+ (y .- mu) ./ w
        Xw = X .* w
        coef = Symmetric(X' * Xw) \ (Xw' * z)
        eta = X * coef
        newdev = _binomial_deviance(y, 1 ./ (1 .+ exp.(-eta)))
        converged = abs(newdev - dev) / (abs(newdev) + 0.1) < tol
        dev = newdev
        converged && break
    end

    mu = 1 ./ (1 .+ exp.(-eta))
    w = clamp.(mu .* (1 .- mu), 1e-10, Inf)
    covm = inv(Symmetric(X' * (X .* w)))
    se = sqrt.(max.(diag(covm), 0.0))
    return (coef=coef, se=se, deviance=dev, fitted=mu)
end

"""
    netlogit(y, xs; intercept=true, nullhyp=:qapspp, reps=1000, mode=:auto,
             rng=Random.default_rng()) -> NetLogitResult

Logistic regression of the binary network `y` on one or more predictor
networks `xs`, with QAP null-hypothesis testing, following R
`sna::netlogit`.

The networks are vectorized over off-diagonal dyads exactly as in
[`netlm`](@ref) and a binomial GLM with logit link is fit by iteratively
reweighted least squares. Coefficient significance is assessed by comparing
the observed z-statistics against a permutation null distribution; the
supported null hypotheses (`:qap`/`:qapspp` — Dekker double
semi-partialing, the default — `:qapy`, `:qapx`, and `:classical`) have the
same semantics as in [`netlm`](@ref). As in R `sna`, the semi-partialing
residuals are computed by *linear* regression of each predictor on the
others, then permuted and refit in the logistic model.

`y` must be dichotomous (all dyad values 0 or 1).

# Example
```julia
flo = load_dataset(:florentine_marriage)
biz = load_dataset(:florentine_business)
fit = netlogit(flo, biz)
```
"""
function netlogit(y, xs::Union{Tuple,AbstractVector}; intercept::Bool=true,
                  nullhyp::Symbol=:qapspp, reps::Int=1000, mode::Symbol=:auto,
                  rng::Random.AbstractRNG=Random.default_rng())
    yv, X, Gx, idx, names, directed = _qap_design(y, xs, intercept, mode,
                                                  "netlogit")
    N, nx = size(X)
    N > nx || throw(ArgumentError("more predictors than dyadic observations"))
    all(v -> v == 0.0 || v == 1.0, yv) ||
        throw(ArgumentError("netlogit requires a dichotomous (0/1) " *
                            "dependent network"))

    base = _logit_fit(X, yv)
    coef, se = base.coef, base.se
    tstat = coef ./ se
    deviance = base.deviance
    # sna calls glm.fit(intercept = FALSE) (the intercept is an explicit
    # ones column), so the null model is mu = 1/2
    null_deviance = 2 * N * log(2)
    aic = deviance + 2 * nx
    bic = deviance + nx * log(N)
    df_residual = N - nx

    nullhyp = _resolve_nullhyp(nullhyp, nx, "netlogit")

    if nullhyp == :classical
        # As in sna::netlogit, classical p-values use the t distribution
        # with the residual degrees of freedom
        tdist = TDist(df_residual)
        pleeq = cdf.(tdist, tstat)
        pgreq = ccdf.(tdist, tstat)
        pgreqabs = 2 .* ccdf.(tdist, abs.(tstat))
        return NetLogitResult(coef, names, se, tstat, pleeq, pgreq, pgreqabs,
                              nothing, deviance, null_deviance, aic, bic, N,
                              df_residual, nullhyp, 0, intercept, directed)
    end

    reps > 0 || throw(ArgumentError("reps must be positive"))
    zvals(Xp, yp) = (f = _logit_fit(Xp, yp); f.coef ./ f.se)

    dist = Matrix{Float64}(undef, reps, nx)
    if nullhyp == :qapy
        Y = _sociomatrix(y)
        for r in 1:reps
            dist[r, :] = zvals(X, _gvectorize(_rmperm(rng, Y), idx))
        end
    elseif nullhyp == :qapx
        for i in 1:nx
            Xp = copy(X)
            for r in 1:reps
                Xp[:, i] = _gvectorize(_rmperm(rng, Gx[i]), idx)
                dist[r, i] = zvals(Xp, yv)[i]
            end
        end
    else  # :qapspp — Dekker double semi-partialing
        for i in 1:nx
            others = [j for j in 1:nx if j != i]
            E = _dsp_residual_matrix(X, Gx, idx, i, others, directed)
            Xp = hcat(X[:, others], zeros(N))
            for r in 1:reps
                Xp[:, end] = _gvectorize(_rmperm(rng, E), idx)
                dist[r, i] = zvals(Xp, yv)[end]
            end
        end
    end

    pleeq, pgreq, pgreqabs = _perm_pvalues(dist, tstat)
    return NetLogitResult(coef, names, se, tstat, pleeq, pgreq, pgreqabs,
                          dist, deviance, null_deviance, aic, bic, N,
                          df_residual, nullhyp, reps, intercept, directed)
end

netlogit(y, x::Union{AbstractNetwork,AbstractMatrix}; kwargs...) =
    netlogit(y, (x,); kwargs...)
