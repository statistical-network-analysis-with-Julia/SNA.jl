# ---------------------------------------------------------------------------
# Conditional uniform graph tests
# ---------------------------------------------------------------------------

"""
    CUGTestResult

Result of a [`cug_test`](@ref), with the fields of R's `cug.test` object.

# Fields
- `obs_stat::Float64`: the observed statistic
- `rep_stat::Vector{Float64}`: its values on the simulated networks
- `cmode::Symbol`: what the simulated networks were conditioned on
- `directed::Bool`: whether the networks are directed
- `plteobs::Float64`: the proportion of simulated values `<=` the observed one
- `pgteobs::Float64`: the proportion of simulated values `>=` the observed one
- `reps::Int`: the number of simulated networks (`n_sim`)

# Example
```julia
using SNA, Random
flo = load_dataset(:florentine_marriage)
ct = cug_test(flo, gtrans; cmode=:edges, n_sim=200, rng=Xoshiro(1))
ct isa CUGTestResult, ct.obs_stat, ct.pgteobs
```
"""
struct CUGTestResult
    obs_stat::Float64
    rep_stat::Vector{Float64}
    cmode::Symbol
    directed::Bool
    plteobs::Float64
    pgteobs::Float64
    reps::Int
end

function Base.show(io::IO, r::CUGTestResult)
    println(io, "Univariate Conditional Uniform Graph Test")
    println(io, "=========================================")
    println(io, "Conditioning: $(r.cmode) ($(r.directed ? "directed" : "undirected")), " *
                "$(r.reps) replications")
    println(io)
    println(io, "Observed value: $(round(r.obs_stat, digits=6))")
    println(io, "  Pr(X >= Obs): $(NetworkCore.format_pvalue(r.pgteobs; floor=1/r.reps))")
    print(io, "  Pr(X <= Obs): $(NetworkCore.format_pvalue(r.plteobs; floor=1/r.reps))")
end

# A uniform draw from the digraphs with exactly M mutual, A asymmetric and
# the rest null dyads (sna's rguman(method="exact")).
function _rguman(rng::AbstractRNG, n::Integer, M::Integer, A::Integer)
    n, M, A = Int(n), Int(M), Int(A)
    net = network(n; directed=true)
    dyads = [(i, j) for i in 1:n for j in (i+1):n]
    order = randperm(rng, length(dyads))
    for t in 1:M
        i, j = dyads[order[t]]
        add_edge!(net, i, j)
        add_edge!(net, j, i)
    end
    for t in (M+1):(M+A)
        i, j = dyads[order[t]]
        rand(rng, Bool) ? add_edge!(net, i, j) : add_edge!(net, j, i)
    end
    return net
end

"""
    cug_test(net, f; cmode=:size, n_sim=1000, threaded=true, missing=:error,
             rng=Random.default_rng(), kwargs...) -> CUGTestResult

Conditional uniform graph (CUG) test of the graph-level statistic `f`, R
`sna::cug.test`. The observed `f(net; kwargs...)` is compared with its
distribution over `n_sim` networks drawn uniformly at random from the
networks that share a feature of `net`:

- `cmode=:size` (default): the same number of vertices (every tie present
  independently with probability 1/2);
- `cmode=:edges`: the same number of vertices and of ties;
- `cmode=:dyad_census`: the same numbers of mutual, asymmetric and null dyads.

The simulated networks have the directedness of `net` and no self-loops.
`f` is any function of a `Network` that returns a number, for instance
[`gtrans`](@ref), [`grecip`](@ref) or `g -> centralization(g, degreecent)`;
keyword arguments other than those listed are passed to it.

The result reports `pgteobs` and `plteobs`, the proportions of simulated
values at least and at most the observed one (as in sna, without a `+1`
correction, so they can be 0). A small `pgteobs` says the statistic is
larger than expected from the conditioning feature alone. Simulated values
that are `NaN` count in neither tail.

Self-loops in `net` are ignored (sna's `diag=FALSE`); conditioning with
loops (`diag=TRUE`), valued networks (`ignore.eval=FALSE`) and two-mode
networks are not implemented and raise an `ArgumentError`.

Replicates use independently seeded random number generators drawn from
`rng`, so the result is reproducible and does not depend on the number of
threads; `threaded=false` runs them serially. `f` must be thread-safe when
threading is on. `missing=:error` refuses a network with masked dyads;
`missing=:face` uses their stored face values (`f` itself is called on
unmasked networks).

# Example
```julia
using SNA, Random
flo = load_dataset(:florentine_marriage)
# Is the marriage network more transitive than its density implies?
ct = cug_test(flo, gtrans; cmode=:edges, n_sim=500, rng=Xoshiro(1))
ct.obs_stat                          # 0.1915
ct.pgteobs                           # about 0.33: no
```
"""
function cug_test(net::AbstractNetwork, f; cmode::Symbol=:size, n_sim::Int=1000,
                  threaded::Bool=true, missing::Symbol=:error,
                  rng::Random.AbstractRNG=Random.default_rng(), kwargs...)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="cug_test")
    cmode in (:size, :edges, :dyad_census) ||
        throw(ArgumentError("cug_test: cmode must be :size, :edges or :dyad_census; got :$cmode"))
    n_sim > 0 || throw(ArgumentError("n_sim must be positive"))
    _qap_partition(net) === nothing ||
        throw(ArgumentError("cug_test: two-mode networks are not implemented"))
    n = Int(nv(net))          # the generators draw Int networks whatever the id type
    directed = is_directed(net)
    # The observed network without self-loops (and without a mask, whose
    # policy has been applied above), as the simulated networks are.
    A = _sociomatrix(net)
    obs_net = network_from_matrix(A; directed)
    dc = dyad_census(obs_net)
    m = directed ? 2 * dc.mutual + dc.asymmetric : dc.mutual

    draw = if cmode == :size
        r -> rgnp(n, 0.5; directed, rng=r)
    elseif cmode == :edges || !directed
        r -> rgnm(n, m; directed, rng=r)
    else
        r -> _rguman(r, n, dc.mutual, dc.asymmetric)
    end
    obs = Float64(f(obs_net; kwargs...))
    rep = Vector{Float64}(undef, n_sim)
    _qap_replicates(rng, n_sim, threaded) do r, rrng
        rep[r] = Float64(f(draw(rrng); kwargs...))
    end
    return CUGTestResult(obs, rep, cmode, directed, count(<=(obs), rep) / n_sim,
                         count(>=(obs), rep) / n_sim, n_sim)
end
