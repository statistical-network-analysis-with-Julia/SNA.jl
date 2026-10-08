"""
Random network generators.

Provides Bernoulli and fixed-edge-count random networks in the style of
R `sna::rgraph` and the classic Erdős–Rényi G(n,p) / G(n,m) models.
"""

using Random

"""
    rgraph(n; m=1, tprob=0.5, mode=:digraph, rng=Random.default_rng())

Generate Bernoulli random networks (R `sna::rgraph`).

# Arguments
- `n::Int`: Number of vertices
- `m::Int=1`: Number of networks to generate
- `tprob::Float64=0.5`: Tie probability
- `mode::Symbol=:digraph`: `:digraph` (directed) or `:graph` (undirected)
- `rng`: Random number generator

# Returns
A single `Network` when `m == 1`, otherwise a `Vector{Network}`.

# Example
```julia
using SNA, Random
net = rgraph(10; tprob=0.2, rng=Xoshiro(1))
nets = rgraph(10; m=3, tprob=0.2, mode=:graph, rng=Xoshiro(2))
```
"""
function rgraph(n::Integer; m::Int=1, tprob::Float64=0.5, mode::Symbol=:digraph,
                rng::Random.AbstractRNG=Random.default_rng())
    mode in (:digraph, :graph) ||
        throw(ArgumentError("mode must be :digraph or :graph"))
    directed = mode == :digraph

    nets = [rgnp(n, tprob; directed=directed, rng=rng) for _ in 1:m]
    return m == 1 ? nets[1] : nets
end

"""
    rgnp(n, p; directed=true, rng=Random.default_rng()) -> Network

Generate an Erdős–Rényi G(n, p) random network: each possible edge is
present independently with probability `p`.

# Example
```julia
using SNA, Random
rgnp(20, 0.1; rng=Xoshiro(1))
```
"""
function rgnp(n::Integer, p::Float64; directed::Bool=true,
              rng::Random.AbstractRNG=Random.default_rng())
    n = Int(n)
    0.0 <= p <= 1.0 || throw(ArgumentError("p must be in [0, 1]"))
    net = network(n; directed=directed)

    for i in 1:n
        j_range = directed ? (1:n) : ((i+1):n)
        for j in j_range
            i == j && continue
            if rand(rng) < p
                add_edge!(net, i, j)
            end
        end
    end

    return net
end

"""
    rgnm(n, m; directed=true, rng=Random.default_rng()) -> Network

Generate an Erdős–Rényi G(n, m) random network: exactly `m` edges placed
uniformly at random among the possible dyads. (sna's `rgnm(n, nv, m)`
takes the number of networks first; here a single network is returned.)

# Example
```julia
using SNA, Random
ne(rgnm(10, 15; directed=false, rng=Xoshiro(1)))    # 15
```
"""
function rgnm(n::Integer, m::Integer; directed::Bool=true,
              rng::Random.AbstractRNG=Random.default_rng())
    n, m = Int(n), Int(m)
    max_edges = directed ? n * (n - 1) : n * (n - 1) ÷ 2
    0 <= m <= max_edges ||
        throw(ArgumentError("m must be between 0 and $max_edges"))

    # Enumerate dyads and sample m of them without replacement
    dyads = Tuple{Int,Int}[]
    for i in 1:n
        j_range = directed ? (1:n) : ((i+1):n)
        for j in j_range
            i == j && continue
            push!(dyads, (i, j))
        end
    end

    net = network(n; directed=directed)
    for idx in randperm(rng, length(dyads))[1:m]
        add_edge!(net, dyads[idx][1], dyads[idx][2])
    end

    return net
end
