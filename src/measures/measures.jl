"""
Network-level measures and indices, with R `sna`'s names: `gden`, `grecip`,
`gtrans`, `dyad_census`, `triad_census`, `hierarchy`, `efficiency`,
`connectedness`, `mutuality`, `reachability`; plus local clustering
coefficients (`transitivity`).
"""

"""
    gden(net; diag=false, discount_bipartite=false, missing=:error) -> Float64

Graph density, the proportion of possible ties that are present: R
`sna::gden`. An undirected tie counts once, so a loop-free undirected network
has density `ne(net) / (n(n−1)/2)`. Loops are excluded unless `diag=true`.
Two-mode density includes all actor pairs in its denominator unless
`discount_bipartite=true`, which restricts it to cross-mode pairs. With
`diag=true` an undirected network has `n(n+1)/2` possible dyads (NetworkCore'
convention), not sna's symmetric-matrix `n²`.

`missing=:error` refuses masked (unobserved) dyads. `missing=:face` counts
each masked dyad at its stored face value: the denominator still contains
every dyad. SNA.jl does not estimate a density over observed dyads only.

`Graphs.density` is Graphs.jl's function and keeps its own semantics on a
`Network` (it counts loops); `gden` is the sna measure.

# Example
```julia
using SNA
net = network(4; directed=false)
add_edges!(net, [(1, 2), (2, 3)])
gden(net)                            # 2 / 6 = 0.3333
```
"""
function gden(net::AbstractNetwork; diag::Bool=false, discount_bipartite::Bool=false,
              missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="gden")
    return network_density(net; missing=policy, diag, discount_bipartite)
end

"""
    grecip(net; measure=:dyadic, missing=:error) -> Float64

Graph reciprocity: R `sna::grecip`. With `M`, `A` and `N` the numbers of
mutual, asymmetric and null dyads:

- `:dyadic` (default): the proportion of dyads that are symmetric,
  `(M + N) / (M + A + N)`;
- `:dyadic_nonnull`: the proportion of non-null dyads that are mutual,
  `M / (M + A)`;
- `:edgewise`: the proportion of ties that are reciprocated, `2M / (2M + A)`;
- `:edgewise_lrr`: the log of edgewise reciprocity over density,
  `log(M (M + A + N) / (M + A/2)²)`;
- `:correlation`: the correlation between `y_ij` and `y_ji` over dyads.

An undefined ratio is `NaN`, as in sna. `missing=:error` refuses masked
dyads; `missing=:face` classifies them from their stored face values.

# Example
```julia
using SNA
net = network(3)
add_edges!(net, [(1, 2), (2, 1), (2, 3)])
grecip(net)                          # (1 + 1) / 3 = 0.6667
grecip(net; measure=:edgewise)       # 2 / 3
```
"""
function grecip(net::AbstractNetwork; measure::Symbol=:dyadic, missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="grecip")
    census = dyad_census(net; missing=policy)
    M, A, N = census.mutual, census.asymmetric, census.null
    if measure == :dyadic
        return (M + N) / (M + A + N)
    elseif measure == :dyadic_nonnull
        return M / (M + A)
    elseif measure == :edgewise
        return 2M / (2M + A)
    elseif measure == :edgewise_lrr
        return log(M * (M + A + N) / (M + A / 2)^2)
    elseif measure == :correlation
        return _grecip_correlation(nv(net), M, A, N)
    end
    throw(ArgumentError("grecip: measure must be :dyadic, :dyadic_nonnull, :edgewise, " *
                        ":edgewise_lrr or :correlation; got :$measure"))
end

# sna's dyadic correlation for binary data: the product-moment correlation
# of (y_ij, y_ji) pooled over unordered dyads, around the graph mean.
function _grecip_correlation(n::Integer, M::Integer, A::Integer, N::Integer)
    n, M, A, N = Int(n), Int(M), Int(A), Int(N)
    n < 2 && return NaN
    n == 2 && return A == 0 ? 1.0 : 0.0
    M + A == 0 && return 1.0
    nd = binomial(n, 2)
    ne = 2nd
    gm = (2M + A) / ne
    gv = ((2M + A) * (1 - gm)^2 + (ne - 2M - A) * gm^2) / (ne - 1)
    gv == 0 && return 1.0
    dsum = M * (1 - gm)^2 - A * gm * (1 - gm)
    return 2 * (dsum + gm^2 * N) / ((2nd - 1) * gv)
end

"""
    gtrans(net; measure=:weak, missing=:error) -> Float64

Graph transitivity: R `sna::gtrans`, over ordered triples of distinct
vertices `(i, j, k)` (loops are ignored).

- `:weak` (default): the proportion of two-paths `i→j→k` that are closed by
  `i→k`. On an undirected network this is `3 × triangles / connected
  triples`, igraph's global transitivity.
- `:strong`: the proportion of ordered triples for which `i→k` is present
  exactly when some two-path `i→j→k` is, counted per intermediary `j`.
- `:weakcensus` / `:strongcensus`: the numerators of the two ratios (counts).
- `:correlation`: the correlation between `y_ik` and the number of two-paths
  from `i` to `k`.

A ratio with no qualifying triples is `1.0`, as in sna. sna's `rank`
measure (for valued data) is not implemented.

`missing=:error` refuses masked dyads; `missing=:face` counts them at their
stored face values in numerator and denominator.

# Example
```julia
using SNA
net = network(3)
add_edges!(net, [(1, 2), (2, 3), (1, 3)])
gtrans(net)                          # 1.0: the only two-path 1→2→3 is closed
```
"""
function gtrans(net::AbstractNetwork; measure::Symbol=:weak, missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="gtrans")
    n = nv(net)
    if measure in (:weak, :weakcensus)
        potential = 0
        realized = 0
        for j in 1:n
            for i in inneighbors(net, j)
                i == j && continue
                for k in outneighbors(net, j)
                    (k == j || k == i) && continue
                    potential += 1
                    has_edge(net, i, k) && (realized += 1)
                end
            end
        end
        measure == :weakcensus && return Float64(realized)
        return potential > 0 ? realized / potential : 1.0
    elseif measure in (:strong, :strongcensus, :correlation)
        A = Float64.(_sociomatrix(net) .!= 0)
        S = A * A
        if measure == :correlation
            x = Float64[A[i, k] for i in 1:n, k in 1:n if i != k]
            y = Float64[S[i, k] for i in 1:n, k in 1:n if i != k]
            length(x) < 2 && return 1.0
            tv = var(x) * var(y)
            return tv == 0 ? Float64(x == y) : cor(x, y)
        end
        hits = 0.0
        for i in 1:n, k in 1:n
            i == k && continue
            hits += A[i, k] == 1 ? S[i, k] : n - 2 - S[i, k]
        end
        measure == :strongcensus && return hits
        total = binomial(n, 3) * 6
        return total > 0 ? hits / total : 1.0
    end
    throw(ArgumentError("gtrans: measure must be :weak, :strong, :weakcensus, :strongcensus " *
                        "or :correlation; got :$measure (sna's :rank is not implemented)"))
end

"""
    transitivity(net; type=:global, cmode=:total, missing=:error) -> Float64 or Vector{Float64}

Clustering coefficients, named after igraph's `transitivity` (sna has no
local clustering coefficient).

- `type=:global` (default): sna's weak transitivity, [`gtrans`](@ref)`(net)`.
  On an undirected network this is igraph's global transitivity. On a
  directed network igraph ignores direction; sna and this function do not.
- `type=:local`: a vector of local clustering coefficients. A vertex with
  fewer than two neighbours (loops excluded) has no coefficient and scores
  `NaN`, as in igraph.
    - Undirected networks: the Watts–Strogatz coefficient, the share of the
      `k(k−1)/2` pairs of a vertex's `k` neighbours that are tied.
    - Directed networks, `cmode=:total` (default): Fagiolo's (2007)
      total-degree coefficient,
      `[(A+Aᵀ)³]ᵢᵢ / (2[dᵢᵗᵒᵗ(dᵢᵗᵒᵗ−1) − 2dᵢ↔])`, the share of all directed
      triangles through `i` that are realized. It equals the Watts–Strogatz
      coefficient when every tie is mutual.
    - Directed networks, `cmode=:weak`: the Watts–Strogatz coefficient of
      the weak symmetrization (an arc in either direction is a tie), which is
      what igraph computes on a directed graph.
- `type=:average`: the mean of the local coefficients over vertices that
  have one (degree ≥ 2), igraph's `transitivity(type="average")`; `NaN` if no
  vertex has one.

`missing=:error` refuses masked dyads; `missing=:face` counts them at their
stored face values.

# Example
```julia
using SNA
net = network(4; directed=false)
add_edges!(net, [(1, 2), (1, 3), (2, 3), (3, 4)])
transitivity(net; type=:local)       # [1.0, 1.0, 0.3333, NaN]
transitivity(net; type=:average)     # 0.7778
```

# References
Fagiolo, G. (2007). Clustering in complex directed networks. *Physical
Review E*, 76, 026107.
"""
function transitivity(net::AbstractNetwork; type::Symbol=:global, cmode::Symbol=:total,
                      missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="transitivity")
    cmode in (:total, :weak) ||
        throw(ArgumentError("transitivity: cmode must be :total (Fagiolo) or :weak; got :$cmode"))
    if type == :global
        return gtrans(net; missing=policy)
    elseif type == :local
        return _local_clustering(net, cmode)
    elseif type == :average
        cc = filter(!isnan, _local_clustering(net, cmode))
        return isempty(cc) ? NaN : mean(cc)
    end
    throw(ArgumentError("type must be :global, :local, or :average"))
end

function _local_clustering(net, cmode::Symbol)
    n = nv(net)
    cc = fill(NaN, n)
    if is_directed(net) && cmode == :total
        A = Float64.(_sociomatrix(net) .!= 0)        # loops removed
        S = A + A'
        S3 = S * S * S
        A2 = A * A
        for i in 1:n
            dtot = sum(@view A[i, :]) + sum(@view A[:, i])
            den = 2 * (dtot * (dtot - 1) - 2 * A2[i, i])
            den > 0 && (cc[i] = S3[i, i] / den)
        end
        return cc
    end
    # Watts–Strogatz on the (weakly symmetrized) neighbourhoods
    nbrs = [_union_neighborhood(net, v) for v in 1:n]
    for v in 1:n
        nb = nbrs[v]
        k = length(nb)
        k < 2 && continue
        t = 0
        for a in 1:k, b in (a+1):k
            insorted(nb[b], nbrs[nb[a]]) && (t += 1)
        end
        cc[v] = t / (k * (k - 1) / 2)
    end
    return cc
end

"""
    dyad_census(net; missing=:error) -> NamedTuple

The dyad census, R `sna::dyad.census`: the numbers of mutual, asymmetric
and null dyads, as `(mutual=…, asymmetric=…, null=…)`. Loops are not dyads.

`missing=:error` refuses masked dyads; `missing=:face` classifies each from
its stored face value (there is no NA cell in this census).

# Example
```julia
using SNA
net = network(3)
add_edges!(net, [(1, 2), (2, 1), (2, 3)])
dyad_census(net)                     # (mutual = 1, asymmetric = 1, null = 1)
```
"""
function dyad_census(net::AbstractNetwork; missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="dyad_census")
    n = nv(net)
    mutual = 0
    asymmetric = 0
    null = 0

    for i in 1:n
        for j in (i+1):n
            has_ij = has_edge(net, i, j)
            has_ji = has_edge(net, j, i)

            if has_ij && has_ji
                mutual += 1
            elseif has_ij || has_ji
                asymmetric += 1
            else
                null += 1
            end
        end
    end

    return (mutual=mutual, asymmetric=asymmetric, null=null)
end

"""
    triad_census(net; missing=:error) -> Vector{Int}

The triad census, R `sna::triad.census`.

Masked (unobserved) dyads are rejected by default (`missing=:error`); pass
`missing=:face` to classify triads from the stored face values of their
masked dyads (see `NetworkCore.require_observed`).

For directed networks, returns the 16-element Davis–Leinhardt census over
the M-A-N (mutual/asymmetric/null) isomorphism classes, in the standard
order used by R `sna::triad.census`:

    003, 012, 102, 021D, 021U, 021C, 111D, 111U, 030T, 030C,
    201, 120D, 120U, 120C, 210, 300

For undirected networks, returns the 4-element census by triad edge count
(0, 1, 2, 3 edges).

Uses the edge-driven Batagelj–Mrvar (2001) algorithm: only triads containing
at least one tie are enumerated (each exactly once, from its lowest-labeled
connected pair), so the cost is `O(Σ_(u,v)∈E (deg(u)+deg(v)))` rather than
`O(n³)`; the empty-triad count is recovered by subtraction from `C(n,3)`.

# Example
```julia
using SNA
net = network(3)
add_edges!(net, [(1, 2), (2, 3), (1, 3)])
triad_census(net)[9]                 # 1: one 030T (transitive) triad
```
"""
function triad_census(net::AbstractNetwork; missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="triad_census")
    is_directed(net) || return _triad_census_undirected(net)

    n = nv(net)
    census = zeros(Int, 16)
    nbrs = [_union_neighborhood(net, v) for v in 1:n]
    buf = Int[]

    for v in 1:n
        for u in nbrs[v]
            u > v || continue
            # Third vertices attached to the dyad {v, u}
            _sorted_union!(buf, nbrs[v], nbrs[u], v, u)
            # Triads whose third vertex is attached to neither v nor u have
            # (v,u) as their only non-null dyad: 102 if mutual, 012 if asym
            dyad = (has_edge(net, v, u) && has_edge(net, u, v)) ? 3 : 2
            census[dyad] += n - length(buf) - 2
            for w in buf
                # Count each connected triad exactly once: from the edge to
                # its lowest-labeled attached pair
                if u < w || (v < w && w < u && !insorted(w, nbrs[v]))
                    census[_TRICODE_CLASS[_tricode(net, v, u, w)+1]] += 1
                end
            end
        end
    end

    census[1] = binomial(n, 3) - sum(@view census[2:16])
    return census
end

# Map from the 6-bit dyad code of a triad (v,u,w) to its 1-based position in
# the Davis–Leinhardt census order (Batagelj & Mrvar 2001, Table 1)
const _TRICODE_CLASS = (
    1, 2, 2, 3, 2, 4, 6, 8, 2, 6, 5, 7, 3, 8, 7, 11,
    2, 6, 4, 8, 5, 9, 9, 13, 6, 10, 9, 14, 7, 14, 12, 15,
    2, 5, 6, 7, 6, 9, 10, 14, 4, 9, 9, 12, 8, 13, 14, 15,
    3, 7, 8, 11, 7, 12, 14, 15, 8, 14, 13, 15, 11, 15, 15, 16)

# 6-bit arc code of the triad (v, u, w)
function _tricode(net, v::Integer, u::Integer, w::Integer)
    code = 0
    has_edge(net, v, u) && (code += 1)
    has_edge(net, u, v) && (code += 2)
    has_edge(net, v, w) && (code += 4)
    has_edge(net, w, v) && (code += 8)
    has_edge(net, u, w) && (code += 16)
    has_edge(net, w, u) && (code += 32)
    return code
end

# Sorted union of in- and out-neighbors of v, excluding v itself
function _union_neighborhood(net, v::Integer)
    if !is_directed(net)
        # Undirected storage is a symmetric digraph: out-neighbors suffice
        return Int[w for w in outneighbors(net, v) if w != v]
    end
    return _sorted_union!(Int[], inneighbors(net, v), outneighbors(net, v),
                          v, v)
end

# Merge two sorted vectors into `buf` (deduplicated), skipping two vertices
function _sorted_union!(buf::Vector{Int}, a::AbstractVector{<:Integer},
                        b::AbstractVector{<:Integer}, skip1::Int, skip2::Int)
    empty!(buf)
    i, j = 1, 1
    na, nb = length(a), length(b)
    while i <= na || j <= nb
        w = if i > na
            x = b[j]; j += 1; x
        elseif j > nb
            x = a[i]; i += 1; x
        elseif a[i] < b[j]
            x = a[i]; i += 1; x
        elseif a[i] > b[j]
            x = b[j]; j += 1; x
        else
            x = a[i]; i += 1; j += 1; x
        end
        (w == skip1 || w == skip2) || push!(buf, w)
    end
    return buf
end

# Edge-driven undirected census by triad edge count (0, 1, 2, 3)
function _triad_census_undirected(net)
    n = nv(net)
    census = zeros(Int, 4)
    nbrs = [_union_neighborhood(net, v) for v in 1:n]
    buf = Int[]

    for v in 1:n
        for u in nbrs[v]
            u > v || continue
            _sorted_union!(buf, nbrs[v], nbrs[u], v, u)
            census[2] += n - length(buf) - 2
            for w in buf
                if u < w || (v < w && w < u && !insorted(w, nbrs[v]))
                    m = 1 + (insorted(w, nbrs[v]) ? 1 : 0) +
                        (insorted(w, nbrs[u]) ? 1 : 0)
                    census[m+1] += 1
                end
            end
        end
    end

    census[1] = binomial(n, 3) - census[2] - census[3] - census[4]
    return census
end

"""
    hierarchy(net; measure=:reciprocity, missing=:error) -> Float64

Graph hierarchy: R `sna::hierarchy`.

- `:reciprocity` (default): `1 −` [`grecip`](@ref)`(net)` (dyadic reciprocity).
- `:krackhardt`: Krackhardt's hierarchy, the proportion of reachable
  (unordered) pairs whose reachability is one-way, i.e. `1 −` the
  `:dyadic_nonnull` reciprocity of the reachability matrix.

A network with no dyads, or (for `:krackhardt`) no reachable pair, gives
`NaN`, as in sna. `missing=:error` refuses masked dyads; `missing=:face`
uses their stored face values.

# Example
```julia
using SNA
chain = network(3)
add_edges!(chain, [(1, 2), (2, 3)])
hierarchy(chain; measure=:krackhardt)    # 1.0: every reachable pair is one-way
```
"""
function hierarchy(net::AbstractNetwork; measure::Symbol=:reciprocity, missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="hierarchy")
    if measure == :reciprocity
        return 1.0 - grecip(net; measure=:dyadic, missing=policy)
    elseif measure != :krackhardt
        throw(ArgumentError("Unknown hierarchy measure: $measure"))
    end
    n = nv(net)
    reach = reachability(net; missing=policy)
    hierarchical = 0
    total_connected = 0
    for i in 1:n, j in (i+1):n
        if reach[i, j] || reach[j, i]
            total_connected += 1
            reach[i, j] != reach[j, i] && (hierarchical += 1)
        end
    end
    return hierarchical / total_connected      # NaN with no reachable pair, as sna
end

"""
    efficiency(net; missing=:error) -> Float64

Krackhardt's efficiency: R `sna::efficiency`.

Efficiency is computed per weak component in arc (digraph) terms, exactly
as in sna: a component of size `nᶜ` requires `nᶜ − 1` arcs for weak
connection and can hold at most `nᶜ(nᶜ − 1)`; arcs beyond the requirement
are "excess", and efficiency is `1 − Σ excess / Σ maximum possible excess`.
An undirected edge counts as two arcs (sna treats symmetric data as a
digraph), and self-loops are ignored (sna's `diag=FALSE`). A network whose
components are all single vertices has no possible excess and gives `NaN`,
as in sna.

`missing=:error` refuses masked dyads; `missing=:face` counts their stored
face values as arcs.

# Example
```julia
using SNA
path = network(3; directed=false)
add_edges!(path, [(1, 2), (2, 3)])
efficiency(path)                     # 1 − (4 − 2) / (6 − 2) = 0.5
```
"""
function efficiency(net::AbstractNetwork; missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="efficiency")
    sizes = component_dist(net; connected=:weak, missing=policy).csize
    required = sum(s - 1 for s in sizes; init=0)
    max_excess = sum(s * (s - 1) - (s - 1) for s in sizes; init=0)
    arcs = 0
    for e in edges(net)
        src(e) == dst(e) && continue
        arcs += is_directed(net) ? 1 : 2
    end
    return 1.0 - (arcs - required) / max_excess
end

"""
    connectedness(net; missing=:error) -> Float64

Krackhardt's connectedness: the proportion of (unordered) vertex pairs
joined by a semipath, i.e. lying in the same weak component (R
`sna::connectedness`). A network of one vertex has connectedness 1, as in sna.

`missing=:error` refuses masked dyads; `missing=:face` takes their stored
face values as ties.

# Example
```julia
using SNA
net = network(4; directed=false)
add_edges!(net, [(1, 2), (3, 4)])
connectedness(net)                   # 2 of 6 pairs = 0.3333
```
"""
function connectedness(net::AbstractNetwork; missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="connectedness")
    n = nv(net)
    n <= 1 && return 1.0
    sizes = component_dist(net; connected=:weak, missing=policy).csize
    connected_pairs = sum(s * (s - 1) ÷ 2 for s in sizes)
    return connected_pairs / (n * (n - 1) ÷ 2)
end

"""
    mutuality(net; missing=:error) -> Int

The number of mutual (reciprocated) dyads, R `sna::mutuality` (a count, not
a proportion).

`missing=:error` refuses masked dyads; `missing=:face` classifies them from
their stored face values.

# Example
```julia
using SNA
net = network(3)
add_edges!(net, [(1, 2), (2, 1), (2, 3)])
mutuality(net)                       # 1
```
"""
function mutuality(net::AbstractNetwork; missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="mutuality")
    return dyad_census(net; missing=policy).mutual
end

"""
    reachability(net; missing=:error) -> Matrix{Bool}

The reachability matrix, R `sna::reachability`: entry `(i, j)` is `true`
when there is a directed path from `i` to `j`. Every vertex reaches itself.

`missing=:error` refuses masked dyads; `missing=:face` traverses them at
their stored face values.

# Example
```julia
using SNA
chain = network(3)
add_edges!(chain, [(1, 2), (2, 3)])
reachability(chain)[1, 3], reachability(chain)[3, 1]     # (true, false)
```
"""
function reachability(net::AbstractNetwork; missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="reachability")
    n = nv(net)
    reach = falses(n, n)

    for i in 1:n
        # BFS from vertex i (index-pointer queue keeps pops O(1))
        visited = falses(n)
        queue = [i]
        visited[i] = true
        head = 1

        while head <= length(queue)
            v = queue[head]
            head += 1
            reach[i, v] = true

            for w in outneighbors(net, v)
                if !visited[w]
                    visited[w] = true
                    push!(queue, w)
                end
            end
        end
    end

    return reach
end
