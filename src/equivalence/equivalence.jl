"""
Structural and regular equivalence: `sedist` (R `sna::sedist`),
`regular_equivalence`, `equiv_clust`, `blockmodel` and
`consensus_clustering`.
"""

using LinearAlgebra
using Statistics

"""
    sedist(net; method=:hamming, diag=false, ignore_eval=true, attr=:weight,
           missing=:error) -> Matrix{Float64}

Structural equivalence between every pair of vertices, R `sna::sedist`.

Vertex `i`'s profile is its row and its column of the adjacency matrix, the
ties it sends followed by the ties it receives. As in sna (and Wasserman &
Faust, eq. 9.4), the comparison of `i` and `j` leaves out the cells that are
not about third parties: `i`'s and `j`'s self cells and the ties between `i`
and `j` (profile positions `i`, `j`, `n+i`, `n+j`). Two actors with identical
ties to and from everyone else are therefore equivalent whether or not they
are tied to each other. `diag=true` keeps every cell.

- `:hamming` (default): the number of differing cells (for valued data, the
  sum of absolute differences), sna's default;
- `:correlation`: the Pearson correlation of the two profiles (a
  similarity; `0` where it is undefined, as in sna, including the diagonal
  entry of a vertex with a constant profile);
- `:euclidean`: the Euclidean distance between the profiles;
- `:gamma`: Goodman–Kruskal gamma, `(agreements − disagreements) /
  (agreements + disagreements)`;
- `:exact`: `1` if the profiles differ anywhere, else `0`.

`missing=:error` refuses masked dyads (sna would drop them pairwise);
`missing=:face` keeps each masked dyad in the profiles at its stored face
value.

# Example
```julia
using SNA
net = network(4)
add_edges!(net, [(1, 2), (2, 1), (1, 3), (2, 3), (4, 1), (4, 2)])
sedist(net)[1, 2]                        # 0.0: 1 and 2 are structurally equivalent
sedist(net; method=:correlation)[1, 2]   # 1.0
```
"""
function sedist(net::AbstractNetwork; method::Symbol=:hamming, diag::Bool=false,
                ignore_eval::Bool=true, attr::Symbol=:weight, missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="sedist")
    method in (:hamming, :correlation, :euclidean, :gamma, :exact) ||
        throw(ArgumentError("sedist: method must be :hamming, :correlation, :euclidean, " *
                            ":gamma or :exact; got :$method"))
    A = _sociomatrix(net; ignore_eval, attr, diag=true)
    n = size(A, 1)
    profile(i) = vcat(A[i, :], A[:, i])
    P = [profile(i) for i in 1:n]
    result = zeros(n, n)
    keep = trues(2n)
    for i in 1:n, j in i:n
        fill!(keep, true)
        if !diag
            keep[i] = keep[n+i] = keep[j] = keep[n+j] = false
        end
        x = P[i][keep]
        y = P[j][keep]
        value = if method == :hamming
            sum(abs.(x .- y); init=0.0)
        elseif method == :euclidean
            sqrt(sum(abs2.(x .- y); init=0.0))
        elseif method == :correlation
            c = length(x) < 2 ? NaN : cor(x, y)
            isnan(c) ? 0.0 : c
        elseif method == :gamma
            agree = count(x .== y)
            (2agree - length(x)) / length(x)
        else
            Float64(x != y)
        end
        result[i, j] = result[j, i] = value
    end
    return result
end

"""
    regular_equivalence(net; maxiter=100, tol=1e-6, missing=:error) -> Matrix{Float64}

Regular equivalence similarity between every pair of vertices, in `[0, 1]`.

Two vertices are regularly equivalent if they have equivalent ties to
equivalent others (a recursive definition). The similarity is computed by
iterative neighbour matching (CATREGE-style best-match averaging over in-
and out-neighbourhoods), starting from all pairs equivalent, so it converges
to the *maximal* regular equivalence: on a network in which every actor
both sends and receives ties it can stay at 1 for every pair. It is in the
spirit of White & Reitz regular equivalence but is **not** the Burt/White
REGE algorithm, so scores are not directly comparable to UCINET's REGE. sna
has no regular-equivalence function.

`maxiter` bounds the iterations (a warning is issued if `tol` is not
reached). Weighted ties are not supported (`ignore_eval=false` raises an
`ArgumentError`).

`missing=:error` refuses masked dyads; `missing=:face` includes them in the
neighbourhoods at their stored face values.

# Example
```julia
using SNA
# a boss, two managers, four workers
tree = network(7)
add_edges!(tree, [(1, 2), (1, 3), (2, 4), (2, 5), (3, 6), (3, 7)])
sim = regular_equivalence(tree)
sim[2, 3], sim[4, 7], sim[1, 4]          # (1.0, 1.0, 0.0)
```
"""
function regular_equivalence(net::AbstractNetwork; maxiter::Int=100, tol::Float64=1e-6,
                             ignore_eval::Bool=true, attr::Symbol=:weight,
                             missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="regular_equivalence")
    ignore_eval || throw(ArgumentError("weighted regular equivalence is not implemented; use binary ties or sedist"))
    maxiter > 0 || throw(ArgumentError("maxiter must be positive"))
    tol > 0 || throw(ArgumentError("tol must be positive"))
    n = nv(net)
    sim = ones(n, n)
    n <= 1 && return sim
    A = _sociomatrix(net)
    out = [findall(!iszero, A[i, :]) for i in 1:n]
    inc = [findall(!iszero, A[:, i]) for i in 1:n]
    function match_neighbors(a, b, old)
        isempty(a) && isempty(b) && return 1.0
        (isempty(a) || isempty(b)) && return 0.0
        # Compare both directions so actor ordering cannot change similarity.
        ab = sum(maximum(old[k, l] for l in b) for k in a) / length(a)
        ba = sum(maximum(old[k, l] for k in a) for l in b) / length(b)
        return (ab + ba) / 2
    end
    converged = false
    for _ in 1:maxiter
        old = copy(sim)
        for i in 1:n, j in (i+1):n
            value = (match_neighbors(out[i], out[j], old) +
                     match_neighbors(inc[i], inc[j], old)) / 2
            sim[i, j] = sim[j, i] = value
        end
        if maximum(abs.(sim - old)) < tol
            converged = true
            break
        end
    end
    converged || @warn "Regular equivalence did not converge; increase maxiter"

    return sim
end

"""
    equiv_clust(net; k=nothing, equiv_fun=sedist, method=:hamming,
                cluster_method=:complete, missing=:error) -> Vector{Int}

Cluster the vertices into `k` positions by equivalence, R
`sna::equiv.clust` followed by `cutree(k)`: agglomerative hierarchical
clustering of an equivalence distance, cut at `k` clusters.

With the defaults this is what `cutree(equiv.clust(dat)\$cluster, k)`
returns in R: the [`sedist`](@ref) Hamming distance, clustered by complete
linkage, with ties between equally close pairs broken as R's `hclust` breaks
them (the pair with the lowest vertex first) and clusters numbered by their
first vertex.

# Arguments
- `k`: the number of clusters. Defaults to `max(2, n ÷ 4)`; `k ≥ n` gives
  one vertex per cluster. (R returns the tree and leaves the cut to
  `cutree`; SNA.jl returns the labels of the cut.)
- `equiv_fun`: the equivalence function, `sedist` (sna's default) or
  [`regular_equivalence`](@ref), whose similarity `s` is clustered as the
  distance `1 − s`.
- `method`: the `sedist` method (sna's `method`): `:hamming` (default),
  `:euclidean`, `:exact`, or the similarities `:correlation` and `:gamma`,
  which are negated into dissimilarities as sna's `code.diss=TRUE` does.
- `cluster_method`: the linkage, `:complete` (sna's default), `:average`
  (UPGMA) or `:single`.
- `missing=:error`: masked dyads are refused; `missing=:face` uses their
  stored face values in the equivalence matrix.

Average linkage on the profile correlation is `equiv_clust(net;
method=:correlation, cluster_method=:average)`, and clustering on regular
equivalence is `equiv_clust(net; equiv_fun=regular_equivalence)`.

# Example
```julia
using SNA
net = network(6)
add_edges!(net, [(1, 3), (1, 4), (2, 3), (2, 4), (3, 5), (4, 6)])
equiv_clust(net; k=3)                # [1, 1, 2, 2, 3, 3]
equiv_clust(net; k=3, method=:correlation, cluster_method=:average)
```
"""
function equiv_clust(net::AbstractNetwork; k::Union{Int,Nothing}=nothing,
                     equiv_fun::Function=sedist, method::Symbol=:hamming,
                     cluster_method::Symbol=:complete,
                     ignore_eval::Bool=true, attr::Symbol=:weight,
                     missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="equiv_clust")
    cluster_method in (:complete, :average, :single) ||
        throw(ArgumentError("equiv_clust: cluster_method must be :complete, :average or " *
                            ":single; got :$cluster_method"))
    dist = if equiv_fun === sedist
        d = sedist(net; method, ignore_eval, attr, missing=policy)
        method in (:correlation, :gamma) ? -d : d      # sna's code.diss=TRUE
    elseif equiv_fun === regular_equivalence
        1 .- regular_equivalence(net; ignore_eval, attr, missing=policy)
    else
        throw(ArgumentError("equiv_clust: equiv_fun must be sedist or regular_equivalence"))
    end

    n = nv(net)
    if isnothing(k)
        k = min(n, max(2, n ÷ 4))  # Heuristic: n/4 clusters
    end
    k = clamp(k, 1, n)

    return _hclust(Matrix{Float64}(dist), k, cluster_method)
end

# Agglomerative hierarchical clustering cut at k clusters, with R hclust's
# conventions: the closest pair is the first minimum in (i, j) order with
# i < j, the merged cluster keeps the lower index, and the Lance–Williams
# update is the maximum (complete), minimum (single) or size-weighted mean
# (average). Labels are numbered by each cluster's first vertex, as cutree.
function _hclust(dist::Matrix{Float64}, k::Int, linkage::Symbol)
    n = size(dist, 1)
    clusters = [[i] for i in 1:n]
    d = copy(dist)  # inter-cluster distances, indexed like `clusters`
    active = trues(n)

    n_active = n
    while n_active > k
        # Find the closest active pair
        best = (0, 0)
        best_d = Inf
        for i in 1:n
            active[i] || continue
            for j in (i+1):n
                active[j] || continue
                if d[i, j] < best_d
                    best_d = d[i, j]
                    best = (i, j)
                end
            end
        end
        a, b = best

        na, nb = length(clusters[a]), length(clusters[b])
        for m in 1:n
            (active[m] && m != a && m != b) || continue
            d[a, m] = d[m, a] = linkage === :complete ? max(d[a, m], d[b, m]) :
                                linkage === :single ? min(d[a, m], d[b, m]) :
                                (na * d[a, m] + nb * d[b, m]) / (na + nb)
        end
        append!(clusters[a], clusters[b])
        active[b] = false
        n_active -= 1
    end

    assignments = zeros(Int, n)
    label = 0
    for i in 1:n
        active[i] || continue
        label += 1
        for v in clusters[i]
            assignments[v] = label
        end
    end
    return assignments
end

"""
    blockmodel(net; k, equiv_fun=sedist, method=:hamming,
               cluster_method=:complete, missing=:error) -> NamedTuple

A blockmodel, R `sna::blockmodel(dat, equiv.clust(dat), k=k)`: the `k`
positions found by [`equiv_clust`](@ref) (whose keywords are passed on) and
the density of ties between them.

Returns `(membership, block_matrix, n_blocks)`: the position of each vertex,
the `n_blocks × n_blocks` matrix of block densities and the number of
positions. As in sna, a diagonal block's density is taken over the ordered
pairs of distinct vertices in it, so a position with a single vertex has an
undefined (`NaN`) diagonal density. (R's `block.membership` lists the same
positions in dendrogram order; `membership` is in vertex order.)

`missing=:error` refuses masked dyads; `missing=:face` lets them enter the
clustering and the densities at their stored face values (densities are not
renormalized over observed dyads).

# Example
```julia
using SNA
net = network(4)
add_edges!(net, [(1, 2), (2, 1), (3, 4), (4, 3)])
bm = blockmodel(net; k=2)
bm.block_matrix
```
"""
function blockmodel(net::AbstractNetwork; k::Int, equiv_fun::Function=sedist,
                    method::Symbol=:hamming, cluster_method::Symbol=:complete,
                    ignore_eval::Bool=true, attr::Symbol=:weight,
                    missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="blockmodel")
    n = nv(net)
    membership = equiv_clust(net; k, equiv_fun, method, cluster_method, ignore_eval, attr,
                             missing=policy)
    A = _sociomatrix(net; ignore_eval, attr)
    # equiv_clust clamps k to at most n; size blocks by actual labels
    k = maximum(membership; init=0)

    block_matrix = zeros(k, k)
    block_counts = zeros(Int, k, k)
    for i in 1:n, j in 1:n
        i == j && continue
        bi, bj = membership[i], membership[j]
        block_matrix[bi, bj] += A[i, j]
        block_counts[bi, bj] += 1
    end
    block_matrix ./= block_counts          # 0/0 = NaN for a singleton's own block, as sna

    return (membership = membership, block_matrix = block_matrix, n_blocks = k)
end

"""
    consensus_clustering(clusterings; threshold=0.5) -> Vector{Int}

Combine several partitions of the same vertices into one. Two vertices are
linked when they share a cluster in at least a fraction `threshold` of the
partitions, and the consensus clusters are the connected components of
that graph (so the result does not depend on the vertex order and is
transitive). Clusters are numbered by their lowest vertex.

This is not R `sna::consensus`, which builds a consensus *network* from
several informants' reports, so the name `consensus` is not used.

# Example
```julia
using SNA
consensus_clustering([[1, 1, 2, 2], [1, 1, 2, 3], [1, 2, 2, 2]])   # [1, 1, 2, 2]
```
"""
function consensus_clustering(clusterings::AbstractVector{<:AbstractVector{<:Integer}};
                              threshold::Real=0.5)
    isempty(clusterings) && throw(ArgumentError("consensus_clustering needs at least one clustering"))
    n = length(first(clusterings))
    all(c -> length(c) == n, clusterings) ||
        throw(ArgumentError("every clustering must label the same number of vertices"))
    0 < threshold <= 1 || throw(ArgumentError("threshold must be in (0, 1]"))
    g = Graphs.SimpleGraph(n)
    for i in 1:n, j in (i+1):n
        shared = count(c -> c[i] == c[j], clusterings)
        shared >= threshold * length(clusterings) && Graphs.add_edge!(g, i, j)
    end
    groups = sort!(Graphs.connected_components(g); by=minimum)
    labels = zeros(Int, n)
    for (c, members) in enumerate(groups), v in members
        labels[v] = c
    end
    return labels
end
