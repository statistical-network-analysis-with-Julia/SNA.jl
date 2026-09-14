"""
Structural and regular equivalence analysis.

Provides functions for analyzing structural positions in networks,
including structural equivalence, regular equivalence, and blockmodeling.
"""

using LinearAlgebra
using Statistics

"""
    structural_equivalence(net; method=:correlation, missing=:error) -> Matrix{Float64}

Compute structural equivalence between all pairs of vertices.

Two vertices are structurally equivalent if they have identical ties
to and from all other vertices.

# Arguments
- `net`: Network object
- `method::Symbol=:correlation`: Distance/similarity measure
    - `:correlation`: Pearson correlation of adjacency rows/columns
    - `:euclidean`: Euclidean distance
    - `:hamming`: Hamming distance (proportion of different ties)
- `missing::Symbol=:error`: Missing-dyad policy (`Networks.require_observed`);
  `:error` rejects a network with masked (unobserved) dyads, `:face` keeps
  each masked dyad in the tie profiles at its stored face value (there is no
  pairwise-deletion variant)

# Returns
Matrix of similarity/distance scores.
"""
function structural_equivalence(net; method::Symbol=:correlation,
                                ignore_eval::Bool=true, attr::Symbol=:weight,
                                missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="structural_equivalence")
    A = _sociomatrix(net; ignore_eval, attr)
    n = nv(net)

    # Create profile for each vertex: concatenate row and column
    profiles = Matrix{Float64}(undef, n, 2*n)
    for i in 1:n
        profiles[i, 1:n] = A[i, :]
        profiles[i, (n+1):end] = A[:, i]
    end

    result = zeros(n, n)

    if method == :correlation
        for i in 1:n
            for j in i:n
                if i == j
                    result[i, j] = 1.0
                else
                    c = cor(profiles[i, :], profiles[j, :])
                    result[i, j] = isnan(c) ? 0.0 : c
                    result[j, i] = result[i, j]
                end
            end
        end
    elseif method == :euclidean
        for i in 1:n
            for j in i:n
                result[i, j] = norm(profiles[i, :] - profiles[j, :])
                result[j, i] = result[i, j]
            end
        end
    elseif method == :hamming
        for i in 1:n
            for j in i:n
                result[i, j] = mean(profiles[i, :] .!= profiles[j, :])
                result[j, i] = result[i, j]
            end
        end
    end

    return result
end

"""
    regular_equivalence(net; max_iter=100, tol=1e-6, missing=:error) -> Matrix{Float64}

Compute regular equivalence between all pairs of vertices.

Two vertices are regularly equivalent if they have equivalent ties to
equivalent others (a recursive definition).

Note: this uses an iterative neighbor-matching similarity (CATREGE-style
best-match averaging over in- and out-neighborhoods). It is in the spirit
of White & Reitz regular equivalence but is **not** the Burt/White REGE
algorithm, so scores are not directly comparable to UCINET/R REGE output.

Masked (unobserved) dyads are rejected by default (`missing=:error`); pass
`missing=:face` to include them in the neighborhoods at their stored face
values (see `Networks.require_observed`).
"""
function regular_equivalence(net; max_iter::Int=100, tol::Float64=1e-6,
                             ignore_eval::Bool=true, attr::Symbol=:weight,
                             missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="regular_equivalence")
    ignore_eval || throw(ArgumentError("weighted regular equivalence is not implemented; use binary ties or structural_equivalence"))
    max_iter > 0 || throw(ArgumentError("max_iter must be positive"))
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
    for _ in 1:max_iter
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
    converged || @warn "Regular equivalence did not converge; increase max_iter"

    return sim
end

"""
    equiv_clust(net; method=:structural, k=nothing, missing=:error) -> Vector{Int}

Cluster vertices by equivalence using agglomerative hierarchical clustering
with average linkage (UPGMA) on the equivalence distance matrix, in the
spirit of R `sna::equiv.clust`.

# Arguments
- `net`: Network object
- `method::Symbol=:structural`: Equivalence type (:structural or :regular)
- `k::Union{Int,Nothing}=nothing`: Number of clusters (defaults to
  `max(2, n ÷ 4)`; values `≥ n` yield the trivial one-vertex-per-cluster
  solution)
- `missing::Symbol=:error`: Missing-dyad policy (`Networks.require_observed`);
  `:error` rejects a network with masked (unobserved) dyads, `:face` uses
  their stored face values in the underlying equivalence matrix

# Returns
Vector of cluster assignments in `1:k` for each vertex.
"""
function equiv_clust(net; method::Symbol=:structural,
                     ignore_eval::Bool=true, attr::Symbol=:weight,
                     k::Union{Int,Nothing}=nothing, missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="equiv_clust")
    # Compute distance matrix
    if method == :structural
        sim = structural_equivalence(net; method=:correlation, ignore_eval, attr, missing=policy)
        dist = 1 .- sim  # Convert similarity to distance
    else
        sim = regular_equivalence(net; ignore_eval, attr, missing=policy)
        dist = 1 .- sim
    end

    # Ensure non-negative and zero diagonal
    dist = max.(dist, 0.0)
    for i in 1:size(dist, 1)
        dist[i, i] = 0.0
    end

    n = nv(net)
    if isnothing(k)
        k = min(n, max(2, n ÷ 4))  # Heuristic: n/4 clusters
    end
    k = clamp(k, 1, n)

    return _hclust_average(dist, k)
end

# Agglomerative hierarchical clustering (average linkage), cut at k
# clusters. Returns assignments renumbered to 1:k.
function _hclust_average(dist::Matrix{Float64}, k::Int)
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

        # Average-linkage update: merge b into a
        na, nb = length(clusters[a]), length(clusters[b])
        for m in 1:n
            (active[m] && m != a && m != b) || continue
            d[a, m] = d[m, a] = (na * d[a, m] + nb * d[b, m]) / (na + nb)
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
    blockmodel(net; k::Int, method=:structural, missing=:error) -> NamedTuple

Create a blockmodel from equivalence-based clustering.

Masked (unobserved) dyads are rejected by default (`missing=:error`); pass
`missing=:face` to let them enter both the clustering and the block densities
at their stored face values (see `Networks.require_observed`). Block densities
are *not* renormalized over observed dyads only.

# Returns
NamedTuple with:
- `membership::Vector{Int}`: Block membership for each vertex
- `block_matrix::Matrix{Float64}`: Density of ties between blocks
- `n_blocks::Int`: Number of blocks
"""
function blockmodel(net; k::Int, method::Symbol=:structural,
                    ignore_eval::Bool=true, attr::Symbol=:weight,
                    missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="blockmodel")
    n = nv(net)
    membership = equiv_clust(net; method=method, k=k, ignore_eval, attr, missing=policy)
    A = _sociomatrix(net; ignore_eval, attr)
    # equiv_clust clamps k to at most n; size blocks by actual labels
    k = maximum(membership)

    # Compute block densities (diagonal blocks average over ordered
    # off-diagonal pairs within the block)
    block_matrix = zeros(k, k)
    block_counts = zeros(Int, k, k)

    for i in 1:n
        for j in 1:n
            if i != j
                bi, bj = membership[i], membership[j]
                block_matrix[bi, bj] += A[i, j]
                block_counts[bi, bj] += 1
            end
        end
    end

    # Normalize by block sizes
    for bi in 1:k
        for bj in 1:k
            if block_counts[bi, bj] > 0
                block_matrix[bi, bj] /= block_counts[bi, bj]
            end
        end
    end

    return (
        membership = membership,
        block_matrix = block_matrix,
        n_blocks = k
    )
end

"""
    consensus(clusterings::Vector{Vector{Int}}) -> Vector{Int}

Compute consensus clustering from multiple clustering solutions.
"""
function consensus(clusterings::Vector{Vector{Int}})
    n = length(clusterings[1])
    n_clusterings = length(clusterings)

    # Create co-occurrence matrix
    cooccur = zeros(n, n)
    for clustering in clusterings
        for i in 1:n
            for j in (i+1):n
                if clustering[i] == clustering[j]
                    cooccur[i, j] += 1
                    cooccur[j, i] += 1
                end
            end
        end
    end
    cooccur ./= n_clusterings

    # Cluster based on co-occurrence
    dist = 1 .- cooccur
    for i in 1:n
        dist[i, i] = 0.0
    end

    # Simple thresholding approach
    threshold = 0.5
    assignments = zeros(Int, n)
    current_cluster = 0

    for i in 1:n
        if assignments[i] == 0
            current_cluster += 1
            assignments[i] = current_cluster
            for j in (i+1):n
                if cooccur[i, j] >= threshold && assignments[j] == 0
                    assignments[j] = current_cluster
                end
            end
        end
    end

    return assignments
end
