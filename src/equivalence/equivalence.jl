"""
Structural and regular equivalence analysis.

Provides functions for analyzing structural positions in networks,
including structural equivalence, regular equivalence, and blockmodeling.
"""

using LinearAlgebra
using Statistics

"""
    structural_equivalence(net; method=:correlation) -> Matrix{Float64}

Compute structural equivalence between all pairs of vertices.

Two vertices are structurally equivalent if they have identical ties
to and from all other vertices.

# Arguments
- `net`: Network object
- `method::Symbol=:correlation`: Distance/similarity measure
    - `:correlation`: Pearson correlation of adjacency rows/columns
    - `:euclidean`: Euclidean distance
    - `:hamming`: Hamming distance (proportion of different ties)

# Returns
Matrix of similarity/distance scores.
"""
function structural_equivalence(net; method::Symbol=:correlation)
    A = as_matrix(net)
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
    regular_equivalence(net; max_iter=100, tol=1e-6) -> Matrix{Float64}

Compute regular equivalence between all pairs of vertices.

Two vertices are regularly equivalent if they have equivalent ties to
equivalent others (a recursive definition).

Note: this uses an iterative neighbor-matching similarity (CATREGE-style
best-match averaging over in- and out-neighborhoods). It is in the spirit
of White & Reitz regular equivalence but is **not** the Burt/White REGE
algorithm, so scores are not directly comparable to UCINET/R REGE output.
"""
function regular_equivalence(net; max_iter::Int=100, tol::Float64=1e-6)
    n = nv(net)
    A = as_matrix(net)

    # Initialize similarity matrix
    sim = ones(n, n)
    for i in 1:n
        sim[i, i] = 1.0
    end

    for _ in 1:max_iter
        sim_old = copy(sim)

        for i in 1:n
            for j in (i+1):n
                # Compute similarity based on neighbor similarities
                out_i = findall(!iszero, A[i, :])
                out_j = findall(!iszero, A[j, :])
                in_i = findall(!iszero, A[:, i])
                in_j = findall(!iszero, A[:, j])

                # Match out-neighbors
                out_sim = 0.0
                if !isempty(out_i) && !isempty(out_j)
                    for k in out_i
                        out_sim += maximum(sim_old[k, l] for l in out_j)
                    end
                    out_sim /= length(out_i)
                elseif isempty(out_i) && isempty(out_j)
                    out_sim = 1.0
                end

                # Match in-neighbors
                in_sim = 0.0
                if !isempty(in_i) && !isempty(in_j)
                    for k in in_i
                        in_sim += maximum(sim_old[k, l] for l in in_j)
                    end
                    in_sim /= length(in_i)
                elseif isempty(in_i) && isempty(in_j)
                    in_sim = 1.0
                end

                sim[i, j] = (out_sim + in_sim) / 2
                sim[j, i] = sim[i, j]
            end
        end

        if maximum(abs.(sim - sim_old)) < tol
            break
        end
    end

    return sim
end

"""
    equiv_clust(net; method=:structural, k=nothing) -> Vector{Int}

Cluster vertices by equivalence using agglomerative hierarchical clustering
with average linkage (UPGMA) on the equivalence distance matrix, in the
spirit of R `sna::equiv.clust`.

# Arguments
- `net`: Network object
- `method::Symbol=:structural`: Equivalence type (:structural or :regular)
- `k::Union{Int,Nothing}=nothing`: Number of clusters (defaults to
  `max(2, n ÷ 4)`; values `≥ n` yield the trivial one-vertex-per-cluster
  solution)

# Returns
Vector of cluster assignments in `1:k` for each vertex.
"""
function equiv_clust(net; method::Symbol=:structural, k::Union{Int,Nothing}=nothing)
    # Compute distance matrix
    if method == :structural
        sim = structural_equivalence(net; method=:correlation)
        dist = 1 .- sim  # Convert similarity to distance
    else
        sim = regular_equivalence(net)
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
    blockmodel(net; k::Int, method=:structural) -> NamedTuple

Create a blockmodel from equivalence-based clustering.

# Returns
NamedTuple with:
- `membership::Vector{Int}`: Block membership for each vertex
- `block_matrix::Matrix{Float64}`: Density of ties between blocks
- `n_blocks::Int`: Number of blocks
"""
function blockmodel(net; k::Int, method::Symbol=:structural)
    n = nv(net)
    membership = equiv_clust(net; method=method, k=k)
    A = as_matrix(net)
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
