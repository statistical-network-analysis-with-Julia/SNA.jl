"""
Centrality measures for network analysis.

Provides vertex-level centrality measures including degree, betweenness,
closeness, eigenvector, and Bonacich power centrality.
"""

using Graphs
using LinearAlgebra

"""
    degree_centrality(net; mode=:total, normalized=false) -> Vector{Float64}

Compute degree centrality for all vertices.

# Arguments
- `net`: Network object
- `mode::Symbol=:total`: Type of degree (:in, :out, or :total)
- `normalized::Bool=false`: Normalize by maximum possible degree

# Returns
- Vector of centrality scores, one per vertex
"""
function degree_centrality(net; mode::Symbol=:total, normalized::Bool=false)
    n = nv(net)
    centrality = zeros(Float64, n)

    for v in vertices(net)
        if mode == :in
            centrality[v] = Float64(length(inneighbors(net, v)))
        elseif mode == :out
            centrality[v] = Float64(length(outneighbors(net, v)))
        else  # :total
            centrality[v] = Float64(length(inneighbors(net, v)) + length(outneighbors(net, v)))
        end
    end

    if normalized && n > 1
        # :in/:out can reach at most n-1; :total (Freeman degree) counts
        # both directions, so its maximum is 2(n-1)
        max_degree = mode == :total ? 2 * (n - 1) : (n - 1)
        centrality ./= max_degree
    end

    return centrality
end

"""
    betweenness_centrality(net; normalized=false) -> Vector{Float64}

Compute betweenness centrality for all vertices.

Betweenness centrality measures the extent to which a vertex lies on paths
between other vertices. The default is the *raw* (unnormalized) score,
matching R `sna::betweenness(rescale=FALSE)`; pass `normalized=true` for
scores scaled to [0, 1].
"""
function betweenness_centrality(net; normalized::Bool=false)
    # Use Graphs.jl implementation
    bc = Graphs.betweenness_centrality(net.graph; normalize=normalized)
    # Graphs counts each undirected path once per direction on the
    # digraph storage; halve to match undirected betweenness
    if !is_directed(net) && !normalized
        bc ./= 2
    end
    return bc
end

"""
    closeness_centrality(net; normalized=true) -> Vector{Float64}

Compute closeness centrality for all vertices: `(n-1)` over the total
geodesic distance to all other vertices (with `normalized=true`, the
Freeman closeness used by R `sna::closeness` for connected graphs).
Unreachable vertices are handled Graphs.jl-style (component-based scaling),
which differs from `sna`'s default of treating the score as undefined.
"""
function closeness_centrality(net; normalized::Bool=true)
    cc = Graphs.closeness_centrality(net.graph; normalize=normalized)
    return cc
end

"""
    eigenvector_centrality(net; max_iter=1000, tol=1e-10) -> Vector{Float64}

Compute eigenvector centrality for all vertices.

A vertex has high eigenvector centrality if it is connected to other
high-centrality vertices. As in R `sna::evcent`, the centrality is the
principal (right) eigenvector of the adjacency matrix — for directed
networks this weights vertices by the centrality of the vertices pointing
*at* their out-neighbors' pattern; symmetrize the network first if you
want the undirected notion. Scores are reported with non-negative
orientation and unit L2 norm.
"""
function eigenvector_centrality(net; max_iter::Int=1000, tol::Float64=1e-10)
    n = nv(net)
    A = as_matrix(net)

    # Power iteration
    x = ones(n) / sqrt(n)
    for _ in 1:max_iter
        x_new = A * x
        norm_x = norm(x_new)
        if norm_x > 0
            x_new ./= norm_x
        end

        if norm(x_new - x) < tol
            x = x_new
            break
        end
        x = x_new
    end

    # Perron orientation: report non-negative scores
    return abs.(x)
end

"""
    bonacich_power(net; exponent=1.0, rescale=false, tol=1e-7) -> Vector{Float64}

Compute Bonacich power centrality, following R `sna::bonpow`:

    c = α (I − β A)⁻¹ A 𝟙,   with α chosen so that Σᵢ cᵢ² = n

# Arguments
- `net`: Network object
- `exponent::Float64=1.0`: The attenuation/decay rate β. Must satisfy
  `|β| < 1/λ_max` for the underlying series to converge; a positive β
  rewards being connected to well-connected others, a negative β rewards
  being connected to poorly-connected others.
- `rescale::Bool=false`: If true, rescale so scores sum to 1 (as in sna)
- `tol::Float64=1e-7`: Solver tolerance for detecting singularity
"""
function bonacich_power(net; exponent::Float64=1.0, rescale::Bool=false,
                        tol::Float64=1e-7)
    n = nv(net)
    A = as_matrix(net)
    I_mat = Matrix{Float64}(I, n, n)

    # c = (I - βA)^(-1) * A * 1, scaled so that Σc² = n (sna's α)
    M = I_mat - exponent * A
    F = lu(M; check=false)
    if !issuccess(F) || abs(det(F)) < tol
        @warn "Bonacich power: (I − βA) is singular or near-singular; " *
              "choose |exponent| < 1/λ_max"
        return fill(NaN, n)
    end
    c = F \ (A * ones(n))

    ssq = sum(abs2, c)
    if ssq > 0
        c .*= sqrt(n / ssq)
    end
    if rescale
        c ./= sum(c)
    end

    return c
end

"""
    katz_centrality(net; α=0.1, β=1.0) -> Vector{Float64}

Compute Katz centrality.

Similar to eigenvector centrality but with damping factor α.
"""
function katz_centrality(net; α::Float64=0.1, β::Float64=1.0)
    return Graphs.katz_centrality(net.graph, α)
end

"""
    pagerank(net; α=0.85, max_iter=100, tol=1e-6) -> Vector{Float64}

Compute PageRank centrality.
"""
function pagerank(net; α::Float64=0.85, max_iter::Int=100, tol::Float64=1e-6)
    return Graphs.pagerank(net.graph, α, max_iter, tol)
end

"""
    flowbet(net) -> Vector{Float64}

Compute Freeman flow betweenness centrality (R `sna::flowbet`):

    f(v) = Σ_{i,j ≠ v} [maxflow(i → j) − maxflow(i → j | v removed)]

using edge capacities from the adjacency matrix. Pairs are ordered for
directed networks and unordered for undirected networks. Raw
(unnormalized) scores are returned, matching sna's default.
"""
function flowbet(net)
    n = nv(net)
    A = as_matrix(net)
    directed = is_directed(net)

    fb = zeros(n)
    for i in 1:n
        j_range = directed ? (1:n) : (i+1:n)
        for j in j_range
            i == j && continue
            base = _maxflow(A, n, i, j, 0)
            base == 0.0 && continue
            for v in 1:n
                (v == i || v == j) && continue
                fb[v] += base - _maxflow(A, n, i, j, v)
            end
        end
    end

    return fb
end

# Edmonds–Karp max flow on a dense capacity matrix, optionally with one
# vertex excluded (0 = none). O(V·E²); fine at research scale.
function _maxflow(cap::Matrix{Float64}, n::Int, s::Int, t::Int, excluded::Int)
    flow = zeros(n, n)
    total = 0.0

    parent = zeros(Int, n)
    while true
        # BFS for an augmenting path in the residual graph
        fill!(parent, 0)
        parent[s] = s
        queue = [s]
        head = 1
        while head <= length(queue) && parent[t] == 0
            u = queue[head]
            head += 1
            for w in 1:n
                if parent[w] == 0 && w != excluded && cap[u, w] - flow[u, w] > 1e-12
                    parent[w] = u
                    push!(queue, w)
                end
            end
        end
        parent[t] == 0 && break

        # Bottleneck capacity along the path
        aug = Inf
        w = t
        while w != s
            u = parent[w]
            aug = min(aug, cap[u, w] - flow[u, w])
            w = u
        end

        # Augment
        w = t
        while w != s
            u = parent[w]
            flow[u, w] += aug
            flow[w, u] -= aug
            w = u
        end
        total += aug
    end

    return total
end
