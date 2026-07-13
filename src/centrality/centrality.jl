"""
Centrality measures for network analysis.

Provides vertex-level centrality measures including degree, betweenness,
closeness, eigenvector, and Bonacich power centrality.
"""

using Graphs
using LinearAlgebra

"""
    degree_centrality(net; mode=:total, normalized=false, missing=:error) -> Vector{Float64}

Compute degree centrality for all vertices.

For undirected networks each edge contributes 1 to the degree of its two
endpoints (single-counted, matching R `sna::degree(gmode="graph")`); the
`mode` argument is ignored since in-, out-, and total degree coincide.

# Arguments
- `net`: Network object
- `mode::Symbol=:total`: Type of degree (:in, :out, or :total; directed only)
- `normalized::Bool=false`: Normalize by maximum possible degree
- `missing::Symbol=:error`: Missing-dyad policy (`Networks.require_observed`);
  `:error` rejects a network with masked (unobserved) dyads, `:face` counts
  each masked dyad at its stored face value

# Returns
- Vector of centrality scores, one per vertex
"""
function degree_centrality(net::AbstractNetwork; mode::Symbol=:total,
                           normalized::Bool=false, missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="degree_centrality")
    n = nv(net)
    centrality = zeros(Float64, n)
    directed = is_directed(net)

    for v in vertices(net)
        if !directed
            # Undirected edges are stored symmetrically in the backing
            # digraph; count each neighbor once (sna gmode="graph")
            centrality[v] = Float64(length(outneighbors(net, v)))
        elseif mode == :in
            centrality[v] = Float64(length(inneighbors(net, v)))
        elseif mode == :out
            centrality[v] = Float64(length(outneighbors(net, v)))
        else  # :total
            centrality[v] = Float64(length(inneighbors(net, v)) + length(outneighbors(net, v)))
        end
    end

    if normalized && n > 1
        # :in/:out and undirected degree can reach at most n-1; directed
        # :total (Freeman degree) counts both directions, max 2(n-1)
        max_degree = (directed && mode == :total) ? 2 * (n - 1) : (n - 1)
        centrality ./= max_degree
    end

    return centrality
end

"""
    betweenness_centrality(net; normalized=false, missing=:error) -> Vector{Float64}

Compute betweenness centrality for all vertices.

Betweenness centrality measures the extent to which a vertex lies on paths
between other vertices. The default is the *raw* (unnormalized) score,
matching R `sna::betweenness(rescale=FALSE)`; pass `normalized=true` for
scores scaled to [0, 1].

Masked (unobserved) dyads are rejected by default (`missing=:error`); pass
`missing=:face` to build the paths from their stored face values
(see `Networks.require_observed`).
"""
function betweenness_centrality(net::AbstractNetwork; normalized::Bool=false,
                                missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="betweenness_centrality")
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
    closeness_centrality(net; normalized=true, missing=:error) -> Vector{Float64}

Compute closeness centrality for all vertices: `(n-1)` over the total
geodesic distance to all other vertices (with `normalized=true`, the
Freeman closeness used by R `sna::closeness` for connected graphs).
Unreachable vertices are handled Graphs.jl-style (component-based scaling),
which differs from `sna`'s default of treating the score as undefined.

Masked (unobserved) dyads are rejected by default (`missing=:error`); pass
`missing=:face` to build the paths from their stored face values (see
`Networks.require_observed`).
"""
function closeness_centrality(net::AbstractNetwork; normalized::Bool=true,
                              missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="closeness_centrality")
    cc = Graphs.closeness_centrality(net.graph; normalize=normalized)
    return cc
end

"""
    eigenvector_centrality(net; max_iter=1000, tol=1e-10, missing=:error) -> Vector{Float64}

Compute eigenvector centrality for all vertices.

A vertex has high eigenvector centrality if it is connected to other
high-centrality vertices. As in R `sna::evcent`, the centrality is the
principal (right) eigenvector of the adjacency matrix — for directed
networks this weights vertices by the centrality of the vertices pointing
*at* their out-neighbors' pattern; symmetrize the network first if you
want the undirected notion. Scores are reported with non-negative
orientation and unit L2 norm.

Masked (unobserved) dyads are rejected by default (`missing=:error`); pass
`missing=:face` to read them from the adjacency matrix at their stored face
value (see `Networks.require_observed`).
"""
function eigenvector_centrality(net::AbstractNetwork; max_iter::Int=1000,
                                tol::Float64=1e-10, missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="eigenvector_centrality")
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
    bonacich_power(net; exponent=1.0, rescale=false, tol=1e-7, missing=:error) -> Vector{Float64}

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
- `missing::Symbol=:error`: Missing-dyad policy (`Networks.require_observed`);
  `:error` rejects a network with masked (unobserved) dyads, `:face` reads
  each masked dyad from the adjacency matrix at its stored face value
"""
function bonacich_power(net; exponent::Float64=1.0, rescale::Bool=false,
                        tol::Float64=1e-7, missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="bonacich_power")
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
    katz_centrality(net; α=0.1, β=1.0, missing=:error) -> Vector{Float64}

Compute Katz centrality.

Similar to eigenvector centrality but with damping factor α.

Masked (unobserved) dyads are rejected by default (`missing=:error`); pass
`missing=:face` to use their stored face values (see
`Networks.require_observed`).
"""
function katz_centrality(net::AbstractNetwork; α::Float64=0.1, β::Float64=1.0,
                         missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="katz_centrality")
    return Graphs.katz_centrality(net.graph, α)
end

"""
    pagerank(net; α=0.85, max_iter=100, tol=1e-6, missing=:error) -> Vector{Float64}

Compute PageRank centrality.

Masked (unobserved) dyads are rejected by default (`missing=:error`); pass
`missing=:face` to use their stored face values (see
`Networks.require_observed`).
"""
function pagerank(net::AbstractNetwork; α::Float64=0.85, max_iter::Int=100,
                  tol::Float64=1e-6, missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="pagerank")
    return Graphs.pagerank(net.graph, α, max_iter, tol)
end

"""
    flowbet(net; missing=:error) -> Vector{Float64}

Compute Freeman flow betweenness centrality (R `sna::flowbet`):

    f(v) = Σ_{i,j ≠ v} [maxflow(i → j) − maxflow(i → j | v removed)]

using edge capacities from the adjacency matrix. Pairs are ordered for
directed networks and unordered for undirected networks. Raw
(unnormalized) scores are returned, matching sna's default.

Masked (unobserved) dyads are rejected by default (`missing=:error`); pass
`missing=:face` to take their stored face values as edge capacities (see
`Networks.require_observed`).
"""
function flowbet(net; missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="flowbet")
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

"""
    centralization(net, measure; mode=:total, normalized=true, missing=:error) -> Float64

Compute Freeman graph centralization for a vertex centrality `measure`,
following R `sna::centralization`:

    C = Σᵢ (c_max − cᵢ) / C_max

where `C_max` is the theoretical maximum deviation sum for a network of the
same size (attained by the star for the classic measures). With
`normalized=false` the raw deviation sum `Σᵢ (c_max − cᵢ)` is returned.

# Arguments
- `net`: Network object
- `measure::Symbol`: Centrality measure to centralize
    - `:degree`: [`degree_centrality`](@ref); `mode` selects `:in`, `:out`,
      or `:total` (Freeman) degree for directed networks
    - `:betweenness`: [`betweenness_centrality`](@ref) (raw scores)
    - `:closeness`: Freeman closeness `(n-1)/Σⱼ d(i,j)`, with the score
      taken as 0 when some vertex is unreachable (sna's convention — a
      disconnected network has closeness centralization 0)
    - `:eigenvector`: [`eigenvector_centrality`](@ref)
- `mode::Symbol=:total`: Degree type for `measure == :degree` (ignored
  otherwise, and for undirected networks)
- `normalized::Bool=true`: Divide by the theoretical maximum (sna
  `normalize=TRUE`)
- `missing::Symbol=:error`: Missing-dyad policy (`Networks.require_observed`).
  Centralization is an inferential summary of the whole structure, so the
  default `:error` refuses a network with masked (unobserved) dyads rather
  than centralizing a partly invented one; `:face` opts in to the stored
  face values and is forwarded to the underlying centrality measure

The theoretical maxima match `sna`'s `tmaxdev` values: e.g. `(n-1)(n-2)`
for undirected degree, `(n-1)²(n-2)` / `(n-1)²(n-2)/2` for directed /
undirected betweenness, `(n-1)(1-1/n)` / `(n-1)(n-2)/(2n-3)` for directed /
undirected closeness, and `n-1` / `√2(n-2)/2` for directed / undirected
eigenvector centrality.
"""
function centralization(net, measure::Symbol; mode::Symbol=:total,
                        normalized::Bool=true, missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="centralization")
    n = nv(net)
    directed = is_directed(net)

    cv, tmax = if measure == :degree
        tmax = if directed
            mode == :total ? Float64((n - 1) * (2 * (n - 1) - 2)) :
                             Float64((n - 1) * (n - 1))
        else
            Float64((n - 1) * (n - 2))
        end
        (degree_centrality(net; mode=mode, missing=policy), tmax)
    elseif measure == :betweenness
        tmax = directed ? Float64((n - 1)^2 * (n - 2)) :
                          (n - 1)^2 * (n - 2) / 2
        (betweenness_centrality(net; missing=policy), tmax)
    elseif measure == :closeness
        tmax = directed ? (n - 1) * (1 - 1 / n) :
                          (n - 2) * (n - 1) / (2 * n - 3)
        (_freeman_closeness(net; missing=policy), tmax)
    elseif measure == :eigenvector
        tmax = directed ? Float64(n - 1) : sqrt(2) / 2 * (n - 2)
        (eigenvector_centrality(net; missing=policy), tmax)
    else
        throw(ArgumentError("Unknown centralization measure: $measure " *
                            "(use :degree, :betweenness, :closeness, or " *
                            ":eigenvector)"))
    end

    isempty(cv) && return 0.0
    cent = sum(maximum(cv) .- cv)
    if normalized
        tmax > 0 || return 0.0
        cent /= tmax
    end
    return cent
end

# Freeman closeness with sna's unreachability convention: (n-1) over the
# total geodesic distance to all other vertices, 0 when any is unreachable
# (this is what sna::closeness computes and what centralization expects;
# closeness_centrality uses Graphs.jl's component-based scaling instead).
function _freeman_closeness(net; missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    n = nv(net)
    dist = geodesic_distance(net; missing=policy)
    clo = zeros(n)
    for i in 1:n
        total = 0.0
        for j in 1:n
            i == j && continue
            total += dist[i, j]
        end
        clo[i] = isfinite(total) && total > 0 ? (n - 1) / total : 0.0
    end
    return clo
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
