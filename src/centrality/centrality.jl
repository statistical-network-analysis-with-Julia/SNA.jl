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

Binary ties are used by default (`ignore_eval=true`), and loops are excluded
(`diag=false`). Set `ignore_eval=false, attr=:weight` to use numeric edge
attributes; every present edge must have that attribute.

"""
function degree_centrality(net::AbstractNetwork; mode::Symbol=:total,
                           normalized::Bool=false, diag::Bool=false,
                           ignore_eval::Bool=true, attr::Symbol=:weight, missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="degree_centrality")
    n = nv(net)
    mode in (:in, :out, :total) || throw(ArgumentError("mode must be :in, :out, or :total"))
    centrality = zeros(Float64, n)
    directed = is_directed(net)
    if !ignore_eval
        A = _sociomatrix(net; ignore_eval, attr, diag)
        centrality = !directed || mode == :out ? vec(sum(A; dims=2)) :
                     mode == :in ? vec(sum(A; dims=1)) :
                     vec(sum(A; dims=1)) + vec(sum(A; dims=2))
        # sna Freeman degree counts a diagonal entry once, not twice.
        directed && mode == :total && diag && (centrality .-= LinearAlgebra.diag(A))
    else
        for v in vertices(net)
            centrality[v] = !directed || mode == :out ? length(outneighbors(net, v)) :
                            mode == :in ? length(inneighbors(net, v)) :
                            length(inneighbors(net, v)) + length(outneighbors(net, v))
            if has_edge(net, v, v)
                if !diag
                    centrality[v] -= directed && mode == :total ? 2 : 1
                elseif directed && mode == :total
                    centrality[v] -= 1
                end
            end
        end
    end
    if normalized && n > 1
        max_degree = (directed && mode == :total ? 2 : 1) * (n - 1) + diag
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
    bc = Graphs.betweenness_centrality(_graph(net); normalize=normalized)
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
With `cmode=:sna` (default), a vertex scores zero when any other vertex is
unreachable, matching R. `cmode=:component` requests Graphs.jl's reachable
component scaling. `normalized=false` returns the reciprocal distance sum.

Masked (unobserved) dyads are rejected by default (`missing=:error`); pass
`missing=:face` to build the paths from their stored face values (see
`Networks.require_observed`).
"""
function closeness_centrality(net::AbstractNetwork; normalized::Bool=true, cmode::Symbol=:sna,
                              missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="closeness_centrality")
    cmode in (:sna, :component) || throw(ArgumentError("cmode must be :sna or :component"))
    if cmode == :component
        return Graphs.closeness_centrality(_graph(net); normalize=normalized)
    end
    cc = _freeman_closeness(net; missing=policy)
    !normalized && nv(net) > 1 && (cc ./= nv(net) - 1)
    return cc
end

"""
    eigenvector_centrality(net; max_iter=1000, tol=1e-10, missing=:error) -> Vector{Float64}

Compute eigenvector centrality for all vertices.

A vertex has high eigenvector centrality if it is connected to other
high-centrality vertices. As in R `sna::evcent`, the centrality is the
principal (right) eigenvector of the adjacency matrix — for directed
networks this weights an actor by the centrality of its out-neighbors.
Katz uses incoming ties; symmetrize first for the undirected notion. Scores are reported with non-negative
orientation and unit L2 norm.

Masked (unobserved) dyads are rejected by default (`missing=:error`); pass
`missing=:face` to read them from the adjacency matrix at their stored face
value (see `Networks.require_observed`).

Binary ties are used by default (`ignore_eval=true`), and loops are excluded
(`diag=false`). Set `ignore_eval=false, attr=:weight` to use numeric edge
attributes; every present edge must have that attribute.

"""
function eigenvector_centrality(net::AbstractNetwork; max_iter::Int=1000,
                                tol::Float64=1e-10, ignore_eval::Bool=true,
                                attr::Symbol=:weight, diag::Bool=false, missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="eigenvector_centrality")
    n = nv(net)
    n == 0 && return Float64[]
    A = _sociomatrix(net; ignore_eval, attr, diag)
    any(<(0), A) && throw(ArgumentError("eigenvector centrality requires nonnegative weights"))
    # Largest real eigenvalue, not largest modulus: bipartite graphs also have
    # -rho, which causes the old power iteration to oscillate without warning.
    F = eigen(A)
    i = argmax(real.(F.values))
    x = abs.(F.vectors[:, i])
    rho = real(F.values[i])
    count(v -> abs(v - rho) <= tol * max(1.0, abs(rho)), F.values) > 1 &&
        @warn "Eigenvector centrality is not unique: the leading eigenvalue is repeated"
    return x / norm(x)
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

Binary ties are used by default (`ignore_eval=true`), and loops are excluded
(`diag=false`). Set `ignore_eval=false, attr=:weight` to use numeric edge
attributes; every present edge must have that attribute.

"""
function bonacich_power(net; exponent::Float64=1.0, rescale::Bool=false,
                        tol::Float64=1e-7, ignore_eval::Bool=true,
                        attr::Symbol=:weight, diag::Bool=false, missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="bonacich_power")
    n = nv(net)
    A = _sociomatrix(net; ignore_eval, attr, diag)
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

Returns `(I - α Aᵀ)⁻¹ β𝟙`: incoming ties carry centrality, `β` is the
baseline score, and `α` must be below the inverse spectral radius. Scores
are unnormalized by default so changing `β` changes their magnitude.
`normalized=true` requests the previous unit-L2 convention.

Masked (unobserved) dyads are rejected by default (`missing=:error`); pass
`missing=:face` to use their stored face values (see
`Networks.require_observed`).

Binary ties are used by default (`ignore_eval=true`), and loops are excluded
(`diag=false`). Set `ignore_eval=false, attr=:weight` to use numeric edge
attributes; every present edge must have that attribute.

"""
function katz_centrality(net::AbstractNetwork; α::Float64=0.1, β::Float64=1.0,
                         normalized::Bool=false, ignore_eval::Bool=true,
                         attr::Symbol=:weight, diag::Bool=false, missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="katz_centrality")
    isfinite(α) && α >= 0 || throw(ArgumentError("α must be finite and nonnegative"))
    isfinite(β) && β >= 0 || throw(ArgumentError("β must be finite and nonnegative"))
    A = _sociomatrix(net; ignore_eval, attr, diag)
    isempty(A) && return Float64[]
    any(<(0), A) && throw(ArgumentError("Katz centrality requires nonnegative weights"))
    rho = maximum(abs, eigvals(A))
    α * rho < 1 || throw(ArgumentError("Katz centrality requires α < 1 / spectral_radius(A)"))
    c = (I - α * A') \ fill(β, nv(net))
    normalized && norm(c) > 0 && (c ./= norm(c))
    return c
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
    return Graphs.pagerank(_graph(net), α, max_iter, tol)
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

Binary ties are used by default (`ignore_eval=true`), and loops are excluded
(`diag=false`). Set `ignore_eval=false, attr=:weight` to use numeric edge
attributes; every present edge must have that attribute.

"""
function flowbet(net; ignore_eval::Bool=true, attr::Symbol=:weight,
                 diag::Bool=false, missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="flowbet")
    n = nv(net)
    A = _sociomatrix(net; ignore_eval, attr, diag)
    any(<(0), A) && throw(ArgumentError("flow capacities must be nonnegative"))
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
# closeness_centrality uses this same convention by default).
function _freeman_closeness(net; missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    n = nv(net)
    dist = geodesic_distance(net; missing=policy)
    clo = zeros(n)
    n == 1 && return [NaN]
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
