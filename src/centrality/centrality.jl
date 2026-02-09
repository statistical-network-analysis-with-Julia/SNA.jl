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
        max_degree = is_directed(net) ? (n - 1) : 2 * (n - 1)
        centrality ./= max_degree
    end

    return centrality
end

"""
    betweenness_centrality(net; normalized=true) -> Vector{Float64}

Compute betweenness centrality for all vertices.

Betweenness centrality measures the extent to which a vertex lies on paths
between other vertices.
"""
function betweenness_centrality(net; normalized::Bool=true)
    # Use Graphs.jl implementation
    bc = Graphs.betweenness_centrality(net.graph; normalize=normalized)
    return bc
end

"""
    closeness_centrality(net; normalized=true) -> Vector{Float64}

Compute closeness centrality for all vertices.

Closeness centrality measures how close a vertex is to all other vertices
(inverse of average shortest path length).
"""
function closeness_centrality(net; normalized::Bool=true)
    cc = Graphs.closeness_centrality(net.graph; normalize=normalized)
    return cc
end

"""
    eigenvector_centrality(net; max_iter=100, tol=1e-6) -> Vector{Float64}

Compute eigenvector centrality for all vertices.

A vertex has high eigenvector centrality if it is connected to other
high-centrality vertices.
"""
function eigenvector_centrality(net; max_iter::Int=100, tol::Float64=1e-6)
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
            return x_new
        end
        x = x_new
    end

    return x
end

"""
    bonacich_power(net; β=0.5, normalized=true) -> Vector{Float64}

Compute Bonacich power centrality.

# Arguments
- `net`: Network object
- `β::Float64=0.5`: Attenuation factor (-1/λ_max < β < 1/λ_max)
- `normalized::Bool=true`: Normalize results

Bonacich power centrality accounts for both the number and quality of
connections, with β controlling how much weight is given to indirect ties.
"""
function bonacich_power(net; β::Float64=0.5, normalized::Bool=true)
    n = nv(net)
    A = as_matrix(net)
    I_mat = Matrix{Float64}(I, n, n)

    # c = (I - βA)^(-1) * A * 1
    try
        M = I_mat - β * A
        c = M \ (A * ones(n))

        if normalized && maximum(abs.(c)) > 0
            c ./= maximum(abs.(c))
        end

        return c
    catch e
        @warn "Bonacich power computation failed: $e"
        return zeros(n)
    end
end

"""
    katz_centrality(net; α=0.1, β=1.0) -> Vector{Float64}

Compute Katz centrality.

Similar to eigenvector centrality but with damping factor α.
"""
function katz_centrality(net; α::Float64=0.1, β::Float64=1.0)
    return Graphs.katz_centrality(net.graph; α=α)
end

"""
    pagerank(net; α=0.85, max_iter=100, tol=1e-6) -> Vector{Float64}

Compute PageRank centrality.
"""
function pagerank(net; α::Float64=0.85, max_iter::Int=100, tol::Float64=1e-6)
    return Graphs.pagerank(net.graph; α=α, n=max_iter, ϵ=tol)
end

"""
    flowbet(net) -> Vector{Float64}

Compute flow betweenness centrality.

Flow betweenness considers all paths (not just shortest) weighted by
their length.
"""
function flowbet(net)
    # Simplified implementation - full flow betweenness requires
    # computing maximum flow between all pairs
    # For now, return standard betweenness as approximation
    @warn "flowbet: returning standard betweenness as approximation"
    return betweenness_centrality(net)
end
