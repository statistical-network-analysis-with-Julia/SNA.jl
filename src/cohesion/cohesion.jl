"""
Cohesion measures for network analysis.

Provides functions for analyzing network cohesion including components,
cliques, k-cores, cutpoints, and bridges.
"""

using Graphs

"""
    components(net; mode=:weak) -> Vector{Vector{Int}}

Find connected components of the network.

# Arguments
- `net`: Network object
- `mode::Symbol=:weak`: Component type
    - `:weak`: Weakly connected components (ignoring direction)
    - `:strong`: Strongly connected components
"""
function components(net; mode::Symbol=:weak)
    if mode == :strong && is_directed(net)
        return Graphs.strongly_connected_components(net.graph)
    else
        return Graphs.connected_components(net.graph)
    end
end

"""
    largest_component(net; mode=:weak) -> Vector{Int}

Return vertices in the largest connected component.
"""
function largest_component(net; mode::Symbol=:weak)
    comps = components(net; mode=mode)
    return comps[argmax(length.(comps))]
end

"""
    cliques(net; min_size=3) -> Vector{Vector{Int}}

Find all cliques (complete subgraphs) of at least the specified size.

Note: Finding all cliques is NP-complete. This returns maximal cliques
for large networks.
"""
function cliques(net; min_size::Int=3)
    # Use maximal cliques from Graphs.jl
    all_cliques = Graphs.maximal_cliques(net.graph)
    return filter(c -> length(c) >= min_size, all_cliques)
end

"""
    kcores(net; k=1) -> Vector{Int}

Find the k-core of the network.

The k-core is the maximal subgraph where every vertex has degree at least k.

# Returns
Vector of vertices in the k-core.
"""
function kcores(net; k::Int=1)
    n = nv(net)
    core = Graphs.core_number(net.graph)

    return findall(c -> c >= k, core)
end

"""
    cutpoints(net) -> Vector{Int}

Find all cutpoints (articulation points) in the network.

A cutpoint is a vertex whose removal disconnects the network.
"""
function cutpoints(net)
    return Graphs.articulation(net.graph)
end

"""
    bridges(net) -> Vector{Tuple{Int,Int}}

Find all bridges in the network.

A bridge is an edge whose removal disconnects the network.
"""
function bridges(net)
    bridge_edges = Graphs.bridges(net.graph)
    return [(src(e), dst(e)) for e in bridge_edges]
end

"""
    bicomponents(net) -> Vector{Vector{Int}}

Find biconnected components of the network.

A biconnected component is a maximal subgraph that has no cutpoints.
"""
function bicomponents(net)
    return Graphs.biconnected_components(net.graph)
end

"""
    geodesic_distance(net) -> Matrix{Float64}

Compute the matrix of geodesic (shortest path) distances.

Returns Inf for unreachable pairs.
"""
function geodesic_distance(net)
    n = nv(net)
    dist = fill(Inf, n, n)

    for i in 1:n
        # BFS from vertex i
        distances = Graphs.gdistances(net.graph, i)
        for j in 1:n
            if distances[j] < typemax(Int)
                dist[i, j] = Float64(distances[j])
            end
        end
    end

    return dist
end

"""
    diameter(net) -> Float64

Compute the diameter of the network (longest shortest path).

Returns Inf if the network is disconnected.
"""
function diameter(net)
    dist = geodesic_distance(net)
    finite_dist = filter(isfinite, dist)
    return isempty(finite_dist) ? Inf : maximum(finite_dist)
end

"""
    average_path_length(net) -> Float64

Compute the average shortest path length over all reachable pairs.
"""
function average_path_length(net)
    dist = geodesic_distance(net)
    finite_dist = filter(d -> isfinite(d) && d > 0, dist)
    return isempty(finite_dist) ? Inf : mean(finite_dist)
end
