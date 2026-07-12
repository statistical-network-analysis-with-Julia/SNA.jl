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

Find all maximal cliques (complete subgraphs) of at least the specified size.

Cliques are an undirected concept; as in R `sna::clique.census`, directed
networks are symmetrized first (weak rule: an undirected tie exists if an
arc exists in either direction).

Note: Finding all cliques is NP-complete. This returns maximal cliques
for large networks.
"""
function cliques(net; min_size::Int=3)
    # The backing store is always a digraph (symmetric for undirected
    # networks); maximal_cliques needs a SimpleGraph. For directed networks
    # this symmetrizes with the weak (either-direction) rule, as sna does.
    g = Graphs.SimpleGraph(net.graph)
    all_cliques = Graphs.maximal_cliques(g)
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
    g = net.graph
    if !is_directed(net)
        g = Graphs.SimpleGraph(g)
    end
    core = Graphs.core_number(g)

    return findall(c -> c >= k, core)
end

"""
    cutpoints(net) -> Vector{Int}

Find all cutpoints (articulation points) in the network.

A cutpoint is a vertex whose removal disconnects the network.
"""
function cutpoints(net)
    g = Graphs.SimpleGraph(net.graph)
    return Graphs.articulation(g)
end

"""
    bridges(net) -> Vector{Tuple{Int,Int}}

Find all bridges in the network.

A bridge is an edge whose removal disconnects the network.

Extends `Graphs.bridges` for `AbstractNetwork` types.
"""
function bridges(net::AbstractNetwork)
    g = Graphs.SimpleGraph(net.graph)
    bridge_edges = Graphs.bridges(g)
    return [(src(e), dst(e)) for e in bridge_edges]
end

"""
    bicomponents(net) -> Vector{Vector{Tuple{Int,Int}}}

Find biconnected components of the network (as lists of edges).

Biconnectivity is an undirected concept; directed networks are treated as
their underlying undirected graph.
"""
function bicomponents(net)
    g = Graphs.SimpleGraph(net.graph)
    comps = Graphs.biconnected_components(g)
    return [[(src(e), dst(e)) for e in comp] for comp in comps]
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

Extends `Graphs.diameter` for `AbstractNetwork` types.
"""
function diameter(net::AbstractNetwork)
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
