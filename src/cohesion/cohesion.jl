"""
Cohesion measures for network analysis.

Provides functions for analyzing network cohesion including components,
cliques, k-cores, cutpoints, and bridges.
"""

using Graphs

"""
    components(net; mode=:weak, missing=:error) -> Vector{Vector{Int}}

Find connected components of the network.

# Arguments
- `net`: Network object
- `mode::Symbol=:weak`: Component type
    - `:weak`: Weakly connected components (ignoring direction)
    - `:strong`: Strongly connected components
- `missing::Symbol=:error`: Missing-dyad policy (`Networks.require_observed`);
  `:error` rejects a network with masked (unobserved) dyads, `:face` treats
  each masked dyad as its stored face value
"""
function components(net; mode::Symbol=:weak, missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="components")
    if mode == :strong && is_directed(net)
        return Graphs.strongly_connected_components(net.graph)
    else
        return Graphs.connected_components(net.graph)
    end
end

"""
    largest_component(net; mode=:weak, missing=:error) -> Vector{Int}

Return vertices in the largest connected component.

Masked (unobserved) dyads are rejected by default (`missing=:error`); pass
`missing=:face` to treat them as their stored face values (see
`Networks.require_observed`).
"""
function largest_component(net; mode::Symbol=:weak, missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="largest_component")
    comps = components(net; mode=mode, missing=policy)
    return comps[argmax(length.(comps))]
end

"""
    cliques(net; min_size=3, missing=:error) -> Vector{Vector{Int}}

Find all maximal cliques (complete subgraphs) of at least the specified size.

Cliques are an undirected concept; as in R `sna::clique.census`, directed
networks are symmetrized first (weak rule: an undirected tie exists if an
arc exists in either direction).

Note: Finding all cliques is NP-complete. This returns maximal cliques
for large networks.

Masked (unobserved) dyads are rejected by default (`missing=:error`); pass
`missing=:face` to let them enter (or break) cliques at their stored face
values (see `Networks.require_observed`).
"""
function cliques(net; min_size::Int=3, missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="cliques")
    # The backing store is always a digraph (symmetric for undirected
    # networks); maximal_cliques needs a SimpleGraph. For directed networks
    # this symmetrizes with the weak (either-direction) rule, as sna does.
    g = Graphs.SimpleGraph(net.graph)
    all_cliques = Graphs.maximal_cliques(g)
    return filter(c -> length(c) >= min_size, all_cliques)
end

"""
    kcores(net; k=1, missing=:error) -> Vector{Int}

Find the k-core of the network.

The k-core is the maximal subgraph where every vertex has degree at least k.

Masked (unobserved) dyads are rejected by default (`missing=:error`); pass
`missing=:face` to take their stored face values as ties (see
`Networks.require_observed`).

# Returns
Vector of vertices in the k-core.
"""
function kcores(net; k::Int=1, missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="kcores")
    g = net.graph
    if !is_directed(net)
        g = Graphs.SimpleGraph(g)
    end
    core = Graphs.core_number(g)

    return findall(c -> c >= k, core)
end

"""
    cutpoints(net; missing=:error) -> Vector{Int}

Find all cutpoints (articulation points) in the network.

A cutpoint is a vertex whose removal disconnects the network.

Masked (unobserved) dyads are rejected by default (`missing=:error`); pass
`missing=:face` to take their stored face values as ties (see
`Networks.require_observed`).
"""
function cutpoints(net; missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="cutpoints")
    g = Graphs.SimpleGraph(net.graph)
    return Graphs.articulation(g)
end

"""
    bridges(net; missing=:error) -> Vector{Tuple{Int,Int}}

Find all bridges in the network.

A bridge is an edge whose removal disconnects the network.

Extends `Graphs.bridges` for `AbstractNetwork` types.

Masked (unobserved) dyads are rejected by default (`missing=:error`); pass
`missing=:face` to take their stored face values as ties (see
`Networks.require_observed`).
"""
function bridges(net::AbstractNetwork; missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="bridges")
    g = Graphs.SimpleGraph(net.graph)
    bridge_edges = Graphs.bridges(g)
    return [(src(e), dst(e)) for e in bridge_edges]
end

"""
    bicomponents(net; missing=:error) -> Vector{Vector{Tuple{Int,Int}}}

Find biconnected components of the network (as lists of edges).

Biconnectivity is an undirected concept; directed networks are treated as
their underlying undirected graph.

Masked (unobserved) dyads are rejected by default (`missing=:error`); pass
`missing=:face` to take their stored face values as ties (see
`Networks.require_observed`).
"""
function bicomponents(net; missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="bicomponents")
    g = Graphs.SimpleGraph(net.graph)
    comps = Graphs.biconnected_components(g)
    return [[(src(e), dst(e)) for e in comp] for comp in comps]
end

"""
    geodesic_distance(net; missing=:error) -> Matrix{Float64}

Compute the matrix of geodesic (shortest path) distances.

Returns Inf for unreachable pairs.

Masked (unobserved) dyads are rejected by default (`missing=:error`); pass
`missing=:face` to traverse them at their stored face values (see
`Networks.require_observed`).
"""
function geodesic_distance(net; missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="geodesic_distance")
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
    diameter(net; missing=:error) -> Float64

Compute the diameter of the network (longest shortest path).

Returns Inf if the network is disconnected.

Extends `Graphs.diameter` for `AbstractNetwork` types.

Masked (unobserved) dyads are rejected by default (`missing=:error`); pass
`missing=:face` to traverse them at their stored face values (see
`Networks.require_observed`).
"""
function diameter(net::AbstractNetwork; missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="diameter")
    dist = geodesic_distance(net; missing=policy)
    finite_dist = filter(isfinite, dist)
    return isempty(finite_dist) ? Inf : maximum(finite_dist)
end

"""
    average_path_length(net; missing=:error) -> Float64

Compute the average shortest path length over all reachable pairs.

Masked (unobserved) dyads are rejected by default (`missing=:error`); pass
`missing=:face` to traverse them at their stored face values (see
`Networks.require_observed`).
"""
function average_path_length(net; missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="average_path_length")
    dist = geodesic_distance(net; missing=policy)
    finite_dist = filter(d -> isfinite(d) && d > 0, dist)
    return isempty(finite_dist) ? Inf : mean(finite_dist)
end
