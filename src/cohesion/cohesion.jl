"""
Cohesion measures for network analysis.

Provides functions for analyzing network cohesion including components,
cliques, k-cores, cutpoints, and bridges.
"""

using Graphs

"""
    components(net; connected=:strong, missing=:error) -> Vector{Vector{Int}}

Find connected components of the network.

# Arguments
- `net`: Network object
- `connected::Symbol=:strong`: Component type (the legacy `mode` keyword is accepted)
    - `:weak`: Weakly connected components (ignoring direction)
    - `:strong`: Strongly connected components
- `missing::Symbol=:error`: Missing-dyad policy (`Networks.require_observed`);
  `:error` rejects a network with masked (unobserved) dyads, `:face` treats
  each masked dyad as its stored face value
"""
function components(net; connected::Symbol=:strong, mode::Union{Nothing,Symbol}=nothing,
                    missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="components")
    connected = something(mode, connected)
    connected in (:strong, :weak) || throw(ArgumentError("connected must be :strong or :weak"))
    if connected == :strong && is_directed(net)
        return Graphs.strongly_connected_components(_graph(net))
    else
        return Graphs.weakly_connected_components(_graph(net))
    end
end

"""
    largest_component(net; connected=:strong, missing=:error) -> Vector{Int}

Return vertices in the largest connected component.

Masked (unobserved) dyads are rejected by default (`missing=:error`); pass
`missing=:face` to treat them as their stored face values (see
`Networks.require_observed`).
"""
function largest_component(net; connected::Symbol=:strong, mode::Union{Nothing,Symbol}=nothing,
                           missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="largest_component")
    comps = components(net; connected, mode, missing=policy)
    isempty(comps) && return Int[]
    return comps[argmax(length.(comps))]
end

"""
    cliques(net; min_size=3, symmetrize=:strong, missing=:error) -> Vector{Vector{Int}}

Find all maximal cliques (complete subgraphs) of at least the specified size.

Cliques are an undirected concept; as in R `sna::clique.census`, directed
networks use mutual arcs (`symmetrize=:strong`). `symmetrize=:weak` accepts
an arc in either direction.

Note: Finding all cliques is NP-complete. This returns maximal cliques
for large networks.

Masked (unobserved) dyads are rejected by default (`missing=:error`); pass
`missing=:face` to let them enter (or break) cliques at their stored face
values (see `Networks.require_observed`).
"""
function cliques(net; min_size::Int=3, symmetrize::Symbol=:strong, missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="cliques")
    g = _symmetrized_graph(net, symmetrize)
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
    g = _graph(net)
    if !is_directed(net)
        g = Graphs.SimpleGraph(g)
    end
    core = Graphs.core_number(g)

    return findall(c -> c >= k, core)
end

"""
    cutpoints(net; connected=:strong, symmetrize=nothing, missing=:error) -> Vector{Int}

Find all cutpoints (articulation points) in the network.

A cutpoint is a vertex whose removal increases the component count. Directed
networks default to strong connectivity, matching `sna::cutpoints`.
`connected=:weak` ignores direction; `connected=:recursive` uses mutual arcs.
An explicit `symmetrize=:strong` or `:weak` also requests articulation of
that undirected graph.

Masked (unobserved) dyads are rejected by default (`missing=:error`); pass
`missing=:face` to take their stored face values as ties (see
`Networks.require_observed`).
"""
function cutpoints(net; connected::Symbol=:strong,
                   symmetrize::Union{Nothing,Symbol}=nothing, missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="cutpoints")
    if symmetrize !== nothing
        return Graphs.articulation(_symmetrized_graph(net, symmetrize))
    end
    connected in (:strong, :weak, :recursive) ||
        throw(ArgumentError("connected must be :strong, :weak, or :recursive"))
    if !is_directed(net) || connected != :strong
        return Graphs.articulation(_symmetrized_graph(net, connected == :recursive ? :strong : :weak))
    end
    # sna's strong cutpoints split a strongly connected component on removal.
    # Mutual-arc symmetrization is the distinct `connected=:recursive` option.
    n = nv(net)
    baseline = length(Graphs.strongly_connected_components(_graph(net)))
    result = Int[]
    for v in 1:n
        g, _ = Graphs.induced_subgraph(_graph(net), [u for u in 1:n if u != v])
        length(Graphs.strongly_connected_components(g)) > baseline && push!(result, v)
    end
    return result
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
    g = Graphs.SimpleGraph(_graph(net))
    bridge_edges = Graphs.bridges(g)
    return [(src(e), dst(e)) for e in bridge_edges]
end

"""
    bicomponents(net; symmetrize=:strong, min_size=3, missing=:error) -> Vector{Vector{Tuple{Int,Int}}}

Find biconnected components of the network (as lists of edges).
The default `min_size=3` excludes bridges, as in R `bicomponent.dist`;
`min_size=2` includes bridge components.

Biconnectivity is an undirected concept; directed networks are treated as
mutual arcs by default (`symmetrize=:strong`), matching R. Use
`symmetrize=:weak` for an arc in either direction.

Masked (unobserved) dyads are rejected by default (`missing=:error`); pass
`missing=:face` to take their stored face values as ties (see
`Networks.require_observed`).
"""
function bicomponents(net; symmetrize::Symbol=:strong, min_size::Int=3, missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="bicomponents")
    g = _symmetrized_graph(net, symmetrize)
    comps = Graphs.biconnected_components(g)
    min_size >= 2 || throw(ArgumentError("min_size must be at least 2"))
    return [[(src(e), dst(e)) for e in comp] for comp in comps
            if length(Set(v for e in comp for v in (src(e), dst(e)))) >= min_size]
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
        distances = Graphs.gdistances(_graph(net), i)
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
    return isempty(dist) ? Inf : maximum(dist)
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
