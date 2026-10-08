"""
Cohesion measures for network analysis: components, cliques, k-cores,
cutpoints, bicomponents and geodesics, with R `sna`'s names where sna has
the measure (`component_dist`, `kcores`, `cutpoints`, `geodist`).
"""

using Graphs

"""
    component_dist(net; connected=:strong, missing=:error) -> NamedTuple

The component structure of the network, R `sna::component.dist`:
`(membership, csize, cdist)`, where

- `membership[v]` is the component of vertex `v`, components being numbered
  in the order of their lowest-numbered vertex;
- `csize[c]` is the number of vertices in component `c`;
- `cdist[s]` is the number of components with `s` vertices, `s = 1, …, n`.

`connected` selects the component type:

- `:strong` (default): strong components (mutual reachability);
- `:weak`: weak components (direction ignored);
- `:recursive`: components of the mutual-tie graph.

On an undirected network the three coincide. sna's `unilateral` components
are not unique and are not implemented.

Before 0.2.0 `component_dist` returned the component sizes sorted in
decreasing order; that is `sort(component_dist(net).csize; rev=true)`.

`missing=:error` refuses masked dyads; `missing=:face` takes their stored
face values as ties.

# Example
```julia
using SNA
net = network(5)
add_edges!(net, [(1, 2), (2, 1), (2, 3), (4, 5)])
cd = component_dist(net; connected=:weak)
cd.membership                        # [1, 1, 1, 2, 2]
cd.csize                             # [3, 2]
component_dist(net).csize            # strong: [2, 1, 1, 1]
```
"""
function component_dist(net::AbstractNetwork; connected::Symbol=:strong,
                        missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="component_dist")
    connected in (:strong, :weak, :recursive) ||
        throw(ArgumentError("component_dist: connected must be :strong, :weak or :recursive " *
                            "(sna's unilateral components are not implemented); got :$connected"))
    n = nv(net)
    groups = if !is_directed(net) || connected == :weak
        Graphs.weakly_connected_components(_graph(net))
    elseif connected == :strong
        Graphs.strongly_connected_components(_graph(net))
    else
        Graphs.connected_components(_symmetrized_graph(net, :strong))
    end
    # Number components by their lowest vertex, as sna does.
    sort!(groups; by=minimum)
    membership = zeros(Int, n)
    for (c, g) in enumerate(groups), v in g
        membership[v] = c
    end
    csize = length.(groups)
    cdist = zeros(Int, n)
    for s in csize
        cdist[s] += 1
    end
    return (membership=membership, csize=csize, cdist=cdist)
end

"""
    largest_component(net; connected=:strong, missing=:error) -> Vector{Int}

The vertices of the largest component (the first one, if several share the
largest size), sorted. `connected` is as in [`component_dist`](@ref).
R's `component.largest` returns the same set as a logical vector.

`missing=:error` refuses masked dyads; `missing=:face` takes their stored
face values as ties.

# Example
```julia
using SNA
net = network(5; directed=false)
add_edges!(net, [(1, 2), (2, 3), (4, 5)])
largest_component(net)               # [1, 2, 3]
```
"""
function largest_component(net::AbstractNetwork; connected::Symbol=:strong,
                           missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="largest_component")
    nv(net) == 0 && return Int[]
    cd = component_dist(net; connected, missing=policy)
    return findall(==(argmax(cd.csize)), cd.membership)
end

"""
    cliques(net; min_size=3, symmetrize=:strong, missing=:error) -> Vector{Vector{Int}}

All maximal cliques (complete subgraphs) with at least `min_size` vertices.

Cliques are an undirected concept; as in R `sna::clique.census`, directed
networks use mutual arcs (`symmetrize=:strong`), and `symmetrize=:weak`
accepts an arc in either direction. (sna's `clique.census` returns clique
counts by size; this function lists the cliques.) Enumerating cliques takes
exponential time in the worst case.

`missing=:error` refuses masked dyads; `missing=:face` lets them enter (or
break) cliques at their stored face values.

# Example
```julia
using SNA
net = network(5; directed=false)
add_edges!(net, [(1, 2), (1, 3), (2, 3), (3, 4), (4, 5)])
cliques(net)                         # [[1, 2, 3]] (in some order)
```
"""
function cliques(net::AbstractNetwork; min_size::Int=3, symmetrize::Symbol=:strong,
                 missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="cliques")
    g = _symmetrized_graph(net, symmetrize)
    all_cliques = Graphs.maximal_cliques(g)
    return filter(c -> length(c) >= min_size, all_cliques)
end

"""
    kcores(net; cmode=:freeman, missing=:error) -> Vector{Int}

The core number of every vertex, R `sna::kcores`: vertex `v` has core
number `k` when it belongs to the `k`-core (the largest subgraph in which
every vertex has degree at least `k`) but not to the `(k+1)`-core.

`cmode` is the degree used on a directed network: `:freeman` (in + out,
sna's default), `:indegree` or `:outdegree`. Undirected networks use the
ordinary degree. Self-loops are ignored, as in sna.

Before 0.2.0 `kcores(net; k)` returned the vertices of the `k`-core; they
are `findall(>=(k), kcores(net))`. Passing `k` raises an `ArgumentError`
that says so.

`missing=:error` refuses masked dyads; `missing=:face` takes their stored
face values as ties.

# Example
```julia
using SNA
net = network(5; directed=false)
add_edges!(net, [(1, 2), (1, 3), (2, 3), (3, 4), (4, 5)])
kcores(net)                          # [2, 2, 2, 1, 1]
findall(>=(2), kcores(net))          # the 2-core: [1, 2, 3]
```
"""
function kcores(net::AbstractNetwork; cmode::Symbol=:freeman,
                k::Union{Nothing,Integer}=nothing, missing::Symbol=:error)
    k === nothing ||
        throw(ArgumentError("kcores now returns the core number of every vertex, as R " *
                            "sna::kcores does; the vertices of the $k-core are " *
                            "findall(>=($k), kcores(net))"))
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="kcores")
    cmode in (:freeman, :indegree, :outdegree) ||
        throw(ArgumentError("kcores: cmode must be :freeman, :indegree or :outdegree; got :$cmode"))
    n = nv(net)
    directed = is_directed(net)
    # w[u, v]: how much u's degree drops when v is peeled away
    outs = [[w for w in outneighbors(net, v) if w != v] for v in 1:n]
    ins = directed ? [[w for w in inneighbors(net, v) if w != v] for v in 1:n] : outs
    deg = if !directed
        length.(outs)
    elseif cmode == :freeman
        length.(outs) .+ length.(ins)
    elseif cmode == :indegree
        length.(ins)
    else
        length.(outs)
    end
    core = zeros(Int, n)
    removed = falses(n)
    level = 0
    for _ in 1:n
        # peel the remaining vertex of smallest current degree
        v, dv = 0, typemax(Int)
        for u in 1:n
            !removed[u] && deg[u] < dv && ((v, dv) = (u, deg[u]))
        end
        level = max(level, dv)
        core[v] = level
        removed[v] = true
        # v's arcs leave the degrees of its remaining neighbours
        if !directed
            for u in outs[v]
                removed[u] || (deg[u] -= 1)
            end
        else
            if cmode != :outdegree          # u's in-degree loses the arc v→u
                for u in outs[v]
                    removed[u] || (deg[u] -= 1)
                end
            end
            if cmode != :indegree           # u's out-degree loses the arc u→v
                for u in ins[v]
                    removed[u] || (deg[u] -= 1)
                end
            end
        end
    end
    return core
end

"""
    cutpoints(net; connected=:strong, symmetrize=nothing, missing=:error) -> Vector{Int}

The cutpoints (articulation points) of the network, R `sna::cutpoints`: the
vertices whose removal increases the number of components.

Directed networks default to strong connectivity (`connected=:strong`, sna's
default): a cutpoint splits a strong component. `connected=:weak` ignores
direction; `connected=:recursive` uses mutual arcs. An explicit
`symmetrize=:strong` or `:weak` also requests articulation of that
undirected graph.

`missing=:error` refuses masked dyads; `missing=:face` takes their stored
face values as ties.

# Example
```julia
using SNA
path = network(4; directed=false)
add_edges!(path, [(1, 2), (2, 3), (3, 4)])
sort(cutpoints(path))                # [2, 3]
```
"""
function cutpoints(net::AbstractNetwork; connected::Symbol=:strong,
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
    bicomponents(net; symmetrize=:strong, min_size=3, missing=:error) -> Vector{Vector{Tuple{Int,Int}}}

The biconnected components of the network, as lists of edges. The default
`min_size=3` excludes bridges, as R `sna::bicomponent.dist` does;
`min_size=2` includes the bridge components.

Biconnectivity is an undirected concept; directed networks use mutual arcs
by default (`symmetrize=:strong`), as in sna. `symmetrize=:weak` accepts an
arc in either direction.

`missing=:error` refuses masked dyads; `missing=:face` takes their stored
face values as ties.

# Example
```julia
using SNA
net = network(5; directed=false)
add_edges!(net, [(1, 2), (2, 3), (1, 3), (3, 4), (4, 5)])
length(bicomponents(net))            # 1: the triangle
```
"""
function bicomponents(net::AbstractNetwork; symmetrize::Symbol=:strong, min_size::Int=3,
                      missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="bicomponents")
    g = _symmetrized_graph(net, symmetrize)
    comps = Graphs.biconnected_components(g)
    min_size >= 2 || throw(ArgumentError("min_size must be at least 2"))
    return [[(src(e), dst(e)) for e in comp] for comp in comps
            if length(Set(v for e in comp for v in (src(e), dst(e)))) >= min_size]
end

"""
    geodist(net; inf_replace=Inf, missing=:error) -> NamedTuple

Geodesic distances and geodesic counts, R `sna::geodist`: `(counts, gdist)`,
where `gdist[i, j]` is the length of a shortest directed path from `i` to `j`
(`inf_replace`, by default `Inf`, when `j` is unreachable) and
`counts[i, j]` is the number of such paths (0 when unreachable). The
diagonal has distance 0 and count 1. Self-loops play no part.

The diameter of a network, which sna does not provide, is
`maximum(geodist(net).gdist)` (`Inf` when some vertex cannot reach another).
`Graphs.diameter` is Graphs.jl's function and keeps Graphs' semantics.

`missing=:error` refuses masked dyads; `missing=:face` traverses them at
their stored face values.

# Example
```julia
using SNA
net = network(4)
add_edges!(net, [(1, 2), (1, 3), (2, 4), (3, 4)])
g = geodist(net)
g.gdist[1, 4], g.counts[1, 4]        # (2.0, 2.0): two geodesics of length 2
g.gdist[4, 1]                        # Inf
```
"""
function geodist(net::AbstractNetwork; inf_replace::Real=Inf, missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="geodist")
    out = _geodist(_graph(net), nv(net))
    isinf(inf_replace) && inf_replace > 0 && return out
    replace!(out.gdist, Inf => Float64(inf_replace))
    return out
end

# Breadth-first geodesic distances and path counts from every source. `n`
# may be any integer type (`nv` of a `Network{Int32}` is an `Int32`).
function _geodist(g, n::Integer)
    n = Int(n)
    gdist = fill(Inf, n, n)
    counts = zeros(n, n)
    dist = fill(-1, n)
    sigma = zeros(n)
    queue = Vector{Int}(undef, n)
    for s in 1:n
        fill!(dist, -1)
        fill!(sigma, 0.0)
        dist[s] = 0
        sigma[s] = 1.0
        queue[1] = s
        head, tail = 1, 1
        while head <= tail
            v = queue[head]
            head += 1
            for w in Graphs.outneighbors(g, v)
                if dist[w] < 0
                    dist[w] = dist[v] + 1
                    tail += 1
                    queue[tail] = w
                end
                dist[w] == dist[v] + 1 && (sigma[w] += sigma[v])
            end
        end
        for t in 1:n
            if dist[t] >= 0
                gdist[s, t] = dist[t]
                counts[s, t] = sigma[t]
            end
        end
    end
    return (counts=counts, gdist=gdist)
end

"""
    average_path_length(net; missing=:error) -> Float64

The mean geodesic distance over ordered pairs of distinct vertices that are
connected by a path (unreachable pairs are left out). `Inf` when no vertex
reaches another. This is not an sna function.

`missing=:error` refuses masked dyads; `missing=:face` traverses them at
their stored face values.

# Example
```julia
using SNA
path = network(3; directed=false)
add_edges!(path, [(1, 2), (2, 3)])
average_path_length(path)            # (1 + 1 + 2) × 2 / 6 = 1.3333
```
"""
function average_path_length(net::AbstractNetwork; missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="average_path_length")
    dist = geodist(net; missing=policy).gdist
    finite_dist = filter(d -> isfinite(d) && d > 0, dist)
    return isempty(finite_dist) ? Inf : mean(finite_dist)
end
