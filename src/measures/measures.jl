"""
Network-level measures and indices.

Provides global network statistics including density, reciprocity,
transitivity, and census functions.
"""

"""
    density(net) -> Float64

Compute network density (proportion of possible edges that exist).
Alias for network_density from Network.jl. Undirected edges are
single-counted, so an undirected network has `ne(net) / (n(n-1)/2)`.

Extends `Graphs.density` for `AbstractNetwork` types.
"""
density(net::AbstractNetwork) = network_density(net)

"""
    gden(net) -> Float64

Alias for [`density`](@ref) (R `sna` compatibility).
"""
gden(net) = density(net)

"""
    reciprocity(net; method=:dyadic) -> Float64

Compute network reciprocity, following R `sna::grecip`.

# Arguments
- `net`: Network object (should be directed)
- `method::Symbol=:dyadic`: Reciprocity measure type
    - `:dyadic`: Proportion of dyads that are symmetric — mutual or null
      count as reciprocated: `(M + N) / (M + A + N)` (sna default)
    - `:dyadic_nonnull`: Proportion of non-null dyads that are mutual:
      `M / (M + A)`
    - `:edgewise`: Proportion of edges that are reciprocated:
      `2M / (2M + A)`
"""
function reciprocity(net; method::Symbol=:dyadic)
    if !is_directed(net)
        return 1.0  # All ties are mutual in undirected networks
    end

    census = dyad_census(net)
    mutual, asymmetric, null = census.mutual, census.asymmetric, census.null

    if method == :dyadic
        # Symmetric dyads (mutual or null) over all dyads
        total = mutual + asymmetric + null
        return total > 0 ? (mutual + null) / total : 0.0
    elseif method == :dyadic_nonnull
        # Proportion of non-null dyads that are mutual
        total_ties = mutual + asymmetric
        return total_ties > 0 ? mutual / total_ties : 0.0
    elseif method == :edgewise
        # Proportion of edges that are reciprocated
        total_edges = 2 * mutual + asymmetric
        return total_edges > 0 ? (2 * mutual) / total_edges : 0.0
    else
        throw(ArgumentError("Unknown reciprocity method: $method"))
    end
end

"""
    grecip(net; kwargs...) -> Float64

Alias for [`reciprocity`](@ref) (R `sna` compatibility).
"""
grecip(net; kwargs...) = reciprocity(net; kwargs...)

"""
    transitivity(net; type=:global) -> Float64 or Vector{Float64}

Compute network transitivity.

# Arguments
- `net`: Network object
- `type::Symbol=:global`: Type of transitivity
    - `:global`: The weak transitivity of R `sna::gtrans(measure="weak")`:
      the fraction of potentially transitive ordered triples (`i→j` and
      `j→k`) that are transitive (`i→k` also present). For undirected
      networks this equals sna's graph-mode transitivity
      (3 × triangles / 2-paths). Returns `Float64`.
    - `:local`: Vector of local clustering coefficients (`Vector{Float64}`)
    - `:average`: Average of local clustering coefficients (`Float64`)
"""
function transitivity(net; type::Symbol=:global)
    if type == :global
        # sna::gtrans weak transitivity; for undirected networks the
        # symmetric adjacency makes this 3·triangles / 2-paths
        n = nv(net)
        potential = 0
        realized = 0
        for j in 1:n
            for i in inneighbors(net, j)
                i == j && continue
                for k in outneighbors(net, j)
                    (k == j || k == i) && continue
                    potential += 1
                    has_edge(net, i, k) && (realized += 1)
                end
            end
        end
        return potential > 0 ? realized / potential : 1.0
    elseif type == :local
        return Graphs.local_clustering_coefficient(_clustering_graph(net))
    else  # :average
        local_cc = Graphs.local_clustering_coefficient(_clustering_graph(net))
        valid = filter(!isnan, local_cc)
        return isempty(valid) ? 0.0 : mean(valid)
    end
end

# Undirected networks are stored as symmetric digraphs; convert to a
# SimpleGraph so local clustering sees single-counted degrees (otherwise
# the k(k-1) denominator is inflated ~4x).
function _clustering_graph(net)
    g = net.graph
    return is_directed(net) ? g : Graphs.SimpleGraph(g)
end

"""
    gtrans(net; kwargs...) -> Float64

Alias for [`transitivity`](@ref) (R `sna` compatibility).
"""
gtrans(net; kwargs...) = transitivity(net; kwargs...)

"""
    dyad_census(net) -> NamedTuple

Compute the dyad census (mutual, asymmetric, null counts).

# Returns
NamedTuple with fields:
- `mutual::Int`: Number of mutual dyads
- `asymmetric::Int`: Number of asymmetric dyads
- `null::Int`: Number of null dyads
"""
function dyad_census(net)
    n = nv(net)
    mutual = 0
    asymmetric = 0
    null = 0

    for i in 1:n
        for j in (i+1):n
            has_ij = has_edge(net, i, j)
            has_ji = has_edge(net, j, i)

            if has_ij && has_ji
                mutual += 1
            elseif has_ij || has_ji
                asymmetric += 1
            else
                null += 1
            end
        end
    end

    return (mutual=mutual, asymmetric=asymmetric, null=null)
end

"""
    triad_census(net) -> Vector{Int}

Compute the triad census.

For directed networks, returns the 16-element Davis–Leinhardt census over
the M-A-N (mutual/asymmetric/null) isomorphism classes, in the standard
order used by R `sna::triad.census`:

    003, 012, 102, 021D, 021U, 021C, 111D, 111U, 030T, 030C,
    201, 120D, 120U, 120C, 210, 300

For undirected networks, returns the 4-element census by triad edge count
(0, 1, 2, 3 edges).

Uses the edge-driven Batagelj–Mrvar (2001) algorithm: only triads containing
at least one tie are enumerated (each exactly once, from its lowest-labeled
connected pair), so the cost is `O(Σ_(u,v)∈E (deg(u)+deg(v)))` rather than
`O(n³)`; the empty-triad count is recovered by subtraction from `C(n,3)`.
"""
function triad_census(net)
    is_directed(net) || return _triad_census_undirected(net)

    n = nv(net)
    census = zeros(Int, 16)
    nbrs = [_union_neighborhood(net, v) for v in 1:n]
    buf = Int[]

    for v in 1:n
        for u in nbrs[v]
            u > v || continue
            # Third vertices attached to the dyad {v, u}
            _sorted_union!(buf, nbrs[v], nbrs[u], v, u)
            # Triads whose third vertex is attached to neither v nor u have
            # (v,u) as their only non-null dyad: 102 if mutual, 012 if asym
            dyad = (has_edge(net, v, u) && has_edge(net, u, v)) ? 3 : 2
            census[dyad] += n - length(buf) - 2
            for w in buf
                # Count each connected triad exactly once: from the edge to
                # its lowest-labeled attached pair
                if u < w || (v < w && w < u && !insorted(w, nbrs[v]))
                    census[_TRICODE_CLASS[_tricode(net, v, u, w)+1]] += 1
                end
            end
        end
    end

    census[1] = binomial(n, 3) - sum(@view census[2:16])
    return census
end

# Map from the 6-bit dyad code of a triad (v,u,w) to its 1-based position in
# the Davis–Leinhardt census order (Batagelj & Mrvar 2001, Table 1)
const _TRICODE_CLASS = (
    1, 2, 2, 3, 2, 4, 6, 8, 2, 6, 5, 7, 3, 8, 7, 11,
    2, 6, 4, 8, 5, 9, 9, 13, 6, 10, 9, 14, 7, 14, 12, 15,
    2, 5, 6, 7, 6, 9, 10, 14, 4, 9, 9, 12, 8, 13, 14, 15,
    3, 7, 8, 11, 7, 12, 14, 15, 8, 14, 13, 15, 11, 15, 15, 16)

# 6-bit arc code of the triad (v, u, w)
function _tricode(net, v::Int, u::Int, w::Int)
    code = 0
    has_edge(net, v, u) && (code += 1)
    has_edge(net, u, v) && (code += 2)
    has_edge(net, v, w) && (code += 4)
    has_edge(net, w, v) && (code += 8)
    has_edge(net, u, w) && (code += 16)
    has_edge(net, w, u) && (code += 32)
    return code
end

# Sorted union of in- and out-neighbors of v, excluding v itself
function _union_neighborhood(net, v::Int)
    if !is_directed(net)
        # Undirected storage is a symmetric digraph: out-neighbors suffice
        return Int[w for w in outneighbors(net, v) if w != v]
    end
    return _sorted_union!(Int[], inneighbors(net, v), outneighbors(net, v),
                          v, v)
end

# Merge two sorted vectors into `buf` (deduplicated), skipping two vertices
function _sorted_union!(buf::Vector{Int}, a::AbstractVector{<:Integer},
                        b::AbstractVector{<:Integer}, skip1::Int, skip2::Int)
    empty!(buf)
    i, j = 1, 1
    na, nb = length(a), length(b)
    while i <= na || j <= nb
        w = if i > na
            x = b[j]; j += 1; x
        elseif j > nb
            x = a[i]; i += 1; x
        elseif a[i] < b[j]
            x = a[i]; i += 1; x
        elseif a[i] > b[j]
            x = b[j]; j += 1; x
        else
            x = a[i]; i += 1; j += 1; x
        end
        (w == skip1 || w == skip2) || push!(buf, w)
    end
    return buf
end

# Edge-driven undirected census by triad edge count (0, 1, 2, 3)
function _triad_census_undirected(net)
    n = nv(net)
    census = zeros(Int, 4)
    nbrs = [_union_neighborhood(net, v) for v in 1:n]
    buf = Int[]

    for v in 1:n
        for u in nbrs[v]
            u > v || continue
            _sorted_union!(buf, nbrs[v], nbrs[u], v, u)
            census[2] += n - length(buf) - 2
            for w in buf
                if u < w || (v < w && w < u && !insorted(w, nbrs[v]))
                    m = 1 + (insorted(w, nbrs[v]) ? 1 : 0) +
                        (insorted(w, nbrs[u]) ? 1 : 0)
                    census[m+1] += 1
                end
            end
        end
    end

    census[1] = binomial(n, 3) - census[2] - census[3] - census[4]
    return census
end

"""
    hierarchy(net; measure=:reciprocity) -> Float64

Compute a graph hierarchy measure, following R `sna::hierarchy`.

# Arguments
- `measure::Symbol=:reciprocity`: Hierarchy measure
    - `:reciprocity`: `1 −` dyadic reciprocity (sna default)
    - `:krackhardt`: Krackhardt's hierarchy — the proportion of reachable
      (unordered) pairs whose reachability is asymmetric (one can reach the
      other but not vice versa)
"""
function hierarchy(net; measure::Symbol=:reciprocity)
    n = nv(net)
    if n <= 1
        return 0.0
    end

    if measure == :reciprocity
        return 1.0 - reciprocity(net; method=:dyadic)
    elseif measure != :krackhardt
        throw(ArgumentError("Unknown hierarchy measure: $measure"))
    end

    # Krackhardt hierarchy from the reachability matrix
    reach = reachability(net)

    hierarchical = 0
    total_connected = 0

    for i in 1:n
        for j in (i+1):n
            i_reaches_j = reach[i, j]
            j_reaches_i = reach[j, i]

            if i_reaches_j || j_reaches_i
                total_connected += 1
                if i_reaches_j != j_reaches_i  # One direction only
                    hierarchical += 1
                end
            end
        end
    end

    return total_connected > 0 ? hierarchical / total_connected : 0.0
end

"""
    efficiency(net) -> Float64

Compute Krackhardt's efficiency, following R `sna::efficiency`.

Efficiency is computed per weak component in arc (digraph) terms, exactly
as in sna: a component of size `nᶜ` requires `nᶜ − 1` arcs for weak
connection and can hold at most `nᶜ(nᶜ − 1)`; arcs beyond the requirement
are "excess" and efficiency is `1 − Σ excess / Σ maximum possible excess`.
Undirected edges count as two arcs (sna treats symmetric data as a
digraph).
"""
function efficiency(net)
    n = nv(net)
    if n <= 1
        return 1.0
    end

    comps = Graphs.weakly_connected_components(net.graph)

    excess = 0.0
    max_excess = 0.0

    for comp in comps
        n_c = length(comp)
        in_comp = falses(n)
        for v in comp
            in_comp[v] = true
        end

        # Arcs within the component (undirected edges count twice)
        m_c = 0
        for e in edges(net)
            if in_comp[src(e)] && in_comp[dst(e)]
                m_c += is_directed(net) ? 1 : 2
            end
        end

        excess += m_c - (n_c - 1)
        max_excess += n_c * (n_c - 1) - (n_c - 1)
    end

    return max_excess > 0 ? 1.0 - excess / max_excess : 1.0
end

"""
    connectedness(net) -> Float64

Compute Krackhardt's connectedness: the proportion of (unordered) vertex
pairs joined by a semipath, i.e. lying in the same weak component
(R `sna::connectedness`).
"""
function connectedness(net)
    n = nv(net)
    if n <= 1
        return 1.0
    end

    comps = Graphs.weakly_connected_components(net.graph)
    connected_pairs = sum(length(c) * (length(c) - 1) ÷ 2 for c in comps)
    total_pairs = n * (n - 1) ÷ 2

    return connected_pairs / total_pairs
end

"""
    mutuality(net) -> Int

Compute network mutuality: the number of mutual (reciprocated) dyads,
matching R `sna::mutuality` which returns a count, not a proportion.
"""
function mutuality(net)
    return dyad_census(net).mutual
end

"""
    component_dist(net) -> Vector{Int}

Return the distribution of component sizes.
"""
function component_dist(net)
    comps = Graphs.connected_components(net.graph)
    return sort([length(c) for c in comps], rev=true)
end

"""
    reachability(net) -> Matrix{Bool}

Compute the reachability matrix.

Entry (i,j) is true if vertex i can reach vertex j via some path.
"""
function reachability(net)
    n = nv(net)
    reach = falses(n, n)

    for i in 1:n
        # BFS from vertex i (index-pointer queue keeps pops O(1))
        visited = falses(n)
        queue = [i]
        visited[i] = true
        head = 1

        while head <= length(queue)
            v = queue[head]
            head += 1
            reach[i, v] = true

            for w in outneighbors(net, v)
                if !visited[w]
                    visited[w] = true
                    push!(queue, w)
                end
            end
        end
    end

    return reach
end
