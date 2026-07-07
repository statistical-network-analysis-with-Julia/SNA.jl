"""
Network-level measures and indices.

Provides global network statistics including density, reciprocity,
transitivity, and census functions.
"""

"""
    density(net) -> Float64

Compute network density (proportion of possible edges that exist).
Alias for network_density from Network.jl.
"""
density(net) = network_density(net)

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
        return Graphs.local_clustering_coefficient(net.graph)
    else  # :average
        local_cc = Graphs.local_clustering_coefficient(net.graph)
        valid = filter(!isnan, local_cc)
        return isempty(valid) ? 0.0 : mean(valid)
    end
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
"""
function triad_census(net)
    n = nv(net)

    if !is_directed(net)
        census = zeros(Int, 4)
        for i in 1:n, j in (i+1):n, k in (j+1):n
            m = (has_edge(net, i, j) ? 1 : 0) +
                (has_edge(net, i, k) ? 1 : 0) +
                (has_edge(net, j, k) ? 1 : 0)
            census[m+1] += 1
        end
        return census
    end

    census = zeros(Int, 16)
    for i in 1:n, j in (i+1):n, k in (j+1):n
        census[_triad_type(net, i, j, k)] += 1
    end
    return census
end

# Classify the directed triad {a, b, c} into one of the 16 Davis–Leinhardt
# M-A-N classes (1-based index into the standard census order).
function _triad_type(net, a::Int, b::Int, c::Int)
    # Dyad states and the set of asymmetric arcs
    mutual = 0
    asym_arcs = Tuple{Int,Int}[]
    mutual_pair = (0, 0)

    for (i, j) in ((a, b), (a, c), (b, c))
        y_ij = has_edge(net, i, j)
        y_ji = has_edge(net, j, i)
        if y_ij && y_ji
            mutual += 1
            mutual_pair = (i, j)
        elseif y_ij
            push!(asym_arcs, (i, j))
        elseif y_ji
            push!(asym_arcs, (j, i))
        end
    end

    A = length(asym_arcs)

    if mutual == 3
        return 16                     # 300
    elseif mutual == 2
        return A == 1 ? 15 : 11       # 210 : 201
    elseif mutual == 1
        if A == 0
            return 3                  # 102
        elseif A == 1
            # 111D: A<->B<-C (arc points into the mutual pair)
            # 111U: A<->B->C (arc points out of the mutual pair)
            s, d = asym_arcs[1]
            return (d == mutual_pair[1] || d == mutual_pair[2]) ? 7 : 8
        else  # A == 2
            s1, d1 = asym_arcs[1]
            s2, d2 = asym_arcs[2]
            s1 == s2 && return 12     # 120D (common source)
            d1 == d2 && return 13     # 120U (common sink)
            return 14                 # 120C (chain)
        end
    else  # mutual == 0
        if A == 0
            return 1                  # 003
        elseif A == 1
            return 2                  # 012
        elseif A == 2
            s1, d1 = asym_arcs[1]
            s2, d2 = asym_arcs[2]
            s1 == s2 && return 4      # 021D (out-star)
            d1 == d2 && return 5      # 021U (in-star)
            return 6                  # 021C (path)
        else  # A == 3
            # Cyclic if every vertex has out-degree 1 among the arcs
            sources = (asym_arcs[1][1], asym_arcs[2][1], asym_arcs[3][1])
            return allunique(sources) ? 10 : 9   # 030C : 030T
        end
    end
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
