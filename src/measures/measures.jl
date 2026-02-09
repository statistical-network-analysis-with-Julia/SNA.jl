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

# Alias for R compatibility
gden(net) = density(net)

"""
    reciprocity(net; method=:dyadic) -> Float64

Compute network reciprocity (proportion of mutual ties).

# Arguments
- `net`: Network object (should be directed)
- `method::Symbol=:dyadic`: Reciprocity measure type
    - `:dyadic`: Proportion of dyads that are mutual
    - `:edgewise`: Proportion of edges that are reciprocated
"""
function reciprocity(net; method::Symbol=:dyadic)
    if !is_directed(net)
        return 1.0  # All ties are mutual in undirected networks
    end

    mutual = 0
    asymmetric = 0
    null = 0

    n = nv(net)
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

    if method == :dyadic
        # Proportion of non-null dyads that are mutual
        total_ties = mutual + asymmetric
        return total_ties > 0 ? mutual / total_ties : 0.0
    else  # :edgewise
        # Proportion of edges that are reciprocated
        total_edges = 2 * mutual + asymmetric
        return total_edges > 0 ? (2 * mutual) / total_edges : 0.0
    end
end

# Alias for R compatibility
grecip(net; kwargs...) = reciprocity(net; kwargs...)

"""
    transitivity(net; type=:global) -> Float64

Compute network transitivity (clustering coefficient).

# Arguments
- `net`: Network object
- `type::Symbol=:global`: Type of transitivity
    - `:global`: Global clustering coefficient
    - `:local`: Vector of local clustering coefficients
    - `:average`: Average of local clustering coefficients
"""
function transitivity(net; type::Symbol=:global)
    if type == :global
        return Graphs.global_clustering_coefficient(net.graph)
    elseif type == :local
        return Graphs.local_clustering_coefficient(net.graph)
    else  # :average
        local_cc = Graphs.local_clustering_coefficient(net.graph)
        valid = filter(!isnan, local_cc)
        return isempty(valid) ? 0.0 : mean(valid)
    end
end

# Alias for R compatibility
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

Compute the triad census (16-element vector for directed networks).

Returns counts of each of the 16 isomorphism classes of triads:
003, 012, 102, 021D, 021U, 021C, 111D, 111U, 030T, 030C,
201, 120D, 120U, 120C, 210, 300
"""
function triad_census(net)
    # Use Graphs.jl triangles as basis, but full triad census is complex
    # Simplified implementation counting triangles
    n = nv(net)
    census = zeros(Int, 16)

    # This is a simplified placeholder - full triad census requires
    # checking all n*(n-1)*(n-2)/6 triads
    @warn "triad_census: simplified implementation"

    for i in 1:n
        for j in (i+1):n
            for k in (j+1):n
                # Count edges in this triad
                edges = 0
                edges += has_edge(net, i, j) ? 1 : 0
                edges += has_edge(net, j, i) ? 1 : 0
                edges += has_edge(net, i, k) ? 1 : 0
                edges += has_edge(net, k, i) ? 1 : 0
                edges += has_edge(net, j, k) ? 1 : 0
                edges += has_edge(net, k, j) ? 1 : 0

                # Map to census class (simplified)
                if edges == 0
                    census[1] += 1  # 003
                elseif edges == 6
                    census[16] += 1  # 300
                else
                    census[edges + 1] += 1  # Approximate
                end
            end
        end
    end

    return census
end

"""
    hierarchy(net) -> Float64

Compute Krackhardt's hierarchy measure.

Measures the proportion of pairs with a strictly hierarchical relationship
(one can reach the other but not vice versa).
"""
function hierarchy(net)
    n = nv(net)
    if n <= 1
        return 0.0
    end

    # Compute reachability matrix
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

Compute network efficiency (1 minus the proportion of edges in excess of a tree).
"""
function efficiency(net)
    n = nv(net)
    m = ne(net)

    if n <= 1
        return 1.0
    end

    min_edges = n - 1  # Minimum for connected graph (tree)
    max_edges = is_directed(net) ? n * (n - 1) : n * (n - 1) ÷ 2

    if max_edges == min_edges
        return 1.0
    end

    return 1.0 - (m - min_edges) / (max_edges - min_edges)
end

"""
    connectedness(net) -> Float64

Compute network connectedness (proportion of pairs where at least one can reach the other).
"""
function connectedness(net)
    n = nv(net)
    if n <= 1
        return 1.0
    end

    reach = reachability(net)
    connected = 0

    for i in 1:n
        for j in (i+1):n
            if reach[i, j] || reach[j, i]
                connected += 1
            end
        end
    end

    total_pairs = n * (n - 1) ÷ 2
    return connected / total_pairs
end

"""
    mutuality(net) -> Float64

Compute network mutuality (proportion of pairs with symmetric relationships).
"""
function mutuality(net)
    census = dyad_census(net)
    total_connected = census.mutual + census.asymmetric
    return total_connected > 0 ? census.mutual / total_connected : 0.0
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
        # BFS from vertex i
        visited = falses(n)
        queue = [i]
        visited[i] = true

        while !isempty(queue)
            v = popfirst!(queue)
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
