"""
Vertex layout algorithms for network visualization.

Each layout function returns an `n × 2` matrix of (x, y) coordinates.
"""

using Random

"""
    layout_circle(net) -> Matrix{Float64}

Place vertices evenly around the unit circle.
"""
function layout_circle(net)
    n = nv(net)
    coords = Matrix{Float64}(undef, n, 2)
    for v in 1:n
        θ = 2π * (v - 1) / max(n, 1)
        coords[v, 1] = cos(θ)
        coords[v, 2] = sin(θ)
    end
    return coords
end

"""
    layout_random(net; rng=Random.default_rng()) -> Matrix{Float64}

Place vertices uniformly at random in the unit square [-1, 1]².
"""
function layout_random(net; rng::Random.AbstractRNG=Random.default_rng())
    n = nv(net)
    return 2 .* rand(rng, n, 2) .- 1
end

"""
    layout_fruchterman_reingold(net; iterations=100, rng=Random.default_rng()) -> Matrix{Float64}

Force-directed layout (Fruchterman & Reingold 1991). Connected vertices
attract, all vertex pairs repel; positions settle over `iterations` steps
with a cooling schedule.
"""
function layout_fruchterman_reingold(net; iterations::Int=100,
                                     rng::Random.AbstractRNG=Random.default_rng())
    n = nv(net)
    n == 0 && return Matrix{Float64}(undef, 0, 2)
    n == 1 && return zeros(1, 2)

    # Symmetric adjacency for attraction
    adj = falses(n, n)
    for e in edges(net)
        adj[src(e), dst(e)] = true
        adj[dst(e), src(e)] = true
    end

    area = 4.0  # layout in [-1, 1]²
    k = sqrt(area / n)  # ideal distance
    pos = 2 .* rand(rng, n, 2) .- 1
    disp = zeros(n, 2)
    t = 0.5  # initial temperature

    for iter in 1:iterations
        fill!(disp, 0.0)

        # Repulsive forces between all pairs
        for i in 1:n
            for j in (i+1):n
                dx = pos[i, 1] - pos[j, 1]
                dy = pos[i, 2] - pos[j, 2]
                dist = max(sqrt(dx^2 + dy^2), 1e-9)
                f = k^2 / dist
                fx, fy = f * dx / dist, f * dy / dist
                disp[i, 1] += fx
                disp[i, 2] += fy
                disp[j, 1] -= fx
                disp[j, 2] -= fy
            end
        end

        # Attractive forces along edges
        for i in 1:n
            for j in (i+1):n
                adj[i, j] || continue
                dx = pos[i, 1] - pos[j, 1]
                dy = pos[i, 2] - pos[j, 2]
                dist = max(sqrt(dx^2 + dy^2), 1e-9)
                f = dist^2 / k
                fx, fy = f * dx / dist, f * dy / dist
                disp[i, 1] -= fx
                disp[i, 2] -= fy
                disp[j, 1] += fx
                disp[j, 2] += fy
            end
        end

        # Limit displacement by temperature and update
        for i in 1:n
            d = max(sqrt(disp[i, 1]^2 + disp[i, 2]^2), 1e-9)
            step = min(d, t)
            pos[i, 1] += disp[i, 1] / d * step
            pos[i, 2] += disp[i, 2] / d * step
        end

        t *= 0.95  # cool
    end

    return pos
end

"""
    layout_kamada_kawai(net) -> Matrix{Float64}

Distance-based layout in the Kamada–Kawai (1989) spirit: classical
multidimensional scaling of the geodesic distance matrix, so that Euclidean
distances in the plane approximate graph distances. Unreachable pairs are
assigned distance `max finite distance + 1`.
"""
function layout_kamada_kawai(net)
    n = nv(net)
    n == 0 && return Matrix{Float64}(undef, 0, 2)
    n == 1 && return zeros(1, 2)

    D = geodesic_distance(net)
    # Symmetrize (directed distances) and cap unreachable pairs
    D = min.(D, transpose(D))
    finite = filter(isfinite, D)
    cap = isempty(finite) ? 1.0 : maximum(finite) + 1.0
    D = map(d -> isfinite(d) ? d : cap, D)

    # Classical MDS: double-center the squared distances
    D2 = D .^ 2
    J = Matrix{Float64}(I, n, n) .- 1.0 / n
    B = -0.5 .* (J * D2 * J)
    B = (B + transpose(B)) ./ 2  # guard symmetry against roundoff

    ev = eigen(Symmetric(B))
    # Two largest eigenvalues/vectors
    order = sortperm(ev.values; rev=true)
    coords = Matrix{Float64}(undef, n, 2)
    for (c, idx) in enumerate(order[1:2])
        λ = max(ev.values[idx], 0.0)
        coords[:, c] = ev.vectors[:, idx] .* sqrt(λ)
    end

    return coords
end
