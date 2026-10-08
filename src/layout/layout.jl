"""
Vertex layout algorithms for network visualization.

Each layout function returns an `n × 2` matrix of (x, y) coordinates.
"""

using LinearAlgebra
using Random

"""
    layout_circle(net) -> Matrix{Float64}

Place vertices evenly around the unit circle.

# Example
```julia
using SNA
layout_circle(network(4))            # the four points (±1, 0), (0, ±1)
```
"""
function layout_circle(net::AbstractNetwork)
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

Place vertices uniformly at random in the square [-1, 1]².

# Example
```julia
using SNA, Random
layout_random(network(5); rng=Xoshiro(1))
```
"""
function layout_random(net::AbstractNetwork; rng::Random.AbstractRNG=Random.default_rng())
    n = nv(net)
    return 2 .* rand(rng, n, 2) .- 1
end

"""
    layout_fruchterman_reingold(net; iterations=100, rng=Random.default_rng()) -> Matrix{Float64}

Force-directed layout (Fruchterman & Reingold 1991). Connected vertices
attract, all vertex pairs repel; positions settle over `iterations` steps
with a cooling schedule. Direction is ignored. The random starting
positions are drawn from `rng`.

# Example
```julia
using SNA, Random
net = network(5; directed=false)
add_edges!(net, [(1, 2), (2, 3), (3, 4), (4, 5)])
layout_fruchterman_reingold(net; rng=Xoshiro(1))
```
"""
function layout_fruchterman_reingold(net::AbstractNetwork; iterations::Int=100,
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
    layout_kamada_kawai(net; maxiter=500, tol=1e-6) -> Matrix{Float64}

Kamada–Kawai (1989) layout: positions that minimize the Kamada–Kawai
energy (the stress)

    Σ_{i<j} (‖xᵢ − xⱼ‖ − dᵢⱼ)² / dᵢⱼ²

where `dᵢⱼ` is the geodesic distance on the network with direction ignored.
The energy is minimized by stress majorization (Gansner, Koren & North
2004), which decreases it monotonically, starting from classical
multidimensional scaling of the distances; the result is deterministic.
Pairs in different components are placed at distance `max finite distance
+ 1`. Iteration stops after `maxiter` steps or when the relative decrease
of the energy falls below `tol`.

Before 0.2.0 this function returned the classical MDS starting point itself.

# Example
```julia
using SNA
path = network(4; directed=false)
add_edges!(path, [(1, 2), (2, 3), (3, 4)])
xy = layout_kamada_kawai(path)
size(xy)                             # (4, 2)
```
"""
function layout_kamada_kawai(net::AbstractNetwork; maxiter::Int=500, tol::Real=1e-6)
    n = nv(net)
    n == 0 && return Matrix{Float64}(undef, 0, 2)
    n == 1 && return zeros(1, 2)

    D = _geodist(_weak_digraph(net), n).gdist
    finite = filter(isfinite, D)
    cap = isempty(finite) || maximum(finite) == 0 ? 1.0 : maximum(finite) + 1.0
    D = map(d -> isfinite(d) ? d : cap, D)

    X = _classical_mds(D)
    W = [i == j ? 0.0 : 1 / D[i, j]^2 for i in 1:n, j in 1:n]
    # Laplacian of the weights; its pseudo-inverse is fixed across iterations
    LW = -copy(W)
    for i in 1:n
        LW[i, i] = sum(@view W[i, :])
    end
    LWp = pinv(LW)
    stress(Y) = sum(W[i, j] * (norm(Y[i, :] - Y[j, :]) - D[i, j])^2
                    for i in 1:n for j in (i+1):n)
    s_old = stress(X)
    for _ in 1:maxiter
        LZ = zeros(n, n)
        for i in 1:n, j in 1:n
            i == j && continue
            dist = norm(X[i, :] - X[j, :])
            LZ[i, j] = dist > 0 ? -W[i, j] * D[i, j] / dist : 0.0
        end
        for i in 1:n
            LZ[i, i] = -sum(@view LZ[i, :])
        end
        X = LWp * (LZ * X)
        s_new = stress(X)
        done = s_old - s_new <= tol * max(s_old, eps())
        s_old = s_new
        done && break
    end
    return X
end

# Classical (Torgerson) multidimensional scaling into two dimensions.
function _classical_mds(D::AbstractMatrix)
    n = size(D, 1)
    J = Matrix{Float64}(I, n, n) .- 1.0 / n
    B = -0.5 .* (J * (D .^ 2) * J)
    B = (B + transpose(B)) ./ 2  # guard symmetry against roundoff
    ev = eigen(Symmetric(B))
    order = sortperm(ev.values; rev=true)
    coords = Matrix{Float64}(undef, n, 2)
    for (c, idx) in enumerate(order[1:2])
        coords[:, c] = ev.vectors[:, idx] .* sqrt(max(ev.values[idx], 0.0))
    end
    return coords
end
