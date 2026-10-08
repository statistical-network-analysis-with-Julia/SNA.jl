"""
Centrality measures for network analysis.

Vertex-level indices with R `sna`'s names and semantics: `degreecent`
(sna's `degree`), `betweenness`, `closeness`, `evcent`, `bonpow`,
`infocent`, `flowbet`, and Freeman graph `centralization`.

None of these functions is a method of a Graphs.jl function: Graphs.jl's
`degree_centrality`, `betweenness_centrality`, `closeness_centrality`, … keep
Graphs.jl's own semantics on a `Network`.
"""

using LinearAlgebra

# ---------------------------------------------------------------------------
# Degree
# ---------------------------------------------------------------------------

"""
    degreecent(net; cmode=:freeman, rescale=false, normalized=false, diag=false,
               ignore_eval=true, attr=:weight, missing=:error) -> Vector{Float64}

Degree centrality of every vertex: R `sna::degree`. (The function is named
`degreecent`, after sna's own naming of centrality indices such as `infocent`,
because `degree` is Graphs.jl's vertex degree, which SNA re-exports unchanged.)

# Arguments
- `cmode`: `:freeman` (in + out degree, sna's default), `:indegree` or
  `:outdegree`. Ignored for undirected networks, where each edge adds 1 to
  both endpoints (sna's `gmode="graph"`).
- `rescale=true`: divide by the sum of the scores, as sna's `rescale=TRUE`.
- `normalized=true`: divide by the largest possible degree, `n−1` (or
  `2(n−1)` for Freeman degree of a directed network), plus 1 with `diag=true`.
  This is not an sna option.
- `diag=false`: self-loops are ignored (sna's `diag=FALSE`); `diag=true`
  counts a loop once.
- `ignore_eval=true`: binary ties. `ignore_eval=false, attr=:weight` sums the
  numeric edge attribute `attr` (sna's default `ignore.eval=FALSE`).
- `missing=:error`: a network with masked (unobserved) dyads is refused;
  `missing=:face` counts each masked dyad at its stored face value (see
  `NetworkCore.require_observed`).

# Example
```julia
using SNA
net = network(4)
add_edges!(net, [(1, 2), (1, 3), (2, 3), (3, 4)])
degreecent(net)                      # [2.0, 2.0, 3.0, 1.0]  (Freeman)
degreecent(net; cmode=:indegree)     # [0.0, 1.0, 2.0, 1.0]
```
"""
function degreecent(net::AbstractNetwork; cmode::Symbol=:freeman, rescale::Bool=false,
                    normalized::Bool=false, diag::Bool=false, ignore_eval::Bool=true,
                    attr::Symbol=:weight, missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="degreecent")
    cmode in (:freeman, :indegree, :outdegree) ||
        throw(ArgumentError("degreecent: cmode must be :freeman, :indegree or :outdegree " *
                            "(R sna::degree); got :$cmode"))
    _check_scaling(rescale, normalized, "degreecent")
    n = nv(net)
    directed = is_directed(net)
    out = !directed || cmode == :outdegree
    total = directed && cmode == :freeman
    centrality = zeros(Float64, n)
    if !ignore_eval
        A = _sociomatrix(net; ignore_eval, attr, diag)
        centrality = out ? vec(sum(A; dims=2)) :
                     !total ? vec(sum(A; dims=1)) :
                     vec(sum(A; dims=1)) + vec(sum(A; dims=2))
        # sna's Freeman degree counts a diagonal entry once, not twice.
        total && diag && (centrality .-= LinearAlgebra.diag(A))
    else
        for v in vertices(net)
            centrality[v] = out ? length(outneighbors(net, v)) :
                            !total ? length(inneighbors(net, v)) :
                            length(inneighbors(net, v)) + length(outneighbors(net, v))
            if has_edge(net, v, v)
                if !diag
                    centrality[v] -= total ? 2 : 1
                elseif total
                    centrality[v] -= 1
                end
            end
        end
    end
    if normalized && n > 1
        centrality ./= (total ? 2 : 1) * (n - 1) + diag
    end
    rescale && (centrality ./= sum(centrality))
    return centrality
end

function _check_scaling(rescale::Bool, normalized::Bool, fname::String)
    rescale && normalized &&
        throw(ArgumentError("$fname: choose rescale=true (scores sum to 1, as in R sna) " *
                            "or normalized=true (divide by the largest possible score), not both"))
    return nothing
end

# ---------------------------------------------------------------------------
# Betweenness
# ---------------------------------------------------------------------------

"""
    betweenness(net; cmode=:directed, rescale=false, normalized=false,
                missing=:error) -> Vector{Float64}

Shortest-path betweenness of every vertex: R `sna::betweenness`. The default
is the raw score, the number of geodesics between other pairs that pass
through the vertex (each unordered pair counted once on an undirected
network).

# Arguments
- `cmode`: `:directed` (sna's default) uses directed geodesics; `:undirected`
  symmetrizes a directed network first (an arc in either direction is a tie).
  sna's other path-weighting modes (`endpoints`, `proximalsrc`, …) are not
  implemented and raise an `ArgumentError`.
- `rescale=true`: divide by the sum of the scores (sna's `rescale=TRUE`).
- `normalized=true`: Freeman's normalization, divide by the number of pairs
  of other vertices, `(n−1)(n−2)` directed or `(n−1)(n−2)/2` undirected,
  giving scores in [0, 1]. This is not an sna option.
- `missing=:error`: masked (unobserved) dyads are refused; `missing=:face`
  builds the paths from their stored face values.

# Example
```julia
using SNA
star = network(4; directed=false)
add_edges!(star, [(1, 2), (1, 3), (1, 4)])
betweenness(star)                    # [3.0, 0.0, 0.0, 0.0]
betweenness(star; normalized=true)   # [1.0, 0.0, 0.0, 0.0]
```
"""
function betweenness(net::AbstractNetwork; cmode::Symbol=:directed, rescale::Bool=false,
                     normalized::Bool=false, missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="betweenness")
    cmode in (:directed, :undirected) ||
        throw(ArgumentError("betweenness: cmode must be :directed or :undirected; sna's " *
                            "other betweenness variants (endpoints, proximal, length-scaled) " *
                            "are not implemented"))
    _check_scaling(rescale, normalized, "betweenness")
    n = nv(net)
    undirected = !is_directed(net) || cmode == :undirected
    g = is_directed(net) && cmode == :undirected ? _weak_digraph(net) : _graph(net)
    # Graphs.jl counts ordered pairs on the digraph storage: an undirected
    # geodesic is counted once per direction, so halve (sna does the same).
    bc = Graphs.betweenness_centrality(g; normalize=false)
    undirected && (bc ./= 2)
    if normalized && n > 2
        bc ./= undirected ? (n - 1) * (n - 2) / 2 : (n - 1) * (n - 2)
    end
    rescale && (bc ./= sum(bc))
    return bc
end

# The weak symmetrization of a network as a digraph that stores both arcs of
# every tie (loops dropped): the storage convention of an undirected Network.
function _weak_digraph(net)
    g = Graphs.SimpleDiGraph(nv(net))
    for e in edges(net)
        i, j = src(e), dst(e)
        i == j && continue
        Graphs.add_edge!(g, i, j)
        Graphs.add_edge!(g, j, i)
    end
    return g
end

# ---------------------------------------------------------------------------
# Closeness
# ---------------------------------------------------------------------------

"""
    closeness(net; cmode=:directed, rescale=false, missing=:error) -> Vector{Float64}

Closeness centrality of every vertex: R `sna::closeness`.

# Arguments
- `cmode`:
    - `:directed` (default): Freeman closeness `(n−1) / Σⱼ d(i,j)` over
      directed geodesics. A vertex that cannot reach every other vertex scores
      0, as in sna, so every score of a disconnected network is 0.
    - `:undirected`: the same on the weak symmetrization.
    - `:suminvdir` / `:suminvundir`: `Σⱼ 1/d(i,j) / (n−1)`, with `1/∞ = 0`,
      directed or symmetrized; defined on disconnected networks.
    - `:gil_schmidt`: Gil–Schmidt power index, `Σⱼ 1/d(i,j)` over the
      vertices `i` reaches divided by their number (0 if it reaches none).

  On an undirected network `:directed` and `:undirected` (and the two
  `suminv` modes) coincide.
- `rescale=true`: divide by the sum of the scores (sna's `rescale=TRUE`).
- `missing=:error`: masked (unobserved) dyads are refused; `missing=:face`
  builds the paths from their stored face values.

Graphs.jl's `closeness_centrality` scales each vertex by the share of the
network it reaches instead; that convention is Graphs.jl's, not sna's.

# Example
```julia
using SNA
path = network(3; directed=false)
add_edges!(path, [(1, 2), (2, 3)])
closeness(path)                       # [0.6667, 1.0, 0.6667]
closeness(path; cmode=:suminvundir)   # [0.75, 1.0, 0.75]
```
"""
function closeness(net::AbstractNetwork; cmode::Symbol=:directed, rescale::Bool=false,
                   missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="closeness")
    cmode = _closeness_cmode(net, cmode)
    n = nv(net)
    g = is_directed(net) && cmode in (:undirected, :suminvundir) ? _weak_digraph(net) :
        _graph(net)
    D = _geodist(g, n).gdist
    clo = zeros(n)
    for i in 1:n
        if cmode in (:directed, :undirected)
            total = 0.0
            for j in 1:n
                i == j || (total += D[i, j])
            end
            clo[i] = (n - 1) / total      # 0 when unreachable; NaN for n = 1, as sna
        elseif cmode in (:suminvdir, :suminvundir)
            s = 0.0
            for j in 1:n
                i == j || (s += 1 / D[i, j])
            end
            clo[i] = s / (n - 1)
        else  # :gil_schmidt
            s, r = 0.0, 0
            for j in 1:n
                (i == j || !isfinite(D[i, j])) && continue
                s += 1 / D[i, j]
                r += 1
            end
            clo[i] = r == 0 ? 0.0 : s / r
        end
    end
    rescale && (clo ./= sum(clo))
    return clo
end

function _closeness_cmode(net, cmode::Symbol)
    cmode in (:directed, :undirected, :suminvdir, :suminvundir, :gil_schmidt) ||
        throw(ArgumentError("closeness: cmode must be :directed, :undirected, :suminvdir, " *
                            ":suminvundir or :gil_schmidt; got :$cmode"))
    is_directed(net) && return cmode
    return cmode == :directed ? :undirected : cmode == :suminvdir ? :suminvundir : cmode
end

# ---------------------------------------------------------------------------
# Eigenvector
# ---------------------------------------------------------------------------

"""
    evcent(net; rescale=false, tol=1e-10, diag=false, ignore_eval=true,
           attr=:weight, missing=:error) -> Vector{Float64}

Eigenvector centrality: R `sna::evcent`. The scores are the principal right
eigenvector of the adjacency matrix (a directed network weights an actor by
the centrality of the actors it sends ties to), oriented non-negative with
unit Euclidean norm.

The eigenvector of the largest *real* eigenvalue is taken by a direct
eigen-decomposition, which is sna's `use.eigen=TRUE` (up to sign). sna's
default power iteration does not converge on bipartite or other periodic
networks and then returns a vector that is not an eigenvector; this function
does not have that failure. When the leading eigenvalue is repeated the
eigenvector is not unique and a warning says so. When the leading eigenvalue
is 0 (an empty network, or a directed acyclic one) there is no positive
eigenvector and every score is `NaN`, as sna returns.

# Arguments
- `rescale=true`: divide by the sum of the scores (sna's `rescale=TRUE`).
- `tol`: relative tolerance of the repeated-eigenvalue check.
- `diag`, `ignore_eval`, `attr`: as in [`degreecent`](@ref); weights must be
  non-negative.
- `missing=:error`: masked dyads are refused; `:face` uses their stored face values.

# Example
```julia
using SNA
path = network(3; directed=false)
add_edges!(path, [(1, 2), (2, 3)])
evcent(path)                         # [0.5, 0.7071, 0.5]
```
"""
function evcent(net::AbstractNetwork; rescale::Bool=false, tol::Float64=1e-10,
                ignore_eval::Bool=true, attr::Symbol=:weight, diag::Bool=false,
                missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="evcent")
    n = nv(net)
    n == 0 && return Float64[]
    A = _sociomatrix(net; ignore_eval, attr, diag)
    any(<(0), A) && throw(ArgumentError("evcent requires nonnegative weights"))
    # Largest real eigenvalue, not largest modulus: bipartite graphs also have
    # -rho, on which an unshifted power iteration oscillates without warning.
    # A non-negative matrix is nilpotent exactly when its digraph has no
    # cycle; its eigenvalues are then all 0, so there is no Perron vector.
    # (Tested on the graph, not the computed eigenvalues: a nilpotent Jordan
    # block of size k perturbs them by about eps^(1/k).)
    _has_cycle(A) || return fill(NaN, n)
    F = eigen(A)
    i = argmax(real.(F.values))
    rho = real(F.values[i])
    x = abs.(real.(F.vectors[:, i]))
    count(v -> abs(v - rho) <= tol * max(1.0, abs(rho)), F.values) > 1 &&
        @warn "evcent: the leading eigenvalue is repeated, so the eigenvector is not unique"
    x ./= norm(x)
    rescale && (x ./= sum(x))
    return x
end

function _has_cycle(A::AbstractMatrix)
    g = Graphs.SimpleDiGraph(size(A, 1))
    for j in axes(A, 2), i in axes(A, 1)
        if A[i, j] != 0
            i == j && return true          # a loop kept by diag=true is a cycle
            Graphs.add_edge!(g, i, j)
        end
    end
    return Graphs.is_cyclic(g)
end

# ---------------------------------------------------------------------------
# Bonacich power
# ---------------------------------------------------------------------------

"""
    bonpow(net; exponent=1.0, rescale=false, tol=1e-7, diag=false,
           ignore_eval=true, attr=:weight, missing=:error) -> Vector{Float64}

Bonacich power centrality: R `sna::bonpow`,

    c = α (I − β A)⁻¹ A 𝟙,   with α chosen so that Σᵢ cᵢ² = n.

# Arguments
- `exponent`: the attenuation β. The power series converges for
  `|β| < 1/λ_max`; a positive β rewards ties to well-connected actors, a
  negative β ties to poorly connected ones.
- `rescale=true`: divide by the sum of the scores (sna's `rescale=TRUE`).
- `tol`: `(I − βA)` is treated as singular when the reciprocal of its
  1-norm condition number is below `tol`, the test R's `solve(tol=)` applies.
  A singular system gives `NaN` scores and a warning. (A determinant test,
  used before 0.2.0, underflows on large well-conditioned systems.)
- `diag`, `ignore_eval`, `attr`: as in [`degreecent`](@ref).
- `missing=:error`: masked dyads are refused; `:face` uses their stored face values.

A network with no ties has no scale `α` and returns `NaN` scores, as sna does.

# Example
```julia
using SNA
star = network(4; directed=false)
add_edges!(star, [(1, 2), (1, 3), (1, 4)])
bonpow(star; exponent=0.2)
```
"""
function bonpow(net::AbstractNetwork; exponent::Real=1.0, rescale::Bool=false,
                tol::Real=1e-7, ignore_eval::Bool=true,
                attr::Symbol=:weight, diag::Bool=false, missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="bonpow")
    n = nv(net)
    A = _sociomatrix(net; ignore_eval, attr, diag)
    M = Matrix{Float64}(I, n, n) - Float64(exponent) * A
    if n > 0
        F = lu(M; check=false)
        rc = issuccess(F) ? _rcond(F, M) : 0.0
        if rc < tol
            @warn "bonpow: (I − βA) is computationally singular (reciprocal condition " *
                  "number $(rc) < tol = $(tol)); choose |exponent| < 1/λ_max"
            return fill(NaN, n)
        end
        c = F \ (A * ones(n))
    else
        c = Float64[]
    end
    c .*= sqrt(n / sum(abs2, c))     # NaN when every score is 0, as in sna
    rescale && (c ./= sum(c))
    return c
end

# Reciprocal 1-norm condition number estimate from an LU factorization: the
# quantity LAPACK's dgecon returns and R's `solve(tol=)` compares with `tol`.
function _rcond(F::LinearAlgebra.LU, M::AbstractMatrix{Float64})
    anorm = opnorm(M, 1)
    anorm == 0 && return 0.0
    return LinearAlgebra.LAPACK.gecon!('1', copy(F.factors), anorm)
end

# ---------------------------------------------------------------------------
# Information centrality
# ---------------------------------------------------------------------------

"""
    infocent(net; cmode=:weak, rescale=false, diag=false, ignore_eval=true,
             attr=:weight, missing=:error) -> Vector{Float64}

Stephenson–Zelen information centrality: R `sna::infocent`. Information
centrality is defined on undirected networks: a directed network is first
symmetrized, with `cmode=:weak` (a tie in either direction, sna's default)
or `cmode=:strong` (mutual ties only). Isolates score 0, and, as in sna, `n`
in the formula counts every vertex, isolates included.

# Example
```julia
using SNA
path = network(4; directed=false)
add_edges!(path, [(1, 2), (2, 3), (3, 4)])
infocent(path)
```
"""
function infocent(net::AbstractNetwork; cmode::Symbol=:weak, rescale::Bool=false,
                  diag::Bool=false, ignore_eval::Bool=true, attr::Symbol=:weight,
                  missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="infocent")
    cmode in (:weak, :strong) ||
        throw(ArgumentError("infocent: cmode must be :weak or :strong; got :$cmode"))
    m = _sociomatrix(net; ignore_eval, attr, diag)
    n = size(m, 1)
    if m != m'
        m = cmode == :weak ? max.(m, m') : min.(m, m')
    end
    # sna sets the diagonal to NA before summing rows unless diag=TRUE.
    rowsum = diag ? vec(sum(m; dims=2)) : vec(sum(m; dims=2)) .- LinearAlgebra.diag(m)
    isolate = [all(iszero, (m[i, j] for j in 1:n if diag || j != i)) &&
               all(iszero, (m[j, i] for j in 1:n if diag || j != i)) for i in 1:n]
    ix = findall(!, isolate)
    cent = zeros(n)
    if !isempty(ix)
        k = length(ix)
        B = Matrix{Float64}(undef, k, k)
        for (a, i) in enumerate(ix), (b, j) in enumerate(ix)
            B[a, b] = a == b ? 1 + rowsum[i] : (m[i, j] == 0 ? 1.0 : 1 - m[i, j])
        end
        F = lu(B; check=false)
        issuccess(F) ||
            throw(ArgumentError("infocent: the information matrix is singular"))
        C = inv(F)
        Tr = tr(C)
        R = vec(sum(C; dims=2))
        cent[ix] = 1 ./ (LinearAlgebra.diag(C) .+ (Tr .- 2 .* R) ./ n)
    end
    rescale && (cent ./= sum(cent))
    return cent
end

# ---------------------------------------------------------------------------
# Flow betweenness
# ---------------------------------------------------------------------------

"""
    flowbet(net; ignore_eval=true, attr=:weight, diag=false, missing=:error) -> Vector{Float64}

Freeman flow betweenness of every vertex: R `sna::flowbet`,

    f(v) = Σ_{i,j ≠ v} [maxflow(i → j) − maxflow(i → j | v removed)]

with edge capacities from the adjacency matrix. Pairs are ordered for
directed networks and unordered for undirected networks, and raw scores are
returned (sna's default). Ties are binary by default (`ignore_eval=true`);
sna's default `ignore.eval=FALSE` uses edge values, which
`ignore_eval=false, attr=:weight` reproduces.

`missing=:error` refuses masked dyads; `missing=:face` takes their stored
face values as capacities.

# Example
```julia
using SNA
path = network(3; directed=false)
add_edges!(path, [(1, 2), (2, 3)])
flowbet(path)                        # [0.0, 1.0, 0.0]
```
"""
function flowbet(net::AbstractNetwork; ignore_eval::Bool=true, attr::Symbol=:weight,
                 diag::Bool=false, missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="flowbet")
    n = nv(net)
    A = _sociomatrix(net; ignore_eval, attr, diag)
    any(<(0), A) && throw(ArgumentError("flow capacities must be nonnegative"))
    directed = is_directed(net)

    fb = zeros(n)
    for i in 1:n
        j_range = directed ? (1:n) : (i+1:n)
        for j in j_range
            i == j && continue
            base = _maxflow(A, n, i, j, 0)
            base == 0.0 && continue
            for v in 1:n
                (v == i || v == j) && continue
                fb[v] += base - _maxflow(A, n, i, j, v)
            end
        end
    end

    return fb
end

# ---------------------------------------------------------------------------
# Centralization
# ---------------------------------------------------------------------------

"""
    centralization(net, f; normalized=true, missing=:error, kwargs...) -> Float64
    centralization(net, measure::Symbol; mode=:total, normalized=true, missing=:error) -> Float64

Freeman graph centralization of a vertex index: R `sna::centralization`,

    C = Σᵢ (c_max − cᵢ) / C_max

where `cᵢ = f(net; kwargs...)` and `C_max` is the largest value of the
numerator over all networks of the same size (sna's `tmaxdev`). Keyword
arguments other than `normalized` and `missing` are passed to `f`, as R
passes `...` to `FUN`; the theoretical maximum follows the same `cmode`.

`f` is one of [`degreecent`](@ref), [`betweenness`](@ref),
[`closeness`](@ref), [`evcent`](@ref) or [`bonpow`](@ref). Any other function
taking `(net; missing, kwargs...)` works with `normalized=false`, which
returns the raw deviation sum. The `Symbol` form accepts `:degree`,
`:betweenness`, `:closeness` and `:eigenvector`, with `mode` = `:total`,
`:in` or `:out` selecting the degree type.

A network too small for the maximum to be positive (for instance two
vertices) gives `NaN`, as sna does.

`missing=:error` refuses a network with masked dyads; `missing=:face` is
forwarded to `f`.

# Example
```julia
using SNA
star = network(5; directed=false)
add_edges!(star, [(1, 2), (1, 3), (1, 4), (1, 5)])
centralization(star, degreecent)     # 1.0
centralization(star, betweenness)    # 1.0
```
"""
function centralization(net::AbstractNetwork, f::Function; normalized::Bool=true,
                        missing::Symbol=:error, kwargs...)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="centralization")
    cv = f(net; missing=policy, kwargs...)
    isempty(cv) && return NaN
    cent = sum(maximum(cv) .- cv)
    normalized || return cent
    return cent / _tmaxdev(f, net; kwargs...)
end

function centralization(net::AbstractNetwork, measure::Symbol; mode::Symbol=:total,
                        normalized::Bool=true, missing::Symbol=:error)
    if measure == :degree
        mode in (:total, :in, :out) ||
            throw(ArgumentError("mode must be :total, :in or :out"))
        cmode = mode == :total ? :freeman : mode == :in ? :indegree : :outdegree
        return centralization(net, degreecent; cmode, normalized, missing)
    end
    f = measure == :betweenness ? betweenness :
        measure == :closeness ? closeness :
        measure == :eigenvector ? evcent :
        throw(ArgumentError("Unknown centralization measure: $measure " *
                            "(use :degree, :betweenness, :closeness, or :eigenvector, " *
                            "or pass the function itself)"))
    return centralization(net, f; normalized, missing)
end

# sna's tmaxdev values: the largest possible deviation sum for a network of
# this size and directedness, under the same cmode.
function _tmaxdev(::typeof(degreecent), net; cmode::Symbol=:freeman, diag::Bool=false,
                  kwargs...)
    n = nv(net)
    if !is_directed(net)
        return Float64((n - 1) * (n - 2 + diag))
    end
    return cmode == :freeman ? Float64((n - 1) * (2 * (n - 1) - 2 + diag)) :
                               Float64((n - 1) * (n - 1 + diag))
end
function _tmaxdev(::typeof(betweenness), net; cmode::Symbol=:directed, kwargs...)
    n = nv(net)
    undirected = !is_directed(net) || cmode == :undirected
    return undirected ? (n - 1)^2 * (n - 2) / 2 : Float64((n - 1)^2 * (n - 2))
end
function _tmaxdev(::typeof(closeness), net; cmode::Symbol=:directed, kwargs...)
    n = nv(net)
    cmode = _closeness_cmode(net, cmode)
    return cmode == :directed ? (n - 1) * (1 - 1 / n) :
           cmode == :undirected ? (n - 2) * (n - 1) / (2n - 3) :
           cmode == :suminvdir ? Float64((n - 1)^2) :
           cmode == :suminvundir ? n - 1 - n / 2 :
           is_directed(net) ? Float64(n - 1) : (n - 2) / 2
end
_tmaxdev(::typeof(evcent), net; kwargs...) =
    is_directed(net) ? Float64(nv(net) - 1) : sqrt(2) / 2 * (nv(net) - 2)
_tmaxdev(::typeof(bonpow), net; kwargs...) =
    is_directed(net) ? sqrt(nv(net)) * (nv(net) - 1) : (nv(net) - 2) * sqrt(nv(net) / 2)
_tmaxdev(f, net; kwargs...) =
    throw(ArgumentError("centralization: no theoretical maximum is known for $f; " *
                        "pass normalized=false for the raw deviation sum"))

# Edmonds–Karp max flow on a dense capacity matrix, optionally with one
# vertex excluded (0 = none). O(V·E²); fine at research scale.
function _maxflow(cap::Matrix{Float64}, n::Integer, s::Integer, t::Integer, excluded::Integer)
    n, s, t, excluded = Int(n), Int(s), Int(t), Int(excluded)
    flow = zeros(n, n)
    total = 0.0

    parent = zeros(Int, n)
    while true
        # BFS for an augmenting path in the residual graph
        fill!(parent, 0)
        parent[s] = s
        queue = [s]
        head = 1
        while head <= length(queue) && parent[t] == 0
            u = queue[head]
            head += 1
            for w in 1:n
                if parent[w] == 0 && w != excluded && cap[u, w] - flow[u, w] > 1e-12
                    parent[w] = u
                    push!(queue, w)
                end
            end
        end
        parent[t] == 0 && break

        # Bottleneck capacity along the path
        aug = Inf
        w = t
        while w != s
            u = parent[w]
            aug = min(aug, cap[u, w] - flow[u, w])
            w = u
        end

        # Augment
        w = t
        while w != s
            u = parent[w]
            flow[u, w] += aug
            flow[w, u] -= aug
            w = u
        end
        total += aug
    end

    return total
end
