"""
    brokerage(net, cl; missing=:error) -> NamedTuple

Gould–Fernandez brokerage, R `sna::brokerage`. Actor `j` brokers the ordered
pair `(i, k)` when `i → j → k` and there is no tie `i → k`. Given the class
of every actor, `cl` (a vector with one entry per vertex, or the name of a
vertex attribute), each brokered pair counts towards one role of `j`:

| Role | sna label | Classes of `i → j → k` |
|---|---|---|
| coordinator | `w_I` | A → A → A |
| itinerant broker | `w_O` | A → B → A |
| representative | `b_IO` | A → A → B |
| gatekeeper | `b_OI` | A → B → B |
| liaison | `b_O` | A → B → C |
| total | `t` | any |

The result has sna's components, with the six roles in the order above
(`roles`):

- `raw_nli`: the `n × 6` matrix of each actor's role counts;
- `exp_nli`, `sd_nli`, `z_nli`: their expectation and standard deviation,
  and the z-scores `(raw − exp) / sd`, under Gould and Fernandez's null
  model, in which ties are independent with the observed density and the
  class sizes are fixed;
- `raw_gli`, `exp_gli`, `sd_gli`, `z_gli`: the same for the network totals;
- `exp_grp`, `sd_grp`: expectation and standard deviation per class (rows
  in the order of `clid`);
- `clid`, `n`: the classes in order of first appearance and their sizes.

A z-score whose standard deviation is 0 is `NaN` or infinite, as in R. The
z-scores refer to a normal approximation that is rough for small groups; a
conditional uniform graph test ([`cug_test`](@ref)) gives a simulation-based
alternative. On an undirected network every tie is read in both directions,
as sna does. Self-loops are ignored (sna counts them in the density).

`missing=:error` refuses masked dyads; `missing=:face` uses their stored
face values.

# Example
```julia
using SNA
net = network(5)
add_edges!(net, [(1, 2), (2, 3), (3, 4), (4, 5)])
b = brokerage(net, [:a, :a, :b, :b, :c])
b.roles                              # ["w_I", "w_O", "b_IO", "b_OI", "b_O", "t"]
b.raw_nli[2, :]                      # actor 2 is a representative: 1 → 2 → 3
b.z_gli                              # z-scores of the network totals
```

# References
Gould, R.V., Fernandez, R.M. (1989). Structures of mediation: a formal
approach to brokerage in transaction networks. *Sociological Methodology*,
19, 89–126.
"""
function brokerage(net::AbstractNetwork, cl; missing::Symbol=:error)
    policy = missing  # local alias; `missing` here is the kwarg, not `Base.missing`
    require_observed(net, policy; context="brokerage")
    N = nv(net)
    classes = cl isa Symbol ? [get_vertex_attribute(net, cl, v) for v in 1:N] : collect(cl)
    length(classes) == N ||
        throw(ArgumentError("brokerage: cl must have one class per vertex ($N), got $(length(classes))"))
    clid = unique(classes)
    icl = [findfirst(isequal(c), clid) for c in classes]
    A = _sociomatrix(net) .!= 0                    # loops removed

    br = zeros(N, 6)
    for j in 1:N, i in 1:N
        (i == j || !A[i, j]) && continue
        for k in 1:N
            (k == i || k == j || !A[j, k] || A[i, k]) && continue
            role = icl[i] == icl[j] ? (icl[k] == icl[j] ? 1 : 3) :
                   icl[i] == icl[k] ? 2 :
                   icl[j] == icl[k] ? 4 : 5
            br[j, role] += 1
            br[j, 6] += 1
        end
    end
    gbr = vec(sum(br; dims=1))

    # Expectations and variances: sna's formulas (Gould & Fernandez 1989)
    d = count(A) / (N * (N - 1))
    n = Float64[count(==(c), icl) for c in eachindex(clid)]
    G = length(clid)
    p2 = d^2 * (1 - d)
    p3 = d^3 * (1 - d)^3
    ch2(x) = x * (x - 1) / 2
    ebr = zeros(G, 6)
    vbr = zeros(G, 6)
    for i in 1:G
        others = [n[g] for g in 1:G if g != i]
        ebr[i, 1] = p2 * (n[i] - 1) * (n[i] - 2)
        vbr[i, 1] = ebr[i, 1] * (1 - p2) + 2 * (n[i] - 1) * (n[i] - 2) * (n[i] - 3) * p3
        ebr[i, 2] = p2 * sum(o * (o - 1) for o in others; init=0.0)
        vbr[i, 2] = ebr[i, 2] * (1 - p2) +
                    2 * sum(o * (o - 1) * (o - 2) for o in others; init=0.0) * p3
        ebr[i, 3] = p2 * (N - n[i]) * (n[i] - 1)
        vbr[i, 3] = ebr[i, 3] * (1 - p2) +
                    2 * ((n[i] - 1) * ch2(N - n[i]) + (N - n[i]) * ch2(n[i] - 1)) * p3
        ebr[i, 4] = ebr[i, 3]
        vbr[i, 4] = vbr[i, 3]
        ebr[i, 5] = p2 * (sum(others; init=0.0)^2 - sum(abs2, others; init=0.0))
        vbr[i, 5] = ebr[i, 5] * (1 - p2) +
                    4 * sum(o * ch2(N - o - n[i]) * p3 for o in others; init=0.0)
        ebr[i, 6] = p2 * (N - 1) * (N - 2)
        vbr[i, 6] = ebr[i, 6] * (1 - p2) + 2 * (N - 1) * (N - 2) * (N - 3) * p3
    end
    exp_nli = ebr[icl, :]
    sd_nli = sqrt.(vbr[icl, :])
    z_nli = (br .- exp_nli) ./ sd_nli

    d4 = d^4 * (1 - d)^2
    d5 = d^5 * (1 - d)
    egbr = zeros(6)
    vgbr = zeros(6)
    egbr[1] = p2 * sum(@. n * (n - 1) * (n - 2))
    vgbr[1] = egbr[1] * (1 - p2) +
              sum(@. n * (n - 1) * (n - 2) * ((4n - 10) * p3 - 4 * (n - 3) * d4 + (n - 3) * d5))
    egbr[2] = p2 * sum(@. n * (N - n) * (n - 1))
    vgbr[2] = egbr[2] * (1 - p2) +
              sum(x * y * (x - 1) * ((2x + 2y - 6) * p3 + (N - x - 1) * d5) for x in n, y in n) -
              sum(@. n * n * (n - 1) * ((4n - 6) * p3 + (N - n - 1) * d5))
    egbr[3] = egbr[2]
    vgbr[3] = egbr[3] * (1 - p2) +
              sum(@. n * (N - n) * (n - 1) * ((N - 3) * p3 + (n - 2) * d5))
    egbr[4] = egbr[3]
    vgbr[4] = vgbr[3]
    egbr[5] = p2 * (sum(x * y * (N - x - y) for x in n, y in n) - sum(@. n * n * (N - 2n)))
    vgbr[5] = egbr[5] * (1 - p2)
    for i in 1:G, j in 1:G, k in 1:G
        (i != j && j != k && i != k) || continue
        vgbr[5] += n[i] * n[j] * n[k] *
                   ((4 * (N - n[j]) - 2 * (n[i] + n[k] + 1)) * p3 -
                    (4 * (N - n[k]) - 2 * (n[i] + n[j] + 1)) * d4 +
                    (N - (n[i] + n[k] + 1)) * d5)
    end
    egbr[6] = p2 * N * (N - 1) * (N - 2)
    vgbr[6] = egbr[6] * (1 - p2) +
              N * (N - 1) * (N - 2) * ((4N - 10) * p3 - 4 * (N - 3) * d4 + (N - 3) * d5)
    sd_gli = sqrt.(vgbr)

    return (roles = ["w_I", "w_O", "b_IO", "b_OI", "b_O", "t"],
            raw_nli = br, exp_nli = exp_nli, sd_nli = sd_nli, z_nli = z_nli,
            raw_gli = gbr, exp_gli = egbr, sd_gli = sd_gli, z_gli = (gbr .- egbr) ./ sd_gli,
            exp_grp = ebr, sd_grp = sqrt.(vbr), clid = clid, n = Int.(n))
end
