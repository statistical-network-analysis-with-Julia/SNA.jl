# Centrality

Centrality indices describe the prominence of each actor. SNA.jl provides R
`sna`'s indices under sna's names, with sna's defaults: `degreecent` (sna's
`degree`), `betweenness`, `closeness`, `evcent`, `bonpow`, `infocent` and
`flowbet`, plus Freeman graph `centralization`. Each returns a
`Vector{Float64}` with one score per vertex, indexed by vertex ID.

!!! note "Graphs.jl's centrality functions are different functions"
    Graphs.jl's `degree_centrality`, `betweenness_centrality`,
    `closeness_centrality`, `eigenvector_centrality`, `katz_centrality` and
    `pagerank` also work on a `Network`, with Graphs.jl's conventions (for
    example normalized betweenness, and closeness scaled to the reachable
    component). SNA.jl does not change them. Before 0.2.0 SNA.jl added
    methods to those functions; its measures now have their own names.

## Example network

```julia
using NetworkCore, SNA

net = network(7; directed=true)
for (i, j) in [(1, 2), (2, 1), (1, 3), (2, 4), (3, 4), (3, 5), (4, 5),
               (4, 6), (5, 6), (5, 7), (6, 7), (7, 4)]
    add_edge!(net, i, j)
end
```

Actors 1 and 2 are tied both ways and both send to the rest of the network;
actors 4 to 7 form a cycle-rich core that 3 feeds into.

## Degree

```julia
degreecent(net)                          # Freeman degree: in + out
degreecent(net; cmode=:indegree)         # popularity
degreecent(net; cmode=:outdegree)        # activity
degreecent(net; normalized=true)         # ÷ the largest possible degree, 2(n−1)
degreecent(net; rescale=true)            # scores sum to 1, as sna's rescale=TRUE
```

| Keyword | Default | Meaning |
|---|---|---|
| `cmode` | `:freeman` | `:freeman`, `:indegree` or `:outdegree` (sna's `cmode`); ignored on undirected networks |
| `rescale` | `false` | divide by the sum of the scores (sna) |
| `normalized` | `false` | divide by the largest possible degree (not an sna option) |
| `diag` | `false` | count self-loops (once) |
| `ignore_eval`, `attr` | `true`, `:weight` | `ignore_eval=false` sums the edge attribute `attr` |

On an undirected network each edge adds 1 to the degree of both endpoints.

## Betweenness

Betweenness counts how many geodesics between other pairs pass through an
actor; high-betweenness actors are brokers.

```julia
betweenness(net)                         # raw scores (the default, as in sna)
betweenness(net; normalized=true)        # ÷ (n−1)(n−2): scores in [0, 1]
betweenness(net; cmode=:undirected)      # on the symmetrized network
```

$$C_B(i) = \sum_{s \neq i \neq t} \frac{\sigma_{st}(i)}{\sigma_{st}}$$

where $\sigma_{st}$ is the number of geodesics from $s$ to $t$ and
$\sigma_{st}(i)$ the number through $i$. An undirected pair counts once.

| Keyword | Default | Meaning |
|---|---|---|
| `cmode` | `:directed` | `:directed` or `:undirected` (sna's other path weightings are not implemented) |
| `normalized` | `false` | divide by the number of pairs of other vertices, `(n−1)(n−2)`, halved for undirected networks |
| `rescale` | `false` | divide by the sum of the scores (sna) |

## Closeness

```julia
closeness(net)                           # (n−1) / Σ distances; 0 if someone is unreachable
closeness(net; cmode=:suminvdir)         # Σ 1/distance / (n−1): defined on disconnected networks
closeness(net; cmode=:gil_schmidt)       # Gil–Schmidt power index
```

The default is Freeman closeness, $C_C(i) = (n-1) / \sum_j d(i,j)$. As in
sna, a vertex that cannot reach every other vertex scores 0, so on a network
that is not strongly connected many or all scores are 0. For disconnected
networks use the harmonic form (`:suminvdir`, or `:suminvundir` with
direction ignored) or Gil–Schmidt's index. A single vertex has undefined
closeness (`NaN`).

## Eigenvector centrality

```julia
evcent(net)
```

`evcent` is the principal eigenvector of the adjacency matrix, $Ax =
\lambda_1 x$, oriented non-negative with unit length. It weights an actor by
the centrality of the actors it sends ties to. SNA.jl solves the eigenproblem
directly (sna's `use.eigen=TRUE`): sna's default power iteration does not
converge on bipartite or other periodic networks. A repeated leading
eigenvalue (for instance two identical components) makes the vector
non-unique and triggers a warning; a network with no directed cycle has no
leading eigenvector and every score is `NaN`, as in sna.

## Bonacich power

```julia
bonpow(net; exponent=0.2)                # power from ties to powerful actors
bonpow(net; exponent=-0.2)               # power from ties to dependent actors
bonpow(net; exponent=0.2, rescale=true)  # scores sum to 1
```

$$c(\beta) = \alpha (I - \beta A)^{-1} A \mathbf{1}, \qquad \sum_i c_i^2 = n$$

The attenuation `exponent` (β) must satisfy $|\beta| < 1/\lambda_1$ for the
walk series to converge. A positive β rewards ties to well-connected actors
(prestige), a negative β ties to poorly connected ones (bargaining power). If
`(I − βA)` is computationally singular (reciprocal condition number below
`tol`, the test R's `solve` applies) the scores are `NaN` and a warning
names the problem.

## Information centrality

```julia
infocent(net)
```

Stephenson and Zelen's information centrality (sna's `infocent`) uses every
path, weighted by its information. It is defined on undirected networks: a
directed network is symmetrized (`cmode=:weak`, a tie in either direction,
the default; or `:strong`, mutual ties only). Isolates score 0.

## Flow betweenness

```julia
flowbet(net)
```

Freeman's flow betweenness: for each pair, the maximum flow that cannot
pass when the actor is removed, summed over pairs. Unlike `betweenness` it
counts every path, not just geodesics. Ties are binary unless
`ignore_eval=false` (sna's default uses edge values).

## Comparing measures

```julia
scores = (degree = degreecent(net; normalized=true),
          betweenness = betweenness(net; normalized=true),
          closeness = closeness(net; cmode=:suminvdir),
          eigenvector = evcent(net))
for i in 1:nv(net)
    println(i, ": ", join((round(s[i]; digits=3) for s in scores), "  "))
end
```

| Question | Measure |
|---|---|
| Who has the most ties? | `degreecent` |
| Who brokers between others? | `betweenness`, `flowbet` |
| Who reaches others quickly? | `closeness` (`:suminvdir` on disconnected networks) |
| Who is tied to central actors? | `evcent`, `bonpow` with β > 0 |
| Who has power over dependents? | `bonpow` with β < 0 |

## Graph centralization

Freeman centralization summarizes how concentrated an index is: 0 when every
vertex scores the same, 1 for the most centralized network of that size (a
star for the classic indices). `centralization(net, f)` follows R's
`centralization(dat, FUN)`: keyword arguments go to `f`, and the theoretical
maximum follows them.

```julia
star = network(5; directed=false)
for v in 2:5
    add_edge!(star, 1, v)
end

centralization(star, degreecent)          # 1.0
centralization(star, betweenness)         # 1.0
centralization(star, closeness)           # 1.0
centralization(net, degreecent; cmode=:indegree)
centralization(net, degreecent; normalized=false)   # raw deviation sum
```

`f` can be `degreecent`, `betweenness`, `closeness`, `evcent` or `bonpow`;
any other function works with `normalized=false`. A network too small for
the maximum to be positive gives `NaN`, as in sna.

## Valued and two-mode inputs

Ties are binary by default. Degree, eigenvector, Bonacich, information and
flow centrality accept `ignore_eval=false, attr=:weight`; every edge must
then carry a finite value of that attribute, and spectral and flow measures
require non-negative values. Self-loops are ignored unless `diag=true`.
Shortest paths are always binary (weighted geodesics are not implemented).

Two-mode networks are analysed on the square adjacency matrix of all actors.
Dense spectral solves need O(n²) memory and O(n³) time; flow betweenness is
much more expensive and suits networks of moderate size.
