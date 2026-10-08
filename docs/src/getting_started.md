# Getting Started

Build a small observed network, compare centrality and cohesion, then test the
association between two bundled relations. All examples are deterministic or
use an explicit random seed.

## Prepare the environment

```@raw html
<p>Use Julia <strong>1.12+</strong> and the <a href="/getting-started/">workspace installation guide</a> for this <strong>unreleased 0.2.0 development version</strong>. From the prepared workspace:</p>
```

```bash
julia --project=SNA.jl
```

```@raw html
<p><a href="/NetworkCore.jl/dev/">NetworkCore.jl</a> supplies network objects, attributes and datasets. The remaining examples perform analysis in the already prepared environment.</p>
```

## Build a network with two groups

```julia
using NetworkCore, SNA

net = network(6; directed=false)
for (i, j) in [(1, 2), (1, 3), (2, 3), (3, 4), (4, 5), (4, 6), (5, 6)]
    add_edge!(net, i, j)
end
(nv(net), ne(net))                  # (6, 7)
```

Actors 1–3 and 4–6 each form a triangle. The tie between 3 and 4 connects the two
groups. `network` defaults to a directed network, so specify `directed=false`
when each observed relation is mutual by definition.

## Compare actor positions

```julia
deg = degreecent(net)
btw = betweenness(net)
[(actor=i, degree=deg[i], betweenness=btw[i]) for i in 1:nv(net)]
```

Degree counts direct ties. Betweenness counts an actor's participation in shortest
paths between other actors; actors 3 and 4 have the largest values here. These are
structural descriptions. Calling a high-scoring actor influential requires further
substantive evidence.

SNA.jl's measures carry R `sna`'s names (`degreecent` is sna's `degree`,
`betweenness`, `closeness`, `evcent`, `bonpow`, …) and follow sna's defaults,
including raw betweenness counts. Graphs.jl's `degree_centrality`,
`betweenness_centrality`, … are different functions with Graphs.jl's
conventions. On directed data choose `cmode=:indegree`, `:outdegree` or
`:freeman` for degree; closeness and eigenvector centrality answer different
questions and have their own direction and connectivity conventions. See
[centrality](guide/centrality.md) and the [R concordance](r_concordance.md).

## Describe groups and connections

```julia
gden(net)                          # 7 / 15
gtrans(net)                        # 0.6
component_dist(net).csize          # [6]: one component
sort(cutpoints(net))               # [3, 4]
cliques(net; min_size=3)           # the two triangles
kcores(net)                        # [2, 2, 2, 2, 2, 2]
```

Density is the proportion of possible ties present. The two triangles give
local closure, while the bridge is the only connection between groups. Removing
that tie, or either cutpoint, disconnects the network. Maximal cliques identify
fully connected groups; k-cores (`kcores` returns each actor's core number)
use a minimum-degree criterion instead.

For directed data, component and cutpoint defaults use strong connectivity.
Weak connectivity ignores edge direction. Cliques and bicomponents default to
mutual-arc symmetrization. Choose these conventions explicitly for the question
at hand; see [cohesion](guide/cohesion.md).

## Compare tie patterns

```julia
distances = sedist(net)            # Hamming distances, as R sna::sedist
blocks = blockmodel(net; k=2)
blocks.membership
blocks.block_matrix
```

Structural equivalence compares ties to the same other actors (the ties
between the two actors compared are left out, as in sna). Blockmodeling
clusters these profiles and summarizes between-block densities; the requested
number of blocks is an analyst choice. Regular equivalence asks about similar
roles and uses a different criterion. Neither is a guarantee of a substantive
community structure. See [equivalence](guide/equivalence.md).

## Test association between two relations

The bundled Florentine marriage and business networks use the same actor order.
QAP compares their observed association with associations obtained by permuting
actor labels in one network. This preserves its internal relational structure.

```julia
using Random

marriage = load_dataset(:florentine_marriage)
business = load_dataset(:florentine_business)
gcor(marriage, business)           # 0.372, as R sna::gcor
qap = qaptest(gcor, marriage, business; n_sim=999, rng=Xoshiro(2026))
qap
```

`gcor` is the graph correlation of R `sna::gcor`, over off-diagonal dyads.
The result reports the observed statistic and empirical permutation tails.
Permutation inference requires an appropriate exchangeability assumption; a seed
provides reproducibility, not that assumption. More replicates improve Monte Carlo
precision. Use `netlm` or `netlogit` for network regression with multiple predictors;
read their [QAP and uncertainty contracts](api/measures.md#QAP-Inference-and-Network-Regression)
before interpreting p-values.

## Analysis checks

| Check | Why it matters |
|---|---|
| Actor alignment | Comparing two networks requires the same actors in the same order; equal dimensions alone are insufficient. |
| Missing ties | Masked dyads throw by default. `missing=:face` uses stored values explicitly and does not estimate the missing ties. |
| Values and paths | Binary ties are the default. Selected routines accept `ignore_eval=false, attr=:weight`; weighted shortest-path centrality is not implemented. |
| Two-mode data | Measures use both actor modes. QAP and network regression use only the cross-mode dyads, and permutations preserve modes. |
| Uncertainty | Regression covariance and Wald intervals assume independent dyads. QAP p-values are a separate permutation result. |
| Logistic identification | A separated observed design has no finite estimate: `netlogit` warns, reports `converged == false`, names the separated coefficients and withholds its statistics, p-values and intervals. The default likelihood-ratio QAP statistic exists on separated permutations; with sna's Wald statistic (`statistic=:wald`) the p-values a separated permutation touches are withheld rather than computed from the remaining replicates. |

For layouts, graph generators and the full missing-data policy, see
[utilities](api/utilities.md). For global network indices and censuses, see
[network measures](guide/measures.md).

```@raw html
<p>SNA describes observed structure; <a href="/ERGM.jl/dev/">ERGM.jl</a> fits probability models for networks and <a href="/TSNA.jl/dev/">TSNA.jl</a> handles temporal network measures.</p>
```
