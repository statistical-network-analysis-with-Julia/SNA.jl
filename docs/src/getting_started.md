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
<p><a href="/Networks.jl/dev/">Networks.jl</a> supplies network objects, attributes and datasets. The remaining examples perform analysis in the already prepared environment.</p>
```

## Build a network with two groups

```julia
using Networks, SNA

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
deg = degree_centrality(net)
btw = betweenness_centrality(net)
[(actor=i, degree=deg[i], betweenness=btw[i]) for i in 1:nv(net)]
```

Degree counts direct ties. Betweenness counts an actor's participation in shortest
paths between other actors; actors 3 and 4 have the largest values here. These are
structural descriptions. Calling a high-scoring actor influential requires further
substantive evidence.

On directed data, choose `mode=:in`, `:out` or `:total` for degree. Closeness,
eigenvector centrality and PageRank answer different questions and have different
direction and connectivity conventions. See [centrality](guide/centrality.md).

## Describe groups and connections

```julia
density(net)
transitivity(net)
components(net)
cutpoints(net)                     # [3, 4]
bridges(net)                       # [(3, 4)]
cliques(net; min_size=3)
```

Density is the proportion of possible ties present. The two triangles give
local closure, while the bridge is the only connection between groups. Removing
that edge, or either cutpoint, disconnects the network. Maximal cliques identify
fully connected groups; k-cores use a minimum-degree criterion instead.

For directed data, component and cutpoint defaults use strong connectivity.
Weak connectivity ignores edge direction. Cliques and bicomponents default to
mutual-arc symmetrization. Choose these conventions explicitly for the question
at hand; see [cohesion](guide/cohesion.md).

## Compare tie patterns

```julia
similarity = structural_equivalence(net; method=:correlation)
blocks = blockmodel(net; k=2)
blocks.membership
blocks.block_matrix
```

Structural equivalence compares ties to the same other actors. Blockmodeling
clusters these profiles and summarizes between-block densities; the requested
number of blocks is an analyst choice. Regular equivalence asks about similar
roles and uses a different criterion. Neither is a guarantee of a substantive
community structure. See [equivalence](guide/equivalence.md).

## Test association between two relations

The bundled Florentine marriage and business networks use the same actor order.
QAP compares their observed association with associations obtained by permuting
actor labels in one network. This preserves its internal relational structure.

```julia
using Random, Statistics

marriage = load_dataset(:florentine_marriage)
business = load_dataset(:florentine_business)
function dyad_correlation(A, B)
    a = [A[i, j] for i in axes(A, 1) for j in axes(A, 2) if i != j]
    b = [B[i, j] for i in axes(B, 1) for j in axes(B, 2) if i != j]
    cor(a, b)
end
qap = qaptest(dyad_correlation, marriage, business;
              reps=999, rng=Xoshiro(2026))
qap
```

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
| Two-mode data | Measures use both actor modes; QAP permutations preserve modes. The interpretation is not an automatic one-mode projection. |
| Uncertainty | Regression covariance and Wald intervals assume independent dyads. QAP p-values are a separate permutation result. |
| Logistic identification | Proved separating directions are rejected, with a nonexhaustive search. A nonidentified QAP replicate fails the requested analysis rather than being discarded. |

For layouts, graph generators and the full missing-data policy, see
[utilities](api/utilities.md). For global network indices and censuses, see
[network measures](guide/measures.md).

```@raw html
<p>SNA describes observed structure; <a href="/ERGM.jl/dev/">ERGM.jl</a> fits probability models for networks and <a href="/TSNA.jl/dev/">TSNA.jl</a> handles temporal network measures.</p>
```
