# SNA.jl

```@raw html
<p>Describe the structure of an observed network: who connects to whom, which actors lie between groups, how cohesive the network is, and which positions have similar tie patterns. SNA.jl provides centrality, cohesion, network measures, equivalence analysis and QAP inference on <a href="/Networks.jl/dev/">Networks.jl</a> objects.</p>
```

**Start here:** [Getting started](getting_started.md) ·
[Centrality](guide/centrality.md) · [Network measures](guide/measures.md) ·
[Cohesion](guide/cohesion.md) · [API](api/measures.md)

## Find the connection between two groups

```@raw html
<p>Use Julia <strong>1.12+</strong> and the <a href="/getting-started/">workspace installation guide</a> for the current <strong>0.2.0 development version, unreleased</strong>. The examples assume that environment is already prepared.</p>
```

```julia
using Networks, SNA

# Two triangles joined by the tie between actors 3 and 4.
net = network(6; directed=false)
for (i, j) in [(1, 2), (1, 3), (2, 3), (3, 4), (4, 5), (4, 6), (5, 6)]
    add_edge!(net, i, j)
end

betweenness_centrality(net)
cutpoints(net)                  # [3, 4]
bridges(net)                    # the tie joining the groups
density(net)
```

Actors 3 and 4 connect the groups: removing either disconnects the network.
Betweenness quantifies their position on shortest paths. Such a structural position
alone does not establish influence or a causal role.

## Choose a measure that matches the question

| Question | Start with |
|---|---|
| Who has connections, reach or brokerage? | [Degree, closeness, betweenness and spectral centrality](guide/centrality.md) |
| How dense, reciprocal or transitive is the network? | [Global measures and dyad/triad censuses](guide/measures.md) |
| Which groups or vulnerable connections exist? | [Components, cores, cliques, bridges and cutpoints](guide/cohesion.md) |
| Who has similar tie patterns or roles? | [Structural/regular equivalence and blockmodels](guide/equivalence.md) |
| Are two relations associated? | [QAP tests and network regression](api/measures.md#QAP-Inference-and-Network-Regression) |

Measures use binary ties by default. Selected routines support valued ties;
shortest-path centralities currently use binary paths. Directed and two-mode
conventions differ by measure, so read its API before comparing results.

Masked dyads are rejected by default. `missing=:face` explicitly uses their stored
values; it is not a missing-data estimator. Regression covariance assumes independent
dyads, while QAP p-values come from a separate permutation procedure. Selected R
`sna` fixtures pin supported conventions; equivalence routines are not a complete
`sedist` port. See the [tutorial's analysis checks](getting_started.md#Analysis-checks).

```@raw html
<p>For statistical models of a single network, see <a href="/ERGM.jl/dev/">ERGM.jl</a>. For temporal measures, see <a href="/TSNA.jl/dev/">TSNA.jl</a>; for network panels, <a href="/Siena.jl/dev/">Siena.jl</a>.</p>
```

Package citation: [CITATION.bib](https://github.com/statistical-network-analysis-with-Julia/SNA.jl/blob/main/CITATION.bib).

## Module

```@docs
SNA
```
