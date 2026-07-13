# Utilities

Layout algorithms, random-graph generators, and the missing-data contract that
every SNA.jl measure honours.

## Missing Data

Every descriptive measure in SNA.jl takes a `missing=` keyword. A masked dyad is
**unobserved**, not absent, so by default (`missing=:error`) a measure refuses to
run on a network with masked dyads rather than silently computing a number from
the stored face value of a tie nobody observed. Pass `missing=:face` to opt in to
face values explicitly.

```julia
using SNA

degree_centrality(net)                 # throws if `net` has masked dyads
degree_centrality(net; missing=:face)  # explicit, auditable opt-in
```

This is the ecosystem missing-data contract. It is *defined in Networks.jl* and
re-exported here, so `require_observed`, `supports_missing` and
`MISSING_POLICIES` are all callable unqualified after `using SNA`.
`require_observed` is the guard the measures call; `supports_missing` is the
trait a routine opts into.

No SNA measure implements a principled missing-data estimator — there is no
listwise deletion, no density-over-observed-dyads, no imputation — so none
declares `supports_missing`. The guard plus a documented `:face` treatment *is*
the contract here. Layouts and random-graph generators are exempt: a layout
draws a picture rather than reporting a statistic, and the generators create
fresh unmasked networks.

The full contract is documented in the Networks.jl manual:
[Ecosystem Contracts](https://Statistical-network-analysis-with-Julia.github.io/Networks.jl/dev/api/contracts/).

## Result Metadata

`netlm` and `netlogit` implement the ecosystem's shared result-metadata
protocol, so `Networks.fit_metadata(result)` reports what the fit actually did.
The per-result-type methods are documented on the
[Measures page](measures.md#Result-Metadata); the protocol itself is in the
Networks.jl manual under
[Result Metadata](https://Statistical-network-analysis-with-Julia.github.io/Networks.jl/dev/api/metadata/).

## Layouts

Vertex-coordinate algorithms for plotting, following R's `sna::gplot.layout.*`.

```@docs
layout_fruchterman_reingold
layout_kamada_kawai
layout_circle
layout_random
```

## Random Graphs

Random-graph generators following R's `sna::rgraph`, `sna::rgnm` and
`sna::rgnp`.

```@docs
rgraph
rgnm
rgnp
```
