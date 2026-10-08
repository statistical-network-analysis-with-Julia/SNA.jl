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

net = Network(3)
add_edge!(net, 1, 2)
degreecent(net)
set_missing_dyad!(net, 1, 2)
try
    degreecent(net)
catch err
    @assert err isa ArgumentError
end
degreecent(net; missing=:face)  # explicit opt-in to stored values
```

This is the ecosystem missing-data contract. It is *defined in NetworkCore.jl* and
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

The full contract is documented in the NetworkCore.jl manual:
[Ecosystem Contracts](https://Statistical-network-analysis-with-Julia.github.io/NetworkCore.jl/dev/api/contracts/).

## Result Metadata

`netlm` and `netlogit` implement the ecosystem's shared result-metadata
protocol, so `NetworkCore.fit_metadata(result)` reports what the fit actually did.
The per-result-type methods are documented on the
[Measures page](measures.md#Result-Metadata); the protocol itself is in the
NetworkCore.jl manual under
[Result Metadata](https://Statistical-network-analysis-with-Julia.github.io/NetworkCore.jl/dev/api/metadata/).

## Layouts

Vertex-coordinate algorithms for plotting, following R's `sna::gplot.layout.*`.
Each returns an `n × 2` coordinate matrix.

```@docs
layout_fruchterman_reingold
layout_kamada_kawai
layout_circle
layout_random
```

## Random Graphs

Random-graph generators following R's `sna::rgraph` and `sna::rgnm`, and
the Erdős–Rényi `rgnp`. All draw from an `rng` keyword.

```@docs
rgraph
rgnm
rgnp
```

## Renamed in 0.2.0

These names were used during development and were never part of a release;
they have no aliases. Use the new names. The full table, with the Graphs.jl
functions SNA no longer extends, is in the
[R concordance](../r_concordance.md) and the README.

| Old name | Use instead |
|---|---|
| `bonacich_power(net)` | [`bonpow`](@ref)`(net)` |
| `geodesic_distance(net)` | [`geodist`](@ref)`(net).gdist` |
| `reciprocity(net; method=m)` | [`grecip`](@ref)`(net; measure=m)` |
| `structural_equivalence(net; method=m)` | [`sedist`](@ref)`(net; method=m)` (default Hamming, a count) |
| `consensus(clusterings)` | [`consensus_clustering`](@ref)`(clusterings)` |
| `SNA.components(net)` | [`component_dist`](@ref)`(net).membership` |
| `equiv_clust(net; method=:structural)` | `equiv_clust(net; method=:correlation, cluster_method=:average)` |
| `equiv_clust(net; method=:regular)` | `equiv_clust(net; equiv_fun=regular_equivalence, cluster_method=:average)` |
| `reps=` (`qaptest`, `netlm`, `netlogit`) | `n_sim=` |
| `max_iter=` (`regular_equivalence`) | `maxiter=` |
