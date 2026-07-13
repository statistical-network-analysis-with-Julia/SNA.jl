# Measures API Reference

This page documents network-level measures, path functions, and census functions available in SNA.jl.

## Network-Level Measures

### Density and Reciprocity

```@docs
density
gden
reciprocity
grecip
mutuality
```

### Transitivity and Hierarchy

```@docs
transitivity
gtrans
hierarchy
efficiency
connectedness
```

## Census Functions

Functions for classifying dyads, triads, and components.

```@docs
dyad_census
triad_census
component_dist
```

## Path Functions

Functions for computing distances, reachability, and path-based summaries.

```@docs
geodesic_distance
reachability
diameter
average_path_length
```

## QAP Inference and Network Regression

Permutation tests and network regression following R `sna::qaptest`,
`sna::netlm`, and `sna::netlogit`.

```@docs
qaptest
netlm
netlogit
```

### Result Types

```@docs
QAPTestResult
NetLMResult
NetLogitResult
```

### Result Metadata

SNA.jl implements the ecosystem's
[result-metadata protocol](https://Statistical-network-analysis-with-Julia.github.io/Networks.jl/dev/api/metadata/)
for its regression results, so what a fit did is inspectable via
`Networks.fit_metadata(result)`.

```@docs
objective(::NetLMResult)
objective(::NetLogitResult)
is_exact(::NetLMResult)
is_exact(::NetLogitResult)
se_method(::NetLMResult)
se_method(::NetLogitResult)
```

