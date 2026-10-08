# Centrality API Reference

Vertex-level centrality indices with R `sna`'s names and semantics. None of
them is a method of a Graphs.jl function; see the
[R concordance](../r_concordance.md) for the mapping from sna and from the
pre-0.2 names.

## Degree

```@docs
degreecent
```

## Path-based centrality

```@docs
betweenness
closeness
flowbet
```

## Spectral and information centrality

```@docs
evcent
bonpow
infocent
```

## Graph centralization

```@docs
centralization
```
