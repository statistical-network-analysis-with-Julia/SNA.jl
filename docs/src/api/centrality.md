# Centrality API Reference

This page documents all centrality functions available in SNA.jl.

## Degree-Based Centrality

Centrality measures based on the number of direct connections.

```@docs
degree_centrality
```

## Path-Based Centrality

Centrality measures based on shortest paths and network flow.

```@docs
betweenness_centrality
closeness_centrality
flowbet
```

## Spectral Centrality

Centrality measures based on eigenvalues and iterative weighting of the adjacency matrix.

```@docs
eigenvector_centrality
bonacich_power
katz_centrality
pagerank
```

## Graph Centralization

Freeman centralization of a vertex centrality measure, following R
`sna::centralization`.

```@docs
centralization
```
