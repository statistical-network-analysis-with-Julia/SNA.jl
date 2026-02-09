# SNA.jl

Social Network Analysis tools for Julia.

## Overview

SNA.jl provides a comprehensive suite of descriptive analysis tools for social networks, including centrality measures, cohesion analysis, structural equivalence, and network-level statistics.

This package is a Julia port of the R `sna` package from the StatNet collection.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/Statistical-network-analysis-with-Julia/SNA.jl")
```

## Features

- **Centrality**: Degree, betweenness, closeness, eigenvector, Bonacich power, PageRank
- **Cohesion**: Components, cliques, k-cores, cutpoints, bridges
- **Equivalence**: Structural equivalence, regular equivalence, blockmodels
- **Measures**: Density, reciprocity, transitivity, dyad/triad census

## Quick Start

```julia
using Network
using SNA

# Create a network
net = Network{Int}(; n=10, directed=false)
for (i, j) in [(1,2), (2,3), (3,4), (4,5), (1,3), (2,4), (3,5)]
    add_edge!(net, i, j)
end

# Centrality measures
deg = degree_centrality(net)
bet = betweenness_centrality(net)
clo = closeness_centrality(net)

# Network-level measures
d = density(net)
t = transitivity(net)
```

## Centrality Measures

```julia
# Degree centrality
degree_centrality(net)
degree_centrality(net; mode=:in)   # In-degree (directed)
degree_centrality(net; mode=:out)  # Out-degree (directed)

# Betweenness centrality
betweenness_centrality(net; normalized=true)

# Closeness centrality
closeness_centrality(net)

# Eigenvector centrality
eigenvector_centrality(net)

# Bonacich power centrality
bonacich_power(net; beta=0.5)

# Katz centrality
katz_centrality(net; alpha=0.1)

# PageRank
pagerank(net; damping=0.85)
```

## Network-Level Measures

```julia
# Density
density(net)

# Reciprocity (directed networks)
reciprocity(net)

# Transitivity (clustering coefficient)
transitivity(net)

# Dyad census
dc = dyad_census(net)  # Returns (mutual, asymmetric, null)

# Triad census (16 isomorphism classes)
tc = triad_census(net)
```

## Cohesion Analysis

```julia
# Connected components
comps = components(net)

# Largest component
lc = largest_component(net)

# K-cores
cores = kcores(net)

# Cliques
cliques(net; min_size=3)

# Cutpoints (articulation points)
cutpoints(net)

# Bridges
bridges(net)
```

## Structural Equivalence

```julia
# Structural equivalence matrix
se = structural_equivalence(net)

# Regular equivalence
re = regular_equivalence(net)

# Blockmodel
bm = blockmodel(net, membership)
```

## Path Analysis

```julia
# Geodesic distances
dists = geodesic_distance(net)

# Diameter
d = diameter(net)

# Average path length
apl = average_path_length(net)

# Reachability matrix
reach = reachability(net)
```

## Random Graphs

```julia
# Erdos-Renyi G(n,m)
net = rgnm(n=100, m=200)

# Erdos-Renyi G(n,p)
net = rgnp(n=100, p=0.1)

# General random graph
net = rgraph(n=100, tprob=0.1)
```

## License

MIT License
