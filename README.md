# SNA.jl

[![Network Analysis](https://img.shields.io/badge/Network-Analysis-orange.svg)](https://github.com/statistical-network-analysis-with-Julia/SNA.jl)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://statistical-network-analysis-with-Julia.github.io/SNA.jl/stable/)
[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://statistical-network-analysis-with-Julia.github.io/SNA.jl/dev/)
[![Julia](https://img.shields.io/badge/Julia-1.12+-purple.svg)](https://julialang.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

<p align="center">
  <img src="docs/src/assets/logo.svg" alt="SNA.jl icon" width="160">
</p>

A Julia implementation of **Social Network Analysis** tools for descriptive analysis of social networks.

## Overview

SNA.jl provides a comprehensive suite of descriptive analysis tools for social networks, accounting for:

- **Centrality**: Vertex-level importance measures including degree, betweenness, closeness, eigenvector, and PageRank
- **Network measures**: Global statistics such as density, reciprocity, transitivity, and dyad/triad census
- **Cohesion**: Substructure detection including components, cliques, k-cores, cutpoints, and bridges
- **Equivalence**: Positional analysis including structural equivalence, regular equivalence, and blockmodeling

SNA.jl is a port of the R [`sna`](https://cran.r-project.org/package=sna) package from the StatNet collection, providing efficient tools for analysing both directed and undirected networks.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/statistical-network-analysis-with-Julia/SNA.jl")
```

## Functions Implemented

### 1. Centrality Measures

Vertex-level measures of structural importance.

#### Degree-Based Centrality

<!-- skip-check -->
```julia
degree_centrality(net; mode=:total, normalized=false)  # In-degree, out-degree, or total
```

#### Path-Based Centrality

<!-- skip-check -->
```julia
betweenness_centrality(net; normalized=true)   # Fraction of shortest paths through vertex
closeness_centrality(net; normalized=true)      # Inverse average distance to all others
flowbet(net)                                    # Flow betweenness (all paths, not just shortest)
```

#### Spectral Centrality

<!-- skip-check -->
```julia
eigenvector_centrality(net; max_iter=100, tol=1e-6)  # Centrality weighted by neighbor centrality
bonacich_power(net; β=0.5, normalized=true)           # Bonacich power (direct and indirect ties)
katz_centrality(net; α=0.1, β=1.0)                    # Katz centrality with damping factor
pagerank(net; α=0.85, max_iter=100, tol=1e-6)         # Google PageRank
```

### 2. Network-Level Measures

Global statistics characterising the network as a whole.

#### Density and Reciprocity

<!-- skip-check -->
```julia
density(net)                          # Proportion of possible edges present
reciprocity(net; method=:dyadic)      # Proportion of mutual ties (:dyadic or :edgewise)
mutuality(net)                        # Proportion of connected pairs that are symmetric
gden(net)                             # Alias for density (R sna compatibility)
grecip(net)                           # Alias for reciprocity (R sna compatibility)
```

#### Transitivity and Hierarchy

<!-- skip-check -->
```julia
transitivity(net; type=:global)       # Clustering coefficient (:global, :local, or :average)
hierarchy(net)                        # Krackhardt's hierarchy (asymmetric reachability)
efficiency(net)                       # 1 minus proportion of excess edges over a tree
connectedness(net)                    # Proportion of pairs where one can reach the other
gtrans(net)                           # Alias for transitivity (R sna compatibility)
```

#### Census Functions

<!-- skip-check -->
```julia
dyad_census(net)                      # Counts of mutual, asymmetric, and null dyads
triad_census(net)                     # 16-element vector of triad isomorphism classes
component_dist(net)                   # Distribution of component sizes
```

### 3. Cohesion Analysis

Functions for detecting substructures and assessing network vulnerability.

#### Components

<!-- skip-check -->
```julia
components(net; mode=:weak)           # Connected components (:weak or :strong)
largest_component(net; mode=:weak)    # Vertices in the largest component
bicomponents(net)                     # Biconnected components (no cutpoints)
```

#### Subgroup Detection

<!-- skip-check -->
```julia
cliques(net; min_size=3)              # Maximal cliques of at least min_size
kcores(net; k=1)                      # Vertices in the k-core
```

#### Vulnerability

<!-- skip-check -->
```julia
cutpoints(net)                        # Vertices whose removal disconnects the network
bridges(net)                          # Edges whose removal disconnects the network
```

### 4. Path Analysis

Functions for computing distances and reachability.

<!-- skip-check -->
```julia
geodesic_distance(net)                # Matrix of shortest path distances (Inf if unreachable)
diameter(net)                         # Longest shortest path
average_path_length(net)              # Mean shortest path over reachable pairs
reachability(net)                     # Boolean matrix of reachability between all pairs
```

### 5. Structural Equivalence

Functions for analysing positional similarity and building blockmodels.

#### Equivalence Measures

<!-- skip-check -->
```julia
structural_equivalence(net; method=:correlation)  # Similarity matrix (:correlation, :euclidean, :hamming)
regular_equivalence(net; max_iter=100, tol=1e-6)  # REGE algorithm for regular equivalence
```

#### Clustering and Blockmodeling

<!-- skip-check -->
```julia
equiv_clust(net; method=:structural, k=nothing)   # Cluster vertices by equivalence
blockmodel(net; k=2, method=:structural)           # Block densities from equivalence clustering
consensus(clusterings)                              # Consensus from multiple clustering solutions
```

### 6. Random Graph Generators

<!-- skip-check -->
```julia
rgraph(10; tprob=0.1)                # Bernoulli random digraph (sna-style)
rgnm(100, 200)                       # Erdős-Rényi G(n,m)
rgnp(100, 0.1)                       # Erdős-Rényi G(n,p)
```

### 7. Layout Algorithms

Functions for computing vertex positions for network visualisation.

<!-- skip-check -->
```julia
layout_fruchterman_reingold(net)      # Force-directed layout
layout_kamada_kawai(net)              # Energy-based layout
layout_circle(net)                    # Circular layout
layout_random(net)                    # Random layout
```

## Usage

### Basic Example

```julia
using Network
using SNA

# Create a network
net = network(5)
add_edge!(net, 1, 2)
add_edge!(net, 1, 3)
add_edge!(net, 2, 3)
add_edge!(net, 3, 4)
add_edge!(net, 4, 5)

# Centrality measures
deg = degree_centrality(net; mode=:out)
bet = betweenness_centrality(net)
clo = closeness_centrality(net)

# Network-level measures
d = density(net)
r = reciprocity(net)
t = transitivity(net)
```

### Dyad and Triad Census

```julia
# Complete directed graph on 3 vertices
net = network(3)
for (i, j) in [(1,2), (2,1), (1,3), (3,1), (2,3), (3,2)]
    add_edge!(net, i, j)
end

census = dyad_census(net)   # (mutual=3, asymmetric=0, null=0)
tc = triad_census(net)      # 16-element vector of triad counts
```

### Cohesion Analysis

```julia
# Network with two components
net = network(6)
add_edge!(net, 1, 2)
add_edge!(net, 2, 3)
add_edge!(net, 4, 5)
add_edge!(net, 5, 6)

comps = components(net; mode=:weak)     # 2 components
largest = largest_component(net)         # Vertices in the larger component

# Undirected network for k-cores and cutpoints
net = network(5; directed=false)
add_edge!(net, 1, 2)
add_edge!(net, 1, 3)
add_edge!(net, 2, 3)
add_edge!(net, 3, 4)
add_edge!(net, 4, 5)

core_2 = kcores(net; k=2)              # [1, 2, 3] (the triangle)
cuts = cutpoints(net)                   # [3, 4] (articulation points)
br = bridges(net)                       # Edges whose removal disconnects
```

### Structural Equivalence and Blockmodeling

```julia
net = network(4)
add_edge!(net, 1, 3)
add_edge!(net, 1, 4)
add_edge!(net, 2, 3)
add_edge!(net, 2, 4)

# Vertices 1 and 2 have identical tie patterns
se = structural_equivalence(net; method=:correlation)
se[1, 2]  # 1.0 (perfectly equivalent)

# Build a blockmodel
bm = blockmodel(net; k=2)
bm.membership      # Block assignments
bm.block_matrix     # Density of ties between blocks
```

### Path Analysis

```julia
net = network(4)
add_edge!(net, 1, 2)
add_edge!(net, 2, 3)
add_edge!(net, 3, 4)

dist = geodesic_distance(net)
dist[1, 4]          # 3.0 (three hops)
dist[4, 1]          # Inf (unreachable in directed network)

diameter(net)               # 3.0
average_path_length(net)    # Mean over reachable pairs
```

## Running Tests

<!-- skip-check -->
```julia
include("test/runtests.jl")
```

## Documentation

For more detailed documentation, see:

- [Stable Documentation](https://statistical-network-analysis-with-Julia.github.io/SNA.jl/stable/)
- [Development Documentation](https://statistical-network-analysis-with-Julia.github.io/SNA.jl/dev/)

## References

1. Wasserman, S., Faust, K. (1994). *Social Network Analysis: Methods and Applications*. Cambridge University Press.

2. Butts, C.T. (2008). Social network analysis with sna. *Journal of Statistical Software*, 24(6), 1-51.

3. Butts, C.T. (2020). sna: Tools for Social Network Analysis. R package. [https://cran.r-project.org/package=sna](https://cran.r-project.org/package=sna)

## License

MIT License - see [LICENSE](LICENSE) for details.
