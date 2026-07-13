# SNA.jl

*Social Network Analysis for Julia*

A Julia package providing descriptive analysis tools for social networks, including centrality measures, cohesion analysis, structural equivalence, and network-level indices.

## Overview

SNA.jl is a comprehensive toolkit for descriptive social network analysis. It provides functions for measuring the structural properties of networks at the vertex, dyad, triad, and whole-network levels. The package covers five core areas: centrality (who is important?), cohesion (how connected is the network?), network measures (what are the global properties?), structural equivalence (which actors occupy similar positions?), and QAP inference (are two relations associated?).

SNA.jl is a port of the R [sna](https://cran.r-project.org/package=sna) package from the [StatNet](https://statnet.org/) collection, adapted for Julia's multiple dispatch and type system. It builds on [Networks.jl](https://github.com/statistical-network-analysis-with-Julia/Networks.jl), which provides the core network data structure implementing the Graphs.jl interface.

### What is Social Network Analysis?

Social network analysis (SNA) is a set of methods for studying the structure of relationships among social actors. A social network consists of:

```text
Vertices (actors) connected by Edges (ties/relations)
```

Examples include:

- Friendship networks among individuals
- Trade relationships between countries
- Collaboration among scientists
- Communication patterns in organizations
- Citation networks among publications

### Key Concepts

| Concept | Description |
|---------|-------------|
| **Centrality** | Measures of vertex importance based on structural position |
| **Cohesion** | Properties related to connectivity and subgroup structure |
| **Density** | Proportion of possible ties that are present |
| **Reciprocity** | Extent to which ties are mutual in directed networks |
| **Transitivity** | Tendency for friends of friends to be friends |
| **Structural Equivalence** | Similarity of actors' tie patterns to all others |

### Applications

SNA.jl supports a wide range of analytical tasks:

- **Identifying key actors**: Centrality measures reveal influential, brokering, or prestigious vertices
- **Assessing cohesion**: Component analysis, k-cores, and clique detection reveal subgroup structure
- **Characterizing networks**: Global indices like density, reciprocity, and transitivity summarize overall structure
- **Position analysis**: Structural equivalence and blockmodeling uncover role structures
- **Network comparison**: Standardized measures allow comparison across networks

## Features

- **8 centrality measures**: Degree, betweenness, closeness, eigenvector, Bonacich power, Katz, PageRank, and flow betweenness — plus Freeman graph [`centralization`](@ref) of any of them
- **13 network-level indices**: Density, reciprocity, transitivity, mutuality, hierarchy, efficiency, connectedness, dyad census, triad census (edge-driven Batagelj–Mrvar algorithm), and more
- **7 cohesion functions**: Components, cliques, k-cores, cutpoints, bridges, bicomponents, and largest component
- **5 equivalence tools**: Structural equivalence, regular equivalence, blockmodeling, equivalence clustering, and consensus
- **4 path functions**: Geodesic distance, reachability, diameter, and average path length
- **4 layout algorithms**: Fruchterman-Reingold, Kamada-Kawai, circle, and random
- **3 random graph generators**: Bernoulli (rgraph/rgnp), fixed edge count (rgnm)
- **QAP inference**: [`qaptest`](@ref) for arbitrary graph-level statistics, and [`netlm`](@ref)/[`netlogit`](@ref) network regression with Dekker double-semi-partialing permutation tests

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/statistical-network-analysis-with-Julia/Networks.jl")
Pkg.add(url="https://github.com/statistical-network-analysis-with-Julia/SNA.jl")
```

Or for development:

```julia
using Pkg
Pkg.develop(path="/path/to/SNA.jl")
```

SNA.jl depends on [Networks.jl](https://github.com/statistical-network-analysis-with-Julia/Networks.jl), which will be installed automatically.

## Quick Start

```julia
using Networks, SNA

# Create a small directed network
net = network(5; directed=true)
add_edge!(net, 1, 2)
add_edge!(net, 1, 3)
add_edge!(net, 2, 3)
add_edge!(net, 3, 1)
add_edge!(net, 3, 4)
add_edge!(net, 4, 5)
add_edge!(net, 5, 4)

# Centrality
dc = degree_centrality(net; mode=:total)
bc = betweenness_centrality(net)
cc = closeness_centrality(net)

println("Degree centrality:      ", round.(dc, digits=3))
println("Betweenness centrality: ", round.(bc, digits=3))
println("Closeness centrality:   ", round.(cc, digits=3))

# Network-level measures
println("Density:      ", round(density(net), digits=3))
println("Reciprocity:  ", round(reciprocity(net), digits=3))
println("Transitivity: ", round(transitivity(net), digits=3))

# Cohesion
comps = components(net)
println("Components: ", length(comps))
println("Cutpoints:  ", cutpoints(net))

# Dyad census
dc_result = dyad_census(net)
println("Mutual: $(dc_result.mutual), Asymmetric: $(dc_result.asymmetric), Null: $(dc_result.null)")
```

## Choosing Measures

| Use Case | Recommended Functions |
|----------|----------------------|
| Who is most connected? | [`degree_centrality`](@ref), [`eigenvector_centrality`](@ref) |
| Who bridges groups? | [`betweenness_centrality`](@ref), [`cutpoints`](@ref) |
| How dense is the network? | [`density`](@ref), [`efficiency`](@ref) |
| Are ties reciprocated? | [`reciprocity`](@ref), [`dyad_census`](@ref) |
| Is there clustering? | [`transitivity`](@ref), [`triad_census`](@ref) |
| What subgroups exist? | [`components`](@ref), [`cliques`](@ref), [`kcores`](@ref) |
| Who occupies similar roles? | [`structural_equivalence`](@ref), [`blockmodel`](@ref) |
| What is the network diameter? | [`geodesic_distance`](@ref), [`diameter`](@ref) |
| How centralized is the network? | [`centralization`](@ref) |
| Are two relations associated? | [`qaptest`](@ref), [`netlm`](@ref), [`netlogit`](@ref) |

## Documentation

```@contents
Pages = [
    "getting_started.md",
    "guide/centrality.md",
    "guide/measures.md",
    "guide/cohesion.md",
    "guide/equivalence.md",
    "api/centrality.md",
    "api/measures.md",
    "api/cohesion.md",
    "api/equivalence.md",
]
Depth = 2
```

## Theoretical Background

### Centrality

Centrality measures quantify the importance of vertices in a network. The most common measures are based on degree (direct connections), paths (shortest routes through the network), and spectral properties (eigenvectors of the adjacency matrix).

For a network with adjacency matrix $A$ and $n$ vertices:

- **Degree centrality**: $C_D(i) = \sum_j A_{ij}$
- **Closeness centrality**: $C_C(i) = \frac{n - 1}{\sum_j d(i, j)}$ where $d(i, j)$ is the geodesic distance
- **Betweenness centrality**: $C_B(i) = \sum_{s \neq i \neq t} \frac{\sigma_{st}(i)}{\sigma_{st}}$ where $\sigma_{st}$ is the number of shortest paths from $s$ to $t$ and $\sigma_{st}(i)$ is the number passing through $i$
- **Eigenvector centrality**: The leading eigenvector of $A$, satisfying $Ax = \lambda x$

### Network-Level Indices

Global network properties summarize the overall structure:

- **Density**: $\Delta = \frac{m}{n(n-1)}$ for directed networks, where $m$ is the number of edges
- **Reciprocity**: Proportion of mutual dyads among all non-null dyads
- **Transitivity**: $T = \frac{3 \times \text{triangles}}{\text{connected triples}}$

### Structural Equivalence

Two vertices $i$ and $j$ are structurally equivalent if they have identical ties to and from all other vertices. In practice, approximate equivalence is measured using correlation, Euclidean distance, or Hamming distance between the vertices' tie profiles.

## References

1. Freeman, L.C. (1978). Centrality in social networks: Conceptual clarification. *Social Networks*, 1(3), 215-239.

2. Wasserman, S., & Faust, K. (1994). *Social Network Analysis: Methods and Applications*. Cambridge University Press.

3. Bonacich, P. (1987). Power and centrality: A family of measures. *American Journal of Sociology*, 92(5), 1170-1182.

4. Burt, R.S. (1976). Positions in networks. *Social Forces*, 55(1), 93-122.

5. Krackhardt, D. (1994). Graph theoretical dimensions of informal organizations. In K.M. Carley & M.J. Prietula (Eds.), *Computational Organization Theory* (pp. 89-111). Lawrence Erlbaum.

6. Seidman, S.B. (1983). Network structure and minimum degree. *Social Networks*, 5(3), 269-287.

7. Page, L., Brin, S., Motwani, R., & Winograd, T. (1999). The PageRank citation ranking: Bringing order to the web. *Stanford InfoLab Technical Report*.

8. Butts, C.T. (2008). Social network analysis with sna. *Journal of Statistical Software*, 24(6), 1-51.

## Module

```@docs
SNA
```
