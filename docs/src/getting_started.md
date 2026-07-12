# Getting Started

This tutorial walks through common use cases for SNA.jl, from basic centrality computation to advanced structural equivalence analysis.

## Installation

Install SNA.jl from GitHub:

```julia
using Pkg
Pkg.add(url="https://github.com/Statistical-network-analysis-with-Julia/SNA.jl")
```

SNA.jl depends on [Network.jl](https://github.com/Statistical-network-analysis-with-Julia/Network.jl), which provides the core network data structure. It will be installed automatically as a dependency.

## Basic Workflow

The typical SNA.jl workflow consists of four steps:

1. **Create or load a network** - Build the network from data
2. **Compute centrality measures** - Identify important actors
3. **Analyze network structure** - Assess density, reciprocity, cohesion
4. **Positional analysis** - Find structurally equivalent actors and build blockmodels

## Step 1: Create a Network

Networks are created using the `Network.jl` package:

```julia
using Network, SNA

# Create a directed network with 6 vertices
net = network(6; directed=true)

# Add edges representing a communication network
add_edge!(net, 1, 2)  # Alice → Bob
add_edge!(net, 1, 3)  # Alice → Carol
add_edge!(net, 2, 1)  # Bob → Alice (reciprocal)
add_edge!(net, 2, 4)  # Bob → David
add_edge!(net, 3, 4)  # Carol → David
add_edge!(net, 3, 5)  # Carol → Eve
add_edge!(net, 4, 5)  # David → Eve
add_edge!(net, 4, 6)  # David → Frank
add_edge!(net, 5, 6)  # Eve → Frank
add_edge!(net, 6, 4)  # Frank → David (reciprocal)

println("Vertices: ", nv(net))   # 6
println("Edges: ", ne(net))      # 10
println("Directed: ", is_directed(net))  # true
```

### Creating Undirected Networks

```julia
# Undirected friendship network
net_undir = network(5; directed=false)
add_edge!(net_undir, 1, 2)
add_edge!(net_undir, 1, 3)
add_edge!(net_undir, 2, 3)
add_edge!(net_undir, 3, 4)
add_edge!(net_undir, 4, 5)
```

### Random Networks

SNA.jl includes random graph generators:

```julia
# Bernoulli random graph (each edge exists with probability p)
net_random = rgnp(20, 0.15)

# Random graph with fixed edge count
net_fixed = rgnm(20, 40)

# Directed random graph with specified density
net_dense = rgraph(10; tprob=0.3)
```

## Step 2: Compute Centrality Measures

Centrality measures identify the most important actors in the network:

```julia
using Network, SNA

net = network(6; directed=true)
add_edge!(net, 1, 2); add_edge!(net, 1, 3)
add_edge!(net, 2, 1); add_edge!(net, 2, 4)
add_edge!(net, 3, 4); add_edge!(net, 3, 5)
add_edge!(net, 4, 5); add_edge!(net, 4, 6)
add_edge!(net, 5, 6); add_edge!(net, 6, 4)

# Degree centrality (number of connections)
dc = degree_centrality(net; mode=:total)
println("Degree centrality: ", round.(dc, digits=2))

# Betweenness centrality (brokerage)
bc = betweenness_centrality(net)
println("Betweenness: ", round.(bc, digits=3))

# Closeness centrality (proximity to all others)
cc = closeness_centrality(net)
println("Closeness: ", round.(cc, digits=3))

# Eigenvector centrality (connection to well-connected others)
ec = eigenvector_centrality(net)
println("Eigenvector: ", round.(ec, digits=3))

# PageRank (prestige via incoming links)
pr = pagerank(net)
println("PageRank: ", round.(pr, digits=3))
```

### Comparing Centrality Measures

Different centrality measures capture different aspects of importance:

| Measure | What It Captures | Best For |
|---------|-----------------|----------|
| Degree | Direct connections | Activity, visibility |
| Betweenness | Brokerage position | Gatekeepers, bridges |
| Closeness | Proximity to all | Efficiency, reach |
| Eigenvector | Quality of connections | Prestige, influence |
| PageRank | Incoming link prestige | Web/citation analysis |
| Bonacich power | Connections + indirect ties | Power, dependence |
| Katz | Weighted path counts | Status, influence |

### Directed Degree Options

For directed networks, degree centrality supports three modes:

```julia
# In-degree: number of incoming ties (popularity)
dc_in = degree_centrality(net; mode=:in)

# Out-degree: number of outgoing ties (activity)
dc_out = degree_centrality(net; mode=:out)

# Total degree: sum of in and out
dc_total = degree_centrality(net; mode=:total)

# Normalized by maximum possible degree
dc_norm = degree_centrality(net; mode=:total, normalized=true)
```

## Step 3: Analyze Network Structure

SNA.jl provides many network-level measures for characterizing overall structure:

```julia
# Global indices
println("Density:       ", round(density(net), digits=3))
println("Reciprocity:   ", round(reciprocity(net), digits=3))
println("Transitivity:  ", round(transitivity(net), digits=3))
println("Connectedness: ", round(connectedness(net), digits=3))
println("Hierarchy:     ", round(hierarchy(net), digits=3))
println("Efficiency:    ", round(efficiency(net), digits=3))

# Dyad census: mutual, asymmetric, null
dc_result = dyad_census(net)
println("Dyad census: M=$(dc_result.mutual), A=$(dc_result.asymmetric), N=$(dc_result.null)")

# Triad census: 16 isomorphism classes
tc = triad_census(net)
println("Triad census: ", tc)
```

### Cohesion Analysis

Assess how well-connected the network is:

```julia
# Connected components
comps = components(net; mode=:weak)
println("Weak components: ", length(comps))

scomps = components(net; mode=:strong)
println("Strong components: ", length(scomps))

# Largest component
lc = largest_component(net)
println("Largest component size: ", length(lc))

# Cutpoints (articulation points)
cp = cutpoints(net)
println("Cutpoints: ", cp)

# Bridges
br = bridges(net)
println("Bridges: ", br)

# K-cores
k2 = kcores(net; k=2)
println("2-core vertices: ", k2)

# Cliques (complete subgraphs)
cl = cliques(net; min_size=3)
println("Cliques (size >= 3): ", cl)
```

### Path Analysis

Understand distances and reachability in the network:

```julia
# Geodesic distance matrix
dist = geodesic_distance(net)
println("Distance matrix:\n", dist)

# Network diameter
d = diameter(net)
println("Diameter: ", d)

# Average path length
apl = average_path_length(net)
println("Average path length: ", round(apl, digits=3))

# Reachability matrix
reach = reachability(net)
println("Reachability (1 can reach 6): ", reach[1, 6])
```

## Step 4: Structural Equivalence and Blockmodels

Identify actors who occupy similar structural positions:

```julia
# Compute structural equivalence (correlation method)
se = structural_equivalence(net; method=:correlation)
println("Structural equivalence matrix:")
for i in 1:nv(net)
    println("  ", round.(se[i, :], digits=2))
end

# Cluster vertices by structural equivalence
assignments = equiv_clust(net; method=:structural, k=3)
println("Cluster assignments: ", assignments)

# Create a blockmodel
bm = blockmodel(net; k=3)
println("Block membership: ", bm.membership)
println("Block density matrix:")
for i in 1:bm.n_blocks
    println("  ", round.(bm.block_matrix[i, :], digits=2))
end
```

### Regular Equivalence

Regular equivalence is a weaker form: actors are equivalent if they have equivalent ties to equivalent others (rather than identical ties to the same others):

```julia
# Regular equivalence similarity
re = regular_equivalence(net)
println("Regular equivalence matrix:")
for i in 1:nv(net)
    println("  ", round.(re[i, :], digits=2))
end

# Cluster by regular equivalence
assignments_re = equiv_clust(net; method=:regular, k=3)
println("Regular equivalence clusters: ", assignments_re)
```

### Consensus Clustering

When you have multiple clustering solutions, compute a consensus:

```julia
# Multiple runs with different methods or parameters
clust1 = equiv_clust(net; method=:structural, k=3)
clust2 = equiv_clust(net; method=:regular, k=3)

# Consensus across clustering solutions
cons = consensus([clust1, clust2])
println("Consensus clustering: ", cons)
```

## Complete Example

```julia
using Network, SNA

# Build a small organizational communication network
net = network(8; directed=true)

# Core group (1-3): dense mutual ties
add_edge!(net, 1, 2); add_edge!(net, 2, 1)
add_edge!(net, 1, 3); add_edge!(net, 3, 1)
add_edge!(net, 2, 3); add_edge!(net, 3, 2)

# Periphery group (6-8): dense mutual ties
add_edge!(net, 6, 7); add_edge!(net, 7, 6)
add_edge!(net, 6, 8); add_edge!(net, 8, 6)
add_edge!(net, 7, 8); add_edge!(net, 8, 7)

# Bridge actors (4-5): connect core to periphery
add_edge!(net, 3, 4); add_edge!(net, 4, 3)
add_edge!(net, 4, 5); add_edge!(net, 5, 4)
add_edge!(net, 5, 6); add_edge!(net, 6, 5)

# --- Centrality Analysis ---
println("=== Centrality Analysis ===")
dc = degree_centrality(net; mode=:total, normalized=true)
bc = betweenness_centrality(net)

for i in 1:nv(net)
    println("Vertex $i: degree=$(round(dc[i], digits=2)), betweenness=$(round(bc[i], digits=3))")
end

# Identify the most central actor
most_central = argmax(bc)
println("\nMost central (betweenness): Vertex $most_central")

# --- Network Structure ---
println("\n=== Network Structure ===")
println("Density:      ", round(density(net), digits=3))
println("Reciprocity:  ", round(reciprocity(net), digits=3))
println("Transitivity: ", round(transitivity(net), digits=3))
println("Diameter:     ", diameter(net))

# --- Cohesion ---
println("\n=== Cohesion ===")
println("Components:   ", length(components(net)))
println("Cutpoints:    ", cutpoints(net))
println("Bridges:      ", bridges(net))

k2 = kcores(net; k=2)
println("2-core:       ", k2)

# --- Blockmodel ---
println("\n=== Blockmodel (3 blocks) ===")
bm = blockmodel(net; k=3)
println("Membership:   ", bm.membership)
println("Block matrix:")
for i in 1:bm.n_blocks
    println("  Block $i: ", round.(bm.block_matrix[i, :], digits=2))
end
```

## Best Practices

1. **Start with descriptives**: Compute density, reciprocity, and transitivity before centrality to understand the overall structure
2. **Use multiple centrality measures**: Different measures capture different aspects of importance; comparing them reveals richer insights
3. **Normalize for comparison**: Use `normalized=true` when comparing centrality across networks of different sizes
4. **Check connectivity**: Use `components()` to verify the network is connected before computing path-based measures like closeness or diameter
5. **Choose equivalence methods carefully**: Structural equivalence is strict (same ties to same alters); regular equivalence is more flexible (similar ties to equivalent alters)
6. **Validate blockmodels**: Try different values of `k` and compare block density matrices to find meaningful role structures

## Next Steps

- Learn about all [Centrality Measures](guide/centrality.md) in detail
- Explore [Network Measures](guide/measures.md) for global indices
- Understand [Cohesion](guide/cohesion.md) for subgroup analysis
- Dive into [Structural Equivalence](guide/equivalence.md) for positional analysis
