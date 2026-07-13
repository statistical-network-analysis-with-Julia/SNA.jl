# Cohesion

Cohesion analysis examines how well-connected a network is and identifies its substructures. SNA.jl provides functions for finding connected components, cliques, k-cores, cutpoints, bridges, biconnected components, and path-based measures like geodesic distance and diameter.

## Example Network

Throughout this guide, we use a network with clear subgroup structure:

```julia
using Networks, SNA

# Network with two dense groups connected by a bridge
net = network(8; directed=false)

# Group 1: vertices 1-4 (nearly complete)
add_edge!(net, 1, 2); add_edge!(net, 1, 3)
add_edge!(net, 1, 4); add_edge!(net, 2, 3)
add_edge!(net, 2, 4); add_edge!(net, 3, 4)

# Bridge: vertex 4 to vertex 5
add_edge!(net, 4, 5)

# Group 2: vertices 5-8 (nearly complete)
add_edge!(net, 5, 6); add_edge!(net, 5, 7)
add_edge!(net, 5, 8); add_edge!(net, 6, 7)
add_edge!(net, 6, 8); add_edge!(net, 7, 8)
```

```text
Network structure:
  Group 1         Bridge      Group 2
  1 - 2            |          6 - 7
  |\ |             |          |\ |
  | \|             |          | \|
  3 - 4 --------- 5 -------- 8
```

## Connected Components

Components are maximal connected subgraphs. Every vertex within a component can reach every other vertex in that component.

```julia
# Weakly connected components (ignoring edge direction)
comps = components(net; mode=:weak)
println("Number of components: ", length(comps))
println("Component sizes: ", [length(c) for c in comps])
println("Component members: ", comps)

# For directed networks: strongly connected components
net_dir = network(5; directed=true)
add_edge!(net_dir, 1, 2); add_edge!(net_dir, 2, 3)
add_edge!(net_dir, 3, 1); add_edge!(net_dir, 3, 4)
add_edge!(net_dir, 4, 5)

scomps = components(net_dir; mode=:strong)
println("Strong components: ", scomps)
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `mode` | `Symbol` | `:weak` | `:weak` (ignore direction) or `:strong` (respect direction) |

### Weak vs. Strong Components

| Type | Definition | Directed Networks |
|------|------------|-------------------|
| **Weak** | $i$ and $j$ are connected ignoring edge direction | Always at least as many members as strong |
| **Strong** | $i$ can reach $j$ AND $j$ can reach $i$ | May have more, smaller components |

### When to Use

Component analysis is the first step in most network analyses. Use it to:

- Determine if the network is connected
- Identify isolated subgroups
- Decide whether to analyze the full network or the largest component

## Largest Component

A convenience function that returns only the largest connected component:

```julia
lc = largest_component(net)
println("Largest component: ", lc)
println("Size: ", length(lc), " out of ", nv(net))

# For directed networks
lc_strong = largest_component(net_dir; mode=:strong)
println("Largest strong component: ", lc_strong)
```

### When to Use

Many network measures (closeness centrality, diameter, average path length) require a connected network. Extracting the largest component is a common preprocessing step.

## Cliques

A clique is a maximal complete subgraph: every vertex in the clique is connected to every other vertex in the clique.

```julia
# Find all cliques of size 3 or larger
cl = cliques(net; min_size=3)
println("Cliques (size >= 3): ", length(cl))
for (i, c) in enumerate(cl)
    println("  Clique $i: ", sort(c))
end

# Find larger cliques only
cl_large = cliques(net; min_size=4)
println("Cliques (size >= 4): ", length(cl_large))
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `min_size` | `Int` | `3` | Minimum clique size to return |

### Interpretation

- Cliques represent the densest possible subgroups
- The distribution of clique sizes reveals how clustered the network is
- Overlapping cliques (sharing members) indicate bridging actors

### Computational Notes

Finding all maximal cliques is NP-complete. SNA.jl uses the Bron-Kerbosch algorithm from Graphs.jl, which is efficient for sparse networks but can be slow for very dense graphs.

### When to Use

Clique analysis is appropriate for identifying cohesive subgroups in undirected networks. It is most useful when you expect well-defined groups with strong internal connections (e.g., friend groups, research teams, criminal organizations).

## K-Cores

The k-core is the maximal subgraph in which every vertex has degree at least $k$. K-cores provide a hierarchical decomposition of the network by connectivity.

```julia
# Find the 2-core (vertices with degree >= 2 in the subgraph)
k2 = kcores(net; k=2)
println("2-core: ", k2)

# Find the 3-core
k3 = kcores(net; k=3)
println("3-core: ", k3)

# Find coreness values for all vertices
using Graphs
core_numbers = Graphs.core_number(net.graph)
println("Core numbers: ", core_numbers)
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `k` | `Int` | `1` | Minimum degree requirement |

### Interpretation

- **Higher k-core**: More tightly connected subgroup
- The **coreness** of a vertex is the maximum $k$ for which it belongs to the k-core
- K-core decomposition reveals a nested hierarchy of increasingly dense subgraphs:
  - 1-core $\supseteq$ 2-core $\supseteq$ 3-core $\supseteq$ ...

### Relationship to Cliques

| Concept | Requirement | Properties |
|---------|------------|------------|
| **Clique** | All pairs connected | Very strict, may be small |
| **K-core** | Minimum degree $k$ | Relaxed, can be larger |

K-cores are less strict than cliques: a 3-core requires each member to have 3 connections within the subgraph, but not necessarily to every other member.

### When to Use

K-cores are useful for:

- Identifying the "core" vs. "periphery" of a network
- Filtering out peripheral vertices before further analysis
- Characterizing network resilience (higher k-core = more robust to vertex removal)
- As a faster alternative to clique detection

## Cutpoints (Articulation Points)

A cutpoint is a vertex whose removal disconnects the network (or increases the number of components).

```julia
cp = cutpoints(net)
println("Cutpoints: ", cp)
println("Number of cutpoints: ", length(cp))
```

### Interpretation

- Cutpoints are **critical vertices** that hold the network together
- Removing a cutpoint creates at least two disconnected components
- Networks with many cutpoints are **fragile** (vulnerable to targeted removal)

### Example

In our example network, vertices 4 and 5 are likely cutpoints because they form the bridge between the two groups. Removing either one would disconnect the network.

### When to Use

Cutpoint analysis is important for:

- Assessing network vulnerability and robustness
- Identifying critical actors whose departure would fragment the network
- Network design (ensuring redundancy around critical points)
- Understanding organizational resilience

## Bridges

A bridge is an edge whose removal disconnects the network.

```julia
br = bridges(net)
println("Bridges: ", br)
println("Number of bridges: ", length(br))
```

### Interpretation

- Bridges are **critical edges** without redundancy
- If a bridge is removed, the network splits into two or more components
- Bridges often connect different communities or groups

### Relationship to Cutpoints

- Every bridge has at least one endpoint that is a cutpoint (unless the bridge connects two isolated vertices)
- A cutpoint is adjacent to at least one bridge
- The bridge in our example is the edge (4, 5) connecting the two groups

### When to Use

Bridge analysis is useful for:

- Identifying bottlenecks in communication or flow networks
- Assessing the redundancy of connections between groups
- Understanding which relationships are critical for network cohesion

## Biconnected Components

A biconnected component is a maximal subgraph with no cutpoints. Every pair of vertices in a biconnected component has at least two vertex-disjoint paths between them.

```julia
bc = bicomponents(net)
println("Number of bicomponents: ", length(bc))
for (i, comp) in enumerate(bc)
    println("  Bicomponent $i: ", sort(comp))
end
```

### Interpretation

- Biconnected components are more robust than plain connected components
- Each bicomponent can tolerate the removal of any single vertex without disconnecting
- The bicomponent decomposition partitions the edge set (not the vertex set -- cutpoints belong to multiple bicomponents)

### Relationship to Other Concepts

| Concept | Tolerance | Definition |
|---------|-----------|------------|
| **Connected component** | None | Removing certain vertices disconnects it |
| **Biconnected component** | 1 vertex | No single vertex removal disconnects it |
| **K-vertex-connected** | k-1 vertices | Removing fewer than k vertices cannot disconnect it |

### When to Use

Bicomponent analysis is useful for:

- Understanding the block structure of a network
- Identifying independently robust subgraphs
- Analyzing the relationship between groups and bridges

## Geodesic Distance

The geodesic distance matrix contains the shortest path length between all pairs of vertices.

```julia
dist = geodesic_distance(net)
println("Distance from 1 to 8: ", dist[1, 8])
println("Distance from 1 to 2: ", dist[1, 2])

# Full distance matrix
println("Distance matrix:")
for i in 1:nv(net)
    println("  ", [dist[i, j] == Inf ? "Inf" : string(Int(dist[i, j])) for j in 1:nv(net)])
end
```

### Properties

- $d(i, i) = 0$ for all $i$
- $d(i, j) = 1$ if there is a direct edge from $i$ to $j$
- $d(i, j) = \infty$ if $j$ is not reachable from $i$
- For undirected networks: $d(i, j) = d(j, i)$
- For directed networks: $d(i, j) \neq d(j, i)$ in general

### When to Use

The geodesic distance matrix is the foundation for many analyses:

- Closeness centrality is based on row sums of the distance matrix
- Diameter is the maximum finite entry
- Average path length is the mean of finite off-diagonal entries
- Distance distributions reveal network structure

## Diameter

The diameter is the length of the longest shortest path in the network.

```julia
d = diameter(net)
println("Diameter: ", d)
```

### Formula

$$\text{Diameter} = \max_{i,j} d(i, j)$$

Where $d(i, j)$ is the geodesic distance, considering only finite distances.

### Interpretation

| Value | Interpretation |
|-------|----------------|
| 1 | Complete graph (every pair directly connected) |
| Small | "Small world" property |
| Large | Network has long chains |
| Inf | Network is disconnected |

### When to Use

Diameter is a simple measure of network "spread." It is commonly reported alongside average path length and is one of the defining properties of small-world networks (small diameter + high transitivity).

## Average Path Length

The average shortest path length over all reachable pairs:

```julia
apl = average_path_length(net)
println("Average path length: ", round(apl, digits=3))
```

### Formula

$$\bar{d} = \frac{1}{|\{(i,j) : d(i,j) < \infty, i \neq j\}|} \sum_{\substack{i \neq j \\ d(i,j) < \infty}} d(i,j)$$

### Interpretation

- **Small average path length**: Information or influence can travel quickly through the network
- Compare with random graphs of the same size and density: $\bar{d}_{random} \approx \frac{\ln n}{\ln \bar{k}}$

### When to Use

Average path length is essential for:

- Characterizing small-world properties
- Assessing the efficiency of information flow
- Comparing network structure across different networks

## Combining Cohesion Measures

A comprehensive cohesion analysis uses multiple measures together:

```julia
using Networks, SNA

net = network(10; directed=false)
# Create a network with interesting structure
add_edge!(net, 1, 2); add_edge!(net, 1, 3); add_edge!(net, 2, 3)
add_edge!(net, 3, 4); add_edge!(net, 4, 5); add_edge!(net, 4, 6)
add_edge!(net, 5, 6); add_edge!(net, 5, 7); add_edge!(net, 6, 7)
add_edge!(net, 7, 8); add_edge!(net, 8, 9); add_edge!(net, 8, 10)
add_edge!(net, 9, 10)

println("=== Cohesion Analysis ===")

# Connectivity
comps = components(net)
println("Connected: ", length(comps) == 1)
println("Components: ", length(comps))

# Vulnerability
cp = cutpoints(net)
br = bridges(net)
println("Cutpoints: ", cp, " (", length(cp), " total)")
println("Bridges: ", br, " (", length(br), " total)")

# Substructure
cl = cliques(net; min_size=3)
println("Cliques (size >= 3): ", length(cl))

k2 = kcores(net; k=2)
println("2-core: ", k2)

# Distances
d = diameter(net)
apl = average_path_length(net)
println("Diameter: ", d)
println("Avg path length: ", round(apl, digits=3))

# Bicomponents
bc = bicomponents(net)
println("Bicomponents: ", length(bc))
```

## Best Practices

1. **Start with components**: Always check connectivity before computing path-based measures. If the network is disconnected, either analyze components separately or focus on the largest component.

2. **Report vulnerability alongside connectivity**: A connected network with many cutpoints and bridges is fragile. Report both to give a complete picture.

3. **Use k-cores for large networks**: Clique detection is expensive for dense graphs. K-core decomposition is much faster and provides a useful alternative for identifying cohesive subgroups.

4. **Compare with random baselines**: Compute the expected diameter and clustering for random graphs of the same size and density to assess whether the observed structure is unusual.

5. **Combine cliques and cutpoints**: Identify cliques as cohesive subgroups and cutpoints as the actors bridging them. This reveals the meso-level structure of the network.

6. **Consider direction**: For directed networks, weak and strong components tell different stories. Strong components identify groups with mutual reachability, which is a stronger form of cohesion.
