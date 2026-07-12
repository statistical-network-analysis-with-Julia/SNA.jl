# Network Measures

Network-level measures summarize the overall structure of a network with single numbers. SNA.jl provides measures for density, reciprocity, transitivity, hierarchy, efficiency, and more, along with census functions that provide detailed dyad- and triad-level breakdowns.

## Example Network

Throughout this guide, we use a small directed network:

```julia
using Network, SNA

net = network(6; directed=true)
add_edge!(net, 1, 2); add_edge!(net, 2, 1)  # Mutual
add_edge!(net, 1, 3)                          # Asymmetric
add_edge!(net, 2, 3); add_edge!(net, 3, 2)   # Mutual
add_edge!(net, 3, 4)                          # Asymmetric
add_edge!(net, 4, 5); add_edge!(net, 5, 4)   # Mutual
add_edge!(net, 4, 6)                          # Asymmetric
add_edge!(net, 5, 6)                          # Asymmetric
```

## Density

Network density is the proportion of possible edges that actually exist. It is the most basic measure of network structure.

```julia
d = density(net)
println("Density: ", round(d, digits=3))

# R-compatible alias
d = gden(net)
```

### Formula

For a directed network with $n$ vertices and $m$ edges:

$$\Delta = \frac{m}{n(n-1)}$$

For an undirected network:

$$\Delta = \frac{2m}{n(n-1)}$$

### Interpretation

| Density Value | Interpretation |
|---------------|----------------|
| 0.0 | No edges (empty graph) |
| 0.0 - 0.1 | Sparse network |
| 0.1 - 0.5 | Moderate density |
| 0.5 - 1.0 | Dense network |
| 1.0 | Complete graph (all possible edges exist) |

### When to Use

Density provides a quick summary of how connected a network is. It is useful for:

- Comparing networks of similar size
- Setting expectations for other measures (many measures are correlated with density)
- Checking whether a network is sparse enough for certain algorithms

!!! note
    Density is inversely related to network size for most real-world networks. Comparing densities across networks of very different sizes can be misleading.

## Reciprocity

Reciprocity measures the extent to which ties are mutual in directed networks. It answers the question: when actor A ties to actor B, does actor B also tie to actor A?

```julia
# Dyadic reciprocity (proportion of non-null dyads that are mutual)
r = reciprocity(net; method=:dyadic)
println("Dyadic reciprocity: ", round(r, digits=3))

# Edgewise reciprocity (proportion of edges that are reciprocated)
r_edge = reciprocity(net; method=:edgewise)
println("Edgewise reciprocity: ", round(r_edge, digits=3))

# R-compatible alias
r = grecip(net)
```

### Formula

**Dyadic reciprocity**:

$$R_{dyadic} = \frac{M}{M + A}$$

Where $M$ is the number of mutual dyads and $A$ is the number of asymmetric dyads.

**Edgewise reciprocity**:

$$R_{edge} = \frac{2M}{2M + A}$$

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `method` | `Symbol` | `:dyadic` | `:dyadic` or `:edgewise` |

### Interpretation

- **High reciprocity** (close to 1.0): Most ties are mutual (common in friendship, collaboration)
- **Low reciprocity** (close to 0.0): Most ties are one-directional (common in advice-seeking, citation)
- **Undirected networks**: Reciprocity is always 1.0 (all ties are inherently mutual)

### Difference Between Methods

| Method | Denominator | Interpretation |
|--------|-------------|----------------|
| `:dyadic` | Non-null dyads ($M + A$) | What proportion of connected pairs are mutual? |
| `:edgewise` | All edges ($2M + A$) | What proportion of edges are reciprocated? |

The edgewise method gives higher values because mutual dyads contribute two edges, while asymmetric dyads contribute only one.

## Transitivity

Transitivity (also called the clustering coefficient) measures the tendency for triadic closure: if A is connected to B and B is connected to C, how often is A also connected to C?

```julia
# Global clustering coefficient
t = transitivity(net; type=:global)
println("Global transitivity: ", round(t, digits=3))

# Local clustering coefficients (per vertex)
t_local = transitivity(net; type=:local)
println("Local transitivity: ", round.(t_local, digits=3))

# Average local clustering coefficient
t_avg = transitivity(net; type=:average)
println("Average transitivity: ", round(t_avg, digits=3))

# R-compatible alias
t = gtrans(net)
```

### Formula

**Global transitivity**:

$$T = \frac{3 \times \text{number of triangles}}{\text{number of connected triples}}$$

**Local clustering coefficient** for vertex $i$:

$$C_i = \frac{\text{triangles involving } i}{\binom{k_i}{2}}$$

Where $k_i$ is the degree of vertex $i$.

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `type` | `Symbol` | `:global` | `:global`, `:local`, or `:average` |

### Interpretation

- **High transitivity** (> 0.3): Strong clustering, "friends of friends are friends"
- **Low transitivity** (< 0.1): Little triadic closure
- **Global vs. average**: Global transitivity weights vertices by degree; average gives equal weight to all vertices

### When to Use

Transitivity is one of the most widely reported network measures. Use it to:

- Characterize the clustering tendency of the network
- Compare with random graphs (Erdos-Renyi graphs have transitivity $\approx p$, the edge probability)
- Identify small-world properties (high transitivity + short path lengths)

## Dyad Census

The dyad census classifies all pairs of vertices into three categories: mutual, asymmetric, and null.

```julia
dc = dyad_census(net)
println("Mutual:     ", dc.mutual)
println("Asymmetric: ", dc.asymmetric)
println("Null:       ", dc.null)
```

### Categories

| Type | Description | Notation |
|------|-------------|----------|
| **Mutual** | Both $i \to j$ and $j \to i$ exist | $M$ |
| **Asymmetric** | Either $i \to j$ or $j \to i$ exists (not both) | $A$ |
| **Null** | Neither $i \to j$ nor $j \to i$ exists | $N$ |

### Relationship to Other Measures

The dyad census provides the building blocks for several other measures:

$$M + A + N = \binom{n}{2}$$

$$\text{Reciprocity (dyadic)} = \frac{M}{M + A}$$

$$\text{Density} = \frac{2M + A}{n(n-1)}$$

### When to Use

The dyad census is the foundation of structural analysis in directed networks. Use it when you need a complete description of pairwise relationships, or when computing derived measures like reciprocity. It is also useful for testing whether a network differs from random expectations (e.g., more mutual ties than expected by chance).

## Triad Census

The triad census extends the dyad census to all triples of vertices. For directed networks, there are 16 isomorphism classes of triads, labeled using the MAN notation:

```julia
tc = triad_census(net)
println("Triad census (16 classes): ", tc)
```

The census is computed with the edge-driven Batagelj–Mrvar (2001)
algorithm: only triads containing at least one tie are enumerated (each
exactly once), and the empty-triad count is recovered by subtraction from
$\binom{n}{3}$. The cost therefore scales with the number of edges rather
than $O(n^3)$, so large sparse networks are censused quickly.

### The 16 Triad Types (MAN Notation)

| Index | Label | M-A-N | Description |
|-------|-------|-------|-------------|
| 1 | 003 | 0-0-3 | Three null dyads (empty triad) |
| 2 | 012 | 0-1-2 | One asymmetric, two null |
| 3 | 102 | 1-0-2 | One mutual, two null |
| 4 | 021D | 0-2-1 | Two asymmetric (both outward from one vertex) |
| 5 | 021U | 0-2-1 | Two asymmetric (both inward to one vertex) |
| 6 | 021C | 0-2-1 | Two asymmetric (chain: $A \to B \to C$) |
| 7 | 111D | 1-1-1 | One mutual, one asymmetric (directed outward) |
| 8 | 111U | 1-1-1 | One mutual, one asymmetric (directed inward) |
| 9 | 030T | 0-3-0 | Three asymmetric (transitive: $A \to B \to C, A \to C$) |
| 10 | 030C | 0-3-0 | Three asymmetric (cyclic: $A \to B \to C \to A$) |
| 11 | 201 | 2-0-1 | Two mutual, one null |
| 12 | 120D | 1-2-0 | One mutual, two asymmetric (outward) |
| 13 | 120U | 1-2-0 | One mutual, two asymmetric (inward) |
| 14 | 120C | 1-2-0 | One mutual, two asymmetric (mixed) |
| 15 | 210 | 2-1-0 | Two mutual, one asymmetric |
| 16 | 300 | 3-0-0 | Three mutual (complete triad) |

### Interpretation

The triad census reveals structural tendencies beyond pairwise analysis:

- **High 030T count**: Transitive closure (hierarchical structure)
- **High 030C count**: Cyclic closure (generalized exchange)
- **High 300 count**: Dense local clustering
- **High 003 count**: Sparse network

### When to Use

The triad census is the standard tool for characterizing local structure in directed networks. It is used in:

- Testing structural balance theory (excess of certain triad types)
- Comparing observed networks to random baselines
- Identifying tendencies toward hierarchy vs. reciprocity vs. transitivity

## Mutuality

Mutuality measures the proportion of connected pairs that have symmetric (mutual) relationships.

```julia
m = mutuality(net)
println("Mutuality: ", round(m, digits=3))
```

### Formula

$$\text{Mutuality} = \frac{M}{M + A}$$

This is equivalent to dyadic reciprocity.

### When to Use

Mutuality is a focused measure of symmetry in relationships. Use it when the theoretical distinction between symmetric and asymmetric ties is important (e.g., mutual friendship vs. one-sided admiration).

## Hierarchy

Krackhardt's hierarchy measure quantifies the extent to which reachability relations are asymmetric. A perfectly hierarchical network has no reciprocal reachability (if $i$ can reach $j$, then $j$ cannot reach $i$).

```julia
h = hierarchy(net)
println("Hierarchy: ", round(h, digits=3))
```

### Formula

$$H = \frac{V}{V + W}$$

Where $V$ is the number of ordered pairs $(i, j)$ such that $i$ can reach $j$ but $j$ cannot reach $i$, and $W$ is the number of pairs where both can reach each other.

### Interpretation

| Value | Interpretation |
|-------|----------------|
| 0.0 | No hierarchy (all reachable pairs are mutually reachable) |
| 0.5 | Mixed hierarchy |
| 1.0 | Perfect hierarchy (strict top-down reachability) |

### When to Use

Hierarchy is appropriate for analyzing organizational structures, power networks, and any directed network where asymmetric reachability is theoretically meaningful. It is one of Krackhardt's four graph theoretical dimensions of informal organizations.

## Efficiency

Network efficiency measures how close the network is to a tree (a minimally connected graph). Trees are "efficient" because they use the minimum number of edges to connect all vertices.

```julia
e = efficiency(net)
println("Efficiency: ", round(e, digits=3))
```

### Formula

$$E = 1 - \frac{m - (n-1)}{m_{max} - (n-1)}$$

Where $m$ is the number of edges, $n-1$ is the minimum for connectivity (a tree), and $m_{max}$ is the maximum possible edges.

### Interpretation

| Value | Interpretation |
|-------|----------------|
| 1.0 | Tree-like structure (minimal edges for connectivity) |
| 0.5 | Moderate redundancy |
| 0.0 | Complete graph (maximum redundancy) |

### When to Use

Efficiency is useful for comparing the "leanness" of network structures. It is one of Krackhardt's four graph theoretical dimensions and is relevant in organizational design (efficient vs. redundant communication structures).

## Connectedness

Connectedness measures the proportion of vertex pairs where at least one can reach the other through some directed path.

```julia
c = connectedness(net)
println("Connectedness: ", round(c, digits=3))
```

### Formula

$$C = \frac{\text{pairs where } i \text{ can reach } j \text{ or } j \text{ can reach } i}{\binom{n}{2}}$$

### Interpretation

| Value | Interpretation |
|-------|----------------|
| 1.0 | Weakly connected (every pair has at least one-directional reachability) |
| 0.0 | Completely disconnected (no pairs are connected) |

### When to Use

Connectedness is the weakest connectivity measure. It assesses the basic reachability structure and is one of Krackhardt's four graph theoretical dimensions. Use it alongside hierarchy, efficiency, and lubness to characterize organizational structure.

## Component Distribution

The component distribution shows the sizes of connected components, revealing fragmentation.

```julia
cd = component_dist(net)
println("Component sizes (descending): ", cd)
```

### Interpretation

- A single component containing all vertices indicates a connected network
- Multiple components indicate fragmentation
- The distribution shape reveals whether fragmentation is balanced (similar-sized components) or skewed (one giant component plus small isolates)

## Reachability

The reachability matrix indicates which vertices can reach which others through directed paths of any length.

```julia
reach = reachability(net)
println("Can vertex 1 reach vertex 6? ", reach[1, 6])
println("Can vertex 6 reach vertex 1? ", reach[6, 1])
```

### Formula

$$R_{ij} = \begin{cases} 1 & \text{if there exists a directed path from } i \text{ to } j \\ 0 & \text{otherwise} \end{cases}$$

The diagonal entries $R_{ii}$ are always 1 (every vertex can reach itself).

### When to Use

Reachability is the foundation for hierarchy, connectedness, and strong components. Use it when you need to know which actors can potentially influence or communicate with which others through any chain of intermediaries.

## Using Measures Together

### Krackhardt's Graph Theoretical Dimensions

Krackhardt (1994) proposed four measures for characterizing informal organizations:

```julia
println("=== Krackhardt's Dimensions ===")
println("Connectedness: ", round(connectedness(net), digits=3))
println("Hierarchy:     ", round(hierarchy(net), digits=3))
println("Efficiency:    ", round(efficiency(net), digits=3))
# Lubness requires additional computation (not yet implemented)
```

These four dimensions together characterize the shape of an organization's informal structure.

### Complete Network Profile

```julia
println("=== Network Profile ===")
println("Size:          ", nv(net), " vertices, ", ne(net), " edges")
println("Density:       ", round(density(net), digits=3))
println("Reciprocity:   ", round(reciprocity(net), digits=3))
println("Transitivity:  ", round(transitivity(net), digits=3))
println("Connectedness: ", round(connectedness(net), digits=3))
println("Hierarchy:     ", round(hierarchy(net), digits=3))
println("Efficiency:    ", round(efficiency(net), digits=3))

dc = dyad_census(net)
println("\nDyad census:")
println("  Mutual:     ", dc.mutual)
println("  Asymmetric: ", dc.asymmetric)
println("  Null:       ", dc.null)

tc = triad_census(net)
println("\nTriad census: ", tc)
```

### Comparing Networks

When comparing two networks, compute a consistent set of measures:

```julia
function network_profile(net, label)
    println("--- $label ---")
    println("  n=$(nv(net)), m=$(ne(net))")
    println("  Density:     $(round(density(net), digits=3))")
    println("  Reciprocity: $(round(reciprocity(net), digits=3))")
    println("  Transitivity:$(round(transitivity(net), digits=3))")
    println("  Components:  $(length(components(net)))")
end

# Compare two networks
net1 = rgnp(20, 0.15)
net2 = rgnp(20, 0.30)

network_profile(net1, "Sparse (p=0.15)")
network_profile(net2, "Dense (p=0.30)")
```

## Best Practices

1. **Report density alongside other measures**: Many structural measures are correlated with density, so reporting density provides context for interpretation.

2. **Use the dyad census for directed networks**: The dyad census gives a complete picture of pairwise relationships and is the basis for many higher-level measures.

3. **Compare with random baselines**: Compute the expected values under a random graph model (e.g., Erdos-Renyi) to assess whether observed values are unusual.

4. **Consider network size**: Most measures are sensitive to network size. Use normalized versions when comparing networks of different sizes.

5. **Use multiple complementary measures**: Density, reciprocity, transitivity, and the census functions capture different aspects of structure. Reporting several together gives a richer picture.

6. **Check assumptions**: Reciprocity and the triad census require directed networks. Transitivity can be computed for both directed and undirected networks but has different interpretations.
