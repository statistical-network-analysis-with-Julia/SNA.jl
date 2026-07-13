# Centrality

Centrality measures quantify the importance or prominence of vertices in a network. SNA.jl provides eight centrality measures, each capturing a different aspect of structural importance. All centrality functions return a `Vector{Float64}` with one score per vertex.

## Centrality Interface

All centrality functions follow a common pattern:

<!-- skip-check -->
```julia
centrality_function(net; keyword_arguments...) -> Vector{Float64}
```

The result vector is indexed by vertex ID, so `result[i]` gives the centrality score for vertex `i`.

## Example Network

Throughout this guide, we use a small directed network to illustrate each measure:

```julia
using Networks, SNA

# Create example network
net = network(7; directed=true)
add_edge!(net, 1, 2); add_edge!(net, 1, 3)
add_edge!(net, 2, 1); add_edge!(net, 2, 4)
add_edge!(net, 3, 4); add_edge!(net, 3, 5)
add_edge!(net, 4, 5); add_edge!(net, 4, 6)
add_edge!(net, 5, 6); add_edge!(net, 5, 7)
add_edge!(net, 6, 7); add_edge!(net, 7, 4)
```

```text
Network structure:
    1 ↔ 2 → 4 ← 7
    ↓       ↓ ↗  ↑
    3 → 4   5 → 6
        ↓   ↑
        5 → 6
```

## Degree Centrality

The simplest centrality measure: the number of direct connections a vertex has.

```julia
# Total degree (in + out)
dc = degree_centrality(net; mode=:total)

# In-degree only (popularity)
dc_in = degree_centrality(net; mode=:in)

# Out-degree only (activity)
dc_out = degree_centrality(net; mode=:out)

# Normalized to [0, 1]
dc_norm = degree_centrality(net; mode=:total, normalized=true)
```

### Formula

For a directed network with $n$ vertices:

- **In-degree**: $C_D^{in}(i) = \sum_j A_{ji}$
- **Out-degree**: $C_D^{out}(i) = \sum_j A_{ij}$
- **Total degree**: $C_D(i) = C_D^{in}(i) + C_D^{out}(i)$
- **Normalized**: $C_D'(i) = \frac{C_D(i)}{n - 1}$ (directed) or $\frac{C_D(i)}{2(n-1)}$ (undirected)

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `mode` | `Symbol` | `:total` | `:in`, `:out`, or `:total` |
| `normalized` | `Bool` | `false` | Normalize by maximum possible degree |

### Interpretation

- **High in-degree**: Popular or prestigious actors (many others choose them)
- **High out-degree**: Active or expansive actors (they choose many others)
- **High total degree**: Well-connected actors overall

### When to Use

Degree centrality is appropriate when direct connections are the primary concern. It is fast to compute ($O(n)$) and easy to interpret. Use it as a baseline before examining more complex measures.

## Betweenness Centrality

Betweenness centrality measures how often a vertex lies on shortest paths between other vertices. High-betweenness actors serve as brokers or bridges.

```julia
# Default: normalized
bc = betweenness_centrality(net)

# Unnormalized (raw counts)
bc_raw = betweenness_centrality(net; normalized=false)
```

### Formula

$$C_B(i) = \sum_{s \neq i \neq t} \frac{\sigma_{st}(i)}{\sigma_{st}}$$

Where $\sigma_{st}$ is the total number of shortest paths from $s$ to $t$, and $\sigma_{st}(i)$ is the number of those paths passing through $i$.

The normalized version divides by $(n-1)(n-2)$ for directed networks or $(n-1)(n-2)/2$ for undirected networks.

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `normalized` | `Bool` | `true` | Normalize by number of vertex pairs |

### Interpretation

- **High betweenness**: The vertex is a gatekeeper or broker who controls information flow
- **Zero betweenness**: The vertex is on no shortest paths between others (peripheral or redundant)

### When to Use

Betweenness is ideal for identifying brokers, gatekeepers, and potential bottlenecks. It is especially informative in networks with clear group structure, where bridges between groups have high betweenness. Computation is $O(nm)$ where $m$ is the number of edges.

## Closeness Centrality

Closeness centrality measures how close a vertex is to all other vertices, based on shortest path distances. Close actors can reach everyone quickly.

```julia
# Default: normalized
cc = closeness_centrality(net)

# Unnormalized
cc_raw = closeness_centrality(net; normalized=false)
```

### Formula

$$C_C(i) = \frac{1}{\sum_j d(i, j)}$$

Where $d(i, j)$ is the geodesic (shortest path) distance from $i$ to $j$.

The normalized version multiplies by $(n-1)$:

$$C_C'(i) = \frac{n - 1}{\sum_j d(i, j)}$$

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `normalized` | `Bool` | `true` | Normalize by `(n-1)` |

### Interpretation

- **High closeness**: The vertex can quickly reach all others (short average distance)
- **Low closeness**: The vertex is on the periphery, far from many others

### Caveats

Closeness centrality is not well-defined for disconnected networks (distances to unreachable vertices are infinite). For disconnected graphs, consider using closeness only within the largest connected component, or use harmonic closeness (the sum of inverse distances).

### When to Use

Closeness is useful when speed of communication or transmission is important. It identifies actors who are well-positioned to disseminate information or coordinate activity efficiently.

## Eigenvector Centrality

A vertex has high eigenvector centrality if it is connected to other high-centrality vertices. It captures the idea that not all connections are equal.

```julia
# Default parameters
ec = eigenvector_centrality(net)

# With custom convergence settings
ec = eigenvector_centrality(net; max_iter=200, tol=1e-8)
```

### Formula

Eigenvector centrality is the leading eigenvector of the adjacency matrix:

$$Ax = \lambda_1 x$$

Where $A$ is the adjacency matrix and $\lambda_1$ is the largest eigenvalue. The centrality vector $x$ satisfies:

$$x_i = \frac{1}{\lambda_1} \sum_j A_{ij} x_j$$

SNA.jl computes this via power iteration.

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `max_iter` | `Int` | `100` | Maximum number of power iterations |
| `tol` | `Float64` | `1e-6` | Convergence tolerance |

### Interpretation

- **High eigenvector centrality**: Connected to many well-connected actors (quality over quantity)
- **Low eigenvector centrality**: Connected to peripheral actors or isolated

### Relationship to Other Measures

Eigenvector centrality generalizes degree centrality: degree counts the number of neighbors, while eigenvector centrality weights neighbors by their own centrality. PageRank and Katz centrality are variants of this idea with different damping mechanisms.

### When to Use

Eigenvector centrality is appropriate when prestige or influence matters more than raw connectivity. It is commonly used in citation networks, social media influence, and organizational hierarchy analysis.

## Bonacich Power Centrality

Bonacich power centrality extends eigenvector centrality with a tunable attenuation parameter $\beta$ that controls how much weight is given to indirect connections.

```julia
# Default: exponent (β) = 1.0
bp = bonacich_power(net)

# Negative beta (power from connections to weak actors)
bp_neg = bonacich_power(net; exponent=-0.3)

# Rescaled so scores sum to 1 (as in R sna)
bp_scaled = bonacich_power(net; rescale=true)
```

### Formula

$$c(\beta) = (I - \beta A)^{-1} A \mathbf{1}$$

Where $\beta$ is the attenuation parameter and must satisfy $|\beta| < 1/\lambda_1$ for convergence.

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `exponent` | `Float64` | `1.0` | Attenuation factor β; requires $|\beta| < 1/\lambda_1$ |
| `rescale` | `Bool` | `false` | If true, rescale scores to sum to 1 |
| `tol` | `Float64` | `1e-7` | Solver tolerance for detecting singularity |

### The Role of Beta

The sign and magnitude of $\beta$ have substantive meaning:

| $\beta$ Value | Interpretation |
|---------------|----------------|
| $\beta > 0$ | Power comes from connections to powerful actors (prestige) |
| $\beta < 0$ | Power comes from connections to weak/dependent actors (dominance) |
| $\beta \to 0$ | Reduces to degree centrality |
| $|\beta|$ large | More weight on indirect ties |

### Interpretation

- **Positive $\beta$**: Actors gain power by being connected to other powerful actors. This captures prestige or influence in cooperative settings.
- **Negative $\beta$**: Actors gain power by being connected to weak or dependent actors. This captures bargaining power in exchange networks.

### When to Use

Bonacich power centrality is ideal for exchange and bargaining networks where the distinction between prestige (positive $\beta$) and dominance (negative $\beta$) is theoretically important. It is widely used in economic sociology and organizational studies.

## Katz Centrality

Katz centrality counts the number of paths from a vertex to all others, with longer paths weighted less via a damping factor $\alpha$.

```julia
# Default parameters
kc = katz_centrality(net)

# Custom damping
kc = katz_centrality(net; α=0.05)
```

### Formula

$$C_K(i) = \sum_{k=1}^{\infty} \sum_j \alpha^k (A^k)_{ij}$$

In matrix form:

$$\mathbf{c} = ((I - \alpha A)^{-1} - I) \mathbf{1}$$

Where $\alpha < 1/\lambda_1$ for convergence.

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `α` | `Float64` | `0.1` | Damping factor (must be less than $1/\lambda_1$) |
| `β` | `Float64` | `1.0` | Weight of direct connections |

### Interpretation

- **High Katz centrality**: The vertex can reach many others through many paths
- Longer paths contribute less due to the $\alpha^k$ damping

### Relationship to Other Measures

Katz centrality is closely related to Bonacich power centrality (with positive $\beta$) and eigenvector centrality. As $\alpha$ approaches $1/\lambda_1$, Katz centrality converges to eigenvector centrality.

### When to Use

Katz centrality is appropriate when influence can spread through indirect connections. It is used in social influence models, information diffusion analysis, and status hierarchies.

## PageRank

PageRank was originally developed for ranking web pages but is widely used in network analysis. It models a random surfer who follows links with probability $\alpha$ and teleports to a random vertex with probability $1 - \alpha$.

```julia
# Default: alpha=0.85
pr = pagerank(net)

# Lower damping (more random teleportation)
pr = pagerank(net; α=0.5)

# Custom convergence
pr = pagerank(net; α=0.85, max_iter=200, tol=1e-8)
```

### Formula

$$PR(i) = \frac{1 - \alpha}{n} + \alpha \sum_{j \to i} \frac{PR(j)}{d_{out}(j)}$$

Where $\alpha$ is the damping factor and $d_{out}(j)$ is the out-degree of vertex $j$.

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `α` | `Float64` | `0.85` | Damping factor (probability of following a link) |
| `max_iter` | `Int` | `100` | Maximum number of iterations |
| `tol` | `Float64` | `1e-6` | Convergence tolerance |

### Interpretation

- **High PageRank**: The vertex receives many incoming links from high-PageRank vertices
- PageRank values sum to 1 (they form a probability distribution)
- The damping factor $\alpha$ controls how much weight is given to the link structure vs. uniform distribution

### When to Use

PageRank is ideal for directed networks where incoming links indicate endorsement, citation, or influence. It handles dangling nodes (vertices with no outgoing links) gracefully via the teleportation mechanism. It is the standard for web page ranking, citation analysis, and any network where prestige flows through directed links.

## Flow Betweenness

Flow betweenness considers all paths (not just shortest paths) when measuring brokerage. It is based on maximum network flow between all pairs of vertices.

```julia
fb = flowbet(net)
```

### Interpretation

- **High flow betweenness**: The vertex carries a large proportion of the maximum flow between other vertex pairs
- More robust than standard betweenness because it accounts for all paths, not just shortest ones

### When to Use

Flow betweenness is appropriate when information or resources can travel along non-shortest paths (e.g., electrical networks, fluid flow). It is computationally more expensive than standard betweenness.

!!! note
    The current implementation returns standard betweenness as an approximation. Full flow betweenness based on maximum flow computation is planned for a future release.

## Comparing Centrality Measures

Different measures can produce very different rankings:

```julia
using Networks, SNA

net = network(7; directed=true)
add_edge!(net, 1, 2); add_edge!(net, 1, 3)
add_edge!(net, 2, 1); add_edge!(net, 2, 4)
add_edge!(net, 3, 4); add_edge!(net, 3, 5)
add_edge!(net, 4, 5); add_edge!(net, 4, 6)
add_edge!(net, 5, 6); add_edge!(net, 5, 7)
add_edge!(net, 6, 7); add_edge!(net, 7, 4)

dc = degree_centrality(net; mode=:total, normalized=true)
bc = betweenness_centrality(net)
cc = closeness_centrality(net)
ec = eigenvector_centrality(net)
pr = pagerank(net)

println("Vertex | Degree | Between. | Closeness | Eigenvec. | PageRank")
println("-------|--------|----------|-----------|-----------|--------")
for i in 1:nv(net)
    println("   $i   | $(round(dc[i], digits=3)) |  $(round(bc[i], digits=3))  |   $(round(cc[i], digits=3))  |   $(round(ec[i], digits=3))  |  $(round(pr[i], digits=3))")
end
```

### Summary Table

| Measure | Based On | Complexity | Directed | Disconnected |
|---------|----------|------------|----------|--------------|
| Degree | Direct ties | $O(n)$ | In/out/total | Yes |
| Betweenness | Shortest paths | $O(nm)$ | Yes | Yes |
| Closeness | Shortest paths | $O(nm)$ | Yes | Problematic |
| Eigenvector | Spectral | $O(n^2 k)$ | Yes | Yes |
| Bonacich | Matrix inverse | $O(n^3)$ | Yes | Yes |
| Katz | Path counting | $O(n^2)$ | Yes | Yes |
| PageRank | Random walk | $O(n^2 k)$ | Yes | Yes |
| Flow betw. | Max flow | $O(n^3 m)$ | Yes | Yes |

### Choosing the Right Measure

| Research Question | Best Measure |
|-------------------|-------------|
| Who has the most connections? | Degree |
| Who brokers between groups? | Betweenness |
| Who can reach others fastest? | Closeness |
| Who is connected to the "right" people? | Eigenvector |
| Who has power in exchange networks? | Bonacich (negative $\beta$) |
| Who has prestige via endorsements? | PageRank |
| Who controls flow in the network? | Flow betweenness |
| Who has high status? | Katz |

## Graph Centralization

Beyond vertex-level scores, Freeman centralization summarizes how
concentrated a centrality measure is in a single number: 0 when every
vertex scores equally, 1 (for the classic measures) when the network is a
perfect star. `centralization` follows R `sna::centralization`, computing
$C = \sum_i (c_{max} - c_i) / C_{max}$ with sna's theoretical maxima:

```julia
using Networks, SNA

star = network(5; directed=false)
for v in 2:5
    add_edge!(star, 1, v)
end

println(centralization(star, :degree))       # 1.0 — the star is maximal
println(centralization(star, :betweenness))  # 1.0
println(centralization(star, :closeness))    # 1.0
```

Supported measures are `:degree` (with `mode=:in/:out/:total` for directed
networks), `:betweenness`, `:closeness`, and `:eigenvector`; pass
`normalized=false` for the raw deviation sum.

## Best Practices

1. **Report multiple measures**: No single centrality measure captures all aspects of importance. Report at least two (e.g., degree and betweenness) to provide a more complete picture.

2. **Normalize for comparison**: When comparing centrality across networks of different sizes, always use normalized versions.

3. **Check for disconnected components**: Closeness centrality is problematic for disconnected networks. Either restrict analysis to the largest component or use harmonic closeness.

4. **Examine the distribution**: Plot the distribution of centrality scores. Highly skewed distributions (e.g., power-law degree distributions) indicate centralization.

5. **Consider direction**: For directed networks, in-degree and out-degree may tell very different stories. Always consider which direction is theoretically meaningful.

6. **Validate substantively**: High centrality scores should be interpretable in the context of the network. Cross-reference with domain knowledge.
