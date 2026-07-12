# Structural Equivalence

Structural equivalence analysis identifies actors who occupy similar positions in the network. SNA.jl provides functions for computing structural and regular equivalence, clustering vertices by position, building blockmodels, and finding consensus across multiple clustering solutions.

## Key Concepts

| Concept | Definition |
|---------|------------|
| **Structural equivalence** | Two actors are equivalent if they have identical ties to and from all other actors |
| **Regular equivalence** | Two actors are equivalent if they have equivalent ties to equivalent others |
| **Blockmodel** | A reduced representation of the network based on equivalence classes |
| **Position** | A group of structurally similar actors |

### Structural vs. Regular Equivalence

The distinction between structural and regular equivalence is fundamental:

**Structural equivalence** (Lorrain & White, 1971): Actors $i$ and $j$ are structurally equivalent if for every other actor $k$:

- $i \to k$ if and only if $j \to k$
- $k \to i$ if and only if $k \to j$

**Regular equivalence** (White & Reitz, 1983): Actors $i$ and $j$ are regularly equivalent if:

- For every $k$ tied to $i$, there exists some $l$ tied to $j$ such that $k$ and $l$ are also regularly equivalent

In practice:

- **Structural equivalence** requires actors to be connected to the *same* others
- **Regular equivalence** requires actors to be connected to *similar* others

```text
Example:
  A → C, A → D        Structural equivalence: A ≡ B (same ties)
  B → C, B → D

  A → C                Regular equivalence: A ≡ B (both send to one other)
  B → D                but NOT structurally equivalent (different targets)
```

## Example Network

Throughout this guide, we use a network with clear positional structure:

```julia
using Network, SNA

# Manager-subordinate network
net = network(8; directed=true)

# Managers (1, 2): send to subordinates
add_edge!(net, 1, 3); add_edge!(net, 1, 4)  # Manager 1 → Workers 3, 4
add_edge!(net, 2, 5); add_edge!(net, 2, 6)  # Manager 2 → Workers 5, 6

# Workers (3-6): send to clients
add_edge!(net, 3, 7); add_edge!(net, 4, 7)  # Workers 3, 4 → Client 7
add_edge!(net, 5, 8); add_edge!(net, 6, 8)  # Workers 5, 6 → Client 8

# Cross-ties
add_edge!(net, 1, 2); add_edge!(net, 2, 1)  # Managers communicate
add_edge!(net, 7, 1); add_edge!(net, 8, 2)  # Clients report to managers
```

```text
Network structure:
  Manager 1 ↔ Manager 2
    ↓   ↓       ↓   ↓
    3   4       5   6
    ↓   ↓       ↓   ↓
  Client 7    Client 8
    ↓               ↓
  Manager 1    Manager 2
```

## Computing Structural Equivalence

The `structural_equivalence` function computes a similarity or distance matrix between all pairs of vertices based on their tie profiles.

```julia
# Correlation-based similarity (default)
se_cor = structural_equivalence(net; method=:correlation)
println("Structural equivalence (correlation):")
for i in 1:nv(net)
    println("  ", round.(se_cor[i, :], digits=2))
end
```

### Methods

Three methods are available for measuring structural equivalence:

```julia
# Pearson correlation of tie profiles
se_cor = structural_equivalence(net; method=:correlation)

# Euclidean distance between tie profiles
se_euc = structural_equivalence(net; method=:euclidean)

# Hamming distance (proportion of different ties)
se_ham = structural_equivalence(net; method=:hamming)
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `method` | `Symbol` | `:correlation` | `:correlation`, `:euclidean`, or `:hamming` |

### How It Works

For each vertex $i$, a **tie profile** is constructed by concatenating:

1. The $i$-th row of the adjacency matrix (outgoing ties)
2. The $i$-th column of the adjacency matrix (incoming ties)

This creates a $2n$-dimensional vector for each vertex. The similarity or distance between two vertices is then computed using the chosen method.

### Method Comparison

| Method | Range | Meaning | Best For |
|--------|-------|---------|----------|
| `:correlation` | $[-1, 1]$ | 1 = identical pattern, 0 = no relationship, -1 = opposite | General use, handles degree differences |
| `:euclidean` | $[0, \infty)$ | 0 = identical, larger = more different | When magnitude matters |
| `:hamming` | $[0, 1]$ | 0 = identical, 1 = completely different | Binary networks |

### Interpreting the Matrix

```julia
# Find the most equivalent pair
se = structural_equivalence(net; method=:correlation)
max_sim = -Inf
best_pair = (0, 0)
for i in 1:nv(net)
    for j in (i+1):nv(net)
        if se[i, j] > max_sim
            max_sim = se[i, j]
            best_pair = (i, j)
        end
    end
end
println("Most equivalent pair: ", best_pair, " (similarity: ", round(max_sim, digits=3), ")")
```

## Computing Regular Equivalence

Regular equivalence uses the REGE algorithm, which iteratively refines similarity estimates based on the similarities of neighbors.

```julia
# Compute regular equivalence
re = regular_equivalence(net)
println("Regular equivalence:")
for i in 1:nv(net)
    println("  ", round.(re[i, :], digits=2))
end

# With custom convergence settings
re = regular_equivalence(net; max_iter=200, tol=1e-8)
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `max_iter` | `Int` | `100` | Maximum iterations for convergence |
| `tol` | `Float64` | `1e-6` | Convergence tolerance |

### How It Works

The REGE algorithm:

1. Initialize all pairwise similarities to 1.0
2. For each pair $(i, j)$:
   - Match $i$'s out-neighbors to $j$'s out-neighbors using current similarity
   - Match $i$'s in-neighbors to $j$'s in-neighbors
   - Update $sim(i, j)$ as the average of best-match similarities
3. Repeat until convergence

### Structural vs. Regular Equivalence Results

In our example network:

- **Structural equivalence**: Vertices 3 and 4 are highly equivalent (same ties to same others). Vertices 3 and 5 are not (different targets).
- **Regular equivalence**: Vertices 3 and 5 may be equivalent (both are workers who send to a client and receive from a manager, even though specific targets differ).

## Equivalence Clustering

The `equiv_clust` function clusters vertices into positions based on equivalence:

```julia
# Cluster by structural equivalence
assign_struct = equiv_clust(net; method=:structural, k=3)
println("Structural clusters: ", assign_struct)

# Cluster by regular equivalence
assign_regular = equiv_clust(net; method=:regular, k=3)
println("Regular clusters: ", assign_regular)

# Auto-detect number of clusters
assign_auto = equiv_clust(net; method=:structural)
println("Auto-detected clusters: ", assign_auto)
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `method` | `Symbol` | `:structural` | `:structural` or `:regular` |
| `k` | `Int` or `nothing` | `nothing` | Number of clusters (auto-detected if `nothing`) |

### How It Works

1. Compute the equivalence matrix (structural or regular)
2. Convert similarities to distances ($d = 1 - sim$)
3. Apply k-medoids clustering to the distance matrix
4. Return cluster assignments for each vertex

### Choosing the Number of Clusters

The choice of $k$ is crucial:

```julia
# Try different values of k
for k in 2:5
    assign = equiv_clust(net; method=:structural, k=k)
    println("k=$k: ", assign)
end
```

Guidelines for choosing $k$:

- Start with theoretical expectations (e.g., manager/worker/client = 3 roles)
- Try multiple values and compare blockmodel fit
- Use the elbow method: plot within-cluster distance vs. $k$

## Blockmodels

A blockmodel reduces the network to a simpler structure by grouping equivalent vertices into blocks and summarizing ties between blocks.

```julia
# Create a 3-block model
bm = blockmodel(net; k=3)
println("Membership: ", bm.membership)
println("Number of blocks: ", bm.n_blocks)
println("\nBlock density matrix:")
for i in 1:bm.n_blocks
    println("  Block $i: ", round.(bm.block_matrix[i, :], digits=3))
end
```

### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `k` | `Int` | Required | Number of blocks |
| `method` | `Symbol` | `:structural` | `:structural` or `:regular` |

### Return Value

A `NamedTuple` with:

| Field | Type | Description |
|-------|------|-------------|
| `membership` | `Vector{Int}` | Block assignment for each vertex |
| `block_matrix` | `Matrix{Float64}` | Density of ties between blocks |
| `n_blocks` | `Int` | Number of blocks |

### Interpreting the Block Matrix

The block matrix shows the density of ties between (and within) blocks:

```julia
bm = blockmodel(net; k=3)

println("Block Density Matrix:")
println("         Block 1  Block 2  Block 3")
for i in 1:bm.n_blocks
    print("Block $i:  ")
    for j in 1:bm.n_blocks
        print("$(round(bm.block_matrix[i, j], digits=2))     ")
    end
    println()
end
```

Common patterns in the block matrix:

| Pattern | Description | Example |
|---------|-------------|---------|
| High diagonal | Cohesive groups | Communities |
| Off-diagonal > diagonal | Between-group ties dominate | Bipartite structure |
| Upper triangle only | Hierarchical flow | Authority structure |
| One row/column high | Core-periphery | Star network |

### Ideal Block Patterns

After computing the block matrix, you can compare with ideal types:

```julia
# Check if a block has density above threshold
threshold = 0.5
println("One-blocks (density > $threshold):")
for i in 1:bm.n_blocks
    for j in 1:bm.n_blocks
        if bm.block_matrix[i, j] > threshold
            println("  Block $i → Block $j: $(round(bm.block_matrix[i, j], digits=2))")
        end
    end
end
```

## Consensus Clustering

When you have multiple clustering solutions (from different methods, parameters, or random starts), consensus clustering finds the most stable grouping.

```julia
# Generate multiple clustering solutions
clust1 = equiv_clust(net; method=:structural, k=3)
clust2 = equiv_clust(net; method=:regular, k=3)
clust3 = equiv_clust(net; method=:structural, k=4)

# Find consensus
cons = consensus([clust1, clust2, clust3])
println("Consensus clustering: ", cons)
```

### How It Works

1. Build a **co-occurrence matrix**: for each pair $(i, j)$, count how often they are in the same cluster across solutions
2. Normalize by the number of solutions
3. Threshold the co-occurrence matrix (default: 0.5)
4. Assign vertices to clusters based on thresholded co-occurrence

### When to Use

Consensus clustering is useful when:

- Results are sensitive to the choice of method or parameters
- You want a robust partitioning that is stable across specifications
- Combining structural and regular equivalence results

## Complete Example

```julia
using Network, SNA

# Create a network with role structure
net = network(10; directed=true)

# Two "senders" (1, 2): send to everyone in group A
add_edge!(net, 1, 3); add_edge!(net, 1, 4); add_edge!(net, 1, 5)
add_edge!(net, 2, 3); add_edge!(net, 2, 4); add_edge!(net, 2, 5)

# Group A (3, 4, 5): send to group B
add_edge!(net, 3, 6); add_edge!(net, 3, 7)
add_edge!(net, 4, 6); add_edge!(net, 4, 7)
add_edge!(net, 5, 8); add_edge!(net, 5, 9)

# Group B (6, 7, 8, 9): send to receiver (10)
add_edge!(net, 6, 10); add_edge!(net, 7, 10)
add_edge!(net, 8, 10); add_edge!(net, 9, 10)

# --- Step 1: Structural Equivalence ---
println("=== Structural Equivalence ===")
se = structural_equivalence(net; method=:correlation)
println("Most similar pairs:")
for i in 1:nv(net)
    for j in (i+1):nv(net)
        if se[i, j] > 0.8
            println("  ($i, $j): $(round(se[i, j], digits=3))")
        end
    end
end

# --- Step 2: Regular Equivalence ---
println("\n=== Regular Equivalence ===")
re = regular_equivalence(net)
println("Regular equivalence clusters:")
assign_re = equiv_clust(net; method=:regular, k=4)
println("  Assignments: ", assign_re)

# --- Step 3: Blockmodel ---
println("\n=== Blockmodel (4 blocks) ===")
bm = blockmodel(net; k=4, method=:structural)
println("Membership: ", bm.membership)
println("\nBlock density matrix:")
for i in 1:bm.n_blocks
    println("  Block $i: ", round.(bm.block_matrix[i, :], digits=2))
end

# --- Step 4: Interpret Roles ---
println("\n=== Role Interpretation ===")
for b in 1:bm.n_blocks
    members = findall(==(b), bm.membership)
    println("Block $b: vertices $members")

    # Characterize the block's ties
    sends_to = findall(x -> x > 0.3, bm.block_matrix[b, :])
    receives_from = findall(x -> x > 0.3, bm.block_matrix[:, b])
    println("  Sends to blocks: $sends_to")
    println("  Receives from blocks: $receives_from")
end
```

## Comparing Equivalence Methods

| Feature | Structural | Regular |
|---------|-----------|---------|
| **Strictness** | Very strict | Relaxed |
| **Requirement** | Same ties to same actors | Similar ties to equivalent actors |
| **Example** | Two people who know exactly the same people | Two managers (different teams, same role) |
| **Computation** | Direct (matrix operations) | Iterative (REGE algorithm) |
| **Result** | Fewer, smaller groups | Fewer, larger groups |
| **Use case** | Exact position analysis | Role analysis |

## Best Practices

1. **Start with structural equivalence**: It is simpler and more interpretable. Move to regular equivalence only if structural analysis produces too many groups or if the theoretical question is about roles rather than positions.

2. **Try multiple methods**: Compute equivalence using correlation, Euclidean, and Hamming distance. If results are consistent, you can be more confident in the groupings.

3. **Validate blockmodels**: Compare the block density matrix with theoretical expectations. A good blockmodel should produce a clear and interpretable pattern.

4. **Use consensus clustering**: When results are sensitive to method or parameters, consensus clustering provides a more robust solution.

5. **Report the full pipeline**: When presenting results, report the equivalence method, the number of clusters (and how it was chosen), and the block density matrix.

6. **Consider network size**: Structural equivalence becomes very strict in large networks (exact position matches are rare). Regular equivalence or relaxed methods may be more appropriate for large networks.

7. **Combine with other analyses**: Use centrality measures to characterize each block (e.g., "the high-betweenness broker block" vs. "the low-degree peripheral block"). This enriches the positional analysis.

## References

- Lorrain, F., & White, H.C. (1971). Structural equivalence of individuals in social networks. *Journal of Mathematical Sociology*, 1(1), 49-80.
- White, D.R., & Reitz, K.P. (1983). Graph and semigroup homomorphisms on networks of relations. *Social Networks*, 5(2), 193-234.
- Burt, R.S. (1976). Positions in networks. *Social Forces*, 55(1), 93-122.
- Wasserman, S., & Faust, K. (1994). *Social Network Analysis: Methods and Applications*, Chapter 9-10. Cambridge University Press.
- Doreian, P., Batagelj, V., & Ferligoj, A. (2005). *Generalized Blockmodeling*. Cambridge University Press.
