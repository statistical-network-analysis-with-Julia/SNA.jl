# Structural Equivalence

Equivalence analysis groups actors who occupy similar positions. SNA.jl
provides structural equivalence distances (`sedist`, R `sna::sedist`),
regular equivalence similarity, hierarchical clustering into positions
(`equiv_clust`), blockmodels, and a consensus of several partitions.

## Two notions of equivalence

**Structural equivalence** (Lorrain and White 1971): actors $i$ and $j$ are
equivalent when they have the same ties to and from every *other* actor $k$:
$i \to k$ exactly when $j \to k$, and $k \to i$ exactly when $k \to j$.
Ties between $i$ and $j$ themselves do not count.

**Regular equivalence** (White and Reitz 1983): $i$ and $j$ are equivalent
when they have ties to and from *equivalent* others: for every actor tied
to $i$ there is an equivalent actor tied to $j$, and the reverse.

Structural equivalence asks "same contacts?"; regular equivalence asks
"same kind of contacts?". Two managers with different subordinates are
regularly but not structurally equivalent.

## Example network

A boss (1), two managers (2, 3) and three workers under each manager:

```julia
using NetworkCore, SNA

net = network(9; directed=true)
for (i, j) in [(1, 2), (1, 3),                    # boss → managers
               (2, 4), (2, 5), (2, 6),            # manager 2 → workers
               (3, 7), (3, 8), (3, 9)]            # manager 3 → workers
    add_edge!(net, i, j)
end
```

```text
              1
            ↙   ↘
          2       3
        ↙ ↓ ↘   ↙ ↓ ↘
       4  5  6  7  8  9
```

The structural positions are {1}, {2}, {3}, {4, 5, 6} and {7, 8, 9}: each
worker has the same tie as their co-workers. The regular positions are
{1}, {2, 3} and {4, …, 9}: both managers receive from a boss and send to
workers.

## Structural equivalence: `sedist`

```julia
sedist(net)                          # Hamming distances (sna's default)
sedist(net; method=:correlation)     # profile correlations (a similarity)
sedist(net)[4, 5]                    # 0: co-workers are equivalent
sedist(net)[4, 7]                    # 2: different managers
```

Each actor's profile is its row (ties sent) and its column (ties received)
of the adjacency matrix. When comparing $i$ and $j$, the cells that are not
about third parties are left out: both self cells and the ties between $i$
and $j$, as in sna and in Wasserman and Faust (1994, eq. 9.4). So two
actors with identical ties to everyone else are at distance 0 whether or not
they are tied to each other (`diag=true` keeps every cell).

| `method` | Value |
|---|---|
| `:hamming` (default) | the number of cells in which the profiles differ |
| `:correlation` | the Pearson correlation of the profiles (1 = equivalent) |
| `:euclidean` | the Euclidean distance between the profiles |
| `:gamma` | Goodman–Kruskal gamma between the profiles |
| `:exact` | 0 if the profiles are identical, else 1 |

## Regular equivalence

```julia
sim = regular_equivalence(net)
sim[2, 3]                            # 1.0: both managers
sim[4, 7]                            # 1.0: all workers
sim[1, 2]                            # ≈ 0: the boss receives no tie
```

`regular_equivalence` returns a similarity in [0, 1] computed by iterative
neighbour matching from a start in which every pair is equivalent; it
converges to the coarsest (maximal) regular equivalence. Because of that
start, on a network where every actor both sends and receives ties it can
leave every pair at 1. It is in the spirit of White and Reitz's definition
but is not UCINET's REGE algorithm, and sna has no counterpart. `maxiter`
and `tol` control the iteration.

## Clustering into positions

```julia
equiv_clust(net; k=5)        # [1, 2, 3, 4, 4, 4, 5, 5, 5]: the structural positions
equiv_clust(net; k=3, equiv_fun=regular_equivalence)   # [1, 2, 2, 3, 3, 3, 3, 3, 3]
equiv_clust(net; k=4)        # [1, 2, 3, 4, 4, 4, 4, 4, 4]
```

`equiv_clust` is R's `equiv.clust` followed by `cutree(k)`: it clusters the
`sedist` Hamming distance by complete linkage, and returns the labels of
the cut at `k` clusters (R returns the tree). `method` selects another
`sedist` method, `cluster_method` the linkage (`:complete`, `:average`,
`:single`), and `equiv_fun=regular_equivalence` clusters regular
equivalence instead.

With five clusters the structural positions are recovered exactly, and
with three the regular ones. Asked for four structural positions, the
clustering merges the two groups of workers, which are closest (distance
2). Clustering is a heuristic; check the result against the distance matrix
and the substance of the network, and try several `k`.

## Blockmodels

```julia
bm = blockmodel(net; k=3, equiv_fun=regular_equivalence)
bm.membership                        # [1, 2, 2, 3, 3, 3, 3, 3, 3]
bm.block_matrix                      # tie density between positions
```

A blockmodel reduces the network to its positions and the density of ties
between them. Here the image is a chain: the boss sends to every manager
(density 1), each manager to half of the workers (density 0.5), and no other
block has ties. A diagonal block's density is taken over the ordered pairs
of distinct actors in it, so the boss's own block, a single actor, is
undefined (`NaN`), as in sna.

## Brokerage

Given a partition of the actors, Gould and Fernandez's brokerage roles
count how often each actor stands between two others who are not directly
tied (`i → j → k` without `i → k`), by the groups of the three:

```julia
b = brokerage(net, [1, 2, 2, 3, 3, 3, 3, 3, 3])    # boss, managers, workers
b.roles                      # ["w_I", "w_O", "b_IO", "b_OI", "b_O", "t"]
b.raw_nli[2, :]              # manager 2 is a liaison (b_O) for 3 pairs
b.raw_gli                    # network totals: 6 liaison paths
b.z_gli                      # z-scores against a random network of the same density
```

`brokerage` is R's `brokerage`: role counts per actor (`raw_nli`) and in
total (`raw_gli`), with their expectations, standard deviations and
z-scores under Gould and Fernandez's density-conditioned null model.

## Consensus of several partitions

```julia
partitions = [equiv_clust(net; k=k) for k in 3:5]
consensus_clustering(partitions)
```

`consensus_clustering` links two actors when they share a cluster in at
least half of the partitions (`threshold`) and returns the connected groups
of that graph, numbered by their lowest actor. The result does not depend on
the order of the actors. (R's `sna::consensus` is a different procedure,
which combines informants' reports into one network; SNA.jl's function was
called `consensus` before 0.2.0.)

## Practice

1. Choose the notion that matches the question: structural equivalence for
   "who could substitute for whom", regular equivalence for roles.
2. Inspect the distance or similarity matrix, not only the clusters.
3. Try several numbers of positions and compare the block images.
4. Missing (masked) ties are refused by default; equivalence on partially
   observed data needs `missing=:face`, which reads the stored values.
