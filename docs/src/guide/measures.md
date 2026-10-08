# Network Measures

Graph-level indices summarize a whole network with one number. SNA.jl
provides R `sna`'s indices under sna's names: `gden` (density), `grecip`
(reciprocity), `gtrans` (transitivity), `dyad_census`, `triad_census`,
`mutuality`, `hierarchy`, `efficiency` and `connectedness`, plus local
clustering coefficients (`transitivity`) and path summaries (`geodist`,
`reachability`, `average_path_length`).

## Example network

```julia
using NetworkCore, SNA

net = network(6; directed=true)
for (i, j) in [(1, 2), (2, 1), (1, 3), (2, 3), (3, 2), (3, 4),
               (4, 5), (5, 4), (4, 6), (5, 6)]
    add_edge!(net, i, j)
end
```

## Density

```julia
gden(net)                    # 10 ties / 30 ordered pairs = 0.3333
```

Density is the proportion of possible ties that are present:
$m / n(n-1)$ for a directed network, $m / (n(n-1)/2)$ for an undirected one
(each undirected tie counts once). Loops are excluded unless `diag=true`.
For a two-mode network, `discount_bipartite=true` restricts the denominator
to cross-mode pairs.

`Graphs.density` is Graphs.jl's function: on a `Network` it keeps Graphs.jl's
definition, which counts self-loops.

## Reciprocity

With $M$, $A$ and $N$ the numbers of mutual, asymmetric and null dyads:

```julia
grecip(net)                              # (M + N) / (M + A + N): sna's default
grecip(net; measure=:dyadic_nonnull)     # M / (M + A)
grecip(net; measure=:edgewise)           # 2M / (2M + A)
grecip(net; measure=:edgewise_lrr)       # log of edgewise reciprocity over density
grecip(net; measure=:correlation)        # correlation of y_ij and y_ji
```

| `measure` | Question |
|---|---|
| `:dyadic` | What share of pairs is symmetric (both tied or both untied)? |
| `:dyadic_nonnull` | What share of tied pairs is mutual? |
| `:edgewise` | What share of ties is returned? |
| `:edgewise_lrr` | How much more often is a tie returned than density predicts? |
| `:correlation` | How correlated are the two directions of a pair? |

Every pair of an undirected network is symmetric, so its reciprocity is 1.
An undefined ratio (for instance `:edgewise` on an empty network) is `NaN`.

## Transitivity and clustering

```julia
gtrans(net)                              # share of two-paths i→j→k closed by i→k
gtrans(net; measure=:strong)             # i→k exactly when a two-path exists
gtrans(net; measure=:weakcensus)         # the number of closed two-paths
```

`gtrans` is sna's graph transitivity. The default weak measure is the share
of two-paths `i→j→k` (distinct actors) that are closed by `i→k`; on an
undirected network this is $3 \times$ triangles / connected triples. With no
two-path at all it is 1, as in sna.

Local clustering coefficients are not in sna; `transitivity` provides them,
named after igraph's function:

```julia
transitivity(net; type=:local)           # one coefficient per actor
transitivity(net; type=:average)         # mean over actors that have one
transitivity(net; type=:global)          # the same as gtrans(net)
```

On an undirected network the local coefficient of actor $i$ with $k_i$
neighbours is the share of the $k_i(k_i-1)/2$ neighbour pairs that are tied
(Watts–Strogatz). On a directed network the default (`cmode=:total`) is
Fagiolo's (2007) total coefficient,

$$C_i = \frac{[(A + A^T)^3]_{ii}}{2\,[d_i^{tot}(d_i^{tot} - 1) - 2 d_i^{\leftrightarrow}]},$$

the share of all possible directed triangles through $i$ that are present
($d_i^{tot}$ is in- plus out-degree, $d_i^{\leftrightarrow}$ the number of
mutual ties). It equals the undirected coefficient when every tie is
mutual. `cmode=:weak` instead ignores direction, as igraph does. An actor
with fewer than two neighbours has no coefficient (`NaN`), and the average
leaves it out; self-loops are never neighbours.

## Dyad and triad census

```julia
dyad_census(net)             # (mutual = 3, asymmetric = 4, null = 8)
triad_census(net)            # 16 counts in the order below
```

The dyad census counts the $\binom{n}{2}$ pairs by type: mutual, asymmetric
or null. Density and reciprocity are functions of it: on a directed network
$m = 2M + A$.

The triad census counts the $\binom{n}{3}$ triples of a directed network in
the 16 Davis–Leinhardt classes, in R `sna::triad.census`'s order (an
undirected network has four classes, by number of ties):

| Index | Class | Description |
|---|---|---|
| 1 | 003 | empty |
| 2 | 012 | one asymmetric tie |
| 3 | 102 | one mutual tie |
| 4 | 021D | two ties out of one actor |
| 5 | 021U | two ties into one actor |
| 6 | 021C | a two-path |
| 7 | 111D | a mutual tie and an asymmetric tie into it |
| 8 | 111U | a mutual tie and an asymmetric tie out of it |
| 9 | 030T | a transitive triple |
| 10 | 030C | a cycle |
| 11 | 201 | two mutual ties |
| 12 | 120D | a mutual tie and two ties out of the third actor |
| 13 | 120U | a mutual tie and two ties into the third actor |
| 14 | 120C | a mutual tie and a two-path |
| 15 | 210 | two mutual ties and an asymmetric one |
| 16 | 300 | complete |

The census is computed by the edge-driven Batagelj–Mrvar (2001) algorithm,
which enumerates only triads that contain a tie, so its cost grows with the
number of ties rather than with $n^3$.

## Mutuality

```julia
mutuality(net)               # 3: the number of mutual dyads
```

As in sna, `mutuality` is a count, not a proportion; for proportions use
`grecip`.

## Krackhardt's dimensions

Krackhardt (1994) characterized informal organizations by connectedness,
hierarchy, efficiency and least-upper-boundedness. SNA.jl provides the first
three, as in sna (`lubness` is not implemented).

```julia
connectedness(net)                       # share of pairs joined by a semipath
hierarchy(net; measure=:krackhardt)      # share of reachable pairs reachable one way only
efficiency(net)                          # 1 − excess ties / maximum possible excess
hierarchy(net)                           # sna's default: 1 − grecip(net)
```

- **Connectedness** is the share of unordered pairs in the same weak
  component (1 for a weakly connected network).
- **Hierarchy** with `measure=:krackhardt` is the share of pairs in which
  one actor can reach the other but not the reverse, among pairs where at
  least one can reach the other. sna's default `measure=:reciprocity` is
  `1 − grecip(net)`. With no reachable pair the ratio is undefined (`NaN`).
- **Efficiency** compares, per weak component of size $n_c$, the ties beyond
  the $n_c - 1$ needed to connect it with the most there could be:
  $E = 1 - \sum (m_c - (n_c - 1)) / \sum (n_c(n_c - 1) - (n_c - 1))$,
  counting an undirected tie as two arcs and ignoring self-loops, as sna
  does. A network of isolates has no possible excess and gives `NaN`.

## Components

```julia
cd = component_dist(net)                 # strong components
cd.membership                            # the component of each actor
cd.csize                                 # the size of each component
cd.cdist                                 # the number of components of each size
component_dist(net; connected=:weak).csize
```

`component_dist` is sna's `component.dist`: components are numbered by
their lowest actor. `connected` is `:strong` (default), `:weak` or
`:recursive` (mutual ties only). See the [cohesion guide](cohesion.md).

## Paths and reachability

```julia
g = geodist(net)
g.gdist[1, 6]                # length of a shortest path from 1 to 6
g.counts[1, 6]               # number of shortest paths
maximum(g.gdist)             # the diameter (Inf: not every pair is reachable)
reachability(net)[6, 1]      # false
average_path_length(net)     # mean distance over reachable ordered pairs
```

`geodist` is sna's `geodist`: distances (`Inf` when unreachable, or
`inf_replace`) and geodesic counts. `reachability` is sna's: entry `(i, j)`
says whether a directed path leads from `i` to `j`, and every actor reaches
itself.

## A network profile

```julia
function profile(net)
    dc = dyad_census(net)
    return (n = nv(net), ties = ne(net), density = gden(net),
            reciprocity = grecip(net; measure=:edgewise),
            transitivity = gtrans(net),
            components = length(component_dist(net; connected=:weak).csize),
            mutual = dc.mutual)
end

using Random
profile(net)
profile(rgnp(20, 0.15; rng=Xoshiro(1)))
```

Most of these indices depend on density and network size, so compare
networks of similar size, or against a random-graph baseline (`rgraph`,
`rgnm`).

## Conventions

- Self-loops are ignored by every index (sna's `diag=FALSE`); `gden` and
  `degreecent` count them with `diag=true`. With `diag=true`, an undirected
  looped network's density has `n(n+1)/2` possible dyads (NetworkCore'
  convention), while sna sums the symmetric matrix over `n²` cells.
- Undefined ratios are `NaN`, as in sna; `gtrans` with no two-path is 1, as
  in sna.
- A network with fewer than three actors has a triad census of zeros (sna
  returns `NaN`).
