# Cohesion

Cohesion analysis asks how well a network holds together and which dense
subgroups it contains. SNA.jl provides components (`component_dist`,
`largest_component`), cliques, k-cores (`kcores`), cutpoints and
bicomponents, with R `sna`'s names where sna has the measure.

## Example network

Two complete groups of four actors joined by a single tie, a pendant actor
and an isolate:

```julia
using NetworkCore, SNA

net = network(10; directed=false)
for (i, j) in [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4),   # group 1
               (4, 5),                                           # the joining tie
               (5, 6), (5, 7), (5, 8), (6, 7), (6, 8), (7, 8),   # group 2
               (8, 9)]                                           # a pendant
    add_edge!(net, i, j)
end
```

```text
  1 ─ 2               6 ─ 7
  │ ╳ │               │ ╳ │
  3 ─ 4 ───────────── 5 ─ 8 ─ 9        10 (isolate)
```

In each group every pair is tied (the `╳` stands for the two diagonals).

## Components

```julia
cd = component_dist(net)
cd.membership                # [1, 1, 1, 1, 1, 1, 1, 1, 1, 2]
cd.csize                     # [9, 1]
cd.cdist                     # one component of size 1, one of size 9
largest_component(net)       # actors 1 to 9
```

`component_dist` is R's `component.dist`: the component of each actor
(numbered by their lowest actor), the size of each component, and the number
of components of each size. On a directed network `connected=:strong` (the
default) requires mutual reachability, `:weak` ignores direction and
`:recursive` uses mutual ties only:

```julia
dnet = network(5; directed=true)
for (i, j) in [(1, 2), (2, 3), (3, 1), (3, 4), (4, 5)]
    add_edge!(dnet, i, j)
end
component_dist(dnet).csize                   # [3, 1, 1]: the cycle and two singletons
component_dist(dnet; connected=:weak).csize  # [5]
```

## Cliques

A clique is a maximal complete subgraph. On a directed network a clique
needs mutual ties (`symmetrize=:strong`, as sna's `clique.census`);
`symmetrize=:weak` accepts a tie in either direction.

```julia
cliques(net)                 # the two groups, {1, 2, 3, 4} and {5, 6, 7, 8}
cliques(net; min_size=2)     # also the dyads {4, 5} and {8, 9}
```

Enumerating cliques takes exponential time in the worst case (the
Bron–Kerbosch algorithm of Graphs.jl); on large or dense networks prefer
k-cores.

## K-cores

The k-core is the largest subgraph in which every actor has at least `k`
ties to other members. `kcores` returns each actor's core number, the
largest `k` whose core contains it, as R `sna::kcores` does:

```julia
kcores(net)                  # [3, 3, 3, 3, 3, 3, 3, 3, 1, 0]
findall(>=(3), kcores(net))  # the 3-core: the two groups
```

On a directed network `cmode` selects the degree: `:freeman` (in + out,
the default), `:indegree` or `:outdegree`. Self-loops are ignored.

Graphs.jl's `core_number` gives the same core numbers when called on the
`Network` itself. Do not call Graphs.jl functions on the field `net.graph`:
an undirected `Network` stores each tie as two arcs there, so every
Graphs.jl algorithm would count each tie twice (core numbers of 6 instead of
3 here).

## Cutpoints and bridges

A cutpoint is an actor whose removal increases the number of components;
here the two ends of the joining tie, and actor 8, which holds the pendant:

```julia
sort(cutpoints(net))         # [4, 5, 8]
```

On a directed network `cutpoints` uses strong components by default, as
sna does (`connected=:weak` ignores direction, `:recursive` uses mutual
ties).

A bridge is a tie whose removal increases the number of components. sna has
no bridge function; Graphs.jl's `bridges` works on an undirected `Network`:

```julia
using Graphs
Graphs.bridges(net)          # the ties 4–5 and 8–9
```

## Bicomponents

A bicomponent is a maximal subgraph that no single actor's removal can
disconnect. `bicomponents` returns each as a list of ties. As in R's
`bicomponent.dist`, a bridge on its own (two actors) is not a bicomponent
unless `min_size=2`:

```julia
length(bicomponents(net))                # 2: the two groups
length(bicomponents(net; min_size=2))    # 4: also the two bridges
```

Directed networks use mutual ties by default (`symmetrize=:strong`).

## Distances

```julia
g = geodist(net)
g.gdist[1, 8]                # 3: 1 → 4 → 5 → 8
g.counts[1, 8]               # 1 geodesic
maximum(g.gdist)             # Inf: the isolate is unreachable
average_path_length(net)     # mean over the pairs that are connected
```

`geodist` is sna's `geodist`: the geodesic distance between every ordered
pair (`Inf` when unreachable) and the number of geodesics. sna has no
diameter function; the diameter is `maximum(geodist(net).gdist)`, which is
`Inf` for a disconnected network. Path-based indices (closeness, diameter)
are often reported for the largest component:

```julia
core = largest_component(net)
sub = NetworkCore.get_induced_subgraph(net, core)
maximum(geodist(sub).gdist)  # 4: from 1, 2 or 3 to 9
```

## Practice

1. Check the components before computing path-based indices; analyse
   components separately or report the largest.
2. Report vulnerability with connectivity: a connected network with many
   cutpoints and bridges is fragile.
3. Use k-cores rather than cliques on large networks.
4. Compare with a random baseline of the same size and density (`rgraph`,
   `rgnm`).
5. On directed networks, say whether components are weak or strong.
