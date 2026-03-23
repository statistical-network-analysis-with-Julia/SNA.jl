# Measures API Reference

This page documents network-level measures, path functions, random graph generators, and layout algorithms available in SNA.jl.

## Network-Level Measures

### Density and Reciprocity

```@docs
density
gden
reciprocity
grecip
mutuality
```

### Transitivity and Hierarchy

```@docs
transitivity
gtrans
hierarchy
efficiency
connectedness
```

## Census Functions

Functions for classifying dyads, triads, and components.

```@docs
dyad_census
triad_census
component_dist
```

## Path Functions

Functions for computing distances, reachability, and path-based summaries.

```@docs
geodesic_distance
reachability
diameter
average_path_length
```

## Random Graph Generators

```@docs
rgraph
rgnm
rgnp
```

## Layout Algorithms

Functions for computing vertex positions for network visualization.

```@docs
layout_fruchterman_reingold
layout_kamada_kawai
layout_circle
layout_random
```
