# Cohesion API Reference

Components, cliques, k-cores, cutpoints, bicomponents and geodesics.

## Components

```@docs
component_dist
largest_component
bicomponents
```

## Subgroups

```@docs
cliques
kcores
```

## Vulnerability

```@docs
cutpoints
```

Bridges (ties whose removal disconnects the network) are Graphs.jl's
`bridges`, which works on an undirected `Network`.

## Paths

```@docs
geodist
reachability
average_path_length
```
