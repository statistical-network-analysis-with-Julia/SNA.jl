# Cohesion API Reference

This page documents cohesion functions available in SNA.jl.

## Components

Functions for identifying connected subgroups.

```@docs
components
largest_component
bicomponents
```

## Subgroup Detection

Functions for finding dense substructures within the network.

```@docs
cliques
kcores
```

## Vulnerability

Functions for identifying vertices and edges whose removal disconnects the network.

```@docs
cutpoints
bridges
```
