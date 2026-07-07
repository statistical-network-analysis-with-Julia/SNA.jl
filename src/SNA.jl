"""
    SNA.jl - Social Network Analysis for Julia

A Julia package providing descriptive analysis tools for social networks,
including centrality measures, structural equivalence, cohesion analysis,
and network visualization layouts.

Port of the R sna package from the StatNet collection.
"""
module SNA

using Graphs
using LinearAlgebra
using Random
using Statistics
using Network

# Centrality measures
export degree_centrality, betweenness_centrality, closeness_centrality
export eigenvector_centrality, bonacich_power, katz_centrality
export pagerank, flowbet

# Network-level measures
export density, reciprocity, transitivity, mutuality
export dyad_census, triad_census, component_dist
export hierarchy, efficiency, connectedness

# Cohesion
export components, cliques, kcores, cutpoints, bridges
export bicomponents, largest_component

# Structural equivalence
export structural_equivalence, regular_equivalence
export blockmodel, equiv_clust, consensus

# Paths
export geodesic_distance, reachability, diameter
export average_path_length

# Graph-level indices
export gden, grecip, gtrans

# Visualization layouts
export layout_fruchterman_reingold, layout_kamada_kawai
export layout_circle, layout_random

# Random graphs
export rgraph, rgnm, rgnp

# Include source files
include("centrality/centrality.jl")
include("measures/measures.jl")
include("cohesion/cohesion.jl")
include("equivalence/equivalence.jl")
include("random/random.jl")
include("layout/layout.jl")

end # module
