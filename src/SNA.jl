"""
    SNA.jl - Social Network Analysis for Julia

A Julia package providing descriptive analysis tools for social networks,
including centrality measures, structural equivalence, cohesion analysis,
and network visualization layouts.

Port of the R sna package from the StatNet collection.
"""
module SNA

using Distributions
using Graphs
using LinearAlgebra
using Random
using Statistics
using Networks

# Extend (rather than shadow) the Graphs.jl generics that share a name with
# SNA functions, so that `using SNA, Graphs` never produces ambiguous
# bindings. The methods below are defined on Network's AbstractNetwork types.
import Graphs: density, diameter, bridges,
    degree_centrality, betweenness_centrality, closeness_centrality,
    eigenvector_centrality, katz_centrality, pagerank

# The shared result-metadata protocol (Networks.jl `src/results.jl`): the
# generic accessors that say what a fit actually did. Imported by name because
# SNA adds methods for `NetLMResult`/`NetLogitResult`; `fit_metadata(fit)`
# collects them.
import Networks: estimand, objective, is_exact, se_method, missing_method,
                 approximations

# Re-export the Networks.jl public API so that `using SNA` alone provides the
# network constructors and accessors, mirroring R's library(sna) working with
# network objects out of the box. The `Network` name itself is skipped: inside
# this module it is bound to the struct, and exporting it would collide with
# the package module binding in downstream namespaces (a plain @reexport fails
# for the same reason).
for _network_export in names(parentmodule(Network))
    _network_export === :Network && continue
    Core.eval(@__MODULE__, Expr(:export, _network_export))
end

# Centrality measures
export degree_centrality, betweenness_centrality, closeness_centrality
export eigenvector_centrality, bonacich_power, katz_centrality
export pagerank, flowbet, centralization

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

# QAP inference and network regression
export qaptest, netlm, netlogit
export QAPTestResult, NetLMResult, NetLogitResult

# Random graphs
export rgraph, rgnm, rgnp

# Include source files
include("centrality/centrality.jl")
include("measures/measures.jl")
include("cohesion/cohesion.jl")
include("equivalence/equivalence.jl")
include("random/random.jl")
include("layout/layout.jl")
include("qap/qap.jl")

end # module
