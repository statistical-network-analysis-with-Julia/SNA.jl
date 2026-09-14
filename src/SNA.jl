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
import StatsAPI: coef, stderror, vcov, confint, loglikelihood, nobs, dof, aic, bic,
                 coeftable

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

# The user-facing network API is deliberate: implementation and fixture
# helpers in Networks must not leak into a session via `using SNA`.
export AbstractNetwork, Network, BipartiteNetwork, network, network_initialize
export nv, ne, vertices, edges, src, dst, has_vertex, has_edge, is_directed
export neighbors, inneighbors, outneighbors, degree, indegree, outdegree
export add_vertex!, add_vertices!, rem_vertex!, add_edge!, add_edges!, rem_edge!
export get_vertex_attribute, set_vertex_attribute!, delete_vertex_attribute!
export get_edge_attribute, set_edge_attribute!, delete_edge_attribute!
export get_network_attribute, set_network_attribute!, delete_network_attribute!
export list_vertex_attributes, list_edge_attributes, list_network_attributes
export vertex_attribute_vector, as_matrix, as_adjacency_matrix, as_edgelist, as_dataframe
export network_from_matrix, network_from_edgelist, network_from_dataframe
export read_pajek, write_pajek, write_graphml, write_edgelist_csv, load_dataset
export set_missing_dyad!, is_missing_dyad, delete_missing_dyad!, clear_missing_dyads!
export missing_dyads, n_missing_dyads, supports_missing, require_observed
export missing_policies, MISSING_POLICIES
export fit_metadata, estimand, objective, is_exact, se_method, missing_method, approximations
export coef, stderror, vcov, confint, loglikelihood, nobs, dof, aic, bic, coeftable

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
include("sociomatrix.jl")
include("centrality/centrality.jl")
include("measures/measures.jl")
include("cohesion/cohesion.jl")
include("equivalence/equivalence.jl")
include("random/random.jl")
include("layout/layout.jl")
include("qap/qap.jl")

# Accepted policies are a capability declaration, not a missing-data estimator.
for f in (degree_centrality, betweenness_centrality, closeness_centrality,
          eigenvector_centrality, bonacich_power, katz_centrality, pagerank,
          flowbet, centralization, density, gden, reciprocity, grecip,
          transitivity, gtrans, mutuality, dyad_census, triad_census,
          component_dist, hierarchy, efficiency, connectedness, components,
          cliques, kcores, cutpoints, bridges, bicomponents, largest_component,
          structural_equivalence, regular_equivalence, equiv_clust, blockmodel,
          geodesic_distance, reachability, diameter, average_path_length,
          qaptest, netlm, netlogit)
    @eval Networks.missing_policies(::typeof($f)) = (:error, :face)
end

end # module
