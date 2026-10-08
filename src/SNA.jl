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
using NetworkCore
using PrecompileTools: @setup_workload, @compile_workload
import StatsAPI: coef, stderror, vcov, confint, loglikelihood, nobs, dof, aic, bic,
                 coeftable, coefnames

# SNA adds no method to any Graphs.jl function. Graphs.jl's
# `degree_centrality`, `betweenness_centrality`, `density`, `diameter`, …
# keep Graphs.jl's own semantics on a `Network` whether or not SNA is loaded;
# SNA's R-semantics measures carry R sna's names (`degreecent`,
# `betweenness`, `closeness`, `evcent`, `gden`, `geodist`, …), none of which
# Graphs.jl, Distributions.jl or NetworkCore.jl exports.

# The shared result-metadata protocol (NetworkCore.jl `src/results.jl`): the
# generic accessors that say what a fit actually did. Imported by name because
# SNA adds methods for `NetLMResult`/`NetLogitResult`; `fit_metadata(fit)`
# collects them.
import NetworkCore: estimand, objective, is_exact, se_method, missing_method,
                 approximations

# The user-facing network API is deliberate: implementation and fixture
# helpers in NetworkCore must not leak into a session via `using SNA`.
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
export coef, stderror, vcov, confint, loglikelihood, nobs, dof, aic, bic, coeftable,
       coefnames

# Centrality (R sna: degree, betweenness, closeness, evcent, bonpow,
# infocent, flowbet, centralization)
export degreecent, betweenness, closeness, evcent, bonpow, infocent
export flowbet, centralization

# Graph-level indices and censuses
export gden, grecip, gtrans, transitivity, mutuality
export dyad_census, triad_census
export hierarchy, efficiency, connectedness, brokerage

# Cohesion and paths
export component_dist, largest_component, cliques, kcores, cutpoints
export bicomponents, geodist, reachability, average_path_length

# Structural equivalence and positions
export sedist, regular_equivalence, equiv_clust, blockmodel, consensus_clustering

# Visualization layouts
export layout_fruchterman_reingold, layout_kamada_kawai
export layout_circle, layout_random

# Graph correlation, QAP inference and network regression
export gcor, gcov, qaptest, netlm, netlogit
export QAPTestResult, NetLMResult, NetLogitResult
export cug_test, CUGTestResult

# Random graphs
export rgraph, rgnm, rgnp

# Include source files
include("sociomatrix.jl")
include("centrality/centrality.jl")
include("measures/measures.jl")
include("measures/brokerage.jl")
include("cohesion/cohesion.jl")
include("equivalence/equivalence.jl")
include("random/random.jl")
include("layout/layout.jl")
include("qap/qap.jl")
include("qap/cug.jl")

# Accepted policies are a capability declaration, not a missing-data estimator.
for f in (degreecent, betweenness, closeness, evcent, bonpow, infocent,
          flowbet, centralization, gden, grecip, gtrans, transitivity,
          mutuality, dyad_census, triad_census, hierarchy, efficiency,
          connectedness, component_dist, largest_component, cliques, kcores,
          cutpoints, bicomponents, geodist, reachability, average_path_length,
          sedist, regular_equivalence, equiv_clust, blockmodel,
          gcor, gcov, qaptest, netlm, netlogit, brokerage, cug_test)
    @eval NetworkCore.missing_policies(::typeof($f)) = (:error, :face)
end

# Precompile the main entry points on small networks, so the first call in
# a session does not pay for compilation (first use of the measures and the
# three QAP routines: about 15 s before, under 1 s after).
@setup_workload begin
    edges_d = [(1, 2), (2, 1), (2, 3), (3, 4), (4, 5), (5, 3), (6, 7), (7, 8), (8, 6), (1, 6)]
    edges_u = [(1, 2), (2, 3), (3, 4), (4, 5), (5, 3), (6, 7), (7, 8), (8, 6), (1, 6), (2, 7)]
    @compile_workload begin
        with_logger(NullLogger()) do
            for directed in (true, false)
                net = network(8; directed)
                add_edges!(net, directed ? edges_d : edges_u)
                degreecent(net); betweenness(net); closeness(net); evcent(net)
                bonpow(net; exponent=0.1); infocent(net); flowbet(net)
                centralization(net, degreecent); centralization(net, betweenness)
                gden(net); grecip(net); gtrans(net); transitivity(net; type=:local)
                dyad_census(net); triad_census(net); hierarchy(net; measure=:krackhardt)
                efficiency(net); connectedness(net); mutuality(net)
                component_dist(net); largest_component(net); cliques(net); kcores(net)
                cutpoints(net); bicomponents(net); geodist(net); reachability(net)
                average_path_length(net); sedist(net); regular_equivalence(net)
                equiv_clust(net; k=3); blockmodel(net; k=3)
                brokerage(net, [1, 1, 1, 2, 2, 2, 3, 3])
                cug_test(net, gtrans; cmode=:edges, n_sim=4, rng=Xoshiro(1), threaded=false)
                layout_kamada_kawai(net); layout_circle(net)
                other = rgnp(8, 0.4; directed, rng=Xoshiro(2))
                x = Float64[abs(i - j) for i in 1:8, j in 1:8]
                gcor(net, other)
                sprint(show, qaptest(gcor, net, other; n_sim=4, rng=Xoshiro(3), threaded=false))
                sprint(show, netlm(net, [other, x]; n_sim=4, rng=Xoshiro(4), threaded=false))
                sprint(show, netlm(net, other; nullhyp=:classical))
                sprint(show, netlogit(net, x; n_sim=4, rng=Xoshiro(5), threaded=false))
                coeftable(netlogit(net, x; nullhyp=:classical))
            end
        end
    end
end

end # module
