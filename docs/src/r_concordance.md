# Coming from R sna

SNA.jl ports the descriptive measures and QAP inference of R's `sna`
package (Butts 2008). A function that ports an sna function has sna's name
(with `.` written as `_`) and sna's default semantics; deliberate differences
are listed below. Agreement with sna 2.8 is pinned by golden fixtures generated
by R (`test/fixtures/r/*.R`), including a corpus of 60 random networks with
self-loops, isolates and disconnected parts on which every measure is
compared with R, and 24 more for `brokerage`, `equiv_clust` and `blockmodel`.

## Functions

| R `sna` | SNA.jl | Notes |
|---|---|---|
| `degree(dat, cmode=)` | `degreecent(net; cmode=)` | `cmode` is `:freeman`, `:indegree` or `:outdegree`. Named `degreecent` because `degree` is Graphs.jl's vertex degree. |
| `betweenness(dat, cmode=)` | `betweenness(net; cmode=)` | `:directed` and `:undirected`; sna's endpoint/proximal/length-scaled variants are not implemented. |
| `closeness(dat, cmode=)` | `closeness(net; cmode=)` | `:directed`, `:undirected`, `:suminvdir`, `:suminvundir`, `:gil_schmidt`. |
| `evcent(dat)` | `evcent(net)` | Direct eigen solve: sna's `use.eigen=TRUE`, oriented non-negative. |
| `bonpow(dat, exponent=)` | `bonpow(net; exponent=)` | Same singularity test (`rcond < tol`); a singular system gives `NaN` with a warning, where R stops. |
| `infocent(dat, cmode=)` | `infocent(net; cmode=)` | `:weak` and `:strong` symmetrization. |
| `flowbet(dat)` | `flowbet(net)` | Binary ties by default (sna uses values). |
| `centralization(dat, FUN, ...)` | `centralization(net, f; kwargs...)` | `f` is `degreecent`, `betweenness`, `closeness`, `evcent` or `bonpow`. |
| `gden(dat)` | `gden(net)` | |
| `grecip(dat, measure=)` | `grecip(net; measure=)` | All five measures; `"dyadic.nonnull"` is `:dyadic_nonnull`. |
| `gtrans(dat, measure=)` | `gtrans(net; measure=)` | `:weak`, `:strong`, `:weakcensus`, `:strongcensus`, `:correlation`; `rank` is not implemented. |
| `dyad.census(dat)` | `dyad_census(net)` | A named tuple. |
| `triad.census(dat)` | `triad_census(net)` | Zeros (not `NaN`) for fewer than three actors. |
| `mutuality(dat)` | `mutuality(net)` | |
| `hierarchy(dat, measure=)` | `hierarchy(net; measure=)` | |
| `efficiency(dat)` | `efficiency(net)` | |
| `connectedness(dat)` | `connectedness(net)` | |
| `component.dist(dat, connected=)` | `component_dist(net; connected=)` | `(membership, csize, cdist)`; `unilateral` is not implemented. |
| `components(dat)` | `length(component_dist(net).csize)` | sna's `components` counts components. |
| `component.largest(dat)` | `largest_component(net)` | Vertex ids rather than a logical vector. |
| `geodist(dat)` | `geodist(net)` | `(counts, gdist)`. |
| `reachability(dat)` | `reachability(net)` | |
| `kcores(dat, cmode=)` | `kcores(net; cmode=)` | Core numbers. |
| `cutpoints(dat, connected=)` | `cutpoints(net; connected=)` | |
| `bicomponent.dist(dat)` | `bicomponents(net)` | Edge lists of the bicomponents. |
| `clique.census(dat)` | `cliques(net)` | Lists the maximal cliques instead of counting them. |
| `sedist(dat, method=)` | `sedist(net; method=)` | `:hamming`, `:correlation`, `:euclidean`, `:gamma`, `:exact`. |
| `cutree(equiv.clust(dat)$cluster, k)` | `equiv_clust(net; k)` | Same defaults (Hamming, complete linkage); returns the labels of the cut rather than the tree. `cluster.method` is `cluster_method`. |
| `blockmodel(dat, ec, k=)` | `blockmodel(net; k)` | Clusters internally with `equiv_clust`; densities only; `membership` is in vertex order. |
| `brokerage(g, cl)` | `brokerage(net, cl)` | Same components (`raw_nli`, `z_nli`, `raw_gli`, …); self-loops ignored. |
| `cug.test(dat, FUN, cmode=, reps=)` | `cug_test(net, f; cmode=, n_sim=)` | `:size`, `:edges`, `:dyad_census`; `diag=TRUE` and `ignore.eval=FALSE` are not implemented. |
| `gcor(g1, g2)`, `gcov(g1, g2)` | `gcor(g1, g2)`, `gcov(g1, g2)` | |
| `qaptest(dat, FUN, reps=)` | `qaptest(f, g1, g2; n_sim=)` | Two-mode data: `f` gets the incidence matrices. |
| `netlm(y, x, nullhyp=, reps=)` | `netlm(y, xs; nullhyp=, n_sim=)` | Two-mode data: cross-mode dyads only. |
| `netlogit(y, x, nullhyp=, reps=)` | `netlogit(y, xs; nullhyp=, n_sim=)` | QAP statistic: signed root likelihood ratio by default (`statistic=:wald` is sna's). |
| `rgraph(n, m, tprob)` | `rgraph(n; m, tprob)` | |
| `rgnm(n, nv, m)` | `rgnm(nv, m)` | One network. |
| `gplot.layout.fruchtermanreingold` | `layout_fruchterman_reingold(net)` | |
| `gplot.layout.kamadakawai` | `layout_kamada_kawai(net)` | Stress majorization from a classical MDS start. |
| `gplot.layout.circle`, `.random` | `layout_circle(net)`, `layout_random(net)` | |

## Deliberate differences

| Quantity | SNA.jl | sna |
|---|---|---|
| Edge values | Binary ties by default (`ignore_eval=true`) in every measure | Several functions use edge values by default |
| Raw matrix in `netlm`/`netlogit` | `mode=:auto` treats symmetric data as undirected | `mode="digraph"`: each pair counted twice |
| `netlogit` QAP statistic | Signed root likelihood ratio, defined on separated permutations | Wald z, meaningless on separated permutations |
| Two-mode QAP | Cross-mode dyads only | Not supported |
| `evcent` | Direct eigen solve | Power iteration (fails on periodic networks) |
| `bonpow` on a singular system | `NaN` and a warning | Error |
| `equiv_clust` | Returns the cut at `k` clusters | Returns the tree |
| `brokerage` density | Self-loops ignored | Self-loops counted |
| `regular_equivalence` | Iterative similarity | No sna counterpart (not REGE) |
| `gden(diag=true)`, undirected | `n(n+1)/2` possible dyads | `n²` matrix cells |
| `triad_census` with `n < 3` | Zeros | `NaN` |
| Monte Carlo p-values | Proportions (can be 0), as in sna | Same |

## Not implemented

These sna functions have no SNA.jl counterpart yet: network
autocorrelation models (`lnam`, `nacf`), canonical correlation
(`netcancor`), structural correlation and covariance (`gscor`, `gscov`,
`structdist`, `hdist`), Bayesian network inference (`bbnam`, `bn`),
informant consensus networks (`consensus`), the prestige, stress, graph, load and Gil–Schmidt indices as
separate functions (`prestige`, `stresscent`, `graphcent`, `loadcent`,
`gilschmidt`; Gil–Schmidt is `closeness(net; cmode=:gil_schmidt)`),
`lubness`, path and cycle censuses (`kpath.census`, `kcycle.census`),
`simmelian`, `ego.extract`, `neighborhood`, the random-graph generators
`rguman`, `rgws`, `rgbn` and `rgnmix`, and plotting (`gplot`). `cug_test`
does not condition on self-loops or edge values. Weighted
shortest paths are not implemented in any measure.

## Graphs.jl functions

SNA.jl adds no method to any Graphs.jl function. Graphs.jl's
`degree_centrality`, `betweenness_centrality`, `closeness_centrality`,
`eigenvector_centrality`, `katz_centrality`, `pagerank`, `density`,
`diameter`, `bridges` and `core_number` work on a `Network` with Graphs.jl's
own conventions (call them on the `Network`, never on its `graph` field,
which stores an undirected tie as two arcs). Before 0.2.0 SNA.jl added
methods to the first nine with R semantics; the [changelog](https://github.com/statistical-network-analysis-with-Julia/SNA.jl/blob/main/CHANGELOG.md)
lists every renamed function.
