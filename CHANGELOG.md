# Changelog

All notable changes to SNA.jl are documented in this file. The format is
based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the
package adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] - Unreleased

The first public release. SNA.jl's measures now carry R `sna`'s names and
semantics, no Graphs.jl function is extended, every measure is checked
against sna 2.8, and QAP inference covers network regression with Dekker
double semi-partialing.

**Dependency renamed:** the foundation package is now `NetworkCore` (developed as `Networks`); write `using NetworkCore` where code said `using Networks`. Types and functions keep their names.

### Breaking

- **SNA.jl no longer adds methods to Graphs.jl functions.** Loading SNA used
  to change what `Graphs.degree_centrality`, `betweenness_centrality`,
  `closeness_centrality`, `eigenvector_centrality`, `katz_centrality`,
  `pagerank`, `density`, `diameter` and `bridges` returned for a `Network`.
  They now keep Graphs.jl's meaning and SNA no longer exports them. SNA's
  measures have sna's names:

  | Before | Now |
  |---|---|
  | `degree_centrality(net; mode=:total/:in/:out)` | `degreecent(net; cmode=:freeman/:indegree/:outdegree)` |
  | `betweenness_centrality(net)` | `betweenness(net)` |
  | `closeness_centrality(net; cmode=:sna)` | `closeness(net)`; `cmode=:component` is `Graphs.closeness_centrality(net)` |
  | `eigenvector_centrality(net)` | `evcent(net)` |
  | `katz_centrality`, `pagerank` | `Graphs.katz_centrality`, `Graphs.pagerank` (Graphs.jl's conventions) |
  | `density(net)` | `gden(net)` |
  | `diameter(net)` | `maximum(geodist(net).gdist)`, or `Graphs.diameter(net)` |
  | `bridges(net)` | `Graphs.bridges(net)` |
  | `bonacich_power(net)` | `bonpow(net)` |
  | `geodesic_distance(net)` | `geodist(net).gdist` |
  | `reciprocity(net; method=m)` | `grecip(net; measure=m)` |
  | `structural_equivalence(net; method=m)` | `sedist(net; method=m)` |
  | `consensus(clusterings)` | `consensus_clustering(clusterings)` |
  | `components(net)` | `component_dist(net).membership` (the old name collided with `Distributions.components`) |

  The old names were never part of a release and have no aliases; the
  table is also in the README and on the documentation's Utilities page.

- `component_dist` returns sna's `(membership, csize, cdist)`; the old
  sorted sizes are `sort(component_dist(net).csize; rev=true)`.
- `kcores` returns each vertex's core number, as sna's `kcores`; the members
  of the k-core are `findall(>=(k), kcores(net))`, and the old `k` keyword
  raises an `ArgumentError` saying so.
- `gtrans(net; measure)` and `grecip(net; measure)` are sna's functions with
  sna's keyword and measures (they were aliases of `transitivity` and
  `reciprocity`).
- `sedist` (formerly `structural_equivalence`) defaults to sna's Hamming
  distance, a count.
- `equiv_clust` and `blockmodel` default to sna's clustering: Hamming
  distance, complete linkage, ties broken as R's `hclust`, so
  `equiv_clust(net; k)` equals `cutree(equiv.clust(dat)$cluster, k)`. The
  previous clustering is `method=:correlation, cluster_method=:average`, and
  the old `method=:regular` is `equiv_fun=regular_equivalence`.
  A single-vertex position has a `NaN` diagonal block density, as in sna.
- QAP keywords: `reps` is now `n_sim`, and `regular_equivalence`'s
  `max_iter` is `maxiter`.
- `netlm`/`netlogit` with raw matrices (`mode=:auto`) treat symmetric data
  as undirected; sna's default `"digraph"` counted every pair twice and made
  classical standard errors √2 too small. `mode=:digraph` restores it.
- `netlogit`'s QAP statistic is the signed root likelihood ratio by default
  (`statistic=:wald` is sna's Wald z). `NetLogitResult` has a `statistic` field.
- A separated `netlogit` design no longer throws. It follows the ecosystem's
  separation policy, through NetworkCore.jl's verdict: a warning,
  `converged == false`, the separated coefficients in the new
  `fit.separated` field (the verdict in `fit.separation`), and `NaN` test
  statistics, p-values and confidence intervals, with no permutation run.
  With `statistic=:wald`, a separated permutation replicate no longer aborts
  the test: the p-values of the coefficients it touches are withheld
  (`NaN`, warned).
- Measures that are undefined return `NaN` as sna does: `hierarchy` with no
  dyad or no reachable pair, `efficiency` of a network of isolates,
  `bonpow` with no ties, `evcent` of a network without a directed cycle, and
  `centralization` when the theoretical maximum is 0.
- Every measure refuses a network with masked (unobserved) dyads by default
  (`missing=:error`); `missing=:face` uses the stored values explicitly.
- Undirected degree counts each edge once; betweenness is raw by default;
  closeness is 0 for an actor that cannot reach everyone; dyadic
  reciprocity counts null dyads as symmetric; `mutuality` is a count;
  `hierarchy` defaults to `measure=:reciprocity`; directed components and
  cutpoints use strong connectivity; cliques and bicomponents use mutual
  ties: all as in sna.
- Requires Julia 1.12.

### Added

- `cug_test` (sna's `cug.test`): conditional uniform graph tests
  conditioning on size, edges or the dyad census, with `n_sim` and `rng`.
- `brokerage` (sna's `brokerage`): Gould–Fernandez roles with sna's
  expected values, standard deviations and z-scores.
- A PrecompileTools workload: the first calls of the measures and QAP
  routines in a session take under 1 s instead of about 15 s.
- `infocent` (Stephenson–Zelen information centrality), `gcor` and `gcov`
  (graph correlation and covariance, diagonal excluded as in sna), and
  sna's other measures and modes: `closeness` `:suminvdir`, `:suminvundir`,
  `:gil_schmidt`; `betweenness` `:undirected`; `grecip` `:edgewise_lrr`,
  `:correlation`; `gtrans` `:strong`, `:weakcensus`, `:strongcensus`,
  `:correlation`; `component_dist` `:recursive`; `kcores` `:indegree` and
  `:outdegree`; `geodist` geodesic counts; `sedist` `:gamma` and `:exact`;
  `rescale` for the centrality indices.
- `centralization(net, f; kwargs...)` takes the measure function, as R's
  `centralization(dat, FUN, ...)`; `bonpow` has a theoretical maximum.
- QAP inference: `qaptest`, `netlm`, `netlogit` (Dekker double
  semi-partialing by default), with StatsAPI accessors (`coef`, `coefnames`,
  `stderror`, `vcov`, `confint`, `loglikelihood`, `nobs`, `dof`, `aic`, `bic`,
  `coeftable`) and the ecosystem's result-metadata protocol. `coefnames(fit)`
  is R's `names(coef(fit))`.
- Two-mode QAP: `qaptest` passes the incidence matrices to the statistic,
  `netlm`/`netlogit` use the cross-mode dyads only, and `p × q` incidence
  matrices are accepted as two-mode data.
- A "Coming from R sna" concordance in the documentation, and an R-generated
  fuzz fixture (60 random networks, every measure) and an inference fixture
  (`brokerage`, `equiv.clust`, `blockmodel`, `cug.test`) beside the Florentine,
  Sampson, two-mode, looped and valued fixtures.
- Random graphs (`rgraph`, `rgnp`, `rgnm`) and layouts (`layout_circle`,
  `layout_random`, `layout_fruchterman_reingold`, `layout_kamada_kawai`).

### Fixed

- `sedist`/structural equivalence left each actor's self cell and the ties
  between the two actors compared in the profile, so two equivalent actors
  tied to each other were not equivalent (correlation 0 instead of 1). They
  are now excluded, as in sna; `equiv_clust` and `blockmodel` inherit the fix.
- Directed local clustering (`transitivity(type=:local)`) used Graphs.jl's
  mixed definition (1/6 on the complete three-actor digraph); it is now
  Fagiolo's (2007) total coefficient (`cmode=:weak` ignores direction, as
  igraph). Loops are not neighbours, vertices of degree below 2 are `NaN`,
  and `:average` leaves them out.
- `bonpow` tested the determinant of `(I − βA)`, which underflows on large
  well-conditioned systems and returned `NaN`; it now tests the reciprocal
  condition number, as R's `solve(tol=)`.
- `efficiency` counted self-loops as arcs (and could go negative); `kcores`
  threw on any self-loop. Both ignore loops, as sna does.
- Two-mode `netlm`/`netlogit` regressed over the impossible within-mode
  pairs as observed zeros (mean slope 0.25 between independent networks).
- `netlogit` QAP aborted whenever any permutation was separated, so the
  canonical marriage ~ business QAP could not run; the likelihood-ratio
  statistic exists on separated permutations.
- `evcent` of an acyclic digraph returned an arbitrary vector; it is `NaN`.
- `consensus_clustering` (formerly `consensus`) depended on vertex order and
  was not transitive; it now takes connected groups of the co-membership graph.
- `layout_kamada_kawai` returned classical MDS; it now minimizes the
  Kamada–Kawai stress.
- Errors inside threaded QAP replicates no longer share a boxed variable.
- `closeness`, `flowbet`, `centralization(net, closeness)`, `geodist`,
  `average_path_length` and `layout_kamada_kawai` threw a `MethodError` on a
  `Network{Int32}`; every measure now accepts any integer vertex type.
- Flow betweenness, the triad census (Batagelj–Mrvar), cliques, undirected
  clustering and Katz centrality were wrong or placeholders before 0.2.0.

### Known limitations

- Not implemented: `lnam`, `nacf`, `netcancor`,
  `gscor`, `gscov`, `structdist`, `hdist`, `bbnam`, `bn`, sna's informant
  `consensus`, `prestige`, `stresscent`, `graphcent`, `loadcent`,
  `gilschmidt` (use `closeness(net; cmode=:gil_schmidt)`), `lubness`, path
  and cycle censuses, `simmelian`, `ego.extract`, the `rguman`, `rgws`,
  `rgbn` and `rgnmix` generators, and plotting.
- sna's endpoint/proximal/length-scaled betweenness, `gtrans`'s `rank`
  measure and `unilateral` components raise an `ArgumentError`.
- Weighted shortest paths and weighted regular equivalence are not
  implemented.
- `equiv_clust` returns the cut at `k` clusters, not the tree;
  `regular_equivalence` is an iterative similarity, not REGE.
- `cug_test` does not condition on self-loops or edge values and refuses
  two-mode networks; `brokerage` ignores self-loops (sna counts them in the
  density).
- Regression covariance and intervals assume independent dyads; inference is
  in the QAP p-values, which are proportions without a `+1` correction (as in
  sna) and can be 0.
- Dense spectral and flow routines suit networks of moderate size.

## [0.1.0] - 2026-02-09

Initial development version (not publicly released): centrality, cohesion,
measures, and equivalence functions ported from R `sna`.
