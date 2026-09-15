# Changelog

All notable changes to SNA.jl are documented in this file. The format is
based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the
package adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] - Unreleased

Release driven by the 2026-07 expert-panel review: three confirmed
correctness bugs (cliques, local clustering, undirected degree) are fixed,
SNA now extends the Graphs.jl generics instead of shadowing them, QAP
inference and Freeman centralization close the main R `sna` feature gaps, and
several measures are realigned to R `sna` semantics (breaking where the old
numbers were wrong or R-divergent).

### Breaking

- Logistic regression now rejects verified quasi-complete separation as well
  as complete separation. Candidate rays from Newton fitting are corrected
  and verified with exact rational arithmetic; tolerances alone cannot label
  near-overlapping observations as separated. If any requested QAP replicate
  fails, the error names its index and no permutation p-value is returned.
  Sparse binary predictors may therefore fail even when the observed fit is
  finite. No failed draws are omitted. These are certificate checks, not an
  exhaustive linear-programming search for all possible separating rays.

- Regular equivalence now averages matching in both directions so actor labels
  cannot change its scores, and warns on non-convergence. Valued regular
  equivalence is explicitly unsupported. Regression rejects redundant predictors
  and complete or certified quasi-complete logistic separation with actionable errors.

- Directed `components`, `largest_component`, and `component_dist` now default
  to strong connectivity; use `connected=:weak` (or legacy `mode=:weak`).
  `cutpoints` uses R's directed strong articulation by default; weak and
  recursive connectivity are explicit options. Cliques/bicomponents default to
  mutual arcs; request `symmetrize=:weak` for the former either-direction rule.
  Bicomponents exclude bridge components by default (`min_size=2` restores them).
- Closeness now uses R's zero score for actors unable to reach every other actor;
  `cmode=:component` restores the prior reachable-component scaling.
- Katz now uses `β` as the baseline in `(I - α Aᵀ)⁻¹ β𝟙`, returns raw scores,
  validates convergence of its walk series, and supports `normalized=true`.
- Degree and SNA density exclude self-loops by default (`diag=true` includes them).
  Undefined reciprocity ratios return `NaN`, matching R. Transitivity with
  no qualifying two-paths stays `1.0`, verified against R `sna` 2.8 source.
  Diameter now returns `Inf` if any ordered pair is unreachable, as documented.
- SNA exports the `Network` constructor through an explicit API list, removing
  accidental fixture, bootstrap, optimizer and presentation-helper re-exports.
- Seeded QAP draws change because each replicate now owns a seeded RNG;
  serial/threaded output is identical and two-mode permutations preserve modes.
  QAP regression results now retain iid-dyad covariance, and OLS `se_method`
  reports `:ols`. Shared stable Newton fitting replaces local logistic IRLS;
  non-converged fits warn and failed permutation fits throw.


- **Every exported measure now refuses a network with masked (unobserved)
  dyads** rather than silently computing from their face values. Measures take
  a `missing::Symbol=:error` keyword and call Networks.jl's `require_observed`;
  `missing=:face` is the explicit opt-in to the previous behaviour and returns
  exactly the old numbers. Covers the centrality, cohesion, measure,
  equivalence, and QAP/network-regression routines (`netlm`/`netlogit` guard
  the response *and* every predictor). *Migration:* pass `missing=:face`, or
  `clear_missing_dyads!(net)`, to analyse a masked network as recorded.
- **`degree_centrality` single-counts undirected edges** (Medici now scores
  the textbook 6, not 12) and normalization uses the correct `n−1` ceiling.
  *Migration:* multiply by 2 if you calibrated against the old doubled
  values.
- **`betweenness_centrality` default changed to `normalized=false`**,
  returning raw R `sna`-style scores (undirected halved). *Migration:* pass
  `normalized=true` for the old scaling.
- **`transitivity(type=:global)` now computes R `sna`'s weak transitivity**
  (directed triple enumeration), which differs numerically from the old
  `Graphs.global_clustering_coefficient` on directed networks. *Migration:*
  call `Graphs.global_clustering_coefficient(net.graph)` if you need the old
  quantity.
- **`reciprocity(method=:dyadic)` matches the R `sna` default** — null dyads
  now count as reciprocated, `(M+N)/(M+A+N)`. The old `M/(M+A)` value is
  available as `method=:dyadic_nonnull`. *Migration:* pass
  `method=:dyadic_nonnull` for the 0.1 number.
- **`mutuality` returns the integer count of mutual dyads** (R `sna`
  semantics; was a `Float64` proportion). *Migration:* use `reciprocity` for
  a proportion.
- **`hierarchy` defaults to `measure=:reciprocity`** (1 − dyadic
  reciprocity, the R `sna` default); the old Krackhardt behavior is
  `measure=:krackhardt`. *Migration:* pass `measure=:krackhardt`.
- **`bonacich_power` keywords renamed** `β` → `exponent`, `normalized` →
  `rescale` (sna `bonpow` semantics, Σc²=n scaling), and it returns `NaN`
  with a warning on a singular `(I−βA)` instead of zeros. *Migration:*
  rename the keywords; handle `NaN`.
- **`bicomponents` returns `Vector{Vector{Tuple{Int,Int}}}`** (edge lists
  per bicomponent) instead of Graphs' component output. *Migration:*
  consume edge tuples.
- **SNA methods now dispatch only on `::AbstractNetwork`.** Because SNA
  extends (rather than shadows) `Graphs.density`, `diameter`, `bridges`,
  and the centrality generics, calling them on a raw Graphs graph resolves
  to the Graphs implementation. *Migration:* pass `Network` objects to get
  SNA/R-`sna` semantics.
- **Minimum Julia raised to 1.12**; package UUID regenerated. *Migration:*
  upgrade Julia and re-resolve environments pinning the old UUID.

### Added

- Square sociomatrices support both two-mode network representations without
  dropping reverse-mode arcs. Density forwards `discount_bipartite` to Networks.
- Valued matrix measures accept `ignore_eval=false, attr=:weight`; eigenvector
  centrality uses a direct solve and warns when the solution is non-unique.
- Provenanced R `sna` fixtures cover full Florentine/Sampson vectors, directed
  cohesion, two-mode density/centrality, loops, valued ties, and regression fits.
- Full StatsAPI methods for QAP regressions and shared coefficient tables with
  permutation p-value display floor `1/reps`. Covariance/intervals are iid-dyad
  quantities; permutation p-values provide QAP inference.


- QAP inference: `qaptest`, `netlm`, `netlogit` with Dekker
  double-semi-partialing (`nullhyp=:qapspp`) as the default, plus
  `QAPTestResult`/`NetLMResult`/`NetLogitResult` result types.
- `centralization(net, measure; mode, normalized)` — Freeman graph
  centralization for degree/betweenness/closeness/eigenvector with R `sna`
  `tmaxdev` theoretical maxima.
- Random network generators implemented: `rgraph` (Bernoulli), `rgnp`,
  `rgnm` (the names existed in 0.1 exports but had no definitions).
- Layout functions implemented: `layout_circle`, `layout_random`,
  `layout_fruchterman_reingold`, `layout_kamada_kawai` (previously exported
  but undefined).
- `reciprocity(method=:dyadic_nonnull)` and `hierarchy(measure=...)`
  options.
- BenchmarkTools suite (`benchmark/`) exercising `triad_census`.

### Changed

- Documentation uses the default Documenter themes, with a new package-specific
  SVG icon and browser favicon in the official Julia logo colors.
- `density`/`gden` forward the caller's `missing=` policy to
  `Networks.network_density(net; missing=policy)`. Networks.jl's
  `network_density` now guards itself with `require_observed` (panel 2026-09,
  item 4), so the bare call `density(masked; missing=:face)` used to make would
  have been refused one layer down. Behaviour is unchanged: `:error` refuses,
  `:face` counts face values over every dyad.

- Composition fix: SNA extends the Graphs.jl generics (`density`,
  `diameter`, `bridges`, and centrality functions) via `import Graphs:`
  instead of defining same-named local functions — `using SNA, Graphs` no
  longer produces ambiguous bindings. `using SNA` also re-exports the
  Networks.jl public API.
- `eigenvector_centrality` defaults tightened (`max_iter` 100 → 1000, `tol`
  1e-6 → 1e-10) and returns the non-negative Perron orientation.

### Fixed

- **`cliques()` no longer throws on every input:** the backing digraph is
  converted to a `SimpleGraph` before `maximal_cliques` (directed networks
  symmetrized under the weak rule).
- **`transitivity(type=:local/:average)` correct on undirected networks:**
  the degree-doubled `k(k−1)` denominator inflated results ~4–6×; local
  clustering now runs on the undirected projection.
- `flowbet` is real Freeman flow betweenness via Edmonds–Karp max-flow
  (was a placeholder returning standard betweenness with a warning).
- `triad_census` implements the actual Davis–Leinhardt 16-class directed /
  4-class undirected census (the old placeholder mis-mapped classes).
- `efficiency` and `connectedness` use the correct per-weak-component
  Krackhardt formulas.
- `katz_centrality` and `pagerank` call Graphs.jl with the correct
  positional-argument API.
- `equiv_clust` performs genuine average-linkage (UPGMA) hierarchical
  clustering; `blockmodel` sizes blocks by realized labels;
  `regular_equivalence` documented as CATREGE-style.

### Performance

- `triad_census` uses the edge-driven Batagelj–Mrvar algorithm — cost
  scales with edges instead of O(n³) triples.
- `reachability` BFS uses an index-pointer queue (O(1) pops).

## [0.1.0] - 2026-02-09

Initial release: centrality, cohesion, measures, and equivalence functions
ported from R `sna`.

### Known limitations

Weighted shortest paths are not implemented. Structural equivalence uses full
incoming/outgoing profiles; regular equivalence is not the full R `sedist`
algorithm catalogue. Layouts/generators do not reproduce R's random stream.
Dense eigen and flow routines suit moderate networks. A non-unique leading
eigenspace produces a warning and one valid unit-norm eigenvector.
