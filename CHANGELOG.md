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
