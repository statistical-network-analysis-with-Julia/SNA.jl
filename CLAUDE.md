# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

SNA.jl is a Julia port of the R `sna` package (StatNet collection). It provides social network analysis tools built on top of the `NetworkCore` package and `Graphs.jl`.

## Common Commands

```bash
# Install dependencies (must resolve NetworkCore first — see note below)
julia --project -e 'using Pkg; Pkg.instantiate()'

# Run tests
julia --project -e 'using Pkg; Pkg.test()'

# Load the package interactively
julia --project -e 'using SNA'

# Build documentation locally (install docs deps first)
julia --project=docs -e 'using Pkg; Pkg.instantiate()'
julia --project=docs docs/make.jl
```

There is no built-in test filter; to run a single `@testset`, comment out the others in `test/runtests.jl`. Run once with `JULIA_NUM_THREADS=4` too: the QAP testsets assert thread-count independence.

Golden fixtures (need R with sna 2.8, network, ergm and igraph):

```bash
Rscript test/fixtures/r/sna_reference.R   # Florentine, Sampson, two-mode, loops, weights, regressions
Rscript test/fixtures/r/sna_fuzz.R        # 60 seeded random networks x every measure, plus edge cases
Rscript test/fixtures/r/sna_inference.R   # brokerage, equiv.clust + cutree, blockmodel, cug.test (20 000 reps)
```

## Architecture

The package is organized into submodules, each in its own directory under `src/`:

- **`centrality/centrality.jl`** — `degreecent` (sna's `degree`), `betweenness`, `closeness`, `evcent`, `bonpow`, `infocent`, `flowbet`, `centralization(net, f; kwargs...)` (sna's `tmaxdev` maxima via `_tmaxdev(f, net; kwargs...)`)
- **`measures/measures.jl`** — `gden`, `grecip`, `gtrans`, `transitivity` (local clustering: Fagiolo total coefficient on directed networks), `dyad_census`, `triad_census` (edge-driven Batagelj–Mrvar), `hierarchy`, `efficiency`, `connectedness`, `mutuality`, `reachability`
- **`cohesion/cohesion.jl`** — `component_dist` (sna's `(membership, csize, cdist)`), `largest_component`, `cliques`, `kcores` (core numbers, generalized peeling for in/out/Freeman degree), `cutpoints`, `bicomponents`, `geodist` (distances and geodesic counts), `average_path_length`
- **`measures/brokerage.jl`** — `brokerage` (Gould–Fernandez roles; sna's expectation/variance formulas ported term by term)
- **`equivalence/equivalence.jl`** — `sedist`, `regular_equivalence`, `equiv_clust` (sna defaults: Hamming, complete linkage; `_hclust` reproduces R `hclust`'s tie-breaking — first minimum in (i, j) order, merged cluster keeps the lower index — and `cutree`'s numbering, pinned against R on 24 graphs), `blockmodel`, `consensus_clustering`
- **`qap/cug.jl`** — `cug_test`/`CUGTestResult` (conditioning on size, edges, dyad census; reuses `_qap_replicates`, so per-replicate seeding and thread independence are shared with QAP; Monte Carlo summaries are compared with a 20 000-replicate R run within 5 MC standard errors)
- **`qap/qap.jl`** — `gcor`, `gcov`, `qaptest`, `netlm`, `netlogit` (Dekker double-semi-partialing default; netlogit's QAP statistic is the signed-root likelihood ratio by default). `NetLMResult`/`NetLogitResult` answer the eleven StatsAPI verbs, `coefnames` (a copy of `fit.names`) included; the testset pins them with `check_statsapi(...; required=(STATSAPI_VERBS..., :coefnames), strict=true)`
- **`random/random.jl`** — `rgraph`, `rgnm`, `rgnp`
- **`layout/layout.jl`** — `layout_fruchterman_reingold`, `layout_kamada_kawai` (stress majorization from a classical-MDS start), `layout_circle`, `layout_random`
- The `@compile_workload` at the end of `src/SNA.jl` runs the measures, `cug_test`, `brokerage` and the three QAP routines on 8-actor networks under a `NullLogger` (first calls: about 15 s → under 1 s). Add new entry points to it
- There are **no deprecation shims** (owner decision 2026-10): the development-time names (`bonacich_power`, `geodesic_distance`, `reciprocity`, `structural_equivalence`, `consensus`, `components`, `equiv_clust(method=:structural/:regular)`, the `reps`/`max_iter` keywords) were never released and were deleted; the rename table is in the README, the CHANGELOG and `docs/src/api/utilities.md`. Tests assert the old names are undefined and the old keywords are `MethodError`s. Do not add aliases for unreleased names
- Every internal helper that takes a vertex count or vertex ids accepts any `Integer` and converts to `Int` (`_geodist(g, n::Integer)`, `_maxflow`, `_tricode`, `_union_neighborhood`): `nv` of a `Network{Int32}` is an `Int32`. The fuzz-corpus testset runs every measure on `Int` and `Int32` networks

## Naming: no type piracy, sna's names

**SNA adds no method to any Graphs.jl function.** Before 0.2.0 it extended `Graphs.degree_centrality`, `betweenness_centrality`, `closeness_centrality`, `eigenvector_centrality`, `katz_centrality`, `pagerank`, `density`, `diameter` and `bridges` for `AbstractNetwork` with R semantics (type piracy: loading SNA, even through TSNA, changed what Graphs' own calls returned). Those functions now keep Graphs' meaning on a `Network`, SNA does not export them, and SNA's measures carry R sna's names (dots → underscores). The rule for new functions: mirror the sna name; a name sna does not have (local clustering `transitivity`, `cliques`, `largest_component`, `average_path_length`, `consensus_clustering`) must not mean something else in R; and no export may clash with Graphs.jl, Distributions.jl, NetworkCore.jl or any ecosystem package (the "No type piracy" testset co-loads Graphs, Distributions, NetworkCore and SNA in a subprocess and asserts every exported name stays defined; Aqua's piracy check runs in the suite). `degree` is Graphs' vertex degree, hence `degreecent` (sna's own naming pattern: `infocent`, `stresscent`). The old→new table is in the README and CHANGELOG.

## Key Dependencies

- **`NetworkCore`** — The core network data structure (`Network{T}`). All SNA functions accept `Network` objects and use its API (`nv`, `vertices`, `inneighbors`, `outneighbors`, `add_edge!`, `is_directed`, `network_density`, `as_matrix`, etc.). **This is an unregistered package** (UUID `027a387e-...`); it must be added via a local path or git URL before anything else will work. The `docs/Project.toml` `[sources]` paths may need updating to match the actual location of the NetworkCore package on the current machine.
- **`Graphs.jl`** — Graph algorithms, called on `_graph(net)` (the stored digraph, through which an undirected tie is two arcs; correct the counts accordingly, as `betweenness` halves) or on symmetrized `SimpleGraph`s. Never extend a Graphs.jl function
- Hierarchical clustering for `equiv_clust`/`blockmodel` is implemented in-package (average-linkage UPGMA in `equivalence.jl`); there is no Clustering.jl dependency

## Missing-dyad policy (the ecosystem missing-data contract)

A masked dyad in a `Network` means the tie status is **unobserved**, not absent — the backing graph still stores a "face value" (edge present or absent) that every structural query reports. Reading that face value silently is how a partially observed network becomes a plausible, wrong number.

Every exported SNA measure therefore takes a `missing::Symbol=:error` keyword and calls NetworkCore.jl's `require_observed(net, policy; context="<function name>")` before touching the graph:

- `:error` (default) — throws an `ArgumentError` naming the routine if the network has any masked dyad.
- `:face` — the explicit, auditable opt-in: masked dyads are used at their stored face values, reproducing exactly the pre-policy answer. Each docstring states what that means for the measure.

Rules for new code:

- **Every new exported measure must take `missing::Symbol=:error` and guard.** Undirected/directed, vertex-level or graph-level, descriptive or inferential — no exceptions; the inferential routines (`qaptest`, `netlm`, `netlogit`, `centralization`) are where it matters most, and `netlm`/`netlogit` guard `y` *and* every predictor (raw matrix arguments carry no mask and are taken as given).
- **Forward the policy to internal calls** (e.g. `centralization` → `degreecent`, `blockmodel` → `equiv_clust` → `sedist`), otherwise `:face` would hit an `:error` default one level down.
- The keyword is *named* `missing` (per the issue) but must be **rebound to a local `policy = missing`** at the top of the body — inside the body the name shadows `Base.missing`.
- `require_observed` / `supports_missing` / `MISSING_POLICIES` are exported by NetworkCore.jl and re-exported by SNA, so they are callable unqualified. (Before the v0.2 `NetworkCore` rename, qualified access was not merely a style choice — `Network.require_observed(...)` resolved to field access on the exported type and failed to precompile. That hazard is gone.)
- **Do not invent statistics.** No SNA measure implements a principled missing-data estimator (no listwise deletion, no density-over-observed-dyads, no imputation), so none declares `NetworkCore.supports_missing`; the guard plus a documented `:face` treatment is the contract. Adding such an estimator means declaring the trait *and* justifying the statistics.
- Layouts (`layout_*`) and random-graph generators (`rgraph`, `rgnm`, `rgnp`) are exempt: layouts draw a picture rather than report a statistic, and the generators create fresh unmasked networks.

## Conventions

- Functions carry R `sna`'s names and semantics where a counterpart exists (see "Naming" above and `docs/src/r_concordance.md`, which lists every sna function and every deliberate difference)
- Full-vector sna 2.8 fixtures (Sampson, Florentine, directed, two-mode, looped and weighted examples, regressions) come from `test/fixtures/r/sna_reference.R`; `test/fixtures/r/sna_fuzz.R` freezes 60 seeded random networks (directed/undirected, loops, isolates) with R's value of every measure (and igraph's local clustering), plus the structural-equivalence, Bonacich and two-mode-regression edge cases. The "Fuzz corpus" testset compares every key (NaN must match NaN). Both are loaded with `NetworkCore.load_golden`
- Undefined ratios are `NaN` as in sna (hierarchy with no reachable pair, efficiency of isolates, bonpow with no ties, evcent of an acyclic digraph, centralization with a zero maximum); `gtrans` with no two-path is 1 (sna). `triad_census` with n < 3 is zeros (sna: NaN), a documented difference
- Keyword vocabulary: `maxiter`, `n_sim`, `rng` (pinned by a `Base.kwarg_decl` testset, which also asserts `max_iter`/`reps` are gone)
- Every export SNA owns has a docstring with a runnable ```` ```julia ```` `# Example`; the "Every export is documented with a runnable example" testset executes each block in a fresh module
- Directed networks are the default; undirected networks are created with `network(n; directed=false)`
- Centrality functions return `Vector{Float64}`; equivalence functions return `Matrix{Float64}`
- Tests are in a single file `test/runtests.jl` using nested `@testset` blocks
- Documentation uses Documenter.jl with source pages in `docs/src/`
- Requires Julia >= 1.12 (NetworkCore.jl cannot load on earlier versions)

## September semantics

- `_sociomatrix` always expands two-mode inputs to square adjacency; `_graph`
  unwraps both `Network` and `BipartiteNetwork`. Never call `as_matrix` directly
  from a matrix-based measure, because its default is a rectangular incidence matrix.
- Components use `connected=:strong`; cutpoints use directed strong articulation
  by default (R's `connected="recursive"` is mutual-arc symmetrization).
  Cliques/bicomponents default to `symmetrize=:strong`, and bicomponents exclude
  bridges unless `min_size=2`. `closeness` takes sna's `cmode` values.
- `sedist` excludes profile positions {i, j, n+i, n+j} for each pair (sna's
  `diag.remove` + pairwise deletion; Wasserman & Faust eq. 9.4); `diag=true`
  keeps them. Never put the self cells or the i–j ties back into the profile.
- Binary ties and `diag=false` are the measure defaults. `ignore_eval=false`
  requests the named edge attribute in matrix measures. Weighted shortest paths
  are not implemented; equivalence routines are not a complete R `sedist` port.
- QAP RNGs are seeded per replicate before threading, and permutations preserve
  two-mode actor partitions. Two-mode QAP uses cross-mode dyads only
  (`_dyad_indices(n, directed, partition; reverse)`; the mode-2 → mode-1 block of a
  directed two-mode design only if some argument has a tie there), `qaptest`
  hands the statistic the p × q incidence matrices, and rectangular matrices are
  two-mode data. Raw square matrices are undirected under `mode=:auto` when y and
  every predictor are symmetric. Regression covariance/intervals assume iid dyads;
  QAP p-values are a separate inference layer. The shared NetworkCore optimizer,
  likelihood derivatives and coefficient table provide numerical/presentation contracts.
- The re-export list is curated. In particular `Network` is exported and fixture,
  optimizer and bootstrap implementation helpers are not.

## Logistic identification (the ecosystem separation policy)

Numerical Newton convergence does not prove a finite logistic MLE. `_logit_fit`
asks NetworkCore.jl's shared verdict, `NetworkCore.logistic_separation(X, y)`
(a linear programme over the design, the separating direction certified in
exact arithmetic; never a coefficient cutoff), and returns it as
`fit.verdict`; `converged` is false on a separated design, and the
optimizer's singular-information warning is silenced there. `netlogit`
applies the ecosystem policy to the OBSERVED design: `warn_separation`, a
result with `converged == false`, `separated` (term names) and `separation`
(the verdict) fields, `NaN` `tstat`/p-values/`confint`, no permutation run
(`dist === nothing`, `reps == 0`), and `separation_caveat` in `show` and
`approximations`. It never throws for separation. Permutation replicates
need not be identified: the default QAP statistic (`statistic=:lr`) is the
signed root likelihood ratio, `sign(β̂ₖ)√(D₋ₖ − D)`, built from
`_logit_deviance` (the shared `newton_fit` run to the deviance infimum,
logging silenced), which exists on separated permutations, so no replicate
is dropped. With `statistic=:wald` a separated replicate's z is `NaN`, and
the p-values of every coefficient with a `NaN` replicate are withheld
(`NaN`, one warning, an `approximations` entry); a replicate that fails for
any other reason still aborts the test through `_qap_replicates`.
