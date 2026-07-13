# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

SNA.jl is a Julia port of the R `sna` package (StatNet collection). It provides social network analysis tools built on top of the `Network` package and `Graphs.jl`.

## Common Commands

```bash
# Install dependencies (must resolve Network first — see note below)
julia --project -e 'using Pkg; Pkg.instantiate()'

# Run tests
julia --project -e 'using Pkg; Pkg.test()'

# Load the package interactively
julia --project -e 'using SNA'

# Build documentation locally (install docs deps first)
julia --project=docs -e 'using Pkg; Pkg.instantiate()'
julia --project=docs docs/make.jl
```

There is no built-in test filter; to run a single `@testset`, comment out the others in `test/runtests.jl`.

## Architecture

The package is organized into submodules, each in its own directory under `src/`:

- **`centrality/centrality.jl`** — Vertex-level centrality measures (degree, betweenness, closeness, eigenvector, Bonacich power, Katz, PageRank, flow betweenness)
- **`measures/measures.jl`** — Network-level statistics (density, reciprocity, transitivity, dyad/triad census, hierarchy, efficiency, connectedness)
- **`cohesion/cohesion.jl`** — Substructure detection (components, cliques, k-cores, cutpoints, bridges, bicomponents)
- **`equivalence/equivalence.jl`** — Position analysis (structural equivalence, regular equivalence, blockmodeling, consensus clustering)
- **`qap/qap.jl`** — QAP inference (`qaptest`, `netlm`, `netlogit` with Dekker double-semi-partialing as the default null)
- **`random/random.jl`** — Random graph generators (`rgraph`, `rgnm`, `rgnp`)
- **`layout/layout.jl`** — Layout algorithms (Fruchterman-Reingold, Kamada-Kawai, circle, random)

All modules are included and all public functions exported from the top-level `src/SNA.jl` module. `centralization` (Freeman centralization with sna's `tmaxdev` maxima) lives in `centrality/centrality.jl`; `triad_census` uses the edge-driven Batagelj–Mrvar algorithm in `measures/measures.jl`.

## Key Dependencies

- **`Network`** — The core network data structure (`Network{T}`). All SNA functions accept `Network` objects and use its API (`nv`, `vertices`, `inneighbors`, `outneighbors`, `add_edge!`, `is_directed`, `network_density`, `as_matrix`, etc.). **This is an unregistered package** (UUID `027a387e-...`); it must be added via a local path or git URL before anything else will work. The `docs/Project.toml` `[sources]` paths may need updating to match the actual location of the Network package on the current machine.
- **`Graphs.jl`** — Graph algorithms accessed via `net.graph` (the underlying `SimpleDiGraph`/`SimpleGraph`)
- Hierarchical clustering for `equiv_clust`/`blockmodel` is implemented in-package (average-linkage UPGMA in `equivalence.jl`); there is no Clustering.jl dependency

## Missing-dyad policy (the ecosystem missing-data contract)

A masked dyad in a `Network` means the tie status is **unobserved**, not absent — the backing graph still stores a "face value" (edge present or absent) that every structural query reports. Reading that face value silently is how a partially observed network becomes a plausible, wrong number.

Every exported SNA measure therefore takes a `missing::Symbol=:error` keyword and calls Networks.jl's `require_observed(net, policy; context="<function name>")` before touching the graph:

- `:error` (default) — throws an `ArgumentError` naming the routine if the network has any masked dyad.
- `:face` — the explicit, auditable opt-in: masked dyads are used at their stored face values, reproducing exactly the pre-policy answer. Each docstring states what that means for the measure.

Rules for new code:

- **Every new exported measure must take `missing::Symbol=:error` and guard.** Undirected/directed, vertex-level or graph-level, descriptive or inferential — no exceptions; the inferential routines (`qaptest`, `netlm`, `netlogit`, `centralization`) are where it matters most, and `netlm`/`netlogit` guard `y` *and* every predictor (raw matrix arguments carry no mask and are taken as given).
- **Forward the policy to internal calls** (e.g. `centralization` → `degree_centrality`, `blockmodel` → `equiv_clust` → `structural_equivalence`), otherwise `:face` would hit an `:error` default one level down.
- The keyword is *named* `missing` (per the issue) but must be **rebound to a local `policy = missing`** at the top of the body — inside the body the name shadows `Base.missing`.
- `require_observed` / `supports_missing` / `MISSING_POLICIES` are exported by Networks.jl and re-exported by SNA, so they are callable unqualified. (Before the v0.2 `Networks` rename, qualified access was not merely a style choice — `Network.require_observed(...)` resolved to field access on the exported type and failed to precompile. That hazard is gone.)
- **Do not invent statistics.** No SNA measure implements a principled missing-data estimator (no listwise deletion, no density-over-observed-dyads, no imputation), so none declares `Networks.supports_missing`; the guard plus a documented `:face` treatment is the contract. Adding such an estimator means declaring the trait *and* justifying the statistics.
- Layouts (`layout_*`) and random-graph generators (`rgraph`, `rgnm`, `rgnp`) are exempt: layouts draw a picture rather than report a statistic, and the generators create fresh unmasked networks.

## Conventions

- Functions follow R `sna` naming where possible (e.g., `gden`, `grecip`, `gtrans` as aliases)
- Measures follow R `sna` conventions exactly where a counterpart exists; golden-master values from sna 2.8 (Sampson samplike, Padgett flomarriage) are baked into the tests
- Directed networks are the default; undirected networks are created with `network(n; directed=false)`
- Centrality functions return `Vector{Float64}`; equivalence functions return `Matrix{Float64}`
- Tests are in a single file `test/runtests.jl` using nested `@testset` blocks
- Documentation uses Documenter.jl with source pages in `docs/src/`
- Requires Julia >= 1.12 (Networks.jl cannot load on earlier versions)
