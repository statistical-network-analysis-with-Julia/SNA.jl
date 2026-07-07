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

The package is organized into four submodules, each in its own directory under `src/`:

- **`centrality/centrality.jl`** — Vertex-level centrality measures (degree, betweenness, closeness, eigenvector, Bonacich power, Katz, PageRank, flow betweenness)
- **`measures/measures.jl`** — Network-level statistics (density, reciprocity, transitivity, dyad/triad census, hierarchy, efficiency, connectedness)
- **`cohesion/cohesion.jl`** — Substructure detection (components, cliques, k-cores, cutpoints, bridges, bicomponents)
- **`equivalence/equivalence.jl`** — Position analysis (structural equivalence, regular equivalence, blockmodeling, consensus clustering)

All modules are included and all public functions exported from the top-level `src/SNA.jl` module.

## Key Dependencies

- **`Network`** — The core network data structure (`Network{T}`). All SNA functions accept `Network` objects and use its API (`nv`, `vertices`, `inneighbors`, `outneighbors`, `add_edge!`, `is_directed`, `network_density`, `as_matrix`, etc.). **This is an unregistered package** (UUID `027a387e-...`); it must be added via a local path or git URL before anything else will work. The `docs/Project.toml` `[sources]` paths may need updating to match the actual location of the Network package on the current machine.
- **`Graphs.jl`** — Graph algorithms accessed via `net.graph` (the underlying `SimpleDiGraph`/`SimpleGraph`)
- Hierarchical clustering for `equiv_clust`/`blockmodel` is implemented in-package (average-linkage UPGMA in `equivalence.jl`); there is no Clustering.jl dependency

## Conventions

- Functions follow R `sna` naming where possible (e.g., `gden`, `grecip`, `gtrans` as aliases)
- Measures follow R `sna` conventions exactly where a counterpart exists; golden-master values from sna 2.8 (Sampson samplike, Padgett flomarriage) are baked into the tests
- Directed networks are the default; undirected networks are created with `Network(n; directed=false)`
- Centrality functions return `Vector{Float64}`; equivalence functions return `Matrix{Float64}`
- Tests are in a single file `test/runtests.jl` using nested `@testset` blocks
- Documentation uses Documenter.jl with source pages in `docs/src/`
- Requires Julia >= 1.12 (Network.jl cannot load on earlier versions)
