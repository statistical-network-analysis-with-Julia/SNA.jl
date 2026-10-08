# SNA.jl

[![Network Analysis](https://img.shields.io/badge/Network-Analysis-orange.svg)](https://github.com/statistical-network-analysis-with-Julia/SNA.jl)
[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://statistical-network-analysis-with-Julia.github.io/SNA.jl/dev/)
[![Julia](https://img.shields.io/badge/Julia-1.12+-purple.svg)](https://julialang.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

<p align="center">
  <img src="docs/src/assets/logo.svg" alt="SNA.jl icon" width="160">
</p>

A Julia port of R's [`sna`](https://cran.r-project.org/package=sna) package
(Butts 2008): descriptive social network analysis and QAP inference on
[NetworkCore.jl](https://github.com/statistical-network-analysis-with-Julia/NetworkCore.jl)
networks.

- **Centrality**: `degreecent` (sna's `degree`), `betweenness`, `closeness`,
  `evcent`, `bonpow`, `infocent`, `flowbet`, and Freeman `centralization`
- **Graph-level indices**: `gden`, `grecip`, `gtrans`, `dyad_census`,
  `triad_census`, `mutuality`, `hierarchy`, `efficiency`, `connectedness`,
  Gould–Fernandez `brokerage`, and local clustering coefficients
  (`transitivity`)
- **Cohesion**: `component_dist`, `largest_component`, `cliques`, `kcores`,
  `cutpoints`, `bicomponents`, `geodist`, `reachability`
- **Positions**: `sedist`, `regular_equivalence`, `equiv_clust`,
  `blockmodel`, `consensus_clustering`
- **Inference**: `gcor`, `gcov`, `qaptest`, `netlm`, `netlogit` (Dekker
  double semi-partialing by default), and conditional uniform graph tests
  (`cug_test`)
- **Utilities**: layouts and random graphs (`rgraph`, `rgnm`, `rgnp`)

A function that ports an sna function has sna's name (with `.` written as
`_`) and sna's semantics. Agreement with sna 2.8 is pinned by R-generated
golden fixtures, including a corpus of 60 random networks (directed and
undirected, with self-loops, isolates and disconnected parts) on which every
measure is compared with R; `brokerage`, `equiv_clust`, `blockmodel` and
`cug_test` have their own R fixture.

## Names changed in 0.2.0

Before 0.2.0 SNA.jl added methods to nine Graphs.jl functions for `Network`,
so loading SNA (even through TSNA) changed what `Graphs.betweenness_centrality`
and the others returned. SNA.jl no longer extends any Graphs.jl function:
those functions keep Graphs.jl's meaning on a `Network`, and SNA's
R-semantics measures carry sna's names.

| Before 0.2.0 | Now | Notes |
|---|---|---|
| `degree_centrality(net; mode=:total/:in/:out)` | `degreecent(net; cmode=:freeman/:indegree/:outdegree)` | Graphs.jl's `degree_centrality` |
| `betweenness_centrality(net)` | `betweenness(net)` | Graphs.jl's function normalizes by default |
| `closeness_centrality(net; cmode=:sna)` | `closeness(net)` | `cmode` now takes sna's values; `cmode=:component` is `Graphs.closeness_centrality(net)` |
| `eigenvector_centrality(net)` | `evcent(net)` | |
| `katz_centrality`, `pagerank` | `Graphs.katz_centrality`, `Graphs.pagerank` | not sna measures; Graphs.jl's conventions |
| `density(net)` | `gden(net)` | `Graphs.density` counts self-loops |
| `diameter(net)` | `maximum(geodist(net).gdist)` | or `Graphs.diameter(net)` |
| `bridges(net)` | `Graphs.bridges(net)` | undirected networks |
| `bonacich_power(net)` | `bonpow(net)` | sna's name |
| `geodesic_distance(net)` | `geodist(net).gdist` | `geodist` also returns the geodesic counts |
| `reciprocity(net; method=m)` | `grecip(net; measure=m)` | sna's name and keyword |
| `structural_equivalence(net; method=m)` | `sedist(net; method=m)` | `sedist` defaults to Hamming, a count |
| `consensus(clusterings)` | `consensus_clustering(clusterings)` | sna's `consensus` is a different procedure |
| `components(net)` | `component_dist(net).membership` | the old name collided with `Distributions.components` |
| `component_dist(net)` (sorted sizes) | `sort(component_dist(net).csize; rev=true)` | `component_dist` now returns sna's `(membership, csize, cdist)` |
| `kcores(net; k)` (members of the k-core) | `findall(>=(k), kcores(net))` | `kcores` now returns core numbers, as sna; passing `k` raises an error that says so |
| `gtrans(net; type=…)`, `grecip(net; method=…)` | `gtrans(net; measure=…)`, `grecip(net; measure=…)` | sna's keywords and measures |
| `equiv_clust(net; method=:structural)` | `equiv_clust(net; method=:correlation, cluster_method=:average)` | the default is now sna's (Hamming distance, complete linkage) |
| `equiv_clust(net; method=:regular)` | `equiv_clust(net; equiv_fun=regular_equivalence, cluster_method=:average)` | |
| `reps=` (QAP), `max_iter=` | `n_sim=`, `maxiter=` | the ecosystem's keyword names |

The old names were never part of a release, so they are not kept as
deprecated aliases: calling one is an `UndefVarError`, and an old keyword is
a `MethodError` naming it.

## Installation

Requires Julia 1.12+. SNA.jl depends on the unregistered
[NetworkCore.jl](https://github.com/statistical-network-analysis-with-Julia/NetworkCore.jl), which must be added first:

```julia
using Pkg
Pkg.add(url="https://github.com/statistical-network-analysis-with-Julia/NetworkCore.jl")
Pkg.add(url="https://github.com/statistical-network-analysis-with-Julia/SNA.jl")
```

The examples below also use the standard library `Random` and, for one
line, `Graphs` (`Pkg.add("Graphs")`). For development of the whole
ecosystem, the website's
[workspace guide](https://statistical-network-analysis-with-julia.github.io/getting-started/)
prepares side-by-side clones with `tools/prepare_workspace.jl`.

## Usage

### Centrality and network indices

```julia
using SNA

net = network(6; directed=false)
for (i, j) in [(1, 2), (1, 3), (2, 3), (3, 4), (4, 5), (4, 6), (5, 6)]
    add_edge!(net, i, j)
end

degreecent(net)              # [2, 2, 3, 3, 2, 2]
betweenness(net)             # [0, 0, 6, 6, 0, 0]: raw counts, as in sna
closeness(net)               # Freeman closeness
evcent(net)                  # unit-length principal eigenvector
centralization(net, degreecent)
gden(net)                    # 7 / 15
gtrans(net)                  # 0.6
transitivity(net; type=:local)
```

### Censuses and cohesion

```julia
using SNA

d = network(4)
for (i, j) in [(1, 2), (2, 1), (2, 3), (3, 4)]
    add_edge!(d, i, j)
end
dyad_census(d)               # (mutual = 1, asymmetric = 2, null = 3)
grecip(d)                    # (1 + 3) / 6
triad_census(d)              # 16 Davis–Leinhardt classes
component_dist(d).csize      # strong components: [2, 1, 1]
kcores(d)                    # core numbers (Freeman degree)
sort(cutpoints(d; connected=:weak))   # [2, 3]
g = geodist(d)
g.gdist[1, 4], g.gdist[4, 1]          # (3.0, Inf)
maximum(g.gdist)                      # Inf: the diameter of a network that is not strongly connected
```

### Structural equivalence and blockmodels

```julia
using SNA

net = network(4)
for (i, j) in [(1, 2), (2, 1), (1, 3), (2, 3), (4, 1), (4, 2)]
    add_edge!(net, i, j)
end
sedist(net)[1, 2]                     # 0: same ties to and from 3 and 4
sedist(net; method=:correlation)[1, 2]   # 1
bm = blockmodel(net; k=3)             # equiv_clust: Hamming, complete linkage, as sna
bm.membership, bm.block_matrix

b = brokerage(net, [:a, :a, :b, :b])  # Gould–Fernandez roles given a partition
b.raw_nli, b.z_gli
```

### QAP test and network regression

```julia
using SNA, Random

flo = load_dataset(:florentine_marriage)
biz = load_dataset(:florentine_business)

gcor(flo, biz)                                      # 0.372, as R sna::gcor
qt = qaptest(gcor, flo, biz; n_sim=1000, rng=Xoshiro(1))
qt.pgreq                                            # about 0.001

# Marriage ties regressed on business ties: logistic, DSP null
fit = netlogit(flo, biz; n_sim=1000, rng=Xoshiro(2))
fit.pgreqabs

# Is the marriage network more transitive than its number of ties implies?
ct = cug_test(flo, gtrans; cmode=:edges, n_sim=1000, rng=Xoshiro(3))
ct.pgteobs                                          # about 0.33: no
```

`netlogit`'s QAP statistic is the signed root likelihood ratio, which exists
on the permutations whose logistic fit is separated (about 6 % here); sna's
Wald z (`statistic=:wald`) does not, and the p-values it would need are then
withheld (`NaN`, with a warning). A separated *observed* design has no
finite estimate: `netlogit` warns, returns `converged == false`, names the
separated coefficients in `fit.separated` and withholds its p-values and
intervals. Regression standard errors and intervals assume independent
dyads; inference comes from the permutation p-values.

## Conventions

- Ties are binary by default (`ignore_eval=true`); degree, eigenvector,
  Bonacich, information and flow centrality and `sedist` accept
  `ignore_eval=false, attr=:weight`. Self-loops are ignored unless `diag=true`.
- A masked (unobserved) dyad is refused by every measure (`missing=:error`);
  `missing=:face` reads stored values explicitly and is not a missing-data
  estimator.
- Undefined ratios are `NaN`, as in sna. `gtrans` with no two-path is 1, as
  in sna.
- Directed components and cutpoints use strong connectivity by default;
  cliques and bicomponents use mutual ties. `connected=:weak` and
  `symmetrize=:weak` ignore direction.
- `evcent` solves the eigenproblem directly (sna's `use.eigen=TRUE`) and is
  `NaN` on a network without a directed cycle.
- `netlm`/`netlogit` with raw matrices treat symmetric data as undirected
  (`mode=:auto`); sna treats every matrix as directed. Two-mode data use the
  cross-mode dyads only. Monte Carlo p-values are proportions without a `+1`
  correction, as in sna.
- `equiv_clust` is sna's `equiv.clust` (Hamming distance, complete linkage)
  cut at `k` clusters, as `cutree` cuts it; it returns the labels, not the tree.
- `regular_equivalence` is an iterative similarity, not UCINET's REGE.
- `layout_kamada_kawai` minimizes the Kamada–Kawai stress by stress
  majorization. Layouts and random graphs draw from an `rng` keyword and do
  not reproduce R's random-number stream.
- Undirected density with `diag=true` has `n(n+1)/2` possible dyads; sna
  uses `n²`.

The documentation's [Coming from R sna](https://statistical-network-analysis-with-Julia.github.io/SNA.jl/dev/r_concordance/)
page maps every sna function.

## Not implemented

- Network autocorrelation models (`lnam`, `nacf`) and canonical correlation
  (`netcancor`).
- `cug_test` conditioning with self-loops (`diag=TRUE`) or on edge values
  (`ignore.eval=FALSE`), and on two-mode networks (an `ArgumentError`).
- Structural correlation and distance (`gscor`, `gscov`, `structdist`,
  `hdist`), Bayesian network inference (`bbnam`, `bn`) and informant
  consensus networks (`consensus`).
- The `prestige`, `stresscent`, `graphcent`, `loadcent` and `gilschmidt`
  indices (Gil–Schmidt is available as `closeness(net; cmode=:gil_schmidt)`),
  `lubness`, path and cycle censuses, `simmelian`, `ego.extract`, and the
  `rguman`, `rgws`, `rgbn` and `rgnmix` generators.
- sna's endpoint, proximal and length-scaled `betweenness` modes, `gtrans`'s
  `rank` measure and `component.dist`'s `unilateral` components (each raises
  an `ArgumentError`).
- Weighted shortest paths; weighted regular equivalence.
- Plotting (`gplot`); the layouts return coordinates only.

## Documentation

- [Development Documentation](https://statistical-network-analysis-with-Julia.github.io/SNA.jl/dev/)

## References

1. Wasserman, S., Faust, K. (1994). *Social Network Analysis: Methods and Applications*. Cambridge University Press.

2. Butts, C.T. (2008). Social network analysis with sna. *Journal of Statistical Software*, 24(6), 1-51.

3. Butts, C.T. (2024). sna: Tools for Social Network Analysis. R package version 2.8. [https://cran.r-project.org/package=sna](https://cran.r-project.org/package=sna)

4. Fagiolo, G. (2007). Clustering in complex directed networks. *Physical Review E*, 76, 026107.

5. Dekker, D., Krackhardt, D., Snijders, T.A.B. (2007). Sensitivity of MRQAP tests to collinearity and autocorrelation conditions. *Psychometrika*, 72(4), 563-581.

## Citation

If you use SNA.jl in your work, please cite it using the entry in
[`CITATION.bib`](CITATION.bib), and please also cite the R package `sna` and
the methods papers it ports (Butts 2008, and the papers of the measures you
use). The per-package list is on the ecosystem's
[How to cite](https://statistical-network-analysis-with-julia.github.io/citing/)
page.

```biblatex
@misc{SNWJSNAJL,
  author = {Santoni, Simone},
  title = {SNA.jl: Social Network Analysis in Julia},
  year = {2026},
  url = {https://github.com/statistical-network-analysis-with-Julia/SNA.jl},
  note = {Homepage: https://statistical-network-analysis-with-Julia.github.io/SNA.jl; GitHub: https://github.com/statistical-network-analysis-with-Julia}
}
```

## License

MIT License - see [LICENSE](LICENSE) for details.
