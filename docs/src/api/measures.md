# Measures API Reference

This page documents network-level measures, path functions, and census functions available in SNA.jl.

## Network-Level Measures

### Density and Reciprocity

```@docs
density
gden
reciprocity
grecip
mutuality
```

### Transitivity and Hierarchy

```@docs
transitivity
gtrans
hierarchy
efficiency
connectedness
```

## Census Functions

Functions for classifying dyads, triads, and components.

```@docs
dyad_census
triad_census
component_dist
```

## Path Functions

Functions for computing distances, reachability, and path-based summaries.

```@docs
geodesic_distance
reachability
diameter
average_path_length
```

## QAP Inference and Network Regression

Permutation tests and network regression following R `sna::qaptest`,
`sna::netlm`, and `sna::netlogit`.

```@docs
qaptest
netlm
netlogit
```

### Result Types

```@docs
QAPTestResult
NetLMResult
NetLogitResult
```

### Result Metadata

SNA.jl implements the ecosystem's
[result-metadata protocol](https://Statistical-network-analysis-with-Julia.github.io/Networks.jl/dev/api/metadata/)
for its regression results, so what a fit did is inspectable via
`Networks.fit_metadata(result)`.

```@docs
objective(::NetLMResult)
objective(::NetLogitResult)
is_exact(::NetLMResult)
is_exact(::NetLogitResult)
se_method(::NetLMResult)
se_method(::NetLogitResult)
```


### Regression accessors and reproducibility

`netlm` and `netlogit` implement `coef`, `stderror`, `vcov`, `confint`,
`loglikelihood`, `nobs`, `dof`, `aic`, `bic`, and `coeftable` from StatsAPI.
The OLS degrees of freedom for information criteria include residual variance.
Covariance and confidence intervals assume independent dyads; the permutation
p-values provide QAP inference. `netlogit` exposes `converged` and warns on a
failed point fit; failed permutation fits raise an error. Complete and
quasi-complete separation are checked using candidate rays from the fit and
its next Newton step, then verified exactly against the input design. This
certificate check is not an exhaustive linear-programming existence test.

A finite observed fit does not ensure that all QAP fits are identified. Sparse
binary predictors can have a permuted stratum containing only successes or only
failures. When any requested permutation is separated or cannot be fitted,
`netlogit` names the failing replicate and returns no p-value; no replicate is
silently omitted. `nullhyp=:classical` only fits the observed design and reports
independent-dyad reference inference, not QAP inference.

```julia
using SNA, Random
flo = load_dataset(:florentine_marriage)
biz = load_dataset(:florentine_business)
fit = netlm(flo, biz; reps=50, rng=Xoshiro(2026))
coef(fit)
stderror(fit)
vcov(fit)
confint(fit)
coeftable(fit)
```

Replicate RNGs are seeded before scheduling, so `threaded=false` and
`threaded=true` produce the same draws for the same initial RNG.
A `qaptest` callback must be thread-safe when parallel execution is enabled.
Two-mode network inputs restrict permutations to actors in the same mode.
Zero permutation tail counts display as less than `1/reps`, preserving the
raw empirical p-values in the result.
