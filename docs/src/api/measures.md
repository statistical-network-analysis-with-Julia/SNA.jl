# Measures API Reference

Graph-level indices, censuses, graph correlation, and QAP inference.

## Graph-level indices

```@docs
gden
grecip
gtrans
transitivity
mutuality
hierarchy
efficiency
connectedness
brokerage
```

## Census functions

```@docs
dyad_census
triad_census
```

## Graph correlation

```@docs
gcor
gcov
```

## QAP Inference and Network Regression

Permutation tests and network regression following R `sna::qaptest`,
`sna::netlm`, and `sna::netlogit`.

```@docs
qaptest
netlm
netlogit
```

### Conditional uniform graph tests

```@docs
cug_test
CUGTestResult
```

### Result Types

```@docs
QAPTestResult
NetLMResult
NetLogitResult
```

### Result Metadata

SNA.jl implements the ecosystem's
[result-metadata protocol](https://Statistical-network-analysis-with-Julia.github.io/NetworkCore.jl/dev/api/metadata/)
for its regression results, so what a fit did is inspectable via
`NetworkCore.fit_metadata(result)`.

```@docs
objective(::NetLMResult)
objective(::NetLogitResult)
is_exact(::NetLMResult)
is_exact(::NetLogitResult)
se_method(::NetLMResult)
se_method(::NetLogitResult)
```

### Inference contracts

`netlm` and `netlogit` implement `coef`, `coefnames`, `stderror`, `vcov`,
`confint`, `loglikelihood`, `nobs`, `dof`, `aic`, `bic`, and `coeftable` from
StatsAPI.
The point estimates are an ordinary OLS or logistic fit that treats the dyads
as independent, and so are the covariance and confidence intervals; dyadic
dependence makes them anti-conservative. Inference about the coefficients
comes from the QAP permutation p-values (`pgreqabs`, printed by `show`).

- The default null hypothesis is Dekker's double semi-partialing
  (`:qapspp`), which keeps its size when predictors are correlated;
  `:qapy` and `:qapx` are also available. `nullhyp=:classical` reports
  parametric tests that assume independent dyads: for reference only.
- A p-value is the proportion of permuted statistics at least as extreme as
  the observed one, without a `+1` correction, as in sna; it can be 0 and is
  printed as `< 1/n_sim`.
- `netlogit`'s QAP statistic is the signed root likelihood-ratio statistic
  by default. Sparse binary data often have permutations whose logistic fit
  is separated (an infinite coefficient); the likelihood-ratio statistic
  exists there, so every replicate counts. `statistic=:wald` uses sna's Wald
  z, which does not exist on a separated permutation: the permutation
  p-values of the coefficients such a replicate touches are then withheld
  (`NaN`, with a warning), never computed from the surviving replicates.
- An observed design that is completely or quasi-completely separated has no
  finite estimate. `netlogit` follows the ecosystem's separation policy: it
  warns, returns `converged == false` with the separated coefficients in
  `fit.separated`, and withholds the test statistics, p-values and
  confidence intervals (`NaN`); no permutation is run. The verdict is
  NetworkCore.jl's `logistic_separation`: a linear programme over the
  design, with the separating direction verified in exact arithmetic.
- Two-mode data use the cross-mode dyads only, and permutations reorder each
  mode separately.
- Raw matrices with `mode=:auto` are undirected when the response and every
  predictor are symmetric (sna treats every matrix as directed).

```julia
using SNA, Random
flo = load_dataset(:florentine_marriage)
biz = load_dataset(:florentine_business)
fit = netlm(flo, biz; n_sim=200, rng=Xoshiro(2026))
coef(fit)
coefnames(fit)                       # ["(intercept)", "x1"], R's names(coef(fit))
stderror(fit)
confint(fit)
coeftable(fit)
fit.pgreqabs                         # the QAP p-values
```

Replicate random number generators are seeded from `rng` before
scheduling, so `threaded=false` and `threaded=true` give the same result. A
`qaptest` statistic must be thread-safe when threading is on.
