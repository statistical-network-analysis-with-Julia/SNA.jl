# Regenerate from SNA.jl: Rscript test/fixtures/r/sna_inference.R
# R sna 2.8 reference values for brokerage, equiv.clust + cutree, blockmodel
# and cug.test. All values are computed by R, never copied from Julia output.
# Keys: g<i>_meta = [n, directed]; g<i>_edges = flattened (from, to) pairs;
# g<i>_cl = class of each vertex. Matrices are column-major (as.vector).
suppressPackageStartupMessages(library(sna))
suppressPackageStartupMessages(library(network))
data(florentine, package = "ergm"); data(sampson, package = "ergm")
seed <- 20261003
set.seed(seed)
n_graphs <- 24
cug_reps <- 20000
values <- list()
put <- function(key, x) values[[key]] <<- unname(as.vector(x))
for (g in 1:n_graphs) {
  n <- sample(5:11, 1)
  directed <- runif(1) < 0.65
  p <- sample(c(0.15, 0.25, 0.4, 0.6), 1)
  A <- matrix(0, n, n)
  for (i in 1:n) for (j in 1:n) {
    if (i == j || (!directed && j < i)) next
    if (runif(1) < p) { A[i, j] <- 1; if (!directed) A[j, i] <- 1 }
  }
  E <- which(A == 1, arr.ind = TRUE)
  if (!directed) E <- E[E[, 1] <= E[, 2], , drop = FALSE]
  cl <- sample(1:sample(2:4, 1), n, replace = TRUE)
  id <- paste0("g", g, "_")
  add <- function(key, x) put(paste0(id, key), x)
  add("meta", c(n, directed)); add("edges", as.vector(t(E))); add("cl", cl)
  b <- brokerage(A, cl)
  for (nm in c("raw.nli", "exp.nli", "sd.nli", "z.nli", "raw.gli", "exp.gli", "sd.gli",
               "z.gli", "exp.grp", "sd.grp", "clid", "n"))
    add(paste0("brokerage_", gsub("\\.", "_", nm)), b[[nm]])
  # equiv.clust defaults (sedist Hamming, complete linkage), cut with cutree
  ec <- equiv.clust(A)
  for (k in 2:4) add(paste0("equiv_clust_k", k), cutree(ec$cluster, k))
  add("blockmodel_k3", blockmodel(A, ec, k = 3)$block.model)
  for (cm in c("average", "single")) {
    e2 <- equiv.clust(A, method = "euclidean", cluster.method = cm)
    add(paste0("equiv_clust_euclidean_", cm, "_k3"), cutree(e2$cluster, 3))
  }
}
# cug.test: the observed statistic is exact; the null distribution is Monte
# Carlo, so its mean, sd and tail proportions are frozen from a long R run
# and compared within Monte Carlo error.
cug <- function(key, A, FUN, mode, cmode, ...) {
  ct <- cug.test(A, FUN, mode = mode, cmode = cmode, reps = cug_reps, ...)
  put(paste0("cug_", key, "_obs"), ct$obs.stat)
  put(paste0("cug_", key, "_mean"), mean(ct$rep.stat))
  put(paste0("cug_", key, "_sd"), sd(ct$rep.stat))
  put(paste0("cug_", key, "_pgteobs"), ct$pgteobs)
  put(paste0("cug_", key, "_plteobs"), ct$plteobs)
}
F <- as.matrix.network.adjacency(flomarriage)
S <- as.matrix.network.adjacency(samplike)
for (cm in c("size", "edges", "dyad.census")) {
  ck <- gsub("\\.", "_", cm)
  cug(paste0("flo_gtrans_", ck), F, gtrans, "graph", cm)
  cug(paste0("samp_gtrans_", ck), S, gtrans, "digraph", cm)
  cug(paste0("samp_grecip_", ck), S, grecip, "digraph", cm, FUN.args = list(measure = "edgewise"))
}
put("cug_reps", cug_reps)
out <- c('name = "sna_inference"', '', '[provenance]',
  sprintf('r_version = "%s"', getRversion()),
  sprintf('sna_version = "%s"', packageVersion("sna")),
  sprintf('seed = %d', seed), 'script = "test/fixtures/r/sna_inference.R"',
  sprintf('date = "%s"', Sys.Date()),
  sprintf('dataset = "%d seeded random networks with random classes (brokerage, equiv.clust, blockmodel); ergm florentine and sampson (cug.test, %d replications)"', n_graphs, cug_reps),
  '', '[tolerance]',
  '# Deterministic counts and closed-form moments. The cug_* Monte Carlo',
  '# summaries are compared in the testset within 5 Monte Carlo standard errors.',
  'default = 1e-9', '', '[values]')
encode <- function(x) {
  if (is.logical(x)) x <- as.numeric(x)
  ifelse(is.nan(x), 'nan', ifelse(is.na(x), 'nan', ifelse(is.infinite(x), ifelse(x > 0, 'inf', '-inf'), sprintf('%.17g', x))))
}
for (key in names(values)) out <- c(out, paste0(key, ' = [', paste(encode(values[[key]]), collapse = ', '), ']'))
writeLines(out, 'test/fixtures/sna_inference.toml')
