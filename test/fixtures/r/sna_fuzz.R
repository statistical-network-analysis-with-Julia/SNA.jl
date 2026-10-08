# Regenerate from SNA.jl: Rscript test/fixtures/r/sna_fuzz.R
# A seeded corpus of small random networks (directed and undirected, with and
# without self-loops, isolates and disconnected parts) and the values R sna
# 2.8 (and igraph, for local clustering) computes for them. All values are
# computed by R, never copied from Julia output. Keys: g<i>_meta = [n,
# directed, loops]; g<i>_edges = flattened (from, to) pairs; g<i>_<measure>.
# Matrices are stored column-major (R's as.vector).
suppressPackageStartupMessages(library(sna))
seed <- 20261002
set.seed(seed)
n_graphs <- 60
values <- list()
put <- function(key, x) values[[key]] <<- unname(as.vector(x))
tryput <- function(key, expr) {
  v <- tryCatch(suppressWarnings(expr), error = function(e) NULL)
  if (!is.null(v)) put(key, v)   # sna refuses: the key is simply absent
}
for (g in 1:n_graphs) {
  n <- sample(3:10, 1)
  directed <- runif(1) < 0.6
  loops <- runif(1) < 0.25
  p <- sample(c(0.1, 0.2, 0.35, 0.5, 0.7), 1)
  A <- matrix(0, n, n)
  for (i in 1:n) for (j in 1:n) {
    if (!directed && j < i) next
    if (i == j && !loops) next
    if (runif(1) < p) { A[i, j] <- 1; if (!directed) A[j, i] <- 1 }
  }
  E <- which(A == 1, arr.ind = TRUE)
  if (!directed) E <- E[E[, 1] <= E[, 2], , drop = FALSE]
  id <- paste0("g", g, "_")
  add <- function(key, x) put(paste0(id, key), x)
  try_add <- function(key, expr) tryput(paste0(id, key), expr)
  add("meta", c(n, directed, loops))
  add("edges", as.vector(t(E)))
  gm <- if (directed) "digraph" else "graph"
  add("degree", degree(A, gmode = gm))
  add("degree_in", degree(A, gmode = gm, cmode = "indegree"))
  add("degree_out", degree(A, gmode = gm, cmode = "outdegree"))
  add("betweenness", betweenness(A, gmode = gm))
  add("betweenness_undirected", betweenness(A, gmode = gm, cmode = "undirected"))
  add("closeness", closeness(A, gmode = gm))
  add("closeness_suminv", closeness(A, gmode = gm, cmode = if (directed) "suminvdir" else "suminvundir"))
  add("closeness_gil_schmidt", closeness(A, gmode = gm, cmode = "gil-schmidt"))
  try_add("bonpow", bonpow(A, gmode = gm, exponent = 0.3))
  try_add("infocent", infocent(A, gmode = gm))
  add("flowbet", flowbet(A, gmode = gm))
  add("centralization_degree", centralization(A, degree, mode = gm))
  add("centralization_betweenness", centralization(A, betweenness, mode = gm))
  add("centralization_closeness", centralization(A, closeness, mode = gm))
  add("gden", gden(A, mode = gm))
  for (m in c("dyadic", "dyadic.nonnull", "edgewise", "edgewise.lrr", "correlation"))
    add(paste0("grecip_", gsub("\\.", "_", m)), grecip(A, measure = m))
  for (m in c("weak", "strong", "weakcensus", "strongcensus", "correlation"))
    add(paste0("gtrans_", m), gtrans(A, mode = gm, measure = m))
  add("dyad_census", dyad.census(A))
  add("triad_census", triad.census(A, mode = gm))
  add("hierarchy", hierarchy(A))
  add("hierarchy_krackhardt", hierarchy(A, measure = "krackhardt"))
  add("efficiency", efficiency(A))
  add("connectedness", connectedness(A))
  for (cn in c("strong", "weak", "recursive")) {
    cd <- component.dist(A, connected = cn)
    add(paste0("components_", cn, "_membership"), cd$membership)
    add(paste0("components_", cn, "_csize"), cd$csize)
    add(paste0("components_", cn, "_cdist"), cd$cdist)
  }
  gd <- geodist(A, inf.replace = Inf)
  add("geodist", gd$gdist)
  add("geodist_counts", gd$counts)
  for (cm in c("freeman", "indegree", "outdegree"))
    add(paste0("kcores_", cm), kcores(A, mode = gm, cmode = cm))
  for (m in c("hamming", "correlation", "euclidean", "gamma", "exact"))
    add(paste0("sedist_", m), sedist(A, method = m))
  # graph correlation/covariance of the network with its transpose
  try_add("gcor_transpose", gcor(A, t(A), mode = gm))
  try_add("gcov_transpose", gcov(A, t(A), mode = gm))
  ig <- igraph::as_undirected(igraph::simplify(igraph::graph_from_adjacency_matrix(
    A, mode = if (directed) "directed" else "undirected")), mode = "collapse")
  add("local_clustering_weak", igraph::transitivity(ig, type = "localundirected", isolates = "NaN"))
  add("average_clustering_weak", igraph::transitivity(ig, type = "average", isolates = "NaN"))
}
# Structural equivalence of tied equivalent actors
T <- matrix(0, 4, 4); T[cbind(c(1, 1, 2), c(2, 3, 3))] <- 1; T <- T + t(T)
for (m in c("hamming", "correlation", "euclidean")) put(paste0("se_triangle_", m), sedist(T, method = m))
D <- matrix(0, 4, 4); D[cbind(c(1, 2, 1, 2, 4, 4), c(2, 1, 3, 3, 1, 2))] <- 1
for (m in c("hamming", "correlation", "euclidean")) put(paste0("se_pair_", m), sedist(D, method = m))
# Bonacich power on a well-conditioned system whose determinant underflows:
# a perfect matching of 200 actors, beta = 0.5.
M <- matrix(0, 200, 200); for (i in seq(1, 199, 2)) { M[i, i + 1] <- 1; M[i + 1, i] <- 1 }
put("matching_bonpow", bonpow(M, gmode = "graph", exponent = 0.5))
# Two-mode regression on the cross-mode dyads only: an
# 8 x 12 incidence response and predictor, R lm on the 96 cells.
Y <- matrix(rbinom(96, 1, 0.3), 8, 12); X <- matrix(rnorm(96), 8, 12)
put("twomode_y", Y); put("twomode_x", X)
put("twomode_lm_coefficients", coef(lm(as.vector(Y) ~ as.vector(X))))
put("twomode_lm_tstat", summary(lm(as.vector(Y) ~ as.vector(X)))$coefficients[, 3])
put("twomode_gcor", cor(as.vector(Y), as.vector(X)))
out <- c('name = "sna_fuzz"', '', '[provenance]',
  sprintf('r_version = "%s"', getRversion()),
  sprintf('sna_version = "%s"', packageVersion("sna")),
  sprintf('igraph_version = "%s"', packageVersion("igraph")),
  sprintf('seed = %d', seed), 'script = "test/fixtures/r/sna_fuzz.R"',
  sprintf('date = "%s"', Sys.Date()),
  sprintf('dataset = "%d seeded random networks (n = 3-10, directed and undirected, self-loops), structural-equivalence and Bonacich edge cases, a two-mode regression"', n_graphs),
  '', '[tolerance]',
  '# Deterministic counts, ratios, geodesics and direct linear algebra.',
  'default = 1e-9', '', '[values]')
encode <- function(x) {
  if (is.logical(x)) x <- as.numeric(x)
  ifelse(is.nan(x), 'nan', ifelse(is.na(x), 'nan', ifelse(is.infinite(x), ifelse(x > 0, 'inf', '-inf'), sprintf('%.17g', x))))
}
for (key in names(values)) {
  x <- values[[key]]
  out <- c(out, paste0(key, ' = [', paste(encode(x), collapse = ', '), ']'))
}
writeLines(out, 'test/fixtures/sna_fuzz.toml')
