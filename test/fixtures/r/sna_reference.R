# Regenerate from SNA.jl: Rscript test/fixtures/r/sna_reference.R
# Requires sna, network, and ergm (only for the published data).
# All values below are computed by R, never copied from Julia output.
suppressPackageStartupMessages(library(sna))
suppressPackageStartupMessages(library(network))
data(florentine, package="ergm")
data(sampson, package="ergm")
set.seed(20260914)
values <- list()
put <- function(key, x) values[[key]] <<- unname(as.vector(x))
mask <- function(x) sum(2^(sort(unique(x))-1))
clique_masks <- function(A, mode) {
  cs <- unlist(clique.census(A, mode=mode, enumerate=TRUE)$cliques, recursive=FALSE)
  sort(vapply(Filter(function(x) length(x)>=3, cs), mask, 0.0))
}
for (id in c("flo", "samp")) {
  net <- if(id=="flo") flomarriage else samplike
  A <- as.matrix.network.adjacency(net)
  mode <- if(id=="flo") "graph" else "digraph"
  add <- function(key,x) put(paste0(id,"_",key),x)
  add("adjacency", A)
  add("density", gden(A, mode=mode))
  add("reciprocity", grecip(A))
  add("reciprocity_edgewise", grecip(A, measure="edgewise"))
  add("reciprocity_nonnull", grecip(A, measure="dyadic.nonnull"))
  add("transitivity", gtrans(A, mode=mode))
  add("dyad_census", dyad.census(A))
  add("triad_census", triad.census(A, mode=mode))
  add("connectedness", connectedness(A))
  add("efficiency", efficiency(A))
  add("hierarchy", hierarchy(A))
  add("hierarchy_krackhardt", hierarchy(A, measure="krackhardt"))
  add("degree", degree(A, gmode=mode))
  add("degree_in", degree(A, gmode=mode, cmode="indegree"))
  add("degree_out", degree(A, gmode=mode, cmode="outdegree"))
  add("betweenness", betweenness(A, gmode=mode))
  add("closeness", closeness(A, gmode=mode))
  add("eigenvector", abs(evcent(A, gmode=mode, use.eigen=TRUE)))
  add("bonacich", bonpow(A, gmode=mode, exponent=0.05))
  add("flowbet", flowbet(A, gmode=mode))
  add("strong_sizes", sort(component.dist(A)$csize, decreasing=TRUE))
  add("weak_sizes", sort(component.dist(A, connected="weak")$csize, decreasing=TRUE))
  add("cutpoints", cutpoints(A, mode=mode))
  add("cutpoints_weak", cutpoints(A, mode=mode, connected="weak"))
  add("cutpoints_recursive", cutpoints(A, mode=mode, connected="recursive"))
  add("bicomponents", sort(vapply(bicomponent.dist(A)$members, mask, 0.0)))
  add("cliques", clique_masks(A, mode))
  for(k in c(1,2,5,7)) add(paste0("kcore_",k), which(kcores(A, mode=mode)>=k))
  for (fun in c("degree", "betweenness", "closeness", "evcent")) {
    add(paste0("centralization_",fun), centralization(A, get(fun), mode=mode))
  }
  add("centralization_in", centralization(A, degree, mode=mode, cmode="indegree"))
  add("centralization_out", centralization(A, degree, mode=mode, cmode="outdegree"))
  add("centralization_raw", centralization(A, degree, mode=mode, normalize=FALSE))
}
A <- matrix(0,4,4); A[cbind(c(1,2,3),c(2,3,1))] <- 1
put("cycle_cutpoints", cutpoints(A))
put("cycle_cutpoints_recursive", cutpoints(A, connected="recursive"))
put("cycle_strong_sizes", sort(component.dist(A)$csize, decreasing=TRUE))
put("cycle_closeness", closeness(A))
put("cycle_cliques", clique_masks(A, "digraph"))
E <- matrix(0,3,3)
put("empty_transitivity", gtrans(E))
put("empty_nonnull_nan", is.nan(grecip(E, measure="dyadic.nonnull")))
put("empty_edgewise_nan", is.nan(grecip(E, measure="edgewise")))
put("singleton_reciprocity_nan", is.nan(grecip(matrix(0,1,1))))
for (directed in c(FALSE, TRUE)) {
  id <- if(directed) "bip_dir" else "bip_undir"
  B <- network.initialize(5, directed=directed, bipartite=2)
  add.edges(B, c(1,1,2), c(3,4,5))
  if(directed) add.edges(B, 3, 2)
  A <- as.matrix.network.adjacency(B, expand.bipartite=TRUE)
  mode <- if(directed) "digraph" else "graph"
  put(paste0(id,"_density"), network.density(B, discount.bipartite=FALSE))
  put(paste0(id,"_discount_density"), network.density(B, discount.bipartite=TRUE))
  put(paste0(id,"_degree"), degree(A, gmode=mode))
  put(paste0(id,"_bonacich"), bonpow(A, exponent=0.05, gmode=mode))
  put(paste0(id,"_flowbet"), flowbet(A,gmode=mode))
}
L <- matrix(0,4,4); L[cbind(c(1,1,2,3),c(1,2,3,4))] <- 1
put("loops_degree", degree(L)); put("loops_degree_diag",degree(L,diag=TRUE))
put("loops_density",gden(L)); put("loops_density_diag",gden(L,diag=TRUE))
W <- matrix(c(0,2,1,1,0,3,4,2,0),3,3,byrow=TRUE)
put("weighted_evcent",abs(evcent(W,use.eigen=TRUE)))
put("weighted_bonacich",bonpow(W,exponent=0.05))
put("weighted_flowbet",flowbet(W))
put("weighted_degree",degree(W,ignore.eval=FALSE))
# Regression fixtures use the data's published wealth attribute, in actor order.
F <- as.matrix.network.adjacency(flomarriage)
B <- as.matrix.network.adjacency(flobusiness)
wealth <- get.vertex.attribute(flomarriage,"wealth")
D <- abs(outer(wealth,wealth,"-"))
put("flo_wealth", wealth)
put("flo_business_adjacency", B)
put("qap_gcor",gcor(F,B))
fit <- netlm(F,list(B,D),mode="graph",nullhyp="classical")
put("lm_r_squared", 1-sum(fit$residuals^2)/sum((F[lower.tri(F)]-mean(F[lower.tri(F)]))^2))
for (key in c("coefficients","tstat","pgreqabs")) put(paste0("lm_",gsub("\\.","_",key)),fit[[key]])
fit <- netlm(as.matrix.network.adjacency(samplike),t(as.matrix.network.adjacency(samplike)),nullhyp="classical")
for (key in c("coefficients","tstat")) put(paste0("lm_directed_",key),fit[[key]])
fit <- netlogit(F,list(B,D),mode="graph",nullhyp="classical")
for (key in c("coefficients","se","tstat","pgreqabs","deviance","null.deviance","aic","bic")) put(paste0("logit_",gsub("\\.","_",key)),fit[[key]])
# Versioned tolerance: deterministic measures round at floating precision;
# sna's iterative centralization(evcent) differs slightly from a direct eigen
# solve. Logistic scores use different convergent optimizers (R glm vs Newton).
out <- c('name = "sna_reference"', '', '[provenance]',
 sprintf('r_version = "%s"',getRversion()),
 sprintf('sna_version = "%s"',packageVersion("sna")),
 sprintf('network_version = "%s"',packageVersion("network")),
 sprintf('ergm_version = "%s"',packageVersion("ergm")),
 'seed = 20260914', 'script = "test/fixtures/r/sna_reference.R"',
 sprintf('date = "%s"',Sys.Date()),
 'dataset = "ergm florentine and sampson; explicit cycle, empty, weighted and two-mode examples"',
 '', '[tolerance]',
 '# Deterministic graph counts, ratios and direct linear algebra; no Monte Carlo tolerance.',
 'default = 1e-9',
 '# sna centralization(evcent) uses iterative evcent; Julia uses a direct eigen solve.',
 'flo_centralization_evcent = 1e-6', 'samp_centralization_evcent = 1e-6',
 '# R glm and Julia shared Newton optimizer terminate by different criteria.',
 'logit_coefficients = 1e-5', 'logit_se = 1e-5', 'logit_tstat = 1e-4',
 'logit_pgreqabs = 1e-6', '', '[values]')
encode <- function(x) {
 if(is.logical(x)) return(ifelse(x,"true","false"))
 ifelse(is.nan(x),'nan',ifelse(is.infinite(x),ifelse(x>0,'inf','-inf'),sprintf('%.17g',x)))
}
for (key in names(values)) {
 x <- values[[key]]
 # Preserve vectors (including length-one) where the API returns a vector.
 vector <- length(x)!=1 || grepl('adjacency|degree|betweenness$|closeness$|eigenvector$|bonacich|flowbet|census|sizes|cutpoints|bicomponents|cliques|kcore|wealth|coefficients|tstat|pgreqabs|_se$|weighted_evcent',key)
 # Centralizations are scalar, despite the centrality name in their keys.
 if(grepl('centralization',key)) vector <- FALSE
 val <- if(vector) paste0('[',paste(encode(x),collapse=', '),']') else encode(x)
 out <- c(out,paste0(key,' = ',val))
}
writeLines(out,'test/fixtures/sna_reference.toml')
