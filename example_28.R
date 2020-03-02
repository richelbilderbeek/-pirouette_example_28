# pirouette example 28:
# multiple DD trees for pirouette article
#
# seed | folder
# -----|---------------
# 314  | example_28_314
# 315  | example_28_315
# 316  | example_28_316
# 317  | example_28_317
# 318  | example_28_318
suppressMessages(library(pirouette))
suppressMessages(library(dplyr))
suppressMessages(library(pryr))
suppressMessages(library(babette))
library(testthat)

# Constants
example_no <- 28
rng_seed <- 314
crown_age <- 10
n_phylogenies <- 5
is_testing <- is_on_travis()
if (is_testing) {
  n_phylogenies <- 2
}

# Create simulation function
sim_dd_tree_fun <- function(crown_age) {
  extinction_rate <- 0.1
  n_taxa <- 6
  n_0 <- 2 # Initial number of species at stem/crown of tree
  diff <- (log(n_taxa) - log(n_0)) / crown_age
  speciation_rate <- 3.0 * (diff + extinction_rate)
  carrying_capacity <- n_taxa # clade-level
  dd_parameters <- c(speciation_rate, extinction_rate, carrying_capacity)
  ddmodel <- 1 # linear dependence in speciation rate with parameter K
  dd_sim_result <- DDD::dd_sim(pars = dd_parameters, age  = crown_age, ddmodel = ddmodel)
  phylogeny <- dd_sim_result$tes # Only extant species
  phylogeny
}
sim_tree_fun <- pryr::partial(
  sim_dd_tree_fun,
  crown_age = crown_age
)

# Create phylogenies
phylogenies <- list()
for (i in seq_len(n_phylogenies)) {
  set.seed(314 - 1 + i)
  phylogenies[[i]] <- sim_tree_fun()
}
expect_equal(length(phylogenies), n_phylogenies)
# Create pirouette parameter sets
pir_paramses <- create_std_pir_paramses(n = length(phylogenies))
expect_equal(length(pir_paramses), n_phylogenies)
if (is_testing) {
  pir_paramses <- shorten_pir_params(pir_paramses)
}
# Save tree to files
for (i in seq_along(pir_paramses)) {
  expect_equal(length(pir_paramses), length(phylogenies))
  rng_seed <- pir_paramses[[i]]$alignment_params$rng_seed
  folder_name <- file.path(paste0("example_", example_no, "_", rng_seed))
  # Create folder, do not warn if it already exists
  dir.create(folder_name, showWarnings = FALSE, recursive = TRUE)
  ape::write.tree(
    phylogenies[[i]],
    file = file.path(folder_name, "true_tree.newick")
  )
}

# Do the runs
pir_outs <- pir_runs(
  phylogenies = phylogenies,
  pir_paramses = pir_paramses
)

# Save
for (i in seq_along(pir_outs)) {
  expect_equal(length(pir_paramses), length(pir_outs))
  rng_seed <- pir_paramses[[i]]$alignment_params$rng_seed
  folder_name <- file.path(paste0("example_", example_no, "_", rng_seed))
  pir_save(
    phylogeny = phylogenies[[i]],
    pir_params = pir_paramses[[i]],
    pir_out = pir_outs[[i]],
    folder_name = folder_name
  )
}

