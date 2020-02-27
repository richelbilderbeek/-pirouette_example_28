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

suppressMessages(library(ggplot2))
suppressMessages(library(pirouette))
suppressMessages(library(dplyr))
suppressMessages(library(pryr))
suppressMessages(library(babette))
library(testthat)
testthat::expect_true(mcbette::can_run_mcbette())

################################################################################
# Constants
################################################################################
root_folder <- getwd()
example_no <- 28
rng_seed <- 314
crown_age <- 10
n_phylogenies <- 5
is_testing <- is_on_travis()

if (is_testing) {
  n_phylogenies <- 2
  root_folder <- tempdir()
}

################################################################################
# Create simulation function
################################################################################
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

################################################################################
# Create phylogenies
################################################################################
phylogenies <- list()
for (i in seq_len(n_phylogenies)) {
  set.seed(314 - 1 + i)
  phylogenies[[i]] <- sim_tree_fun()
}
expect_equal(length(phylogenies), n_phylogenies)
################################################################################
# Create pirouette parameter sets
################################################################################
pir_paramses <- list()
for (i in seq_along(phylogenies)) {

  alignment_params <- create_alignment_params(
    sim_tral_fun = get_sim_tral_with_std_nsm_fun(
      mutation_rate = 1.0 / crown_age
    ),
    root_sequence = create_blocked_dna(length = 1000)
  )

  # Hand-pick a generating model
  # By default, this is JC69, strict, Yule
  generative_experiment <- create_gen_experiment()
  # Create the set of candidate birth-death experiments
  candidate_experiments <- create_all_bd_experiments(
    exclude_model = generative_experiment$inference_model
  )
  # Combine all experiments
  experiments <- c(list(generative_experiment), candidate_experiments)

  twinning_params <- create_twinning_params(
    sim_twin_tree_fun = get_sim_bd_twin_tree_fun(),
    sim_twal_fun = get_sim_twal_same_n_muts_fun(
      mutation_rate = 1.0 / crown_age,
      max_n_tries = 10000
    ),
    twin_evidence_filename = get_temp_evidence_filename()
  )

  pir_params <- create_pir_params(
    alignment_params = alignment_params,
    experiments = experiments,
    twinning_params = twinning_params,
    evidence_filename = get_temp_evidence_filename()
  )

  pir_paramses[[i]] <- pir_params
}
expect_equal(length(pir_paramses), n_phylogenies)
################################################################################
# Shorter run on Travis
################################################################################
if (is_testing) {
  for (i in seq_along(pir_paramses)) {
    pir_paramses[[i]]$experiments <- shorten_experiments(
      pir_paramses[[i]]$experiments
    )
  }
}

################################################################################
# Set the RNG seeds
################################################################################
pir_paramses <- renum_rng_seeds(
  pir_paramses = pir_paramses,
  rng_seeds = seq(314, 314 - 1 + length(pir_paramses))
)

################################################################################
# Rename filenames
################################################################################
for (i in seq_along(pir_paramses)) {
  rng_seed <- pir_paramses[[i]]$alignment_params$rng_seed
  pir_paramses[[i]] <- pir_rename_to_std(
    pir_params = pir_paramses[[i]],
    folder_name = file.path(root_folder, paste0("example_", example_no, "_", rng_seed))
  )
}

################################################################################
# Save tree to files
################################################################################
for (i in seq_along(pir_paramses)) {
  testthat::expect_equal(length(pir_paramses), length(phylogenies))
  rng_seed <- pir_paramses[[i]]$alignment_params$rng_seed
  folder_name <- file.path(root_folder, paste0("example_", example_no, "_", rng_seed))

  # Create folder, do not warn if it already exists
  dir.create(folder_name, showWarnings = FALSE, recursive = TRUE)
  ape::write.tree(
    phylogenies[[i]],
    file = file.path(folder_name, "true_tree.newick")
  )
}

################################################################################
# Delete previous files
################################################################################
for (pir_params in pir_paramses) {
  check_pir_params(pir_params)
  rm_pir_param_files(pir_params)
}

################################################################################
# Do the runs
################################################################################
pir_outs <- pir_runs(
  phylogenies = phylogenies,
  pir_paramses = pir_paramses
)

################################################################################
# Save
################################################################################
for (i in seq_along(pir_outs)) {
  testthat::expect_equal(length(pir_paramses), length(pir_outs))
  rng_seed <- pir_paramses[[i]]$alignment_params$rng_seed
  folder_name <- file.path(root_folder, paste0("example_", example_no, "_", rng_seed))

  utils::write.csv(
    x = pir_outs[[i]],
    file = file.path(folder_name, "errors.csv"),
    row.names = FALSE
  )

  pir_plots(pir_outs[[i]]) +
    ggsave(file.path(folder_name, "errors.png"))
}

