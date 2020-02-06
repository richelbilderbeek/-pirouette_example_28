# pirouette example 28:
# multiple DD trees for pirouette article
suppressMessages(library(ggplot2))
library(pirouette)
library(babette)


is_testing <- FALSE
if (is_on_travis() || TRUE) {
  is_testing <- TRUE
}

if (1 == 2) {
  setwd("~/GitHubs/pirouette_example_28")
}
root_folder <- getwd()
example_no <- 28
rng_seed <- 314
example_folder <- file.path(root_folder, paste0("example_", example_no, "_", rng_seed))
dir.create(example_folder, showWarnings = FALSE, recursive = TRUE)
setwd(example_folder)
set.seed(rng_seed)
testit::assert(is_beast2_installed())

n_phylogenies <- 1
phylogenies <- list()

################################################################################
# Creates phylogenies
################################################################################
crown_age <- 10
for (i in seq_len(n_phylogenies)) {
  print(paste(i, "/", n_phylogenies))
  speciation_rate <- 0.8 # lambda
  extinction_rate <- 0.1 # mu
  carrying_capacity <- 40 # clade-level
  if (is_testing) {
    carrying_capacity <- 4
  }
  dd_parameters <- c(speciation_rate, extinction_rate, carrying_capacity)
  ddmodel <- 1 # linear dependence in speciation rate with parameter K
  set.seed(i)
  dd_sim_result <- DDD::dd_sim(pars = dd_parameters, age  = crown_age, ddmodel = ddmodel)
  phylogeny <- dd_sim_result$tes # Only extant species
  phylogenies[[i]] <- phylogeny
}

################################################################################
# Create piriouette parameter sets
################################################################################
pir_paramses <- list()
for (i in seq_len(n_phylogenies)) {

  print(paste(i, "/", n_phylogenies))

  alignment_params <- create_alignment_params(
    sim_tral_fun = get_sim_tral_with_std_nsm_fun(
      mutation_rate = 1.0 / crown_age,
      site_model = beautier::create_jc69_site_model()
    ),
    root_sequence = create_blocked_dna(length = 1000),
    rng_seed = rng_seed,
    fasta_filename = paste0("true_alignment_", i, ".fas")
  )

  # JC69, strict, Yule
  generative_experiment <- create_gen_experiment()
  generative_experiment$beast2_options$input_filename <- paste0("true_alignment_gen_", i ,".xml")
  generative_experiment$beast2_options$output_state_filename <- paste0("true_alignment_gen_", i ,".xml.state")
  generative_experiment$inference_model$mcmc$tracelog$filename <- paste0("true_alignment_gen_", i ,".log")
  generative_experiment$inference_model$mcmc$treelog$filename <- paste0("true_alignment_gen_", i ,".trees")
  generative_experiment$inference_model$mcmc$screenlog$filename <- paste0("true_alignment_gen_", i ,".csv")
  generative_experiment$errors_filename <- paste0("true_errors_gen_", i ,".csv")

  # Only pure-birth and birth-death models
  nee_tree_priors <- list()
  nee_tree_priors[[1]] <- create_yule_tree_prior()
  nee_tree_priors[[2]] <- create_bd_tree_prior()
  strict_clock_models <- list()
  strict_clock_models[[1]] <- create_strict_clock_model()

  candidate_experiments <- create_all_experiments(
    clock_models = strict_clock_models,
    tree_priors = nee_tree_priors,
    exclude_model = generative_experiment$inference_model
  )

  for (j in seq_along(candidate_experiments)) {
    candidate_experiments[[j]]$beast2_options$input_filename <- paste0("true_alignment_best_", i ,".xml")
    candidate_experiments[[j]]$beast2_options$output_state_filename <- paste0("true_alignment_best_", i ,".xml.state")
    candidate_experiments[[j]]$inference_model$mcmc$tracelog$filename <- paste0("true_alignment_best_", i ,".log")
    candidate_experiments[[j]]$inference_model$mcmc$treelog$filename <- paste0("true_alignment_best_", i ,".trees")
    candidate_experiments[[j]]$inference_model$mcmc$screenlog$filename <- paste0("true_alignment_best_", i ,".csv")
    candidate_experiments[[j]]$errors_filename <- paste0("true_errors_best_", i ,".csv")
  }
  experiments <- c(list(generative_experiment), candidate_experiments)

  # Set the RNG seed
  for (j in seq_along(experiments)) {
    experiments[[j]]$beast2_options$rng_seed <- rng_seed
  }

  # Setup estimation of evidences (aka marginal likelihoods)
  for (j in seq_along(experiments)) {
    experiments[[j]]$est_evidence_mcmc <- create_mcmc_nested_sampling(
      chain_length = 1e7,
      store_every = 1e3,
      epsilon = 1e-12
    )
  }

  # Shorter on Travis
  if (is_testing) {
    for (j in seq_along(experiments)) {
      experiments[[j]]$inference_model$mcmc$chain_length <- 3000
      experiments[[j]]$inference_model$mcmc$store_every <- 1000
      experiments[[j]]$est_evidence_mcmc$chain_length <- 3000
      experiments[[j]]$est_evidence_mcmc$store_every <- 1000
      experiments[[j]]$est_evidence_mcmc$epsilon <- 100.0
    }
  }

  twinning_params <- create_twinning_params(
    rng_seed_twin_tree = rng_seed,
    sim_twin_tree_fun = get_sim_bd_twin_tree_fun(),
    rng_seed_twin_alignment = rng_seed,
    sim_twal_fun = get_sim_twal_with_std_nsm_fun(
      mutation_rate = pirouette::create_standard_mutation_rate(
        phylogeny
      )
    ),
    twin_tree_filename = paste0("twin_tree_", i ,".newick"),
    twin_alignment_filename = paste0("twin_alignment_", i ,".fas"),
    twin_evidence_filename = paste0("twin_evidence_", i ,".csv")
  )

  error_measure_params <- pirouette::create_error_measure_params(
    error_fun = pirouette::get_gamma_error_fun()
  )

  pir_params <- create_pir_params(
    alignment_params = alignment_params,
    experiments = experiments,
    twinning_params = twinning_params,
    error_measure_params = error_measure_params
  )

  pir_paramses[[i]] <- pir_params
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

utils::write.csv(
  x = pir_outs,
  file = file.path(example_folder, "errors.csv"),
  row.names = FALSE
)

pir_plots(pir_outs) +# pirouette example 28:
# multiple DD trees for pirouette article
suppressMessages(library(ggplot2))
library(pirouette)
library(babette)

root_folder <- getwd()
example_no <- 28
rng_seed <- 314
example_folder <- file.path(root_folder, paste0("example_", example_no, "_", rng_seed))
dir.create(example_folder, showWarnings = FALSE, recursive = TRUE)
setwd(example_folder)
set.seed(rng_seed)
testit::assert(is_beast2_installed())

n_phylogenies <- 1
phylogenies <- list()

################################################################################
# Creates phylogenies
################################################################################
for (i in seq_len(n_phylogenies)) {
  print(paste(i, "/", n_phylogenies))
  speciation_rate <- 0.8 # lambda
  extinction_rate <- 0.1 # mu
  carrying_capacity <- 40 # clade-level
  crown_age <- 10
  dd_parameters <- c(speciation_rate, extinction_rate, carrying_capacity)
  ddmodel <- 1 # linear dependence in speciation rate with parameter K
  set.seed(i)
  dd_sim_result <- DDD::dd_sim(pars = dd_parameters, age  = crown_age, ddmodel = ddmodel)
  phylogeny <- dd_sim_result$tes # Only extant species
  phylogenies[[i]] <- phylogeny

}

################################################################################
# Create piriouette parameter sets
################################################################################
pir_paramses <- list()
for (i in seq_len(n_phylogenies)) {

  print(paste(i, "/", n_phylogenies))

  alignment_params <- create_alignment_params(
    sim_tral_fun = get_sim_tral_with_std_nsm_fun(
      mutation_rate = 0.1,
      site_model = beautier::create_jc69_site_model()
    ),
    root_sequence = create_blocked_dna(length = 1000),
    rng_seed = rng_seed,
    fasta_filename = paste0("true_alignment_", i, ".fas")
  )

  # JC69, strict, Yule
  generative_experiment <- create_gen_experiment()
  generative_experiment$beast2_options$input_filename <- paste0("true_alignment_gen_", i ,".xml")
  generative_experiment$beast2_options$output_state_filename <- paste0("true_alignment_gen_", i ,".xml.state")
  generative_experiment$inference_model$mcmc$tracelog$filename <- paste0("true_alignment_gen_", i ,".log")
  generative_experiment$inference_model$mcmc$treelog$filename <- paste0("true_alignment_gen_", i ,".trees")
  generative_experiment$inference_model$mcmc$screenlog$filename <- paste0("true_alignment_gen_", i ,".csv")
  generative_experiment$errors_filename <- paste0("true_errors_gen_", i ,".csv")

  site_models <- create_site_models()
  clock_models <- create_clock_models()
  if (is_testing) {
    site_models <- list()
    site_models[[1]] <- create_hky_site_model()
    clock_models <- list()
    clock_models[[1]] <- create_strict_clock_model()
  }
  # Only pure-birth and birth-death models
  tree_priors <- list()
  tree_priors[[1]] <- create_yule_tree_prior()
  tree_priors[[2]] <- create_bd_tree_prior()

  candidate_experiments <- create_all_experiments(
    site_models = site_models,
    clock_models = clock_models,
    tree_priors = tree_priors,
    exclude_model = generative_experiment$inference_model
  )

  for (j in seq_along(candidate_experiments)) {
    candidate_experiments[[j]]$beast2_options$input_filename <- paste0("true_alignment_best_", i ,".xml")
    candidate_experiments[[j]]$beast2_options$output_state_filename <- paste0("true_alignment_best_", i ,".xml.state")
    candidate_experiments[[j]]$inference_model$mcmc$tracelog$filename <- paste0("true_alignment_best_", i ,".log")
    candidate_experiments[[j]]$inference_model$mcmc$treelog$filename <- paste0("true_alignment_best_", i ,".trees")
    candidate_experiments[[j]]$inference_model$mcmc$screenlog$filename <- paste0("true_alignment_best_", i ,".csv")
    candidate_experiments[[j]]$errors_filename <- paste0("true_errors_best_", i ,".csv")
  }
  experiments <- c(list(generative_experiment), candidate_experiments)

  # Set the RNG seed
  for (j in seq_along(experiments)) {
    experiments[[j]]$beast2_options$rng_seed <- rng_seed
  }

  # Setup estimation of evidences (aka marginal likelihoods)
  for (j in seq_along(experiments)) {
    experiments[[j]]$est_evidence_mcmc <- create_mcmc_nested_sampling(
      chain_length = 1e7,
      store_every = 1e3,
      epsilon = 1e-12
    )
  }

  # Shorter on Travis
  if (is_testing) {
    for (j in seq_along(experiments)) {
      experiments[[j]]$inference_model$mcmc$chain_length <- 3000
      experiments[[j]]$inference_model$mcmc$store_every <- 1000
      experiments[[j]]$est_evidence_mcmc$chain_length <- 3000
      experiments[[j]]$est_evidence_mcmc$store_every <- 1000
      experiments[[j]]$est_evidence_mcmc$epsilon <- 100.0
    }
  }

  twinning_params <- create_twinning_params(
    rng_seed_twin_tree = rng_seed,
    sim_twin_tree_fun = get_sim_bd_twin_tree_fun(),
    rng_seed_twin_alignment = rng_seed,
    sim_twal_fun = get_sim_twal_with_std_nsm_fun(
      mutation_rate = pirouette::create_standard_mutation_rate(
        phylogeny
      )
    ),
    twin_tree_filename = paste0("twin_tree_", i ,".newick"),
    twin_alignment_filename = paste0("twin_alignment_", i ,".fas"),
    twin_evidence_filename = paste0("twin_evidence_", i ,".csv")
  )

  error_measure_params <- pirouette::create_error_measure_params(
    error_fun = pirouette::get_gamma_error_fun()
  )

  pir_params <- create_pir_params(
    alignment_params = alignment_params,
    experiments = experiments,
    twinning_params = twinning_params,
    error_measure_params = error_measure_params
  )

  pir_paramses[[i]] <- pir_params
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

utils::write.csv(
  x = pir_outs,
  file = file.path(example_folder, "errors.csv"),
  row.names = FALSE
)

pir_plots(pir_outs) +
  ggsave(file.path(example_folder, "errors.png"))

if (1 == 2) {
  pir_to_pics(
    phylogeny = phylogeny,
    pir_params = pir_params,
    folder = example_folder
  )

  pir_to_tables(
    pir_params = pir_params,
    folder = example_folder
  )
}
