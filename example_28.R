# pirouette example 28:
# multiple DD trees for pirouette article
suppressMessages(library(ggplot2))
suppressMessages(library(pirouette))
library(babette)

is_testing <- FALSE
if (is_on_travis()) {
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

n_phylogenies <- 5
phylogenies <- list()

################################################################################
# Creates phylogenies
################################################################################
crown_age <- 10
for (i in seq_len(n_phylogenies)) {
  print(paste(i, "/", n_phylogenies))
  # Use same parameters as create_dd_tree
  extinction_rate <- 0.1
  n_taxa <- 6
  n_0 <- 2 # Initial number of species at stem/crown of tree
  diff <- (log(n_taxa) - log(n_0)) / crown_age
  speciation_rate <- 3.0 * (diff + extinction_rate)
  carrying_capacity <- n_taxa # clade-level
  dd_parameters <- c(speciation_rate, extinction_rate, carrying_capacity)
  ddmodel <- 1 # linear dependence in speciation rate with parameter K
  set.seed(i)
  dd_sim_result <- DDD::dd_sim(pars = dd_parameters, age  = crown_age, ddmodel = ddmodel)
  phylogeny <- dd_sim_result$tes # Only extant species

  # Save tree to files
  ape::write.tree(phylogeny, file = file.path(example_folder, "true_tree.newick"))
  png(filename = file.path(example_folder, "true_tree.png"), width = 7, height = 7)
  ape::plot.phylo(phylogeny)
  dev.off()

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
      site_model = create_jc69_site_model()
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

  # Create the set of candidate birth-death experiments
  candidate_experiments <- create_all_bd_experiments(
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
    error_fun = pirouette::get_nltt_error_fun()
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

