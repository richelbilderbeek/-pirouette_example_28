suppressMessages(library(ggplot2))
library(pirouette)
library(babette)

root_folder <- getwd()
example_no <- 2
rng_seed <- 314
example_folder <- file.path(root_folder, paste0("example_", example_no, "_", rng_seed))
dir.create(example_folder, showWarnings = FALSE, recursive = TRUE)
setwd(example_folder)
set.seed(rng_seed)
testit::assert(is_beast2_installed())

phylogeny <- create_yule_tree(n_taxa = 6, crown_age = 10)
ape::write.tree(phylogeny, file = "true_tree.newick")

alignment_params <- create_alignment_params(
  sim_tral_fun = get_sim_tral_with_std_nsm_fun(
    mutation_rate = 0.1,
    site_model = beautier::create_jc69_site_model()
  ),
  root_sequence = create_blocked_dna(length = 1000),
  rng_seed = rng_seed,
  fasta_filename = "true_alignment.fas"
)

# JC69, strict, Yule
generative_experiment <- create_gen_experiment()
generative_experiment$beast2_options$input_filename <- "true_alignment_gen.xml"
generative_experiment$beast2_options$output_state_filename <- "true_alignment_gen.xml.state"
generative_experiment$inference_model$mcmc$tracelog$filename <- "true_alignment_gen.log"
generative_experiment$inference_model$mcmc$treelog$filename <- "true_alignment_gen.trees"
generative_experiment$inference_model$mcmc$screenlog$filename <- "true_alignment_gen.csv"
generative_experiment$errors_filename <- "true_errors_gen.csv"
check_experiment(generative_experiment)

# All non-Yule tree priors
candidate_experiments <- create_all_experiments(
  exclude_model = generative_experiment$inference_model
)
for (i in seq_along(candidate_experiments)) {
  candidate_experiments[[i]]$beast2_options$input_filename <- "true_alignment_best.xml"
  candidate_experiments[[i]]$beast2_options$output_state_filename <- "true_alignment_best.xml.state"
  candidate_experiments[[i]]$inference_model$mcmc$tracelog$filename <- "true_alignment_best.log"
  candidate_experiments[[i]]$inference_model$mcmc$treelog$filename <- "true_alignment_best.trees"
  candidate_experiments[[i]]$inference_model$mcmc$screenlog$filename <- "true_alignment_best.csv"
  candidate_experiments[[i]]$errors_filename <- "true_errors_best.csv"
}
check_experiments(candidate_experiments)

experiments <- c(list(generative_experiment), candidate_experiments)

# Set the RNG seed
for (i in seq_along(experiments)) {
  experiments[[i]]$beast2_options$rng_seed <- rng_seed
}

# Setup estimation of evidences (aka marginal likelihoods)
for (i in seq_along(experiments)) {
  experiments[[i]]$est_evidence_mcmc <- create_mcmc_nested_sampling(
    chain_length = 1e7,
    store_every = 1e3,
    epsilon = 1e-12
  )
}

check_experiments(experiments)

# Shorter on Travis
if (is_on_travis()) {
  for (i in seq_along(experiments)) {
    experiments[[i]]$inference_model$mcmc$chain_length <- 3000
    experiments[[i]]$inference_model$mcmc$store_every <- 1000
    experiments[[i]]$est_evidence_mcmc$chain_length <- 3000
    experiments[[i]]$est_evidence_mcmc$store_every <- 1000
    experiments[[i]]$est_evidence_mcmc$epsilon <- 100.0
  }
}

pir_params <- create_pir_params(
  alignment_params = alignment_params,
  experiments = experiments
)

rm_pir_param_files(pir_params)

errors <- pir_run(
  phylogeny,
  pir_params = pir_params
)

utils::write.csv(
  x = errors,
  file = file.path(example_folder, "errors.csv"),
  row.names = FALSE
)

pir_plot(errors) +
  ggsave(file.path(example_folder, "errors.png"))

pir_to_pics(
  phylogeny = phylogeny,
  pir_params = pir_params,
  folder = example_folder
)

pir_to_tables(
  pir_params = pir_params,
  folder = example_folder
)
