# pirouette example 28:
# multiple DD trees for pirouette article
#
library(pirouette)
library(beautier)
library(beastier)
library(testthat)
library(ggplot2)

# Constants
example_no <- 28
rng_seed <- 314
crown_age <- 10
n_phylogenies <- 100
folder_name <- paste0("example_", example_no)
is_testing <- is_on_ci()
if (is_testing) {
  n_phylogenies <- 2
}

# Create phylogenies
phylogenies <- list()
for (i in seq_len(n_phylogenies)) {
  set.seed(314 - 1 + i)
  phylogenies[[i]] <- create_exemplary_dd_tree(n_taxa = n_taxa, crown_age = crown_age)
}
expect_equal(length(phylogenies), n_phylogenies)

# Create pirouette parameter sets
pir_paramses <- create_std_pir_paramses(
  n = length(phylogenies),
  folder_name = folder_name
)
expect_equal(length(pir_paramses), n_phylogenies)
if (is_testing) {
    pir_paramses <- shorten_pir_paramses(pir_paramses)
}

# Do the runs
pir_outs <- pir_runs(
  phylogenies = phylogenies,
  pir_paramses = pir_paramses
)

# Save summary
pir_plots(pir_outs) +
  ggtitle(paste("Number of replicates: ", n_phylogenies)) +
  ggsave("errors.png", width = 7, height = 7)

# Save individual runs
expect_equal(length(pir_paramses), length(pir_outs))
expect_equal(length(pir_paramses), length(phylogenies))
for (i in seq_along(pir_outs)) {
  pir_save(
    phylogeny = phylogenies[[i]],
    pir_params = pir_paramses[[i]],
    pir_out = pir_outs[[i]],
    folder_name = dirname(pir_paramses[[i]]$alignment_params$fasta_filename)
  )
}
