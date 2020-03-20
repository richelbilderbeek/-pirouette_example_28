#!/bin/bash
#
# Re-run the code locally, to re-create the data and figure.
#
# Usage:
#
#   ./scripts/rerun.sh
#
#SBATCH --partition=gelifes
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=1
#SBATCH --mem=10G
#SBATCH --job-name=pirex28
#SBATCH --output=example_28.log
#
rm -rf example_28
rm *.png
time Rscript example_28.R
zip -r pirouette_example_28.zip example_28 example_28.R scripts *.png

