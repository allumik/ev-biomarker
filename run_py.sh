#!/bin/bash
#SBATCH -J render_document
#SBATCH --partition=research
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alvinmeltsov@gmail.com
#SBATCH --mem=128G
#SBATCH --cpus-per-task=16

# load virtualisation environment
module load singularity

singularity exec --nv $HOME/.containers/base_env.sif \
  ~/.local/bin/micromamba run -n sc-base quarto render ./analysis/sc_trajectory.qmd --to html
# singularity exec --nv $HOME/.containers/base_env.sif \
  # ~/.local/bin/micromamba run -n sc-base python3 ../preproc_scripts/raw_data_load.py
