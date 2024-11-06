#!/bin/bash
#SBATCH -J render_document_gpu
#SBATCH --partition=gpu
#SBATCH --time=18:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alvinmeltsov@gmail.com
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8

## Use a Apptainer container environment for newer dependencies.
module load singularity

module load cuda12.0
# micromamba activate sc-env
export CUDA_VISIBLE_DEVICES=0

# micromamba run -n sc-base python3 ../preproc_scripts/scvi_integrate.py
singularity exec --nv \
  --bind /ifs/data/research/CGM/Projects/sc_atlas/:/ifs/data/research/CGM/Projects/sc_atlas/ \
  $HOME/.containers/base_env.sif \
  ~/.local/bin/micromamba run -n sc-base quarto render analysis/ev_correlation.qmd --to html
