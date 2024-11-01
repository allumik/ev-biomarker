#!/bin/bash
#SBATCH -J run_ev_analysis
#SBATCH --partition=research
# #SBATCH --time=18:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alvinmeltsov@gmail.com
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8


## Some examples
# singularity exec ../images/r-bioverse-dev-fat.sif r -e "rmarkdown::render('./huter_analysis.rmd', 'html_document')"
# ANNDATA_FOLDER="~/proj/huter_analysis/working_data/anndatas/" micromamba run -n sc-base Rscript ../preproc_scripts/stacas_run.R
# micromamba run -n sc-base Rscript ../preproc_scripts/loom_stacas.R

module load singularity/3.8.5
singularity exec --nv $HOME/.containers/base_env.sif \
  ~/.local/bin/micromamba run -n renv-ev \
  R -e "rmarkdown::render('./analysis/ev_clinical.rmd', knit_root_dir = '~/proj/ev-biomarker/')"
# singularity exec --nv $HOME/.containers/base_env.sif \
  # ~/.local/bin/micromamba run -n sc-base quarto render ./analysis/ev_clinical.qmd --to html

