# Endometrial receptivity profiling of Extracellular Vesicules

This study builds on previous work done by Vigano et al 2021 and by Meltsov et al. 2023 to move towards a non-invasive method for receptivity testing. Best efforts were made to try to organise things in the following way:

* Scripts to run nf-core/rnaseq preprocessing pipeline on the raw samples in the SLURM HPC.

* R scripts to take the read count matrices emitted by nf-core in `raw_data/` and output phenotype files after formatting to `data/`.

* RMarkdown files to generate interactive analysis reports in `analysis/`.

Environment manager of choice was `micromamba`, rebuild `renv-ev` environment with (TODO: the micromamba command).

# Citation
Parts of this work were published in:

Meltsov, A., Giacomini, E., Vigano, P., Zarovni, N., Salumets, A., Aleksejeva, E., 2023. Letter to the Editor - Pilot proof for RNA biomarker-based minimally invasive endometrial receptivity testing using uterine fluid extracellular vesicles. European Journal of Obstetrics and Gynecology and Reproductive Biology 287, 237–238. https://doi.org/10.1016/j.ejogrb.2023.06.006
