#### This script is to gather together the clinical samples for EV evaluation.
##

source("./scripts/load_deps.R")

## read in the count matrices
clin_mat <-
  read_tsv("./data_raw/rsem_matrices/rsem.merged.gene_counts_add1.tsv") %>%
  full_join(read_tsv("./data_raw/rsem_matrices/rsem.merged.gene_counts_add2.tsv")) %>%
  full_join(read_tsv("./data_raw/rsem_matrices/rsem.merged.gene_counts_add3.tsv")) %>%
    # remove a line from the end of some smaples
    rename_with(~ str_replace_all(.x, "_L001", "")) %>%
    # set NA's as 0's by default and enforce integers (for more compact data format)
    mutate(across(where(is.numeric), ~ as.integer(replace_na(.x, 0)))) %>%
    # remove the corresponding transcript ID's and the Undetermined sample
    select(-c("transcript_id(s)", "Undetermined")) %>%
    ## remove HUT testing set samples
    select(-contains("HUT")) %>%
    ## remove some extra samples with ambigious outcome
    select(-matches("4ZAM|2ZAM|3SLC|1SLC"))

clin_tpm <-
  read_tsv("./data_raw/rsem_matrices/rsem.merged.gene_tpm_add1.tsv") %>%
  full_join(read_tsv("./data_raw/rsem_matrices/rsem.merged.gene_tpm_add2.tsv")) %>%
  full_join(read_tsv("./data_raw/rsem_matrices/rsem.merged.gene_tpm_add3.tsv")) %>%
    # remove a line from the end of some smaples
    rename_with(~ str_replace_all(.x, "_L001", "")) %>%
    # set NA's as 0's by default
    mutate(across(where(is.numeric), ~ replace_na(.x, 0))) %>%
    # remove the corresponding transcript ID's and the Undetermined sample
    select(-c("transcript_id(s)", "Undetermined")) %>%
    ## remove HUT testing set samples
    select(-contains("HUT")) %>%
    ## remove some extra samples with ambigious outcome
    select(-matches("4ZAM|2ZAM|3SLC|1SLC"))

## read in the metadata and consolidate samplenames
pheno_clin <-
  read_tsv("./data_raw/pheno_clin.tsv") %>%
  # apparently some samples in this table are not present in the expression data
  { .[match(colnames(clin_mat)[-1], .$samplename),] }


## write out
data_subfolder <- "./data/"

clin_mat %>% write_feather(paste0(data_subfolder, "clin_counts_raw.feather"))
clin_tpm %>% write_feather(paste0(data_subfolder, "clin_tpm_raw.feather"))
pheno_clin %>% write_tsv(paste0(data_subfolder, "phenotype_clin.tsv"))

