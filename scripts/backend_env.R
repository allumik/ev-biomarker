#### Backend scripts for the cERtain project
##

## load the dependencies
source("./scripts/load_deps.R")
args <- commandArgs(trailingOnly = TRUE)

#### Load the functions and settings
source("./scripts/backend_env_functions.R")

#### Load raw data and format
source("./scripts/format_raw_data.R")


## load in TED genes
ted_genes <-
  read_tsv("./genesets/ted_diaz-gimeno_2021.tsv") %>%
  filter(TED == 1) %$%
  `Gene Name` %>%
  str_trim

## load in the ERA geneset
era_genes <-
  read_tsv("./genesets/era_238genes.tsv") %$%
  `Gene symbol` %>%
  str_trim

## load in the full beready (Altmäe) geneset
beready_genes <-
  tacseqApp::targets_data$`All READY 72 targets` %>%
  filter(type == "biomarker") %$%
  target %>%
  str_trim

finders <-
  tacseqApp::targets_data$`All READY 72 targets` %>%
  filter(type == "housekeeper") %$%
  target %>%
  str_trim

## the Vigano set housekeepers
finders_alt <- c("SNRPG", "OST4", "TOMM7", "NOP10")

## load in the core set of beready genes signf different in
## pre and post groups when compared to rec
beready_genes_core <-
  tacseqApp::targets_data$`Dynamic filtering geneset` %>%
  ## remove some weird columns)
  select(-c("pro", "prepost")) %>%
  filter(if_all(where(is_double), ~ . == 1)) %$%
  target %>%
  str_trim


## Bind the genesets into a frequency table
all_genes <-
  reduce(
    list(
      ted_genes,
      era_genes,
      beready_genes,
      beready_genes_core
      ),
    union
  )

sig_set <-
  bind_cols(
    "genenames" = all_genes,
    "TED" = all_genes %in% ted_genes,
    "ERA" = all_genes %in% era_genes,
    "beREADY" = all_genes %in% beready_genes,
    "beREADY core" = all_genes %in% beready_genes_core
  )


#### Bigger analysis steps that need to be prerun
global_expr_pca <-
  expr_tpm %>%
    bind_rows(expr_tpm_test) %>%
    recipe() %>%
    update_role(samplename, new_role = "id") %>%
    update_role(c("group", "cyclephase"), new_role = "other") %>%
    ## eliminate some genes with zero expression
    step_zv(all_numeric()) %>%
    step_normalize(all_numeric()) %>%
    step_pca(all_numeric(), num_comp = 3) %>%
    prep(strings_as_factors = F) %>%
    bake(new_data = expr_tpm %>% bind_rows(expr_tpm_test)) 


global_expr_umap <-
  expr_tpm %>%
    bind_rows(expr_tpm_test) %>%
    recipe() %>%
    update_role(samplename, new_role = "id") %>%
    update_role(c("group", "cyclephase"), new_role = "other") %>%
    ## eliminate some genes with zero expression
    step_zv(all_numeric()) %>%
    step_normalize(all_numeric()) %>%
    step_umap(all_numeric(), num_comp = 2) %>%
    prep(strings_as_factors = F) %>%
    bake(new_data = expr_tpm %>% bind_rows(expr_tpm_test))
