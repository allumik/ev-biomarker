#### Script to format the raw data for the analysis
# load dependencies
source("./scripts/load_deps.R")
args <- commandArgs(trailingOnly = TRUE)
print(args)

#### Load the functions and settings
source("./scripts/backend_env_functions.R")

## read in the count matrices and fix naming order
raw_mat <-
  read_tsv("./data_raw/rsem_matrices/rsem.merged.gene_counts.tsv") %>%
  full_join(read_tsv("./data_raw/rsem_matrices/rsem.merged.gene_counts_add1.tsv")) %>%
    # set NA's as 0's by default and enforce integers (for compact data format)
    mutate(across(where(is.numeric), ~ as.integer(replace_na(.x, 0)))) %>%
    # remove the corresponding transcript ID's and the Undetermined sample
    select(-c("transcript_id(s)")) %>%
    # fix the sample names
    rename_with(
      ~ .x %>%
      str_replace("HU10", "HUT10") %>%
      str_remove("(?<=\\HUT)[0]+(?=\\d)")
    ) %>%
    # fix missing sample names from the added samples
    rename_with( ~ str_c(.x, "_UF"), matches("HUT9|HUT25|HUT29|HUT32|HUT27"))

raw_tpm <-
  read_tsv("./data_raw/rsem_matrices/rsem.merged.gene_tpm.tsv") %>%
  full_join(read_tsv("./data_raw/rsem_matrices/rsem.merged.gene_tpm_add1.tsv")) %>%
    # set NA's as 0's by default
    mutate(across(where(is.numeric), ~ replace_na(.x, 0))) %>%
    # remove the corresponding transcript ID's and the Undetermined sample
    select(-c("transcript_id(s)")) %>%
    rename_with(
      ~ .x %>%
      str_replace("HU10", "HUT10") %>%
      str_remove("(?<=\\HUT)[0]+(?=\\d)")
    ) %>%
    # fix missing sample names from the added samples
    rename_with( ~ str_c(.x, "_UF"), matches("HUT9|HUT25|HUT29|HUT32|HUT27"))

## read in the metadata and consolidate samplenames
pheno_all <-
  read_tsv("./data_raw/mens_phase.tsv") %>%
  # apparently some samples in this table are not present in the expression data
  mutate(
    # remove leading zeros from samplenames; already did it in data
    # samplename = str_remove(samplename, "0+(?!$)"),
    # add a group column based on eithe UF or biopsy
    group = ifelse(str_detect(samplename, "UF"), "UF", "biopsy"),
    # rename groups to short group standards in this analysis
    cyclephase = case_when(
      cyclephase == "Proliferative" ~ "pro",
      cyclephase == "LH2_3" ~ "pre",
      cyclephase == "LH7_8" ~ "rec",
      cyclephase == "LH11_13" ~ "post",
      T ~ cyclephase
    ),
    dataset = ifelse(
      str_detect(samplename, "^HUT") &
        !str_detect(samplename, "HUT9|HUT25|HUT29|HUT32|HUT27"),
      "train",
      "test")
  ) %>%
  {.[match(colnames(raw_mat)[-1], .$samplename),]}


# switch and filter some samples if needed
if(any(str_detect(names(args), "switch"))) {
  raw_mat %<>%
    rename(
      ## this takes too much time to figure out how to switch programmatically
      HUT23_biopsy = HUT23_UF,
      HUT23_UF = HUT23_biopsy
    )
  raw_tpm %<>%
    rename(
      HUT23_biopsy = HUT23_UF,
      HUT23_UF = HUT23_biopsy
    )
}

if(any(str_detect(names(args), "remove"))) {
  raw_mat %<>%
    select(-args$remove_sample[args$remove_sample %in% colnames(raw_mat)])
  raw_tpm %<>%
    select(-args$remove_sample[args$remove_sample %in% colnames(raw_tpm)])
}

pheno <-
  pheno_all %>%
  # but remove an samples actually not in the raw_mat (supposedly 8 of them)
  # and also quarantee that the samplenames are in the same order
  # as raw_mat colnames (needed for DESeq2)
  {.[match(colnames(raw_mat)[-1], .$samplename),]} %>%
  # make the training dataset the default one for pheno
  filter(dataset == "train")

pheno_test <-
  pheno_all %>%
  {.[match(colnames(raw_mat)[-1], .$samplename),]} %>%
  filter(dataset == "test")


## join count matrix with pheno for visualization or downstream analysis
expr_mat <- raw_mat %>% transform_count_mat(pheno)
expr_mat_test <- raw_mat %>% transform_count_mat(pheno_test)
expr_tpm <- raw_tpm %>% transform_count_mat(pheno)
expr_tpm_test <- raw_tpm %>% transform_count_mat(pheno_test)


## divide the dataset into train and test
raw_mat_test <-
  raw_mat %>%
    select(c("gene_id", pheno_test$samplename[pheno_test$samplename %in% colnames(raw_mat)]))
raw_tpm_test <-
  raw_tpm %>%
    select(c("gene_id", pheno_test$samplename[pheno_test$samplename %in% colnames(raw_tpm)]))
raw_mat <-
  raw_mat %>%
    select(c("gene_id", pheno$samplename[pheno$samplename %in% colnames(raw_mat)]))
raw_tpm <-
  raw_tpm %>%
    select(c("gene_id", pheno$samplename[pheno$samplename %in% colnames(raw_tpm)]))



#### Emit formatted data files
if(!any(str_detect(names(args), "subfolder")))
  data_subfolder <- "./data/"

raw_mat %>% write_feather(paste0(data_subfolder, "counts_raw.feather"))
raw_mat_test %>% write_feather(paste0(data_subfolder, "counts_raw_test.feather"))
raw_tpm %>% write_feather(paste0(data_subfolder, "tpm_raw.feather"))
raw_tpm_test %>% write_feather(paste0(data_subfolder, "tpm_raw_test.feather"))

expr_mat %>% write_feather(paste0(data_subfolder, "counts_format.feather"))
expr_mat_test %>% write_feather(paste0(data_subfolder, "counts_format_test.feather"))
expr_tpm %>% write_feather(paste0(data_subfolder, "tpm_format.feather"))
expr_tpm_test %>% write_feather(paste0(data_subfolder, "tpm_format_test.feather"))

pheno %>% write_tsv(paste0(data_subfolder, "phenotype.tsv"))
pheno_test %>% write_tsv(paste0(data_subfolder, "phenotype_test.tsv"))


## also write out raw file with switched annotation
raw_mat %>%
  rename(ensembl_gene_id = gene_id) %>%
  geneid_converter(annot) %>%
  rename(gene_id = external_gene_name) %>%
  write_feather(paste0(data_subfolder, "annot_raw.feather"))

