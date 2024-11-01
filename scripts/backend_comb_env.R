#### R script for reading in Vigano/NOTED/HUT datasets
## and harmonising the dataset for further downstream analysis

## 1. Load datasets and unify metadata
## 2. Join them and apply batch normalisation with sva-combatch
## 3. Look at the clustering of the joined dataset before and after sva

source("./scripts/load_deps.R")
# load the filtered set of samples and previous functions from backend
# this datafile is generated with diff_run.R in ev_biomarker.rmd
source("./scripts/backend_env_functions.R")
raw_mat <- read_feather("./data/filtered/counts_raw.feather")
raw_tpm <- read_feather("./data/filtered/tpm_raw.feather")
pheno <- read_tsv("./data/filtered/phenotype.tsv")
pheno_test <- read_tsv("./data/filtered/phenotype_test.tsv")

hut_mat <-
  raw_mat %>%
    rename(ensembl_gene_id = gene_id) %>%
    geneid_converter(annot) %>%
    # set NA's as 0's by default and enforce integers (for compact data format)
    mutate(across(where(is.numeric), ~ replace_na(.x, 0))) %>%
    filter(!is.na(external_gene_name))
rm(raw_mat)

hut_tpm <-
  raw_tpm %>%
    rename(ensembl_gene_id = gene_id) %>%
    geneid_converter(annot) %>%
    # set NA's as 0's by default and enforce integers (for compact data format)
    mutate(across(where(is.numeric), ~ replace_na(.x, 0))) %>%
    filter(!is.na(external_gene_name))
rm(raw_tpm)

# and rename the pheno for clear marking of the data origin
hut_pheno <- pheno
rm(pheno)
hut_test_pheno <- pheno_test
rm(pheno_test)



## load in the Vigano paired dataset
# ev_comb <-
  # read_excel("./data/endometrium.xlsx", range = "A5:U56475")
# the colnames need some fixing
# colnames(ev_comb) <- str_trim(colnames(ev_comb))

# Also load in the Vigano ev dataset
ev_comb <-
  read_excel("./data_raw/EV_RNAseq_Vigano.xlsx", range = "D3:AB54784") %>%
  # inner_join(ev_comb, by = "GeneID") %>%
  rename(external_gene_name = GeneID)

# combine the phenotypes
ev_comb_pheno <-
  read_excel("./data_raw/EV_RNAseq_Vigano.xlsx", range = "E3:AB3") %>% {
    tibble(
      samplename = colnames(.),
      cyclephase = c(rep("pre", 12), rep("rec", 12)),
      grouping = "UF"
    )
  }

# ev_comb_pheno <-
  # read_excel("./data/endometrium.xlsx", range = "B5:U5") %>%
  # colnames %>%
  # str_trim %>%
  # tibble(samplename = .) %>%
    # mutate(
      # cyclephase = case_when(
        # str_detect(samplename, "26|27|32|33|35") ~ "pro",
        # str_detect(samplename, "22|25|30|31|34") ~ "rec",
      # ),
      # grouping = ifelse(str_detect(samplename, "UF"), "UF", "biopsy")
    # ) %>%
  # bind_rows(ev_comb_pheno)



# ## load in the NOTED dataset
# # from GSE98386
# noted_set <-
  # read_tsv("./data/GSE98386_Raw_count_matrix.tsv.gz") %>%
  # rename(ensembl_gene_id = ensemblID) %>%
  # geneid_converter(annot) %>%
  # mutate(across(where(is.numeric), ~ replace_na(.x, 0))) %>%
  # filter(!is.na(external_gene_name)) %>%
  # distinct(external_gene_name, .keep_all = T)
# 
# batch1 <- c(
  # "ERTK11", "NOTNV09", "NOTNV16", "NOTNV18", "NOTNV24",
  # "NOTNV25", "NOTNV40", "NOTNV47", "NOTNV48", "NOTNV75"
# )
# batch2 <- c(
  # "NOT1210", "NOT1213", "NOTNV01", "NOTNV02", "NOTNV03",
  # "NOTNV04", "NOTNV11", "NOTNV12", "NOTNV13", "NOTNV15"
# )
# 
# noted_pheno <-
  # colnames(noted_set)[-1] %>% {
    # tibble(
      # samplename = .,
      # cyclephase = case_when(
        # str_detect(., "LH2") ~ "pre",
        # str_detect(., "LH8") ~ "rec"
      # ),
      # batch = case_when(
        # str_detect(., paste(batch1, collapse = "|")) ~ "batch_1",
        # str_detect(., paste(batch2, collapse = "|")) ~ "batch_2"
      # )
    # )
  # }
# 
# ## noted samples need some batch correction from the get-go
# ## apply poisson-gauss sva before tpm normalisation, so that the same data could be used
# noted_set %<>%
  # column_to_rownames("external_gene_name") %>%
  # as.matrix %>%
  # sva::ComBat_seq(
    # batch = noted_pheno$batch, # the batch grouping
    # group = noted_pheno$cyclephase # the biological grouping
  # ) %>%
  # as_tibble(rownames = "external_gene_name")



#### Join the datasets
# join the dataset itself
comb_mat <-
  list(
    "HUT" = hut_mat,
    "Vigano" = ev_comb
    # "NOTED" = noted_set
  ) %>%
  reduce(inner_join, by = "external_gene_name")

# harmonise the pheno file and add batch parameter
# grouping says the sample type (EV/Biopsy)
# cyclephase the grouping of the cycle timing (pre/rec/post)
comb_pheno <-
  list(
    "HUT" = hut_pheno,
    "Vigano" = ev_comb_pheno %>% rename(group = grouping)
    # "NOTED" = noted_pheno %>% mutate(group = "biopsy")
  ) %>%
  bind_rows(.id = "dataset") %>%
  select(c("samplename", "group", "cyclephase", "dataset")) %>%
  # align the sample order with the matrix
  arrange(match(samplename, colnames(comb_mat)[-1]))

# apply batch normalisation on the dataframes
comb_batch <-
  comb_mat %>%
  column_to_rownames("external_gene_name") %>%
  as.matrix %>%
  sva::ComBat_seq(
    batch = comb_pheno$dataset,
    group = paste(comb_pheno$group, comb_pheno$cyclephase)
  ) %>%
  as_tibble(rownames = "gene_id")


## Combine the TPM datasets
## apply TPM norm function
apply_tpm <- function(non_tpm, gene_len)
  non_tpm %>%
  inner_join(gene_len, by = c("external_gene_name" = "gene_id")) %$%
  tpm_norm(
    as_tibble(.) %>%
      select(-len) %>%
      column_to_rownames("external_gene_name"),
    len
  ) %>%
  as_tibble(rownames = "external_gene_name")


# apply batch normalisation on the TPM dataframes
comb_batch_tpm <-
  list(
    "HUT" = hut_tpm,
    "Vigano" = ev_comb %>% apply_tpm(gene_len)
    # "NOTED" = noted_set %>% apply_tpm(gene_len)
  ) %>%
  reduce(inner_join, by = "external_gene_name") %>%
  column_to_rownames("external_gene_name") %>%
  as.matrix %>%
  sva::ComBat(
    batch = comb_pheno$dataset,
    mod =
      comb_pheno %>%
      select(-dataset) %>%
      column_to_rownames("samplename") %>%
      mutate(across(everything(), as.factor)) %>%
      model.matrix(~ group + cyclephase, data = .)
  ) %>%
  as_tibble(rownames = "gene_id")


#### write the expresiion matrix and the phenotype data out
write_feather(comb_batch, "./data/combined/comb_uf.feather")
write_feather(comb_batch_tpm, "./data/combined/comb_uf_tpm.feather")
write_tsv(comb_pheno, "./data/combined/comb_uf_pheno.tsv")
