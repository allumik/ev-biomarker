#### Backend script for xCell cell type enrichment tool
##
## this ofcourse requires dviraran/xCell
## with `devtools::install_github("dviraran/xCell")`
##
source("./scripts/load_deps.R")
## load the old RData, as it has all the samples without filtering out
## informative to apply xCell in this context as it helps to explain some
## outliers
load("./data/orig/ev_data.RData")

## Load endometrial reference data for xCell
xCell.data <- xCell::xCell.data
# xCell "genes" (actually dont know what theyre for)
xcell_genes <- xCell.data$genes
# unique genes for cell type signatures in xCell
xcell_sigs <- sapply(
  names(xCell.data$signatures),
  function(x)
    xCell.data$signatures[[x]]@geneIds
) %>% unlist() %>% unique()



# Combine the NOTED, Vigano and HUT samples together in another script
# and load the resulting batch-corrected dataset here
ev_comb <- read_feather("./data/combined/comb_uf.feather")
ev_comb_tpm <- read_feather("./data/combined/comb_uf_tpm.feather")
ev_comb_pheno <- read_tsv("./data/combined/comb_uf_pheno.tsv")



## Additional set of EV samples undergoing ART treatment
ev_art <-
  read_excel("./data_raw/Vigano_raw counts_ART.xlsx", range = "A3:AX56612")

ev_art_pheno <-
  read_excel("./data_raw/Vigano_raw counts_ART.xlsx", range = "A3:AX3") %>% {
    tibble(samplename = colnames(.))
  } %>%
    slice(-1) %>%
    mutate(
      grouping = "lh7",
      excluded = str_detect(samplename,"EUF(10|19|29|30|21|23|25)$"),
      preg = str_detect(samplename,
                        "EUF(2|5|8|9|11|16|17|21|22|23|25|28|33|34|36|45|47|48|49|50|51|52)$")
    )

ev_art_tpm <-
  inner_join(ev_art, gene_len, by = c("GeneID" = "gene_id")) %$%
    tpm_norm(
      ## filter out the low-quality samples
      select(
        .,
        -c("len", ev_art_pheno %>% filter(excluded) %$% samplename)
      ) %>%
        column_to_rownames("GeneID"),
      len
    ) %>%
    as_tibble(rownames = "GeneID") %>%
    mutate(across(where(is.numeric), ~ replace_na(.x, 0)))



## from SCRATCH2
scratch_set <-
  read_tsv("./data_raw/TPM_allsamples_matrix.zip")



## from cERtain
certain_tpm <-
  read_tsv("./data_raw/rsem.gene_tpm.certain.tsv") %>%
  mutate(across(where(is.numeric), ~ replace_na(.x, 0))) %>%
  select(-HUT26_S8)

certain_pheno <-
  readxl::read_xlsx("./data_raw/Sequenced brushes final - alt groups.xlsx") %>%
  select(-Sample) %>%
  rename("samplename" = "Name on the tube") %>%
  inner_join(
    tibble(
      "fullname" = colnames(certain_tpm),
      "samplename" = str_remove(colnames(certain_tpm), "_S\\d+")
    ),
    by = "samplename"
  ) %>%
  ## make groupnames safer for deseq etc
  mutate(
    daytype = ifelse(Group == "HRT", "hrt", "lh"),
    group = case_when(
      `Alt Group` == "Proliferative"  ~ "pro",
      `Alt Group` == "Pre-receptive"  ~ "pre",
      `Alt Group` == "Receptive"      ~ "rec",
      `Alt Group` == "Post-receptive" ~ "post",
      `Alt Group` == "HRT"            ~ "hrt",
      T ~ `Alt Group`
    ),
    pooled = case_when(
      Group == "Proliferative"  ~ "pooledpre",
      Group == "Pre-receptive"  ~ "pooledpre",
      Group == "Receptive"      ~ "pooledlate",
      Group == "Post-receptive" ~ "pooledlate",
      Group == "HRT"            ~ "hrt",
      T ~ Group
    ))

## from blood analysis
exoeasy1 <-
  read_tsv("./data_raw/GSE133684_exp_TPM-all.txt.gz") %>%
    rename(ensembl_gene_id = `...1`) %>%
    select(c("ensembl_gene_id", starts_with("N_"))) %>%
    geneid_converter(annot)

exoeasy2 <-
  read_tsv("./data_raw/GSE159657_Matrix_file.txt.gz") %>%
    rename(ensembl_gene_id = AccID) %>%
    select(c("ensembl_gene_id", contains("control-"))) %>%
    geneid_converter(annot)



## Run xCell
## Join all the datasets together for xCell
joined_expr <-
  ev_comb_tpm %>%
  rename(external_gene_name = "gene_id") %>%
  # add the ART samples not in the ev_comb dataset
  inner_join(
    ev_art_tpm %>%
      select(c("GeneID", ev_art_pheno %>% filter(!excluded) %$% samplename)),
    by = c("external_gene_name" = "GeneID")
  ) %>%
  # and the SCRATCH dataset
  inner_join(
    scratch_set %>%
      rename(ensembl_gene_id = mRNA) %>%
      geneid_converter(annot),
    by = "external_gene_name") %>%
  # and finally cERtain samples too
  inner_join(certain_tpm, by = "external_gene_name") %>%
  inner_join(exoeasy1, by = "external_gene_name") %>%
  inner_join(exoeasy2, by = "external_gene_name") %>%
  # fix the dataframe if there are NA's in gene names
  ## Do a outer join with the necessary signature genes from xcell
  ## all the missing genes shold have NA's, and next we turn them to 0
  full_join(tibble(xcell_sigs), by = c("external_gene_name" = "xcell_sigs")) %>%
  mutate(across(where(is.numeric), ~ replace_na(.x, 0))) %>%
  drop_na(external_gene_name) %>%
  # remove some samples per Elina's wish
  select(-matches("HUT1_|HUT01|HUT71"))

## prepare datasets to run sequencially with xCell
xcell_runs <- list(
  "HUT" =
    raw_tpm %>%
      rename(ensembl_gene_id = gene_id) %>%
      geneid_converter(annot) %>%
      full_join(tibble(xcell_sigs), by = c("external_gene_name" = "xcell_sigs")) %>%
      mutate(across(where(is.numeric), ~ replace_na(.x, 0))) %>%
      drop_na(external_gene_name) %>%
      column_to_rownames("external_gene_name"),
  "Combined EV" = ev_comb_tpm %>% column_to_rownames("gene_id"),
  "Joined set" = joined_expr %>% column_to_rownames("external_gene_name")
)

## phenotype document for the Joined set of all possible samples in analysis
xcell_objs_pheno <-
  list(
    "EV ART IS" =
      ev_art_pheno %>% 
        filter(preg & !excluded) %>%
        mutate(grouping = "UF rec"),
    "EV ART IF" =
      ev_art_pheno %>%
        filter(!preg & !excluded) %>%
        mutate(grouping = "UF rec"),
    "SCRATCH2" =
      tibble(
        samplename = colnames(scratch_set)[-1],
        grouping = "biopsy rec"
      ),
    "Combined" =
      ev_comb_pheno %>%
        mutate(grouping = paste(group, cyclephase)) %>%
        filter(!(str_detect(samplename, "HUT1_|HUT01|HUT71"))),
    "cERtain" =
      certain_pheno %>%
        select(c("fullname", "group")) %>%
        mutate(
          samplename = fullname,
          grouping = paste("cervix", group)
          ),
    "exoeasy1" =
      tibble(
        samplename = colnames(exoeasy1)[-1],
        grouping = "blood panc"
      ),
    "exoeasy2" =
      tibble(
        samplename = colnames(exoeasy2)[-1],
        grouping = "blood cardiovasc"
      )
  ) %>%
    bind_rows(.id = "sampleset") %>%
    mutate(
      sampleset = ifelse(
        is.na(dataset),
        sampleset,
        paste(sampleset, dataset))
        ) %>%
    select(sampleset, samplename, grouping)



## it requires the "xCell.data" object to be imported to namespace
run_xcell <- function(rna_df) {
  ## need to load it for the xCell to work
  xCell.data <- xCell::xCell.data
  partype <- ifelse(.Platform$OS.type == "windows", "SOCK", "FORK")
  ## parallel backend is handled by biocparallel
  xcell_obj <- list(
    "scores" = xCell::xCellAnalysis(rna_df, parallel.type = partype)
  )
  xcell_obj %>% append(list(
    "p_values" =
      xCell::xCellSignifcanceBetaDist(xcell_obj$scores)
  ))
}

## preflight heck: how well are signatures covered?
table(xcell_genes %in% rownames(xcell_runs[[1]]))
table(xcell_sigs %in% rownames(xcell_runs[[1]]))


## parallel run
## soooo it does not work...
## or it does, just not in some machines and only with micromamba envs
xcell_objs <-
  foreach(
    i = names(xcell_runs),
    .final = function(i) setNames(i, names(xcell_runs)),
    .packages = c("xCell", "tidyverse")
  ) %dopar%
    run_xcell(xcell_runs[[i]])

saveRDS(xcell_objs, file = "./data/xcell_results.rds")
saveRDS(xcell_objs_pheno, file = "./data/xcell_results_pheno.rds")
