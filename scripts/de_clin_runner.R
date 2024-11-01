#### Differential analysis runner for the ev-biomarker studies, combined dataset
## Depends on back.R and back_diff.R

source("./scripts/load_deps.R")
source("./scripts/backend_env_functions.R")
## source the combined data generator
## instead of sourcing the computationally heavy script and all datasets, just load the datafiles, if available
source("./scripts/de_functions.R")

if(!file.exists("./data/clin_counts_raw.feather"))
  callr::rscript("./scripts/format_clinical_data.R")


## read in the clinical samples and convert genenames to HUGO
clin_pheno <-
  read_tsv("./data/phenotype_clin.tsv")

clin_raw_mat <-
  read_feather("./data/clin_counts_raw.feather") %>%
  rename(ensembl_gene_id = gene_id) %>%
  geneid_converter(annot) %>%
  rename(gene_id = external_gene_name) %>%
  filter(!is.na(gene_id))

clin_tpm_mat <-
  read_feather("./data/clin_tpm_raw.feather") %>%
  rename(ensembl_gene_id = gene_id) %>%
  geneid_converter(annot) %>%
  rename(gene_id = external_gene_name) %>%
  filter(!is.na(gene_id))

clin_expr <-
  clin_raw_mat %>%
  transform_count_mat(clin_pheno, id_columns = c("samplename", "implantation")) %>%
  mutate(across(where(is.numeric), ~ replace_na(.x, 0)))

clin_expr_tpm <-
  clin_tpm_mat %>%
  transform_count_mat(clin_pheno, id_columns = c("samplename", "implantation")) %>%
  mutate(across(where(is.numeric), ~ replace_na(.x, 0)))



## Diff analysis comparisons within `group`
comps <- list(
  "implantation_vs_nonimplantation" = c("implantation", "nonimplantation")
)

pval <- .1

diffterm <-
  clin_pheno %>%
  filter(!is.na(implantation)) %>%
  mutate( 
    implantation = factor(ifelse(implantation, "implantation", "nonimplantation")),
    n_transfer_embryo = factor(n_transfer_embryo)
    ) %>%
  column_to_rownames("samplename")

## Filter out counts with edgeR's cpm method
filter_ids <-
  cpm_filter_grouped(
    clin_raw_mat %>%
      select(c("gene_id", rownames(diffterm))),
    setNames(rownames(diffterm), diffterm$implantation)
    )

# also run the differential analysis for between the UF and Biopsy
# but same implantation group
named_groups <-
  setNames(
    as.character(unique(diffterm$implantation)),
    as.character(unique(diffterm$implantation))
    )

samples_in_group <- lapply(
  named_groups,
  function(implantation_if)
    diffterm %>%
      filter(implantation == implantation_if) %>%
      rownames
)

deseq_all <-
  dds_maker(
    clin_raw_mat %>% select(c("gene_id", rownames(diffterm))),
    diffterm,
    ~ implantation,
    filter_ids = filter_ids
    ) %>%
  deseq_results(comps, annot, groupterm = "implantation", annotation_join_on = "external_gene_name", pval = 1)

deseq_emb_one_all <-
  dds_maker(
    clin_raw_mat %>% select(c("gene_id", rownames(diffterm[diffterm$n_transfer_embryo == 1, ]))),
    diffterm[diffterm$n_transfer_embryo == 1, ],
    ~ implantation,
    filter_ids = filter_ids
    ) %>%
  deseq_results(comps, annot, groupterm = "implantation", annotation_join_on = "external_gene_name", pval = 1)

deseq_emb_two_all <-
  dds_maker(
    clin_raw_mat %>% select(c("gene_id", rownames(diffterm[diffterm$n_transfer_embryo == 2, ]))),
    diffterm[diffterm$n_transfer_embryo == 2, ],
    ~ implantation,
    filter_ids = filter_ids
    ) %>%
  deseq_results(comps, annot, groupterm = "implantation", annotation_join_on = "external_gene_name", pval = 1)

deseq_top <- list( "implantation_vs_nonimplantation" = deseq_all$implantation_vs_nonimplantation %>% filter(padj < pval) )
deseq_emb_one_top <- list( "implantation_vs_nonimplantation" = deseq_emb_one_all$implantation_vs_nonimplantation %>% filter(padj < pval) )
deseq_emb_two_top <- list( "implantation_vs_nonimplantation" = deseq_emb_two_all$implantation_vs_nonimplantation %>% filter(padj < pval) )

deseq_res_group <- list(
  "implantation" = c(deseq_top, extract_top(deseq_top), deseq_all),
  "implantation_one_emb" = c(deseq_emb_one_top, extract_top(deseq_emb_one_top), deseq_emb_one_all),
  "implantation_two_emb" = c(deseq_emb_two_top, extract_top(deseq_emb_two_top), deseq_emb_two_all)
  )


edger_obj <-
  edger_maker(
    clin_raw_mat %>% select(c("gene_id", rownames(diffterm))),
    diffterm,
    "implantation",
    filter_ids
  )
edger_all <-
  edger_obj %>% 
  edger_results(comps, annot, annotation_join_on = "external_gene_name", pval = 1)

edger_emb_one_obj <- 
  edger_maker(
    clin_raw_mat %>% select(c("gene_id", rownames(diffterm[diffterm$n_transfer_embryo == 1, ]))),
    diffterm[diffterm$n_transfer_embryo == 1, ],
    "implantation",
    filter_ids
  )
edger_emb_one_all <-
  edger_emb_one_obj %>%
  edger_results(comps, annot, annotation_join_on = "external_gene_name", pval = 1)

edger_emb_two_obj <-
  edger_maker(
    clin_raw_mat %>% select(c("gene_id", rownames(diffterm[diffterm$n_transfer_embryo == 2, ]))),
    diffterm[diffterm$n_transfer_embryo == 2, ],
    "implantation",
    filter_ids
  )
edger_emb_two_all <-
  edger_emb_two_obj %>%
  edger_results(comps, annot, annotation_join_on = "external_gene_name", pval = 1)

edger_top <- list( "implantation_vs_nonimplantation" = edger_all$implantation_vs_nonimplantation %>% filter(FDR < pval) )
edger_emb_one_top <- list( "implantation_vs_nonimplantation" = edger_emb_one_all$implantation_vs_nonimplantation %>% filter(FDR < pval) )
edger_emb_two_top <- list( "implantation_vs_nonimplantation" = edger_emb_two_all$implantation_vs_nonimplantation %>% filter(FDR < pval) )

edger_res_group <- list(
  "implantation" = c(edger_top, extract_top(edger_top), edger_all),
  "implantation_one_emb" = c(edger_emb_one_top, extract_top(edger_emb_one_top), edger_emb_one_all),
  "implantation_two_emb" = c(edger_emb_two_top, extract_top(edger_emb_two_top), edger_emb_two_all)
  )


voom_all <-
  edger_obj %>%
  voom_maker(comps) %>%
  voom_results(comps, annot, annotation_join_on = "external_gene_name", pval = 1)

voom_emb_one_all <-
  edger_emb_one_obj %>%
  voom_maker(comps) %>%
  voom_results(comps, annot, annotation_join_on = "external_gene_name", pval = 1)

voom_emb_two_all <-
  edger_emb_two_obj %>%
  voom_maker(comps) %>%
  voom_results(comps, annot, annotation_join_on = "external_gene_name", pval = 1)

voom_top <- list( "implantation_vs_nonimplantation" = voom_all$implantation_vs_nonimplantation %>% filter(adj.P.Val < pval) )
voom_emb_one_top <- list( "implantation_vs_nonimplantation" = voom_emb_one_all$implantation_vs_nonimplantation %>% filter(adj.P.Val < pval) )
voom_emb_two_top <- list( "implantation_vs_nonimplantation" = voom_emb_two_all$implantation_vs_nonimplantation %>% filter(adj.P.Val < pval) )

voom_res_group <- list(
  "implantation" = c(voom_top, extract_top(voom_top), voom_all),
  "implantation_one_emb" = c(voom_emb_one_top, extract_top(voom_emb_one_top), voom_emb_one_all),
  "implantation_two_emb" = c(voom_emb_two_top, extract_top(voom_emb_two_top), voom_emb_two_all)
  )




#### Formatting some data for the raport

## re/precalculate some large nasty bits
global_expr_pca <-
  clin_expr_tpm %>%
    select(-implantation) %>%
    recipe() %>%
    update_role(samplename, new_role = "id") %>%
    ## eliminate some genes with zero expression
    step_nzv(all_numeric()) %>%
    step_normalize(all_numeric()) %>%
    step_pca(all_numeric(), num_comp = 3) %>%
    prep(strings_as_factors = F) %>%
    bake(new_data = clin_expr_tpm %>% select(-implantation)) %>%
    left_join(clin_expr_tpm %>% select(c("samplename", "implantation")), by = "samplename")

global_expr_umap <-
  clin_expr_tpm %>%
    select(-implantation) %>%
    recipe() %>%
    update_role(samplename, new_role = "id") %>%
    ## eliminate some genes with zero expression
    step_nzv(all_numeric()) %>%
    step_normalize(all_numeric()) %>%
    step_umap(all_numeric(), num_comp = 2) %>%
    prep(strings_as_factors = F) %>%
    bake(new_data = clin_expr_tpm %>% select(-implantation)) %>%
    left_join(clin_expr_tpm %>% select(c("samplename", "implantation")), by = "samplename")


out_join <-
  unique(c(
    deseq_res_group$implantation$implantation_vs_nonimplantation$locus,
    voom_res_group$implantation$implantation_vs_nonimplantation$locus,
    edger_res_group$implantation$implantation_vs_nonimplantation$locus
  )) %>% {
    tibble(
      locus = .,
      deseq_implantation = . %in% deseq_res_group$implantation$implantation_vs_nonimplantation$locus,
      voom_implantation = . %in% voom_res_group$implantation$implantation_vs_nonimplantation$locus,
      edger_implantation = . %in% edger_res_group$implantation$implantation_vs_nonimplantation$locus
      )
    }

out_join_fdrless <-
  unique(c(
    filter(deseq_res_group$implantation[4]$implantation_vs_nonimplantation, pvalue < 0.05 & abs(log2FoldChange > 1.5))$locus,
    filter(voom_res_group$implantation[4]$implantation_vs_nonimplantation, P.Value < 0.05 & abs(logFC > 1.5))$locus,
    filter(edger_res_group$implantation[4]$implantation_vs_nonimplantation, PValue < 0.05 & abs(logFC > 1.5))$locus
  )) %>% {
    tibble(
      locus = .,
      deseq_implantation = . %in% deseq_res_group$implantation[4]$implantation_vs_nonimplantation$locus,
      voom_implantation = . %in% voom_res_group$implantation[4]$implantation_vs_nonimplantation$locus,
      edger_implantation = . %in% edger_res_group$implantation[4]$implantation_vs_nonimplantation$locus
      )
    }
