#### Functions for model metrics
## This file is for separating the functions used in
## beREADY model evaluation.

silent_beready_model <- function(...) {
  suppressMessages(
    capture.output(x <- tacseqApp::beready_model(...)))
  return(x)
}

## data formatter
## format the geneset as this
# tacseqApp::targets_data$`Dynamic filtering geneset`
format_for_beready <- function(
  input_dat,
  input_genes,
  pheno_dat,
  hs_genes
  ) {
  input_dat %>%
    ## select only genes of interest
    filter(gene_id %in% input_genes) %>%
    ## houskeeper normalise the integers
    mutate_if(is.numeric, function(x) {
      x / exp(mean(log(x[.$gene_id %in% hs_genes])))
    }) %>%
    ## pivor longer for the phenotype join
    pivot_longer(where(is.numeric),
      names_to = "samplename",
      values_to = "norm_count"
    ) %>%
    inner_join(pheno_dat, by = "samplename") %>%
    mutate(
      sample = samplename,
      batch = dataset,
      group = case_when(
        cyclephase == "pro" ~ "PE",
        cyclephase == "pre" ~ "LH+2",
        cyclephase == "rec" ~ "LH+7",
        cyclephase == "post" ~ "LH+10",
        T ~ NA
      ),
      .keep = "unused"
    ) %>%
    # pivot genes wider for the beready model format
    pivot_wider(
      id_cols = where(is.character), 
      names_from = "gene_id",
      values_from = "norm_count"
      )
}

# for the accuracy results
beready_results_formatter <- function(score_table, real_groups)
  score_table %>%
    left_join(real_groups, by = c("Sample Name" = "sample")) %>%
    # manually select those groups that match with the prediction
    mutate(match =
      `Predicted Group` == "Pre-receptive" & group == "LH+2" |
      `Predicted Group` == "Receptive" & group == "LH+7" |
      `Predicted Group` == "Post-receptive" & group == "LH+10" |
      `Predicted Group` == "Early-receptive" & group == "LH+7" |
      `Predicted Group` == "Late-receptive" & group == "LH+7"
    )

# depends on `beready_results_formatter` output
beready_summary_stats <- function(formatted_score_table)
  # duplicate results for the summary stats
  bind_rows(
    formatted_score_table %>% mutate(grouping_var =
      ifelse(str_detect(`Sample Name`, "UF"), "UF", "Biopsy")),
    formatted_score_table %>% mutate(grouping_var = "All")
  ) %>%
  group_by(grouping_var) %>%
  summarise(
    "Samples" = length(match),
    "Hits" = sum(match),
    "Misses" = sum(!match),
    "No Results" = sum(`Predicted Group` == "No result"),
    "Accuracy" = sum(match) / length(match),
    "Pre Acc" = sum(match[`Predicted Group` == "Pre-receptive"]) / 
                      length(match[`Predicted Group` == "Pre-receptive"]),
    "Rec Acc" = sum(match[str_detect(`Predicted Group`, "Early|Receptive|Late")]) / 
                      length(match[str_detect(`Predicted Group`, "Early|Receptive|Late")]),
    "Post Acc" = sum(match[`Predicted Group` == "Post-receptive"]) / 
                      length(match[`Predicted Group` == "Post-receptive"]),
  ) %>%
  # next two steps basically just transposes the table...
  pivot_longer(
    where(is_bare_numeric),
    names_to = "Metrics",
    values_to = "Results"
  ) %>%
  pivot_wider(
    names_from = "grouping_var",
    values_from = "Results"
  )

# function to perform CV and emit a results table
beready_cv_fold <- function(split_object, dynamic_genes) {
  train_dat <- training(split_object)
  test_dat <- testing(split_object)

  # # as safety - remove reference samples with near zero variance across genes
  # filter_samples <-
  #   train_dat %>%
  #     pivot_longer(where(is_bare_numeric), values_to = "norm_count", names_to = "gene_id") %>%
  #     group_by(sample) %>%
  #     summarise(var_counts = var(norm_count)) %>%
  #     filter(var_counts > .01) %$%
  #     sample
  # train_dat %<>% filter(sample %in% filter_samples)

  tryCatch(
    silent_beready_model(train_dat, dynamic_genes, test_dat) %$%
      score_table %>%
      beready_results_formatter(
        test_dat %>% select(c("sample", "group"))
      ) %>%
      beready_summary_stats,
    # output empty tibble so that it would work with row_bind() in foreach loops
    error = function(err) tibble()
  )
}



