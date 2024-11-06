#### This script is purely for recreating the base RData files for the reports
source("./scripts/load_deps.R")


#### First the original dataset
## check if theres a R data object already, othewise run the data muncher
# and save the preprocessed data into the RData
data_subfolder <- "./data/orig/"
dataobj <- paste0(data_subfolder, "ev_data.RData")
if (!file.exists(dataobj)) {
  ## load in the backend data muncher
  commandArgs <- function(...)
    list(data_subfolder = data_subfolder)
  dir.create(data_subfolder, showWarnings = F)
  invisible(source("./scripts/de_runner.R"))
  save.image(file = dataobj)
}


#### Then the filtered dataset
## check if theres a R data object already, othewise run the data muncher
# and save the preprocessed data into the RData
## remove previous data from the clean one.
rm(list = ls())
data_subfolder <- "./data/filtered/"
dataobj <- paste0(data_subfolder, "ev_data_filt.RData")
if (!file.exists(dataobj)) {
  ## load in the backend data muncher with arguments
  commandArgs <- function(...)
    list(
      switch_sample = "HUT23",
      remove_sample = c(
        "HUT10_UF", "HUT1_UF", "HUT1_biopsy", "HUT35_UF", "HUT42_UF", "HUT71_biopsy_2", "HUT71_biopsy_3"
        ),
      data_subfolder = data_subfolder
      ) #  c("HUT71_UF", "HUT1_UF", "HUT01_UF", "HUT17_UF", "HUT10_UF"))
  dir.create(data_subfolder, showWarnings = F)
  invisible(source("./scripts/de_runner.R"))
  save.image(file = dataobj)
}


#### And then the combined one
## check if theres a R data object already, othewise run the data muncher
# and save the preprocessed data into the RData
# use "comb_dataobj" to avoid conflict with "dataobj"
## don't remove previous data, otherwise have to load it again...
comb_data_subfolder <- "./data/combined/"
dataobj <- paste0(comb_data_subfolder, "ev_data_comb.RData")
if (!file.exists(dataobj)) {
  ## load in the backend data muncher with arguments
  commandArgs <- function(...) list(comb_data_subfolder = comb_data_subfolder)
  dir.create(comb_data_subfolder, showWarnings = F)
  invisible(source("./scripts/de_comb_runner.R"))
  save.image(file = dataobj)
}

