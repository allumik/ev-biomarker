## Dependency mapping

library("tidyverse")
library("tidymodels")
library("embed")
library("feather")
library("magrittr")
library("stringr")
library("plotly")
library("reactable")
library("foreach")
library("ggvenn")
library("heatmaply")
library("glue")
library("doParallel")
library("BiocParallel")
library("patchwork")
library("tictoc")
library("future")
library("embed")
library("readxl")

## this is for some reason needed now for the venns
library("grid")

## for the fonts
# extrafont::font_import(paths = "~/.fonts/", prompt = F)


## register BiocParallel cores
if(.Platform$OS.type == "windows") {
  register(SnowParam(4))
  makeCluster(4) %>% doParallel::registerDoParallel(4)
} else {
  register(MulticoreParam(4))
  doParallel::registerDoParallel(cores = 4)
}