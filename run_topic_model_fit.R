print("start")
args = commandArgs(trailingOnly=TRUE)

input_path <- args[1]
output_path <- args[2]
K <- as.integer(args[3])
if (length(args) >= 4) {
  NC <- as.integer(args[4])
} else {
  NC <- 1
}
if (length(args) >= 5) {
  do_transpose <- as.integer(args[5])
} else {
  do_transpose <- 1
}
## libraries 
print("load library")
library(Matrix)
library(fastTopics)
library(ggplot2)
library(cowplot)
library(SeuratDisk)
library(SeuratData)
set.seed(1)
print("load object")
obj <- LoadH5Seurat(input_path)

################################################################################
## fetch the data: 
################################################################################
counts <- obj@assays$RNA@counts
if (do_transpose == 1) {
   counts <- t(counts)
} 


################################################################################
## fit the topic model 
################################################################################
print("Starting the fit")
fit <- fit_topic_model(counts, 
                       k = K,
                       control.init = list(nc = NC),
                       control.main = list(nc = NC),
                       control.refine = list(nc = NC)
)
print("Finished the fit")

saveRDS(fit,  paste0(output_path, "fitted_topic_model_k_", K, ".rds"))
print("complite")

