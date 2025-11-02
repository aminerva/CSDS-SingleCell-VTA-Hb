#!/usr/bin/env Rscript

message("Loading packages")
suppressPackageStartupMessages(library("languageserver"))
suppressPackageStartupMessages(library("ashr"))
suppressPackageStartupMessages(library("caret"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("Seurat"))
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("sda"))
source("sda_functions_cv_by_sample.R") # Custom SDA functions

set.seed(123)

#########################################

message("Reading in arguments")
args <- commandArgs(trailingOnly = TRUE)
cluster <- args[1]

output_dir <- "../analysis_outputs/LDA_cv_by_sample"
output_prefix <- paste0(output_dir, "/1D_sda_", cluster)

#########################################

message("Reading Seurat object")
neurons <- readRDS("../analysis_outputs/seurat_objects/neuron_subset.RDS")
# Remove SERT cluster
neurons <- subset(neurons, idents=c("DA","Mixed","Glut","GABA-1","GABA-2",
                                    "GABA-3","GABA-4","GABA-5","GABA-6"))
Idents(neurons) <- factor(Idents(neurons), levels = c("DA","Mixed","Glut","GABA-1","GABA-2",
                                                      "GABA-3","GABA-4","GABA-5","GABA-6"))

if (cluster != "all") {
  message("Subsetting to cluster: ", cluster)
  if (!(cluster %in% levels(neurons))) stop("Cluster '", cluster, "' not found in object")
  obj <- subset(neurons, idents = cluster)
} else {
  message("Using full Seurat object")
  obj <- neurons
}

#########################################

message("Preparing input data for SDA")
data <- as.matrix(GetAssayData(object = obj, layer = "data", assay = "RNA"))
labels <- obj$outcome  
names(labels) <- colnames(obj)
sample_ids <- obj$sampleID
data_t <- as.data.frame(t(data))

# Find top variable genes
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000) 
variable_genes <- VariableFeatures(obj)

#########################################

n_folds <- length(unique(sample_ids[labels %in% c("Control", "Resilient")]))

message("Running ", n_folds, "-fold cross validation SDA with true labels...")
unshuffled_result <- cross_validate_1D_sda(
  data_t = data_t,
  labels = labels,
  genes = variable_genes,
  sample_ids = sample_ids,
  shuffle_bool = FALSE,
  train_conditions = c("Control", "Resilient"))

message("Running ", n_folds, "-fold cross validation SDA with shuffled labels...")
shuffled_result <- cross_validate_1D_sda(
  data_t = data_t,
  labels = labels,
  genes = variable_genes,
  sample_ids = sample_ids,
  shuffle_bool = TRUE,
  train_conditions = c("Control", "Resilient"))

#########################################

message("Saving results...")
saveRDS(unshuffled_result, paste0(output_prefix, "_unshuffled_cv.rds"))
saveRDS(shuffled_result, paste0(output_prefix, "_shuffled_cv.rds"))

# Also save accuracies in a summary csv
accuracy_df <- data.frame(
  Fold = 1:n_folds,
  Unshuffled = unshuffled_result$accuracies,
  Shuffled = shuffled_result$accuracies)
write.csv(accuracy_df, file = paste0(output_prefix, "_accuracy_summary.csv"), 
            row.names = FALSE)

message("Done with cross validation")
message("")

message("Now training on all Control + Resilient data and projecting Susceptible onto LD1")

sda_projection <- run_1D_sda_projection(
  data_t = data_t,
  labels = labels,
  genes = variable_genes,
  shuffle = FALSE,
  train_conditions = c("Control", "Resilient"),
  project_condition = "Susceptible")
saveRDS(sda_projection$sda_df, paste0(output_prefix, "_susceptible_projection.rds"))

message("Classifying susceptible cells based on SDA model")
susceptible_cells <- rownames(data_t)[labels == "Susceptible"]
X_susceptible <- as.matrix(data_t[susceptible_cells, colnames(sda_projection$sda_model$beta), drop = FALSE])

predicted_labels <- predict(sda_projection$sda_model, X_susceptible)$class
classification_summary <- table(predicted_labels)

message("Susceptible cell classification summary:")
print(classification_summary)

# Save classification result
classification_df <- data.frame(
  Cell = susceptible_cells,
  PredictedLabel = predicted_labels
)
write.csv(classification_df, paste0(output_prefix, "_susceptible_classification.csv"), row.names = FALSE)

#########################################

message("Extracting gene contributions to LD1")
gene_contributions <- data.frame(
  Gene = colnames(sda_projection$sda_model$beta),
  LD1 = sda_projection$sda_model$beta[1, ])

write.csv(gene_contributions, paste0(output_prefix, "_gene_contributions.csv"), row.names = FALSE)


message("All done!")

