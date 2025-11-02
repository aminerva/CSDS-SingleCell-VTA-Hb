calc_DEGs <- function(obj, obj_name, cluster, groups, alph) {
                       
    # Function to perform DESeq2 analysis on a specified cluster of the Seurat object
    # Calculates DEGs against control and between conditions (susceptible vs resilient)

    DefaultAssay(obj) <- "RNA"

    comparisons <- switch(groups,
                            SC = "susceptible vs control",
                            RC = "resilient vs control",
                            SR = "susceptible vs resilient")

    message(paste0('Running DESeq2 for ', cluster, " (", comparisons, ")"))
    
    sample_groups <- list(
        SC = c("FCVTA1", "FCVTA2", "FCVTA3", "FSVTA1", "FSVTA2",
            "MCVTA1", "MCVTA2", "MCVTA3", "MSVTA1", "MSVTA2"),
        RC = c("FCVTA1", "FCVTA2", "FCVTA3", "FRVTA1", "FRVTA2", "FRVTA3",
            "MCVTA1", "MCVTA2", "MCVTA3", "MRVTA1", "MRVTA2", "MRVTA3"),
        SR = c("FSVTA1", "FSVTA2", "FRVTA1", "FRVTA2", "FRVTA3",
            "MSVTA1", "MSVTA2", "MRVTA1", "MRVTA2", "MRVTA3"))

    samples <- sample_groups[[groups]]

    # Subset Seurat object for the specified cluster -------------------------------------------------------
    cluster_field <- switch(obj_name,
                        all_samples = "cluster.names",
                        neurons = "subcluster.names",
                        ieg_neurons = "subcluster.names",
                        da_subclusters = "da.subclusters")
    #obj_subset <- subset(obj, obj@meta.data[[cluster_field]] == cluster & sampleID %in% samples)
    keep_cells_logical <- obj@meta.data[[cluster_field]] == cluster & obj@meta.data$sampleID %in% samples
    cells_to_keep <- rownames(obj@meta.data)[keep_cells_logical]
    obj_subset <- subset(obj, cells = cells_to_keep)

    # Extract the psuedobulked count matrix for the specified cluster
    counts_matrix <- AggregateExpression(obj_subset, assays="RNA", slot="counts", group.by=c("sampleID"))$RNA
    counts_matrix <- counts_matrix[!rownames(counts_matrix) %in% c("Xist", "Tsix"), ] # Remove Xist and Tsix

    # Build sample metadata
    sample_metadata <- data.frame(
        sampleID = colnames(counts_matrix),
        condition = factor(substr(colnames(counts_matrix), 2, 2)),
        sex = factor(substr(colnames(counts_matrix), 1, 1))
    )
    
    stopifnot(all(colnames(counts_matrix) == sample_metadata$sampleID))

    # Create DESeqDataSet ----------------------------------------------------------------------------
    dds <- DESeqDataSetFromMatrix(
        countData = counts_matrix,
        colData = sample_metadata,
        design = ~ sex + condition) # Test for effect of condition, while regressing out variation due to sex
    
    dds$condition <- relevel(dds$condition, ref = ifelse(groups == "SR", "R", "C"))
    
    message('Successfully created DESeqDataSet, normalizing counts...')

    # Normalize counts
    dds <- estimateSizeFactors(dds)
    normalized_counts <- counts(dds, normalized=TRUE)

    # Find 5th percentile of all normalized counts for this cell type
    threshold <- ifelse(quantile(as.vector(normalized_counts), probs = 0.05) > 0, 
                        quantile(as.vector(normalized_counts), probs = 0.05), 25)

    # Keep genes with normalized counts > threshold in at least 2 samples
    keep_genes <- rowSums(normalized_counts > threshold) >= 2

    # Subset the DESeqDataSet
    dds <- dds[keep_genes, ]

    # Run DESeq
    dds <- suppressMessages(DESeq(dds))

    # Extract results table from DESeq2 analysis -------------------------------------------------------
    contrast_vector <- switch(groups,
                             SC = c("condition", "S", "C"),
                             RC = c("condition", "R", "C"),
                             SR = c("condition", "S", "R"))

    res <- results(dds, contrast = contrast_vector, cooksCutoff = TRUE, independentFiltering = FALSE, alpha = alph)
  
    # Format results
    results_df <- as.data.frame(res)
    results_df <- rownames_to_column(results_df, "gene")
    results_df$comparison <- comparisons
    results_df$celltype <- cluster

    message(paste("Completed DESeq2 for", cluster))
    cat("\n")
    return(results_df)
}

#########################################################

count_DEGs_per_cluster <- function(df, groups, cluster.names = NULL, adjp.cutoff = 0.1, lfc.cutoff = 0) {
  
    # Function to count the number of DEGs per cluster and direction (up, down, both)

    # If cluster.names not provided, extract from dataframe
    if (is.null(cluster.names)) {
    cluster.names <- unique(df$celltype)
    }

    result <- data.frame(celltype = character(0), 
                        comparison = character(0), 
                        direction = character(0), 
                        count = numeric(0))

    for (g in groups) {
        for (c in cluster.names) {
      
            sub_df <- filter(df, celltype == c, comparison == g)
            
            ndown <- sub_df %>%
                filter(!is.na(padj), padj < adjp.cutoff, log2FoldChange < -abs(lfc.cutoff)) %>%
                nrow()
            
            nup <- sub_df %>%
                filter(!is.na(padj), padj < adjp.cutoff, log2FoldChange > abs(lfc.cutoff)) %>%
                nrow()
            
            result <- rbind(result, 
                            data.frame(celltype = c, comparison = g, direction = "down", count = ndown),
                            data.frame(celltype = c, comparison = g, direction = "up", count = nup),
                            data.frame(celltype = c, comparison = g, direction = "both", count = ndown + nup))
            
            message(paste("in", c, "there are", ndown+nup, "total DEGs between", g))
        }
    }

  return(result)
}

#########################################################

# Run DESeq2 across celltypes for each comparison separately
groups <- c("SC","RC","SR")
cluster.names <- c("neuron1","neuron2","neuron3","astro","micro","endo","OPCs","premyelinating_oligo","myelinating_oligo")

all_cells_deseq <- data.frame()
for (g in groups) { 
    for (c in cluster.names) {
        all_cells_deseq <- rbind(all_cells_deseq, deseq_wald(obj=all_samples, obj_name="all_samples", cluster=c, groups=g, alph=0.1))
    }
}
write.csv(all_cells_deseq, file="../analysis_outputs/DEGs/deseq/all_cells_deseq_results.csv", row.names=FALSE)

# Calculate how many DEGs there are per cluster for each comparison
all_cells_deg_counts <- count_DEGs_per_cluster(all_cells_deseq, groups=c("resilient vs control", "susceptible vs control", "susceptible vs resilient"))
write.csv(all_cells_deg_counts, file="../analysis_outputs/DEGs/deseq/all_cells_deg_counts.csv", row.names=FALSE)

#########################################################

# Run DESeq2 across neuron subclusters for each comparison separately
groups <- c("SC","RC","RS")
subcluster.names <- c("DA","Mixed","Glut","GABA-1","GABA-2","GABA-3","GABA-4","GABA-5","GABA-6")

neurons_deseq <- data.frame()

for (g in groups) { 
    for (c in subcluster.names) {
        neurons_deseq <- rbind(neurons_deseq, deseq_wald(obj=neurons2, obj_name="neurons", cluster=c, groups=g, alph=0.1))
        }
    }
write.csv(neurons_deseq, file="../analysis_outputs/DEGs/deseq/neurons_deseq_results.csv", row.names=FALSE)

# Calculate how many DEGs there are per cluster for each comparison
neurons_deg_counts <- count_DEGs_per_cluster(neurons_deseq, groups=c("resilient vs control", "susceptible vs control", "susceptible vs resilient"))
write.csv(neurons_deg_counts, file="../analysis_outputs/DEGs/deseq/neurons_deg_counts.csv", row.names=FALSE)









