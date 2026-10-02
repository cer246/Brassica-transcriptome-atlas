##### Module - Tissue Assignments ####

library(tidyverse)
library(WGCNA)
library(ggplot2)
library(viridis)

##### Read in the necessary files ####
vsd_clean_tiss_at <- read_csv("/mnt/Chikorita/Brassica_rapa/coexpression/AT/all_ens_lncRNA/avg_vsd_clean_at_tiss_atlas_at_all_ens_lncRNA.csv")
vsd_clean_tiss_br <- read_csv("/mnt/Chikorita/Brassica_rapa/coexpression/core_files/avg_vsd_counts_br_tiss_atlas.csv")

br_topology <- read_delim("/mnt/Chikorita/Brassica_rapa/coexpression/co_expression_comp/br_topology.tsv", 
                          delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE)
at_topology <- read_delim("/mnt/Chikorita/Brassica_rapa/coexpression/co_expression_comp/at_topology.tsv", 
                          delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE) 

#### Assign Tissues to Modules using Eigengenes ####

##### Taking the top 80% of variable genes to assign to assign modules #####
filter_by_variance <- function(expr_matrix, percentile = 0.8) {
  
  print(paste("Original dimensions:", paste(dim(expr_matrix), collapse="x")))
  
  # 1. Calculate variance for every gene (row)
  # (Assumes rows are Genes, columns are Samples)
  gene_variances <- apply(expr_matrix, 1, var, na.rm = TRUE)
  
  # 2. Find the cutoff value
  # To keep top 80%, we cut at the 20th percentile (bottom 20%)
  cutoff_value <- quantile(gene_variances, probs = (1 - percentile))
  
  # 3. Filter the matrix
  # Keep genes with variance strictly greater than cutoff
  expr_filtered <- expr_matrix[gene_variances > cutoff_value, ]
  
  print(paste("Filtered dimensions:", paste(dim(expr_filtered), collapse="x")))
  print(paste("Removed", nrow(expr_matrix) - nrow(expr_filtered), "low-variance genes."))
  
  return(expr_filtered)
}

#### Robust Tissue Assignment ####
run_tissue_assignment <- function(expr_data, mod_df, species_prefix) {
  
  #Keep the genes that are in the top 80% of expressed
  common_genes <- intersect(colnames(expr_data), mod_df$Gene)
  #common_genes <- intersect(expr_data[1], mod_df[1])
  #common_genes <- intersect(expr_data$Gene, mod_df$Gene)
  
  if(length(common_genes) == 0) {
    stop("Error: No common gene IDs found. Check your gene ID formats (e.g., 'AT1G...' vs 'gene:AT1G...').")
  }
  
  #Subset to common genes
  expr_subset <- expr_data[, common_genes]
 
   mod_subset <- mod_df[match(common_genes, mod_df$Gene), ]
  
  #Extract module assignments
  #We assume modules are named like "At_M1", "At_M2" -> we want "1", "2"
  clean_colors <- gsub(".*_M", "", mod_subset_test$Module)
  
  #Calculate Eigengenes
  # This returns a matrix with columns named "ME1", "ME2", etc.
  #expr_subset %>%column_to_rownames("Gene") -> t(expr_subset_test)
  MEs_obj <- moduleEigengenes(t(expr_subset_test), colors = clean_colors)
  MEs <- MEs_obj$eigengenes
  
  #Just check what Eigengenes look like
  print(paste("Calculated", ncol(MEs), "module eigengenes."))
  # print(colnames(MEs)[1:5]) # Optional: verify names are ME0, ME1...
  
  #Identify Tissues from Sample Names
  #Assumes sample names are like "Tissue_more_info" -> splits by "_" and takes first part
  #get_tissue <- function(x) { strsplit(x, "_")[[1]][1] }
  get_tissue <- function(x) {x}
  tissue_labels <- sapply(rownames(expr_subset), get_tissue)
  
  #Verify we found tissues
  print(paste("Identified", length(unique(tissue_labels)), "unique tissues:"))
  print(unique(tissue_labels))
  
  #Assign Dominant Tissue
  #Initialize results DF
  results <- data.frame(Module = character(), 
                        Dominant_Tissue = character(), 
                        stringsAsFactors = FALSE)
  
  for(i in 1:ncol(MEs)) {
    me_vec <- MEs[, i]     # The vector of expression for this module
    mod_col_name <- names(MEs)[i] # e.g., "ME1"
    
    # Clean the name: "ME1" -> "1"
    mod_num <- gsub("ME", "", mod_col_name)
    
    # Construct the final name: "At_M1"
    final_mod_name <- paste0(species_prefix, "_M", mod_num)
    
    # Calculate average expression per tissue
    tissue_means <- tapply(me_vec, tissue_labels, mean, na.rm = TRUE)
    
    # Find the tissue with the highest value
    best_tissue <- names(tissue_means)[which.max(tissue_means)]
    
    # Append to results
    results <- rbind(results, data.frame(Module = final_mod_name, 
                                         Dominant_Tissue = best_tissue))
  }
  
  return(results)
}

#### Put it all together ####
br_topology %>%
  select(Gene, Module) -> modules_Br #42,830 genes

at_topology %>%
  select(Gene, Module) -> modules_At #33,063 genes

at_ordered <- filter_by_variance(vsd_clean_tiss_at, percentile = 0.8) #Removed 9769 low-variance genes
br_ordered <- filter_by_variance(vsd_clean_tiss_br, percentile = 0.8) #Removed 11330 low-variance genes

at_ordered %>% column_to_rownames("Gene") -> at_ordered
br_ordered %>% column_to_rownames("Gene") -> br_ordered

tissue_map_At <- run_tissue_assignment(t(at_ordered), modules_At, "At")
print(head(tissue_map_At))

# Do the same for Brassica
tissue_map_Br <- run_tissue_assignment(t(br_ordered), modules_Br, "Br")
print(head(tissue_map_Br))

#### Update tissue comparison table ####
# 1. Rename columns in the maps to match the merging strategy
map_at_clean <- tissue_map_At %>% 
  rename(At_Module = Module, 
         At_Dominant_Tissue = Dominant_Tissue)

map_br_clean <- tissue_map_Br %>% 
  rename(Br_Module = Module, 
         Br_Dominant_Tissue = Dominant_Tissue)

write.csv(map_at_clean, "tiss_map_at_clean.csv", row.names = FALSE)
write.csv(map_br_clean, "tiss_map_br_clean.csv", row.names = FALSE)

#read these back in and make tables
read.csv("/mnt/Chikorita/Brassica_rapa/manu_sup_files/fig4_sup_files/tiss_map_at_clean.csv") -> map_at_clean
read.csv("/mnt/Chikorita/Brassica_rapa/manu_sup_files/fig4_sup_files/tiss_map_br_clean.csv") -> map_br_clean

at_br_mod_tiss_map_clean_long <- rbind(map_at_clean, map_br_clean) %>%
  separate(Module,c("Species", "Module"), sep = "_" )

at_br_mod_tiss_map_clean_wide <- rbind(map_at_clean, map_br_clean) %>%
  separate(Module,c("Species", "Module"), sep = "_" ) %>%
  pivot_wider(
    names_from  = Species,
    values_from = Module,
    values_fn   = ~ paste(unique(.x), collapse = ", ")
  )