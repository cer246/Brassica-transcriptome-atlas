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

#Prep
br_topology %>%
  select(Gene, Module) -> modules_Br #42,830 genes

at_topology %>%
  select(Gene, Module) -> modules_At #33,063 genes

# Read the provided ortholog file
orthologs <- read_tsv("/mnt/Chikorita/Brassica_rapa/coexpression/core_files/br_at_homologs.txt") %>%
  dplyr::select(
    br_gene = `Gene stable ID`, 
    at_gene = `Arabidopsis thaliana gene stable ID`
  ) %>%
  # Remove rows where there is no Arabidopsis ortholog
  filter(!is.na(at_gene) & at_gene != "") %>%
  distinct() # Remove duplicates if any - 28,620 orthologous pairs

# Clean up ortholog table names for easier merging
colnames(orthologs)[1] <- "Br_GeneID"
colnames(orthologs)[2] <- "At_GeneID"

# Filter orthologs to keep only those present in our filtered WGCNA datasets
valid_orthologs <- orthologs %>%
  filter(Br_GeneID %in% modules_Br$Gene, 
         At_GeneID %in% modules_At$Gene) %>%
  select(Br_GeneID, At_GeneID) %>%
  distinct() # Remove duplicates if any - 25,008 orthologous pairs are present in the co-expression networks

# Merge Module Assignments into the Ortholog Table
comparison_df <- valid_orthologs %>%
  left_join(modules_Br, by = c("Br_GeneID" = "Gene")) %>%
  rename(Br_Module = Module) %>%
  left_join(modules_At, by = c("At_GeneID" = "Gene")) %>%
  rename(At_Module = Module) #now lets assign module placement for Brassica and Arabidopsis

## save what you have so far
write_tsv(comparison_df, "ortho_pairs_in_networks.tsv")

#### Assign modules to tissues used ortholog presence ####
## Statistical Overlap (Hypergeometric Test)
## looking to see if significant place of orthologs in similiar modules###

# Get unique module labels
at_mods <- sort(unique(comparison_df$At_Module)) #vector of the modules in Arabidopsis, 18
br_mods <- sort(unique(comparison_df$Br_Module)) #vector of the modules in Brassica, 13

# Initialize matrices for counts and p-values
overlap_counts <- matrix(0, nrow = length(at_mods), ncol = length(br_mods))
overlap_pvals <- matrix(1, nrow = length(at_mods), ncol = length(br_mods))
rownames(overlap_counts) <- at_mods; colnames(overlap_counts) <- br_mods
rownames(overlap_pvals) <- at_mods; colnames(overlap_pvals) <- br_mods

# Total universe size (number of ortholog pairs used)
n_total <- nrow(comparison_df) #25,008 orthologous pairs

# Loop to calculate overlaps - determine where these orthologous pairs lie in the coexpression networks
for (at_m in at_mods) {
  for (br_m in br_mods) {
    
    # Set A: Ortholog pairs where the At gene is in module 'at_m'
    set_A <- comparison_df %>% filter(At_Module == at_m) %>% nrow()
    
    # Set B: Ortholog pairs where the Br gene is in module 'br_m'
    set_B <- comparison_df %>% filter(Br_Module == br_m) %>% nrow()
    
    # Intersection: Pairs in both
    intersection <- comparison_df %>% filter(At_Module == at_m, Br_Module == br_m) %>% nrow()
    
    overlap_counts[at_m, br_m] <- intersection
    
    # Fisher's Exact Test (Hypergeometric)
    # phyper(q, m, n, k, lower.tail = FALSE)
    # q = intersection - 1
    # m = size of Set A
    # n = total - size of Set A
    # k = size of Set B
    if (set_A > 0 & set_B > 0) {
      pval <- phyper(intersection - 1, set_A, n_total - set_A, set_B, lower.tail = FALSE)
      overlap_pvals[at_m, br_m] <- pval
    }
  }
}

#### Visualization of module placement of orthologs ####
library(gplots)
# Create a heatmap of -log10(p-value) to show significant overlaps
log_pvals <- -log10(overlap_pvals + 1e-300) # avoid log(0)

# Generate PDF of the correspondence
#pdf("At_Br_Module_Correspondence.pdf", width = 12, height = 12)
heatmap.2(log_pvals,
          main = "At vs Br Module Correspondence\n(-log10 P-val)",
          trace = "none",
          col = colorRampPalette(c("white", "#9e2a2b"))(50),
          margins = c(12, 12),
          dendrogram = "both",
          xlab = "Brassica Modules",
          ylab = "Arabidopsis Modules",
          cellnote = overlap_counts, # Show the count of shared orthologs in cells
          notecol = "black",
          key.xlab = "-log10(P-value)",
          keysize = 1.0) 
#dev.off() 

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

write.csv(at_br_mod_tiss_map_clean_wide, "/mnt/Chikorita/Brassica_rapa/manu_sup_files/at_br_mod_tiss_map_clean.csv")
