#### General Exploration of Brassica Transcriptomic Data ####

## Figure 1, Supplemental Figures 2,3, and 4

##### Load packages ####
library(ggdendro)
library(cluster)
library(pheatmap)
library(stringr)
library(RColorBrewer)
library(viridis)
library(matrixStats)
library(ggpubr)
library(cowplot)
library(GenomicFeatures)
library(tidyverse)
library(patchwork)

##### Assign gene types ####
read_csv("/mnt/Chikorita/Brassica_rapa/coexpression/core_files/br_all_genes_of_int_type.csv") %>% dplyr::select(-1) -> br_all_genes_of_int

br_all_genes_of_int %>% 
  group_by(`Gene type`) %>%
  summarize(n = n())

# Number of genes expressed in each tissue
clean_tiss_tpm <- read_csv("/mnt/Chikorita/Brassica_rapa/coexpression/core_files/avg_br_tiss_atlas_tpm.csv")%>% 
  tibble::column_to_rownames(var = "Gene")

#7,737 genes are 0 tpm in all tissues
#12293 genes are less than 0.3 TPM in all tissues
clean_tiss_tpm[rowSums(clean_tiss_tpm <= 0) == ncol(clean_tiss_tpm), ] -> clean_0tpm

#so lets do at least 1 TPM in at least 1 Tissue
keep_clean_tiss_tpm <- rowSums((clean_tiss_tpm) >= 1) >= 1 #1 TPM in at least 1 tissues #starting with 56,609
clean_tiss_tpm_filt <- clean_tiss_tpm[keep_clean_tiss_tpm,] #Now we have 41,074

#### Tissue clustering ####
datExpr <- as.data.frame(t(clean_tiss_tpm_filt[order(apply(clean_tiss_tpm_filt,1,mad), decreasing = T)[1:41074],]))

hc       <- hclust(dist(datExpr), "ave") # heirarchal clustering
dendr    <- dendro_data(hc, type="rectangle") # convert for ggplot
clust    <- cutree(hc,k=10)                    # find 2 clusters
clust.df <- data.frame(label=names(clust), cluster=factor(clust))
# dendr[["labels"]] has the labels, merge with clust.df based on label column
dendr[["labels"]] <- merge(dendr[["labels"]],clust.df, by="label")
# plot the dendrogram; note use of color=cluster in geom_text(...)
sample_clust <- ggplot() + 
  geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) + 
  geom_text(data=dendr$labels, aes(x, y, label=label, hjust=0, color=clust), 
            size=3) +
  labs(x = "Tissue") +
  coord_flip() + 
  scale_y_reverse(expand=c(0.2, 0)) + 
  theme(axis.text.y = element_blank(),
        axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.line.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x = element_blank(),
        panel.background=element_rect(fill="white"),
        panel.grid=element_blank(),
        legend.position="none")

#### Get expressed genes ####
clean_tiss_tpm_filt %>%
  rownames_to_column("Gene")%>%
  pivot_longer(!Gene)%>% 
  left_join(., br_all_genes_of_int,
            by = c("Gene" = "GENEID")) %>%
  group_by(name,`Gene type`)%>%
  summarise(n_expressed = sum(value > 1)) -> n_expressed_tiss_gene_type #ensuring that we are still filtering for all genes that are expressed at 1 TPM in at least 1 Tissue

n_expressed_tiss_gene_type %>% na.omit() -> n_expressed_tiss_gene_type_omit

##### Save them ####
write_tsv(n_expressed_tiss_gene_type_omit, "/mnt/Chikorita/Brassica_rapa/coexpression/core_files/n_expressed_tiss_gene_type.tsv")

### Set the order for the bar graph using the clustering from above ####
order_vector <- c("dry_seed",
                  "mature_seed_24_hour_imbibed_A",
                  "mature_seed_72_hour_imbibed",
                  "stamens",
                  "petals",
                  "internode",
                  "young_roots",
                  "hypocotyls",
                  "mature_roots",
                  "unopened_flower_buds",
                  "carpels",
                  "5_day_seed",
                  "ovules",
                  "cotyledons",
                  "whole_seedlings",
                  "pedicel",
                  "late_adult_leaf",
                  "leaf_3_and_4",
                  "reproductive_leaf",
                  "sepals",
                  "adult_leaf",
                  "leaf_1_and_2")# Specify the desired order


ungroup(n_expressed_tiss_gene_type_omit) -> n_expressed_tiss_gene_type_omit

#Reorder the 'name' column using factor
n_expressed_tiss_gene_type_omit$name <- factor(n_expressed_tiss_gene_type_omit$name, levels = order_vector)

#Arrange the dataframe by the reordered 'name' column
df_ordered <- n_expressed_tiss_gene_type_omit %>% arrange(name)

bar_1TPM_expressed <- ggplot(data=df_ordered, aes(x=name, y=n_expressed, fill = `Gene type`)) +
  geom_bar(stat="identity", color = "#6C6C6C")+
  theme_minimal()+
  scale_fill_manual(values = c("#f94144","#f3722c","#f8961e","#f9c74f","#90be6d","#43aa8b","#4d908e","#577590")) +
  coord_flip()+
  theme(
    legend.position = "none",
    #axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 45)) +
  labs(x = "Tissue",
       y = "Number of Genes expressed at or above 1 TPM")+
  theme(axis.text.x = (element_text(size = 14, color = "black")),
        axis.text.y = (element_text(size = 14, color = "black")),
        axis.title.x = (element_text(size = 16, color = "black")),
        axis.title.y = (element_text(size = 16, color = "black")))


n_expressed_tiss_gene_type_omit %>%
  group_by(name) %>%
  summarise(sum = sum(n_expressed))%>%
  mutate(avg = mean(sum))

#### Put sample dendrogram and bar plot together ####
sample_clust | bar_1TPM_expressed

#### Calculate TAU scores for each gene type ####
#Calculate the TAU value for each gene
#Formula: Sum(1 - (x_i / max(x))) / (N - 1)
calc_tau <- function(x) {
  # Handle cases with all zeros to avoid division by zero
  if (max(x) == 0) return(NA)
  
  # Normalize by the maximal component value
  x_hat <- x / max(x)
  
  # Calculate Tau
  tau <- sum(1 - x_hat) / (length(x) - 1)
  return(tau)
}

# Identify tissue columns (remove ID columns)
tissue_cols <- colnames(clean_tiss_tpm_filt) 

# Prepare matrix for calculation
# We apply log2(TPM + 1) transformation
expression_matrix <- as.matrix(clean_tiss_tpm_filt[, tissue_cols])
log_expression <- log2(expression_matrix + 1)

# Calculate Tau for each gene (row)
# apply(matrix, 1, function) applies the function to each row
tau_values <- apply(log_expression, 1, calc_tau)

# Create a results dataframe
clean_tiss_tpm_filt %>%
  rownames_to_column("Gene") -> clean_tiss_tpm_filt_gene_id

tau_df <- data.frame(
    Gene = clean_tiss_tpm_filt_gene_id$Gene, # Ensure this matches your Gene ID column name
  Tau = tau_values
)

# Clean the type dataframe if necessary (remove index column)
final_df <- inner_join(tau_df, br_all_genes_of_int, by = c("Gene" = "GENEID"))

##### TAU summary stats ####
summary_stats <- final_df %>%
  group_by(`Gene type`) %>%
  summarise(
  Mean_Tau = mean(Tau, na.rm = TRUE),
  Median_Tau = median(Tau, na.rm = TRUE),
  Count = n()
)

print(summary_stats)

ggplot(final_df, aes(x = `Gene type`, y = Tau, fill = `Gene type`)) +
  # Create the violin plot
  geom_violin(trim = TRUE, alpha = 0.6) +
  # Overlay a narrow boxplot inside
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  scale_fill_brewer(palette = "Set2") + # Optional: Use a nice color palette
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    title = "Tissue Specificity (Tau) Distribution by Gene Class",
    y = "Tau Index (0=Broad, 1=Specific)",
    x = "Gene Class",
    fill = "Gene Class"
  )

#order the gene types
final_df %>%
  mutate(`Gene type order` = factor(`Gene type`,
                                    levels = c("antisense lncRNA", "intergenic lncRNA", "TE-associated lncRNA", "Housekeeping ncRNA precursor", "Polycistronic transcript", "Putative novel PCG", "protein_coding"))) -> final_df

ggplot(final_df, aes(x = `Gene type order`, y = Tau, fill = `Gene type order`)) + # fill=name allow to automatically dedicate a color for each group
  geom_violin(alpha = 1)+
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, alpha  = .75)+
  ylim(c(0,1))+
  #geom_jitter()+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45),
        legend.position = "none")+
  scale_fill_manual(values = c("#f97461", "#f8bc48", "#597e90", "#f39f4f", "#f9dd6d", "#55a6aa", "#85be7b"))+
  theme(axis.text.x = (element_text(size = 14, color = "black", angle = 45, hjust = 1)),
        axis.text.y = (element_text(size = 14, color = "black")),
        axis.title.x = (element_text(size = 16, color = "black")),
        axis.title.y = (element_text(size = 16, color = "black"))) -> updated_tau_gene_type

write.csv(final_df, "/mnt/Chikorita/Brassica_rapa/coexpression/core_files/gene_specificity_tau_scores.csv", row.names = FALSE)

##### Additional stats ####
#test for difference in distribution?
kruskal_result <- kruskal.test(Tau ~ `Gene type`, data = final_df)
kruskal_result #yes there is a significant difference in the distributions of tau scores

final_df %>%
  dplyr::select(-1,-4) -> final_df_stats

final_df_stats %>%
  setNames(c("Tau", "Gene_type")) -> final_df_stats

final_df_stats$Gene_type <- as.factor(final_df_stats$Gene_type)

final_df_stats$Gene_type <- as.factor(final_df_stats$Gene_type)
final_df_stats$Tau <- as.numeric(final_df_stats$Tau)

tau.kt <- kruskal.test(Tau ~ Gene_type, data = final_df_stats)#yes there is a significant difference in the distributions of tau scores
tau.aov <- aov(Tau ~ Gene_type,final_df_stats)
summary(tau.aov)
tau.tukey <- glht(tau.aov, linfct=mcp(Gene_type = "Tukey")) 

cld(tau.tukey)

tau_gene_type_stats %>%
  group_by(Gene_type)%>%
  summarize(n= n())


#### Supplemental Exploration ####

library(gapminder)
library(dplyr)
library(tidyr)
library(stringr)
library(ComplexUpset)
library(purrr)
library(gt)
library(mmtable2)
library(gridExtra)
library(gprofiler2)
library(tidyverse)
library(biomaRt)
library(topGO)
library(tidyverse)
#install.packages("curl")

#Load Data
tpm_df <- read.csv("/mnt/Chikorita/Brassica_rapa/coexpression/core_files/avg_br_tiss_atlas_tpm_only_expressed.csv", check.names = FALSE)
type_df <- read.csv("/mnt/Chikorita/Brassica_rapa/coexpression/core_files/br_all_genes_of_int_type.csv", check.names = FALSE)

#Get Tissue Order (from clustering)
data_matrix <- tpm_df[, -c(1, 2)]
binary_matrix <- as.matrix(data_matrix > 1) + 0
tissue_dist <- dist(t(binary_matrix), method = "binary")
tissue_hclust <- hclust(tissue_dist, method = "average")
ordered_tissues <- tissue_hclust$labels[tissue_hclust$order]

#Process Data
#Find Max Tissue for each Gene
max_counts_df <- inner_join(tpm_df, type_df, by = c("Gene" = "GENEID")) %>%
  pivot_longer(
    cols = colnames(data_matrix),
    names_to = "Tissue",
    values_to = "TPM"
  ) %>%
  group_by(Gene) %>%
  slice_max(TPM, n = 1, with_ties = FALSE) %>% # Keep only the max tissue
  ungroup() %>%
  # Filter out genes with 0 max expression if needed
  filter(TPM > 0) %>%
  # Count
  group_by(`Gene type`, Tissue) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  # Pivot to Wide Format
  pivot_wider(names_from = Tissue, values_from = Count, values_fill = 0)

#Reorder Columns
#Ensure columns follow the clustered order
max_counts_ordered <- max_counts_df %>%
  dplyr::select(`Gene type`, all_of(ordered_tissues))

#View and Save
print(max_counts_ordered)
write.csv(max_counts_ordered, "/mnt/Chikorita/Brassica_rapa/manu_sup_files/supp_fig2_files/max_expression_counts_table.csv", row.names = FALSE)


### Determining what genes are expressed, create a supplemental table ###
tpm_all <- read.csv("/mnt/Chikorita/Brassica_rapa/coexpression/core_files/avg_br_tiss_atlas_tpm.csv", check.names = FALSE)
type_df <- read.csv("/mnt/Chikorita/Brassica_rapa/coexpression/core_files/br_all_genes_of_int_type.csv", check.names = FALSE)
#its already here but do it again

#Remove the first column (index) from type_df
type_df <- type_df[, -1] 

##### Determine Expression Status ####
#Identify tissue columns (all columns except 'Gene')
tissue_cols <- colnames(tpm_all)[colnames(tpm_all) != "Gene"]

status_df <- tpm_all %>%
  rowwise() %>%
  mutate(Max_TPM = max(c_across(all_of(tissue_cols)))) %>%
  ungroup() %>%
  mutate(Status = case_when(
    Max_TPM > 1 ~ "Expressed > 1 TPM",
    Max_TPM > 0 ~ "Expressed < 1 TPM",
    TRUE ~ "Not Detected"
  )) %>%
  dplyr::select(Gene, Status)

#Merge and Summarize
summary_table <- inner_join(status_df, type_df, by = c("Gene" = "GENEID")) %>%
  group_by(`Gene type`, Status) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  pivot_wider(names_from = Status, values_from = Count, values_fill = 0) %>%
  mutate(Total = rowSums(across(where(is.numeric))))

#View and Save
print(summary_table)
write.csv(summary_table, "/mnt/Chikorita/Brassica_rapa/manu_sup_files/supp_fig2_files/detailed_expression_status.csv", row.names = FALSE)

style_list <- list(cell_borders(sides = "top",color = "grey"))
summary_table %>%
  pivot_longer(!`Gene type`, names_to = "express_status", values_to = "num_genes") -> summary_table_long

gm_table <- summary_table_long %>% 
  mmtable(cells = num_genes) +
  #header_left(`Gene type`) +
  header_top(express_status)+ 
  header_left_top(`Gene type`)  +
  header_format(`Gene type`, scope = "table", style = style_list)

gridExtra::grid.table(gm_table)
grid.table(summary_table)
dev.off()



##### GO Enrichement of Not Detected Transcripts ####
#Filter for Protein Coding Genes
#Select only genes annotated as 'protein_coding'
pc_genes <- type_df %>%
  filter(`Gene type` == "protein_coding") %>%
  pull(GENEID)

#Identify Unexpressed Genes
#Identify tissue columns (exclude 'Gene')
tissue_cols <- colnames(tpm_all)[colnames(tpm_all) != "Gene"]

#Filter the TPM table for genes where the maximum expression is 0
unexpressed_genes <- tpm_all %>%
  rowwise() %>%
  mutate(Max_TPM = max(c_across(all_of(tissue_cols)))) %>%
  ungroup() %>%
  filter(Max_TPM == 0) %>%
  pull(Gene)

#Find genes that are BOTH protein_coding AND unexpressed
target_genes <- intersect(pc_genes, unexpressed_genes)

#Create Output Table with Cleaned IDs
output_df <- data.frame(Original_ID = target_genes) %>%
  mutate(Cleaned_ID = gsub("_BraROA$", "", Original_ID))

#Save to CSV
write.csv(output_df, "/mnt/Chikorita/Brassica_rapa/manu_sup_files/unexpressed_protein_coding_genes.csv", row.names = FALSE)


###Attempting with the biomart from ensembl plants ###
# Connect to the plants mart
plants_mart <- useMart(biomart = "plants_mart", 
                       host = "https://plants.ensembl.org")

datasets <- listDatasets(plants_mart)
brapa_dataset <- datasets[grep("ro18", datasets$dataset), ]
print(brapa_dataset)

my_mart <- useDataset("brro18_eg_gene", mart = plants_mart)

#Retrieve GO Annotations
print("Retrieving GO annotations...")
go_mapping <- getBM(
  attributes = c("ensembl_gene_id", "go_id"),
  mart = my_mart
)

#Filter empty GO terms and create list
gene_to_go <- go_mapping %>%
  filter(go_id != "") %>%
  dplyr::select(ensembl_gene_id, go_id)

gene2GO_list <- split(gene_to_go$go_id, gene_to_go$ensembl_gene_id)

#Run Enrichment
#Load unexpressed genes
unexpressed_df <- read.csv("/mnt/Chikorita/Brassica_rapa/manu_sup_files/unexpressed_protein_coding_genes.csv")
interesting_genes <- as.character(unexpressed_df$Original_ID)
gene_universe <- names(gene2GO_list)
geneList <- factor(as.integer(gene_universe %in% interesting_genes))
names(geneList) <- gene_universe

#Create topGO object
GOdata <- new("topGOdata",
              description = "Unexpressed Genes in R-o-18", 
              ontology = "BP", 
              allGenes = geneList, 
              annot = annFUN.gene2GO, 
              gene2GO = gene2GO_list,
              nodeSize = 10)

#Test
result_fisher <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")

#Results
all_res <- GenTable(GOdata, 
                    classicFisher = result_fisher, 
                    orderBy = "classicFisher", 
                    topNodes = 50)

#Filter p < 0.05
all_res$classicFisher_numeric <- as.numeric(gsub("[< ]", "", all_res$classicFisher))
sig_results <- all_res %>% filter(classicFisher_numeric < 0.05)

print(head(sig_results))
write.csv(sig_results, "/mnt/Chikorita/Brassica_rapa/manu_sup_files/unexpressed_genes_topGO_enrichment.csv", row.names = FALSE)

###Visualization###
#Prepare Data
#Ensure p-values are numeric and create a -log10 score for plotting
#We limit to top 20 terms
plot_data <- sig_results %>%
  mutate(logP = -log10(classicFisher_numeric)) %>%
  mutate(Term = factor(Term, levels = rev(Term))) %>% # Reorder for plotting
  head(11)

#Create Dot Plot
ggplot(plot_data, aes(x = logP, y = Term, size = Significant, color = classicFisher_numeric)) +
  geom_point(alpha = 1) +
  scale_color_gradient(low = "#d16266" ,high = "#3b4f7d"  , trans = "reverse") + # Lower p-value = Red
  labs(
    title = "GO Enrichment (Biological Process)",
    subtitle = "Top 11 Terms for Unexpressed Genes",
    x = "-log10(P-value)",
    y = "GO Term",
    size = "Gene Count",
    color = "P-value"
  ) +
  theme_minimal() +
  theme(axis.text.x = (element_text(size = 14, color = "black")),
        axis.text.y = (element_text(size = 14, color = "black")),
        axis.title.x = (element_text(size = 16, color = "black")),
        axis.title.y = (element_text(size = 16, color = "black")))

##### Broad vs. Specific PCG Functional Association ####
# Merge expression data with gene types
merged_df <- inner_join(tpm_df, type_df, by = c("Gene" = "GENEID"))

# Filter for only protein coding genes
pc_df <- merged_df %>%
  dplyr::select(-1) %>%
  dplyr::filter(`Gene type` == "protein_coding")

#Re-calculate Tau Specificity Index for PCGs
# Identify tissue columns (exclude Gene, Gene type, and index columns)
tissue_cols #you already have this 
expression_matrix <- as.matrix(pc_df[, tissue_cols])

# Log-transform (log2(TPM + 1))
log_expression <- log2(expression_matrix + 1)

# Function to calculate Tau
calc_tau <- function(x) {
  if (max(x) == 0) return(NA)
  x_hat <- x / max(x)
  tau <- sum(1 - x_hat) / (length(x) - 1)
  return(tau)
}

# Apply function to every gene
tau_values <- apply(log_expression, 1, calc_tau)

# Add Tau to the dataframe
pc_df$Tau <- tau_values

# Remove any genes with NA Tau (e.g., 0 expression everywhere)
pc_df <- pc_df %>% filter(!is.na(Tau))

#Divide PCGs into expression quartiles
# Create quartiles based on Tau
pc_df$Quartile <- ntile(pc_df$Tau, 4)

# Map integers to quartile labels
pc_df <- pc_df %>%
  mutate(Quartile_Label = case_when(
    Quartile == 1 ~ "Q1_Broad",
    Quartile == 2 ~ "Q2",
    Quartile == 3 ~ "Q3",
    Quartile == 4 ~ "Q4_Specific"
  ))

# Extract and Save Lists
# Q1: Broadly Expressed (Low Tau)
broad_genes <- pc_df %>%
  filter(Quartile == 1) %>%
  dplyr::select(Gene) %>%
  mutate(Cleaned_ID = gsub("_BraROA$", "", Gene)) # Clean for GO tools

# Q4: Tissue Specific (High Tau)
specific_genes <- pc_df %>%
  filter(Quartile == 4) %>%
  dplyr::select(Gene) %>%
  mutate(Cleaned_ID = gsub("_BraROA$", "", Gene))

# Save to CSV
write.csv(broad_genes, "/mnt/Chikorita/Brassica_rapa/manu_sup_files/tau_quartile1_broad_genes.csv", row.names = FALSE)
write.csv(specific_genes, "/mnt/Chikorita/Brassica_rapa/manu_sup_files/tau_quartile4_specific_genes.csv", row.names = FALSE)

# Print Summary
print(paste("Saved", nrow(broad_genes), "broad genes and", nrow(specific_genes), "specific genes."))
print(summary(pc_df$Tau))

###### GO Enrichment on these two classes (i.e Broad and Specific) ####
plants_mart <- useMart(biomart = "plants_mart", host = "https://plants.ensembl.org")
my_mart <- useDataset("brro18_eg_gene", mart = plants_mart)

# Retrieve GO mappings for the whole genome (Universe)
go_mapping <- getBM(attributes = c("ensembl_gene_id", "go_id"), mart = my_mart)
gene_to_go <- go_mapping %>% filter(go_id != "")
gene2GO_list <- split(gene_to_go$go_id, gene_to_go$ensembl_gene_id)
gene_universe <- names(gene2GO_list)

run_enrichment <- function(file_name, label) {
  
  # Load gene list
  df <- read.csv(file_name)
  
  # Remove version suffix if present (e.g., .1)
  interesting_genes <- gsub("\\.\\d+$", "", df$Gene)
  
  # Define gene list factor
  geneList <- factor(as.integer(gene_universe %in% interesting_genes), levels = c(0, 1))
  names(geneList) <- gene_universe
  
  # Check overlap
  if (sum(geneList == 1) == 0) {
    warning(paste("No matching genes found for", label))
    return(NULL)
  }
  
  # Create topGO object (Biological Process)
  GOdata <- new("topGOdata",
                description = label, ontology = "BP",
                allGenes = geneList, annot = annFUN.gene2GO,
                gene2GO = gene2GO_list, nodeSize = 10)
  
  # Run Test
  result_fisher <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
  
  # Get Table
  res_table <- GenTable(GOdata, classicFisher = result_fisher, 
                        orderBy = "classicFisher", topNodes = 50)
  
  # Clean p-values
  res_table$classicFisher_numeric <- as.numeric(gsub("[< ]", "", res_table$classicFisher))
  res_table$Group <- label
  
  return(res_table)
}

###### Run for Broadly Expressed (Q1) ####
res_broad <- run_enrichment("/mnt/Chikorita/Brassica_rapa/manu_sup_files/tau_quartile1_broad_genes.csv", "Broadly Expressed (Q1)")

###### Run for Tissue Specific (Q4) ####
res_specific <- run_enrichment("/mnt/Chikorita/Brassica_rapa/manu_sup_files/tau_quartile4_specific_genes.csv", "Tissue Specific (Q4)")

#Combine results for comparison
combined_results <- bind_rows(
  res_broad %>% mutate(Type = "Broad"),
  res_specific %>% mutate(Type = "Specific")
) %>% filter(classicFisher_numeric < 0.05)

print(head(combined_results))
write.csv(combined_results, "/mnt/Chikorita/Brassica_rapa/manu_sup_files/quartile_enrichment_results.csv", row.names = FALSE)

###### Visualize Top Terms #####
top_terms <- combined_results %>%
  group_by(Type) %>%
  slice_min(classicFisher_numeric, n = 10)
#Create Bar Plot
ggplot(top_terms, aes(x = reorder(Term, -classicFisher_numeric), 
                      y = -log10(classicFisher_numeric), 
                      fill = Type)) +
  geom_col(color = "black") +
  coord_flip() +
  facet_wrap(~Type, scales = "free") +
  labs(title = "Top GO Terms: Broad vs Specific", x = "GO Term", y = "-log10(P-value)") +
  theme_minimal()+
  scale_fill_manual(values = c("#92969a", "#a297a9"))+
  theme(axis.text.x = (element_text(size = 14, color = "black")),
        axis.text.y = (element_text(size = 14, color = "black")),
        axis.title.x = (element_text(size = 16, color = "black")),
        axis.title.y = (element_text(size = 16, color = "black")))

##### Clustered and Stacked Barplot - Just LncRNAs ####

tpm_df <- read.csv("/mnt/Chikorita/Brassica_rapa/coexpression/core_files/avg_br_tiss_atlas_tpm_only_expressed.csv", check.names = FALSE)
type_df <- read.csv("/mnt/Chikorita/Brassica_rapa/coexpression/core_files/br_all_genes_of_int_type.csv", check.names = FALSE)

# Get Tissue Order (Same as before)
data_matrix <- tpm_df[, -c(1, 2)]
binary_matrix <- as.matrix(data_matrix > 1) + 0
tissue_dist <- dist(t(binary_matrix), method = "binary")
tissue_hclust <- hclust(tissue_dist, method = "average")
ordered_tissues <- tissue_hclust$labels[tissue_hclust$order]

# Process Data for Specific lncRNAs
merged_df <- inner_join(tpm_df, type_df, by = c("Gene" = "GENEID"))

lncRNA_data <- merged_df %>%
  # Filter for only the two types of interest
  filter(`Gene type` %in% c("antisense lncRNA", "intergenic lncRNA", "TE-associated lncRNA")) %>%
  pivot_longer(
    cols = colnames(data_matrix), # Pivot tissue columns
    names_to = "Tissue",
    values_to = "TPM"
  ) %>%
  filter(TPM > 1) %>%
  group_by(Tissue, `Gene type`) %>%
  summarise(Count = n(), .groups = 'drop')

#Apply Order
lncRNA_data$Tissue <- factor(lncRNA_data$Tissue, levels = ordered_tissues)

###### Visualize ####
ggplot(lncRNA_data, aes(x = Tissue, y = Count, fill = `Gene type`)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#f37462", "#f7bc48", "#5a7e90"))+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(
    title = "lncRNAs Expressed > 1 TPM",
    x = "Tissue",
    y = "Number of Transcripts",
    fill = "Gene Class"
  )+
  theme(axis.text.x = (element_text(size = 14, color = "black")),
        axis.text.y = (element_text(size = 14, color = "black")),
        axis.title.x = (element_text(size = 16, color = "black")),
        axis.title.y = (element_text(size = 16, color = "black")),
        legend.position = "left")

### Overlap of gene expression between tissues####
tpm_df <- read.csv("/mnt/Chikorita/Brassica_rapa/coexpression/core_files/avg_br_tiss_atlas_tpm_only_expressed.csv", check.names = FALSE)

# Load the mapping table
mapping_df <- read.csv("/mnt/Chikorita/Brassica_rapa/coexpression/core_files/tissue_collapsing_cer.csv", header = FALSE, 
                       col.names = c("Original", "Group"), 
                       stringsAsFactors = FALSE)

if (colnames(mapping_df)[1] == "") {
  mapping_df <- mapping_df[, -1]
}

if (colnames(tpm_df)[1] == "") {
  tpm_df <- tpm_df[, -1]
}

# Process and Collapse Tissues
tpm_long <- tpm_df %>%
  pivot_longer(
    cols = -c(1), 
    names_to = "Original",
    values_to = "TPM"
  )

# Join with the mapping table and perform the collapse
tpm_collapsed <- tpm_long %>%
  left_join(mapping_df, by = "Original") %>%
  # 1. Filter out tissues marked as "Omit" or those not in the mapping file
  filter(Group != "Omit", !is.na(Group)) %>%
  # 2. Group by Gene and the NEW Group Name
  group_by(Gene, Group) %>%
  # 3. Calculate the Mean TPM for the group
  summarise(Mean_TPM = mean(TPM), .groups = 'drop') %>%
  # Reshape back to 'Wide' format (Matrix-like)
  pivot_wider(
    names_from = Group,
    values_from = Mean_TPM
  )

# Prepare Data for Clustering
data_matrix <- as.matrix(tpm_collapsed[, -1])

# Set row names to Gene IDs (optional, but good for tracking)
rownames(data_matrix) <- tpm_collapsed$Gene

# Create Binary Matrix
binary_matrix <- (data_matrix > 1) + 0
tissue_dist <- dist(t(binary_matrix), method = "binary")
tissue_hclust <- hclust(tissue_dist, method = "average")

#Plot Dendrogram
par(mar = c(4, 4, 4, 4)) 
plot(tissue_hclust, 
     main = "Cluster Dendrogram (Collapsed Tissues)",
     xlab = "", 
     sub = "",
     ylab = "Jaccard Distance",
     hang = -1) # Aligns labels at the bottom


# Create Binary Matrix
binary_data <- tpm_collapsed %>%
  mutate(across(-Gene, ~ as.integer(. > 1))) # Convert numeric TPM to 0/1

# Merge with Gene Type info
# Ensure we keep only genes that match our type file
upset_df <- inner_join(binary_data, type_df, by = c("Gene" = "GENEID")) %>% select(-8)

# Define the list of sets (the tissue group names)
tissue_sets <- colnames(tpm_collapsed)[-1] # Exclude 'Gene'

# FILTER: Remove genes "outside of known sets" (All Zeros)
# We calculate the row sum for the tissue columns only. 
# If sum == 0, the gene is not expressed in any of the kept groups.
upset_df_filtered <- upset_df %>%
  filter(rowSums(select(., all_of(tissue_sets))) > 0)

# 4. Generate Stacked UpSet Plot
gene_type_colors <- c("antisense lncRNA" = "#f37462", 
                      "Housekeeping ncRNA precursor" = "#f39f50", 
                      "intergenic lncRNA" = "#f7bc48",
                      "Polycistronic transcript" = "#f9dd6d" ,
                      "protein_coding" = "#85be7b",
                      "Putative novel PCG" = "#55a6aa" ,
                      "TE-associated lncRNA" = "#5a7e90")
##### Plot #####
upset(
  upset_df_filtered,
  tissue_sets,
  name = "Tissue Groups",
  width_ratio = 0.1,
  min_size = 350,
  
  base_annotations = list(
    'Intersection size' = intersection_size(
      counts = TRUE,
      mapping = aes(fill = `Gene type`)) +
      scale_fill_manual(values = gene_type_colors) +
      theme_minimal() +
      theme(
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 16, color = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "none"
      )
  ),
  
  themes = upset_modify_themes(
    list(
      'intersections_matrix' = theme(
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 16, color = "black")
      )
    )
  )
)

#install.packages("remotes")
#remotes::install_version("ggplot2", version = "3.5.1", repos = "http://cran.us.r-project.org")


#### Just lncRNAs #####
upset_df_filtered_lncRNA <- upset_df %>%
  filter(rowSums(select(., all_of(tissue_sets))) > 0) %>%
  filter(`Gene type` %in% c("antisense lncRNA","TE-associated lncRNA","intergenic lncRNA"))

##### Plot ####
upset(
  upset_df_filtered_lncRNA, # Use the filtered dataset
  tissue_sets,
  name = "Tissue Groups",
  width_ratio = 0.1,
  min_size = 350,
  
  base_annotations = list(
    'Intersection size' = intersection_size(
      counts = TRUE,
      mapping = aes(fill = `Gene type`)) +
      scale_fill_manual(values = gene_type_colors) +
      theme_minimal() +
      theme(
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 16, color = "black"),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "none"
      )
  ),
  
  themes = upset_modify_themes(
    list(
      'intersections_matrix' = theme(
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.title.y = element_text(size = 16, color = "black")
      )
    )
  )
)

#### Tissue of Maximal Expression by Gene type ####
tpm_df <- read.csv("/mnt/Chikorita/Brassica_rapa/coexpression/core_files/avg_br_tiss_atlas_tpm_only_expressed.csv", check.names = FALSE)
type_df <- read.csv("/mnt/Chikorita/Brassica_rapa/coexpression/core_files/br_all_genes_of_int_type.csv", check.names = FALSE)

#Tissue Order
data_matrix <- tpm_df[, -c(1, 2)]
binary_matrix <- as.matrix(data_matrix > 1) + 0
tissue_dist <- dist(t(binary_matrix), method = "binary")
tissue_hclust <- hclust(tissue_dist, method = "average")
ordered_tissues <- tissue_hclust$labels[tissue_hclust$order]

#Process Data
plot_data <- inner_join(tpm_df, type_df, by = c("Gene" = "GENEID")) %>%
  pivot_longer(
    cols = colnames(data_matrix),
    names_to = "Tissue",
    values_to = "TPM"
  ) %>%
  group_by(Gene) %>%
  slice_max(TPM, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  # Count
  group_by(`Gene type`, Tissue) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  # Percentage
  group_by(`Gene type`) %>%
  mutate(Percentage = Count / sum(Count) * 100) %>%
  ungroup()

#Apply Order
plot_data$Tissue <- factor(plot_data$Tissue, levels = ordered_tissues)

my_colors6 <- c(
  "#5667B3", "#5C7EB1", "#5F93B7", "#66A6B9", "#6FB6B4",
  "#7EC3A5", "#96CE96", "#AED58A", "#C6D880", "#DCD97A",
  "#EFE07C", "#F0CF82", "#EBB689", "#E29A8E", "#D58190",
  "#C46F97", "#AD6EA2", "#9373AD", "#787BB4", "#6486B8",
  "#5A92BB", "#5FA0B8", "#6AACB2", "#7BB5AA", "#8DBDA2"
)

# Alternatively, extend a Brewer palette:
# nb_colors <- length(unique(plot_data$Tissue))
# my_colors <- colorRampPalette(brewer.pal(12, "Paired"))(nb_colors)

plot_data %>%
  mutate(`Gene type order` = factor(`Gene type`,
                                    levels = c("antisense lncRNA", "intergenic lncRNA", "TE-associated lncRNA", "Housekeeping ncRNA precursor", "Polycistronic transcript", "Putative novel PCG", "protein_coding"))) -> plot_data_ordered_df


##### Plot ####
ggplot(plot_data_ordered_df, aes(x = `Gene type order`, y = Percentage, fill = Tissue)) +
  geom_bar(stat = "identity", color = "black", linewidth = .2) +
  scale_fill_manual(values = rev(my_colors6)) + # Use the custom palette
  theme_minimal() +
  labs(
    title = "Percentage of Transcripts with Maximal Expression in Each Tissue",
    x = "Gene Class",
    y = "Percentage of Transcripts (%)",
    fill = "Tissue of Max Expression")+
  theme(axis.text.x = (element_text(size = 14, color = "black",angle = 45, hjust = 1)),
        axis.text.y = (element_text(size = 14, color = "black")),
        axis.title.x = (element_text(size = 16, color = "black")),
        axis.title.y = (element_text(size = 16, color = "black")),
        legend.position = "none")+ 
  guides(fill = guide_legend(nrow = 3))