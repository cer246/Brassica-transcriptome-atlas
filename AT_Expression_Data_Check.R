#### Arabidopsis Expression Data Preparation ####

#### Expression Data from Mergner et al. 2020 - doi: 10.1038/s41586-020-2094-2

##### Load packages #####
library(tidyverse)
library(DESeq2)
library(txdbmaker)
library(GenomicFeatures)

#### Generate counts/expression dataframes ####
# load annotation file as TxDb
AT_TxDb <- txdbmaker::makeTxDbFromGFF("/mnt/Chikorita/Brassica_rapa/coexpression/core_files/at_oct_annotation_gffread_strand.gff3")
# get annotaion ready for tximport, here we are mapping transcript names to the gene

k <- AnnotationDbi::keys(AT_TxDb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(AT_TxDb, k, "GENEID", "TXNAME")

#Pull in our genetypes that we are interested in - the pcgs, lncRNAs, and kyles hc LncRNAs
#kyles lncRNAs hc 
at_hc_lncRNA <- read_csv("/mnt/Chikorita/Brassica_rapa/coexpression/core_files/at_hc_lncRNA.txt", 
                         col_names = FALSE) %>%
  setNames(c("name"))

Araport11_gene_type <- read_delim("/mnt/Chikorita/Brassica_rapa/coexpression/core_files/araport11_gene_type.txt", 
                                  delim = "\t", escape_double = FALSE, 
                                  col_names = TRUE, trim_ws = TRUE)

Araport11_gene_type %>%
  select(1) -> all_araport_gene_ids

all_araport_gene_ids -> at_araport_retain

#filter tx2gene for the genes of interest
tx2gene %>%
  left_join(.,
            Araport11_gene_type,
            by = c("GENEID" = "name")) -> at_araport_retain_tx2gene

at_araport_retain_tx2gene %>%
  mutate(
    gene_model_type = case_when(
      is.na(gene_model_type) & str_detect(GENEID, "TE") ~ "transposable element",
      is.na(gene_model_type) & str_detect(GENEID, "XLOC") ~ "kyle_hc_intergenic_lncRNA",
      is.na(gene_model_type) & str_detect(GENEID, "YLOC") ~ "kyle_hc_intergenic_lncRNA",
      is.na(gene_model_type) & str_detect(GENEID, "AT5G08185") ~ "miRNA_primary_transcript",
      is.na(gene_model_type) & str_detect(GENEID, "AT5G45105") ~ "pseudogene",
      TRUE ~ gene_model_type
    )
  ) -> at_all_int_gene_type

#save this 
write_csv(at_all_int_gene_type, "/mnt/Chikorita/Brassica_rapa/coexpression/core_files/at_all_genes_of_int_type.csv")

tx2gene %>%
  inner_join(.,
             at_hc_lncRNA,
             by = c("TXNAME" = "name")) -> at_hc_lncRNA_retain_tx2gene #only keeping high confidence HC lncRNAs

rbind(at_araport_retain_tx2gene,at_hc_lncRNA_retain_tx2gene) -> tx2gene_retain


#read in metadata
tiss_atlas_info <- read.csv("/mnt/Chikorita/Brassica_rapa/coexpression/core_files/At_tiss_atlas_infoUpdated.csv")

#Create an object with your file names
tiss_atlas_dir <- "/mnt/Chikorita/Brassica_rapa/coexpression/AT/all_ens_lncRNA/quant_output"
tiss_atlas_Qfiles<- file.path(tiss_atlas_dir, tiss_atlas_info$Run, "quant.sf")
head(tiss_atlas_Qfiles)

file.exists(tiss_atlas_Qfiles)

# Combine your quant files and annotation file to generate transcript/gene-level estimation
txi_tiss_atlas <- tximport(files = tiss_atlas_Qfiles, type = "salmon", tx2gene = tx2gene_retain, countsFromAbundance = "lengthScaledTPM")

colnames(txi_tiss_atlas$counts) <- tiss_atlas_info$Run
colnames(txi_tiss_atlas$abundance) <- tiss_atlas_info$Run
colnames(txi_tiss_atlas$length) <- tiss_atlas_info$Run

ddsTxi_tiss_atlas <- DESeqDataSetFromTximport(txi_tiss_atlas,
                                              colData = tiss_atlas_info,
                                              design = ~Condition) #here condition is tissue


colnames(txi_tiss_atlas$counts) <- tiss_atlas_info$Run
colnames(txi_tiss_atlas$abundance) <- tiss_atlas_info$Run

txi_tiss_atlas_counts <- txi_tiss_atlas$counts
txi_tiss_atlas_tpm <- txi_tiss_atlas$abundance

#### Save raw counts and TPM ####
write.csv(txi_tiss_atlas_counts, "/mnt/Chikorita/Brassica_rapa/coexpression/core_files/at_tiss_atlas_raw_counts.csv")
write.csv(txi_tiss_atlas_tpm, "/mnt/Chikorita/Brassica_rapa/coexpression/core_files/at_tiss_atlas_tpm.csv")


#Do the differential expression analysis
dds_at_tiss <- DESeq(ddsTxi_at_tiss)

vsd_tiss <- vst(dds_tiss_atlas, blind = F) 

dds_counts_tiss_atlas <- counts(dds_tiss_atlas, normalized=TRUE)
vsd_counts_tiss_atlas <- data.frame(assay(vsd_tiss, normalized=TRUE)) #focus on this one first

#normalized counts
dds <- estimateSizeFactors(dds_at_tiss)
at_tiss_norm_counts <- counts(dds, normalized=TRUE)


#### Save the normalized deseq counts and vsd counts ##### 
dds_counts_tiss_atlas %>% 
  as.data.frame() %>%
  write_csv(., "/mnt/Chikorita/Brassica_rapa/coexpression/AT/all_ens_lncRNA/dds_counts_tiss_atlas.csv")
vsd_counts_tiss_atlas %>%
  tibble::rownames_to_column(., var = "Gene") %>%
  as.data.frame() %>%
  write_csv(., "/mnt/Chikorita/Brassica_rapa/coexpression/AT/all_ens_lncRNA/vsd_counts_tiss_atlas.csv")

#### Plot PCA ####
vsd <- vst(dds_at_tiss, blind = F)
plotPCA <- function (object, intgroup=c("~Condition"), ntop = 3000, returnData=TRUE) 
{
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, 
                                               drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3], PC4 = pca$x[, 4], group = group, 
                  intgroup.df, name = colnames(object))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:4]
    return(d)
  }
  ggplot(data = d, aes_string(x = "PC3", y = "PC4", color = "Condition")) + 
    geom_point(size = 3) + xlab(paste0("PC3: ", round(percentVar[1] * 
                                                        100), "% variance")) + ylab(paste0("PC4: ", round(percentVar[2] * 
                                                                                                            100), "% variance")) + coord_fixed()
}


### extract the PC1, PC2, PC3, and PC4 data
pcaData <- plotPCA(vsd, intgroup=c("Condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
percentVar

nb.cols <- 28
mycolors <- colorRampPalette(brewer.pal(8, "Spectral"))(nb.cols)

unique(pcaData$Condition)

condition_colors <- c(
  "dry_seed" = "#eae2b7",
  "imbibed_seed" = "#e9c46a",
  "cotyledon" = "#a7c957",
  "hypocotyl" = "#6a994e",
  "seedling" = "#386641",
  "root" = "#dda15e",
  "root_tip" = "#bc6c25",
  "stem_internode" = "#1e6091",
  "stem_node" = "#0077b6",
  "rosette_leaf" = "#168aad",
  "adult_vascular_leaf" = "#34a0a4",
  "cauline_leaf" = "#52b69a",
  "petiole" = "#d9ed92",
  "flower_pedicel" = "#eeef20",
  "unopened_flower_bud" = "#ffc9b9",
  "budding_flower" = "#ffb3c6",
  "flower" = "#fb6f92",
  "sepal" = "#ff5c8a",
  "petal" = "#b388eb",
  "carpel" = "#973aa8",
  "stamen" = "#fc2f00",
  "pollen" = "#ea7317" ,
  "5_day_seed" = "#ff6b35",
  "fruit_septum" = "#dec0f1" ,
  "fruit_valve" = "#b79ced",
  "silique" = "#957fef" ,
  "plant_callus" = "#ef7a85"
)

DESeq2::plotPCA(vsd, intgroup = c("Condition")) + theme_linedraw()+ scale_color_manual(values = condition_colors)+ 
  theme(axis.text.x = (element_text(size = 14, color = "black")),
        axis.text.y = (element_text(size = 14, color = "black")),
        axis.title.x = (element_text(size = 16, color = "black")),
        axis.title.y = (element_text(size = 16, color = "black"))) -> at_illumina_pca

#### Plot Hierarchial clustering ####
vsd <- vst(dds_at_tiss, blind = F)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$Run)
colnames(sampleDistMatrix) <- paste(vsd$Condition)
colors <- colorRampPalette( rev(brewer.pal(9, "YlGnBu")) )(100)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col=colors,
         border_color = "gray",
         show_rownames = TRUE, 
         show_colnames = TRUE,
         angle_col = 45) -> at_illumina_distance_matrix


#### Generate Expression Dataframes ####
#pull in metadata
at_tissue_metadata <- read_csv("/mnt/Chikorita/Brassica_rapa/coexpression/core_files/At_tiss_atlas_infoUpdated.csv")

#define each tissue
at_tissue_metadata %>%
  dplyr::select(Condition) %>%
  unique()%>%
  as.vector() -> tissues

# Reshape the data into a longer format
df_long <- vsd_counts_tiss_atlas %>%
  rownames_to_column("Gene")%>%
  pivot_longer(!Gene, names_to = "sample", values_to = "vsd_counts")%>%
  dplyr::select(1,3,2)%>%
  left_join(.,
            at_tissue_metadata,
            by = c("sample" = "Run"))%>%
  dplyr::select(Gene,vsd_counts,sample,Condition) # Extract tissue names

# Calculate the mean vsd_counts for each tissue
mean_at_tiss_vsd_counts <- df_long %>%
  group_by(Gene,Condition) %>%
  summarise(mean_vsd_counts = mean(vsd_counts)) %>%
  ungroup() %>%
  pivot_wider(names_from = Condition, values_from = mean_vsd_counts)

#alright how many genes have 0 tpm in in all tissues
mean_at_tiss_vsd_counts %>%
  column_to_rownames("Gene") -> mean_at_tiss_vsd_counts_test 


#cool none of the genes are 0 tpm in all tissues
mean_at_tiss_vsd_counts_test[rowSums(mean_at_tiss_vsd_counts_test == 0) == ncol(mean_at_tiss_vsd_counts_test), ]

#what about have 1 TPM in at least 1 or 2 tissues
mean_at_tiss_vsd_counts_test[rowSums(mean_at_tiss_vsd_counts_test >= 1) >= 2, ]

##### Save avg VSD ####
write_csv(mean_at_tiss_vsd_counts,"/mnt/Chikorita/Brassica_rapa/coexpression/AT/all_ens_lncRNA/avg_vsd_clean_at_tiss_atlas_at_all_ens_lncRNA.csv")


#Generate Average TPM 
txi_tiss_atlas$abundance -> at_txi_tiss_atlas_tpm

# Reshape the data into a longer format
df_long <- at_txi_tiss_atlas_tpm %>%
  as.data.frame()%>%
  rownames_to_column("Gene")%>%
  pivot_longer(!Gene, names_to = "sample", values_to = "tpm")%>%
  dplyr::select(1,3,2)%>%
  left_join(.,
            at_tissue_metadata,
            by = c("sample" = "Run"))%>%
  dplyr::select(Gene,tpm,sample,Condition)# Extract tissue names

# Calculate the mean tpm for each tissue
mean_at_tiss_tpm <- df_long %>%
  group_by(Gene,Condition) %>%
  summarise(mean_tpm = mean(tpm)) %>%
  ungroup() %>%
  pivot_wider(names_from = Condition, values_from = mean_tpm)

##### Save avg TPM ####
write_csv(mean_at_tiss_tpm,"/mnt/Chikorita/Brassica_rapa/coexpression/AT/all_ens_lncRNA/avg_tpm_at_tiss_atlas_at_all_ens_lncRNA.csv")

#### LncRNA Expression Exploration ####
avg_tpm_at_tiss_atlas_at_all_ens_lncRNA <- read_csv("/mnt/Chikorita/Brassica_rapa/coexpression/AT/all_ens_lncRNA/avg_tpm_at_tiss_atlas_at_all_ens_lncRNA.csv")

#High confidence lncRNAs 
at_hc_lncRNA <- read_csv("/mnt/Chikorita/Brassica_rapa/core_files/at_hc_lncRNA.txt", 
                         col_names = FALSE) %>% setNames(c("name"))

Araport11_gene_type <- read_delim("/mnt/Chikorita/Brassica_rapa/coexpression/core_files/araport11_gene_type.txt", 
                                  delim = "\t", escape_double = FALSE, 
                                  col_names = TRUE, trim_ws = TRUE)

Araport11_gene_type %>%
  dplyr::select(1) -> all_araport_gene_ids

all_araport_gene_ids -> at_araport_retain

Araport11_gene_type %>%
  filter(gene_model_type %in% c("antisense_long_noncoding_rna", "long_noncoding_rna")) -> araport_lncRNA_meta

Araport11_gene_type %>%
  filter(gene_model_type %in% c("long_noncoding_rna")) -> araport_lincRNA_meta

#filter tx2gene for the genes of interest
tx2gene %>%
  inner_join(.,
             at_araport_retain,
             by = c("GENEID" = "name")) -> at_araport_retain_tx2gene

araport_lincRNA_meta %>%
  dplyr::select(1) -> at_araport_lincRNA_retain

tx2gene %>%
  inner_join(.,
             at_araport_lincRNA_retain,
             by = c("GENEID" = "name")) -> at_araport_lincRNA_retain_tx2gene

tx2gene %>%
  inner_join(.,
             at_hc_lncRNA,
             by = c("TXNAME" = "name")) -> at_hc_lncRNA_retain_tx2gene

rbind(at_araport_retain_tx2gene,at_hc_lncRNA_retain_tx2gene) -> tx2gene_retain
rbind(at_araport_lincRNA_retain_tx2gene,at_hc_lncRNA_retain_tx2gene) -> tx2gene_lincRNA_retain


write_csv(tx2gene_lincRNA_retain, "/mnt/Chikorita/Brassica_rapa/coexpression/core_files/tx2gene_lincRNA.csv")

tx2gene_lincRNA_retain %>%
  group_by(GENEID)

tx2gene_lincRNA_retain %>%
  dplyr::select(2)%>%
  inner_join(.,
             avg_tpm_at_tiss_atlas_at_all_ens_lncRNA,
             by = c("GENEID" = "Gene"))%>%
  distinct()-> avg_tpm_at_tiss_atlas_at_all_ens_lncRNA_lncRNA_only


write_csv(avg_tpm_at_tiss_atlas_at_all_ens_lncRNA_lncRNA_only,"/mnt/Chikorita/Brassica_rapa/coexpression/AT/all_ens_lncRNAavg_tpm_at_tiss_atlas_at_all_ens_lncRNA_lncRNA_only.csv")

read_csv("/mnt/Chikorita/Brassica_rapa/coexpression/AT/all_ens_lncRNA/avg_tpm_at_tiss_atlas_at_all_ens_lncRNA_lncRNA_only.csv") -> at_hc_lncRNA_tpm

#expressed?
avg_tpm_at_tiss_atlas_at_all_ens_lncRNA_lncRNA_only %>%
  column_to_rownames("GENEID") -> avg_tpm_at_tiss_atlas_at_all_ens_lncRNA_lncRNA_only_mat


#cool none of the genes are 0 tpm in all tissues
avg_tpm_at_tiss_atlas_at_all_ens_lncRNA_lncRNA_only_mat[rowSums(avg_tpm_at_tiss_atlas_at_all_ens_lncRNA_lncRNA_only_mat == 0) == ncol(avg_tpm_at_tiss_atlas_at_all_ens_lncRNA_lncRNA_only_mat), ] #6,565 lincRNAs are not expressed at all in this tissue atlas

#what about have 1 TPM in at least 1 or 2 tissues
avg_tpm_at_tiss_atlas_at_all_ens_lncRNA_lncRNA_only_mat[rowSums(avg_tpm_at_tiss_atlas_at_all_ens_lncRNA_lncRNA_only_mat >= 1) >= 2, ] #536 are expressed at 1 TPM in at least 2 tissues

avg_tpm_at_tiss_atlas_at_all_ens_lncRNA_lncRNA_only_mat[rowSums(avg_tpm_at_tiss_atlas_at_all_ens_lncRNA_lncRNA_only_mat >= 1) >= 1, ] #920 are expressed at 1 TPM in at least 1 tissue

avg_tpm_at_tiss_atlas_at_all_ens_lncRNA_lncRNA_only_mat[rowSums(avg_tpm_at_tiss_atlas_at_all_ens_lncRNA_lncRNA_only_mat >= .5) >= 1, ] #1,561 are expressed at .5 TPM in at least 1 tissues
