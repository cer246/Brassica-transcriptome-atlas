#### Brassica Expression Data Preparation ####

#### Load packages ####
library(tidyverse)
library(DESeq2)
library(txdbmaker)


#### Generate the raw counts dataframe ####
 
##### Read in featureCounts data #####
ftc.out <- read.table("/mnt/Chikorita/Brassica_rapa/featureCounts/br_all_tiss.ftc", header=TRUE, row.names=1)
ftc.out <- ftc.out[ ,6:ncol(ftc.out)]
colnames(ftc.out) <- gsub("\\.[sb]am$", "", colnames(ftc.out))
ftc.out.mat <- as.matrix(ftc.out)
counts <- ftc.out.mat

##### Read in sample metadata ####
br_tissue_metadata <- read_csv("/mnt/Chikorita/Brassica_rapa/coexpression/core_files/br_tissue_metadata_clean.csv") %>% select(Sample_ID,Tissue,Rep,Sample_name) %>% column_to_rownames("Sample_ID")

all(rownames(br_tissue_metadata) %in% colnames(ftc.out)) #True, that is good
all(rownames(br_tissue_metadata) == colnames(ftc.out)) #False this bad, we want the rownames and colnames to match

counts <- ftc.out[, rownames(br_tissue_metadata)]
all(rownames(br_tissue_metadata) == colnames(counts)) #Okay now our counts dataframe is in the correct order

br_tissue_metadata %>% 
  rownames_to_column("Sample_ID") -> br_tissue_metadata_id

colnames(counts) <- br_tissue_metadata_id$Sample_name[match(colnames(counts), br_tissue_metadata_id$Sample_ID)]

br_tissue_metadata_id %>%
  filter(!Sample_name %in% c("mature_seed_24_hour_imbibed_1","mature_seed_24_hour_imbibed_2" )) -> br_tissue_metadata_id_clean #we need to drop these

#you should save raw counts 

###### Save raw counts #####
write_tsv(counts, "/mnt/Chikorita/Brassica_rapa/coexpression/core_files/br_tissue_atlas_raw_counts.tsv")

#dropping one set of the mature_seed_imbibed_24hr
#counts %>%
#select(-mature_seed_24_hour_imbibed_1,-mature_seed_24_hour_imbibed_2) -> counts_clean

#change back to th sample_id
#colnames(counts_clean) <- br_tissue_metadata_id$Sample_ID[match(colnames(counts), br_tissue_metadata_id$Sample_name)]

counts_clean <- counts #we are using the clean metadata so we dont need to drop the mature seed 24hr second set


#### Generate the vsd counts dataframe ####
dds_br_tissue <- DESeqDataSetFromMatrix(countData = counts_clean,
                                        colData = br_tissue_metadata_id,
                                        design = ~Tissue)
#Filter for un-expressed genes
#keep_br_tissue <- rowSums(counts(dds_br_tissue) >= 10) >= 3 No filtering for now
#dds_br_tissue_filt <- dds_br_tissue[keep_br_tissue,]

#### Run DESeq ####
dds_tiss_atlas <- DESeq(dds_br_tissue) #101 genes do not converge with these settings, thats okay we can proceed

#### Transform and normalize counts data ####
vsd_br_tissue <- vst(dds_tiss_atlas, blind = F) 

dds_counts_tiss_atlas <- counts(dds_tiss_atlas, normalized=TRUE)
dds_counts_tiss_atlas %>% as.data.frame() -> dds_counts_tiss_atlas_df

######  Save normalized counts data ####
write.csv(dds_counts_tiss_atlas_df,"/mnt/Chikorita/Brassica_rapa/coexpression/core_files/br_tissue_atlas_norm_counts.csv")


vsd_counts_tiss_atlas <- data.frame(assay(vsd_br_tissue, normalized=TRUE)) #focus on this one first

#rename the embryonic seed samples
vsd_counts_tiss_atlas <- vsd_counts_tiss_atlas %>% 
  rename_at('X5_day_seed_1', ~'5_day_seed_1') %>% 
  rename_at('X5_day_seed_2', ~'5_day_seed_2')

#save the normalized deseq counts and vsd counts 

######  Save the transform and normalize counts data as a dataframe ####
dds_counts_tiss_atlas %>% 
  as.data.frame() %>%
  write_csv(., "/mnt/Chikorita/Brassica_rapa/coexpression/BR/dds_counts_tiss_atlas.csv")

######  Save the transformed and normalized counts data as a dataframe ####
vsd_counts_tiss_atlas %>%
  tibble::rownames_to_column(., var = "Gene") %>%
  as.data.frame() %>%
  write_csv(., "/mnt/Chikorita/Brassica_rapa/coexpression/BR/vsd_counts_tiss_atlas.csv")

##### Save the clean br_tissue_metadata ####
br_tissue_metadata_id_clean %>%
  as.data.frame() %>%
  write_csv(., "/mnt/Chikorita/Brassica_rapa/coexpression/core_files/br_tissue_metadata_clean.csv")


#### Generate raw tpm ####

##### Read in featureCounts data #####
ftc.out <- read.table("/mnt/Chikorita/Brassica_rapa/featureCounts/br_all_tiss.ftc", header=TRUE, row.names=1)
ftc.out <- ftc.out[ ,6:ncol(ftc.out)]
colnames(ftc.out) <- gsub("\\.[sb]am$", "", colnames(ftc.out))

ftc.out.mat <- as.matrix(ftc.out)

#Standard procedure for caluclating tpm
counts <- ftc.out.mat
featureLength <- ftc.out$Length


x <- counts / featureLength 
br_tiss_tpm.mat <- t( t(x) * 1e6 / colSums(x) )
br_tiss_tpm <- as.data.frame(br_tiss_tpm.mat) %>%rownames_to_column("Gene")%>%
  mutate(across('Gene', str_replace, 'gene:', ''))%>% select(-Br17b, -Br17c) 

##### Save the raw tpm ####
br_tiss_tpm %>% write.csv(., "/mnt/Chikorita/Brassica_rapa/coexpression/core_files/br_tiss_atlas_tpm.csv") 


#### Clean up the sample/column names ####

#vsd_counts_tiss_atlas <- read_csv("/mnt/Kanta/br_tissue/br_tissue_atlas/expression/co_expression/vsd_counts_tiss_atlas.csv")
#5-day seed looks weird - rename it
vsd_counts_tiss_atlas <- vsd_counts_tiss_atlas %>% 
  rename_at('X5_day_seed_1', ~'5_day_seed_1') %>% 
  rename_at('X5_day_seed_2', ~'5_day_seed_2')

#Pull in metadata
br_tissue_metadata <- read_csv("/mnt/Chikorita/Brassica_rapa/coexpression/core_files/br_tissue_metadata_clean.csv")

#define each tissue
br_tissue_metadata_id_clean %>%
  select(Tissue) %>%
  unique()%>%
  as.vector() -> tissues

#Reshape the data into a longer format
df_long <- vsd_counts_tiss_atlas %>%
  rownames_to_column("Gene")%>%
  pivot_longer(!Gene, names_to = "sample", values_to = "vsd_counts")%>%
  select(1,3,2)%>%
  left_join(.,
            br_tissue_metadata_id_clean,
            by = c("sample" = "Sample_name"))%>%
  select(Gene,vsd_counts,sample,Tissue)# Extract tissue names

#Calculate the mean vsd_counts for each tissue
mean_br_tiss_vsd_counts <- df_long %>%
  group_by(Gene,Tissue) %>%
  summarise(mean_vsd_counts = mean(vsd_counts)) %>%
  ungroup() %>%
  pivot_wider(names_from = Tissue, values_from = mean_vsd_counts)

#Alright how many genes have 0 tpm in in all tissues
mean_br_tiss_vsd_counts %>%
  column_to_rownames("Gene") -> mean_br_tiss_vsd_counts_test 

#Cool, none of the genes are 0 tpm in all tissues
mean_br_tiss_vsd_counts_test[rowSums(mean_br_tiss_vsd_counts_test == 0) == ncol(mean_br_tiss_vsd_counts_test), ]

#What about have 1 TPM in at least 1 or 2 tissues
mean_br_tiss_vsd_counts_test[rowSums(mean_br_tiss_vsd_counts_test >= 1) >= 2, ]

#Clean up the gene names and make the gene names rownames
mean_br_tiss_vsd_counts$Gene <- gsub("gene:", "", mean_br_tiss_vsd_counts$Gene)


##### Save averaged transformed and normalized counts ####
write_csv(mean_br_tiss_vsd_counts,"/mnt/Chikorita/Brassica_rapa/coexpression/BR/avg_vsd_clean_tiss_atlas.csv")
write_csv(mean_br_tiss_vsd_counts,"/mnt/Chikorita/Brassica_rapa/coexpression/core_files/avg_vsd_counts_br_tiss_atlas.csv")


#### Generate average tpm ####

#Read in raw tpm
br_tiss_tpm <- read_csv("/mnt/Chikorita/Brassica_rapa/coexpression/core_files/br_tiss_atlas_tpm.csv") %>% select(-`...1`)

#Pull in metadata
br_tissue_metadata <- read_csv("/mnt/Chikorita/Brassica_rapa/coexpression/core_files/br_tissue_metadata_clean.csv")

#colnames(br_tiss_tpm.mat) <- clean_br_metadata$Sample_name

#Define each tissue
br_tissue_metadata %>%
  dplyr::select(Tissue) %>%
  unique()%>%
  as.vector() -> tissues

#Reshape the data into a longer format
df_long <- br_tiss_tpm %>%
  pivot_longer(!Gene, names_to = "sample", values_to = "tpm")%>%
  dplyr::select(1,3,2)%>%
  #mutate(across('Gene', str_replace, 'gene:', ''))%>%
  left_join(.,
            br_tissue_metadata,
            by = c("sample" = "Sample_ID"))%>%
  dplyr::select(Gene,tpm,sample,Tissue)# Extract tissue names

#Calculate the mean tpm for each tissue
mean_br_tiss_tpm <- df_long %>%
  group_by(Gene,Tissue) %>%
  summarise(mean_tpm = mean(tpm)) %>%
  ungroup() %>%
  pivot_wider(names_from = Tissue, values_from = mean_tpm)

###### Save the averaged tpm #####
write_csv(mean_br_tiss_tpm,"/mnt/Chikorita/Brassica_rapa/coexpression/core_files/avg_br_tiss_atlas_tpm.csv")

#lets clean up the Gene names
#mean_at_tiss_vsd_counts %>%
#mutate(across('Gene', str_replace, 'gene:', '')) -> mean_at_tiss_vsd_counts

#write expressed 1 TPM in at least 2 tissues
#what about have 1 TPM in at least 1 or 2 tissues 

#### Generate expressed only tpm dataframe ####
#meaning at 1 TPM in at least 1 Tissue

#Read in averaged tpm
read_csv("/home/crailey/Brassica_rapa/coexpression/core_files/avg_br_tiss_atlas_tpm.csv") -> mean_br_tiss_tpm

mean_br_tiss_tpm %>% 
  column_to_rownames("Gene") -> mean_br_tiss_tpm

#mean_br_tiss_tpm[rowSums(mean_br_tiss_tpm >= 1) >= 2, ] %>% rownames_to_column("Gene") -> mean_br_tiss_tpm_exp this is 1 TPM in at least 2 Tissues
mean_br_tiss_tpm[rowSums(mean_br_tiss_tpm >= 1) >= 1, ] %>% rownames_to_column("Gene") -> mean_br_tiss_tpm_exp
mean_br_tiss_tpm[rowSums(mean_br_tiss_tpm[ , -1] >= 1) > 0, ] -> mean_br_tiss_tpm_exp

##### Save expressed only TPM #####
mean_br_tiss_tpm_exp %>% write.csv(., "/mnt/Chikorita/Brassica_rapa/coexpression/core_files/avg_br_tiss_atlas_tpm_only_expressed.csv")

# how many of the br_tiss_genes are PCGs and how many meet the expression threshold
mean_br_tiss_tpm %>%
  filter(Gene %in% br_gene_type_pcg$GENEID)%>% 
  column_to_rownames("Gene") -> br_anno_pcg_mean_br_tiss_tpm 

br_anno_pcg_mean_br_tiss_tpm[rowSums(br_anno_pcg_mean_br_tiss_tpm >= 1) >= 1, ] %>% rownames_to_column("Gene") #1 TPM in at least 1 Tissue
br_anno_pcg_mean_br_tiss_tpm[rowSums(br_anno_pcg_mean_br_tiss_tpm >= 1) > 0, ] 

br_anno_pcg_mean_br_tiss_tpm[rowSums(br_anno_pcg_mean_br_tiss_tpm >= 1) == 0, ]#less than 1 TPM in all tissues

br_anno_pcg_mean_br_tiss_tpm[rowSums(br_anno_pcg_mean_br_tiss_tpm != 0) == 0, ] %>%
  rownames_to_column("Gene") -> br_illumnina_not_detected

br_nano_tiss_norm_counts %>%
  as.data.frame() -> br_nano_tiss_norm_counts_df

br_nano_tiss_norm_counts_df[rowSums(br_nano_tiss_norm_counts_df != 0) == 0, ] %>%
  rownames_to_column("Gene") -> br_nano_not_detected

br_nano_not_detected %>%
  mutate(Gene = str_replace_all(.$Gene, "gene:", "")) %>%
  filter(Gene %in% br_illumnina_not_detected$Gene) -> br_not_detected_in_nano_illumina$Gene


br_not_detected_in_nano_illumina %>%
  dplyr::select(Gene) %>%
  distinct() -> br_not_detected_in_nano_illumina_gene_id
##### We will use these not detected pcgs later ####

#### Subset and save the expressed genes in the transformed and normalized counts data ####
mean_br_tiss_tpm_exp %>%
  select("Gene")%>%
  left_join(.,
            mean_br_tiss_vsd_counts,
            by = c("Gene")) %>%
  write.csv(., "/mnt/Kanta/br_tissue/core_files/avg_vsd_counts_br_tiss_atlas_only_expressed.csv")

#### Generate average raw counts ####

#read in raw counts
ftc.out <- read.table("/mnt/Chikorita/Brassica_rapa/featureCounts/br_all_tiss.ftc", header=TRUE, row.names=1)
ftc.out <- ftc.out[ ,6:ncol(ftc.out)]
colnames(ftc.out) <- gsub("\\.[sb]am$", "", colnames(ftc.out))

ftc.out.mat <- as.matrix(ftc.out)

counts <- ftc.out.mat

#pull in the metadata
br_tissue_metadata <- read_csv("/mnt/Kanta/br_tissue/br_tissue_atlas/metadata/br_tissue_metadata.csv") %>% select(Sample_ID,Tissue,Rep,Sample_name) %>% column_to_rownames("Sample_ID")

all(rownames(br_tissue_metadata) %in% colnames(ftc.out)) #True, that is good
all(rownames(br_tissue_metadata) == colnames(ftc.out)) #False this bad, we want the rownames and colnames to match

counts <- ftc.out[, rownames(br_tissue_metadata)]
all(rownames(br_tissue_metadata) == colnames(counts)) #Okay now our counts dataframe is in the correct order

br_tissue_metadata %>% 
  rownames_to_column("Sample_ID") -> br_tissue_metadata_id

colnames(counts) <- br_tissue_metadata_id$Sample_name[match(colnames(counts), br_tissue_metadata_id$Sample_ID)]
br_tissue_metadata_id %>%
  filter(!Sample_name %in% c("mature_seed_24_hour_imbibed_1","mature_seed_24_hour_imbibed_2" )) -> br_tissue_metadata_id_clean

#dropping one set of the mature_seed_imbibed_24hr
counts %>%
  select(-mature_seed_24_hour_imbibed_1,-mature_seed_24_hour_imbibed_2) -> counts_clean

counts_clean

#define each tissue
br_tissue_metadata_id_clean %>%
  select(Tissue) %>%
  unique()%>%
  as.vector() -> tissues

#Reshape the data into a longer format
df_long <- counts_clean %>%
  rownames_to_column("Gene")%>%
  pivot_longer(!Gene, names_to = "sample", values_to = "raw_counts")%>%
  select(1,3,2)%>%
  left_join(.,
            br_tissue_metadata_id_clean,
            by = c("sample" = "Sample_name"))%>%
  select(Gene,raw_counts,sample,Tissue)# Extract tissue names

#Calculate the mean vsd_counts for each tissue
mean_br_tiss_raw_counts <- df_long %>%
  group_by(Gene,Tissue) %>%
  summarise(mean_raw_counts = mean(raw_counts)) %>%
  ungroup() %>%
  pivot_wider(names_from = Tissue, values_from = mean_raw_counts)

##### Save raw counts with updated gene names ####
counts_clean %>%
  rownames_to_column("Gene")%>%
  mutate(across('Gene', str_replace, 'gene:', ''))  %>%
  write.csv(., "/mnt/Kanta/br_tissue/core_files/br_tiss_atlas_raw_counts.csv")

##### Save averaged counts with updated gene names ####
mean_br_tiss_raw_counts %>%
  mutate(across('Gene', str_replace, 'gene:', '')) %>%
  write.csv(., "/mnt/Kanta/br_tissue/core_files/avg_br_tiss_atlas_raw_counts.csv")

#### Long-read data ####


br_nano_metadata <- read_csv("/mnt/Chikorita/Brassica_rapa/coexpression/core_files/br_nano_metadata.csv")

# set working directory to where the quant folders are
nano_quant_dir <- "/home/crailey/Brassica_rapa/br_longread/oarfish/"

# set the file names
files <- paste(nano_quant_dir,br_nano_metadata$ShortID, ".quant", sep = "")

library(tximport)

##### Read in oarfish files ######
# set txOut to F if you want isoform level counts
txi <- tximport(files, type = "oarfish", tx2gene = tx2gene_br, txOut = F)

ddsTxi_br_nano_tiss <- DESeqDataSetFromTximport(txi,
                                                colData = br_nano_metadata,
                                                design = ~Tissue)

dds_br_nano_tiss <- DESeq(ddsTxi_br_nano_tiss)

### Get normalized counts ###
dds <- estimateSizeFactors(dds_br_nano_tiss)
br_nano_tiss_norm_counts <- counts(dds, normalized=TRUE)

#### Generate PCA ####
vsd <- vst(dds_br_nano_tiss, blind = F)
plotPCA <- function (object, intgroup=c("~Tissue"), ntop = 3000, returnData=TRUE) 
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

#### Plot PCA ####
pcaData <- plotPCA(vsd, intgroup=c("Tissue"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
percentVar

nb.cols <- 28
mycolors <- colorRampPalette(brewer.pal(8, "BrBG"))(nb.cols)

DESeq2::plotPCA(vsd, intgroup = c("Tissue")) + theme_linedraw()+ scale_color_manual(values = tissue_colors)+ 
  theme(axis.text.x = (element_text(size = 14, color = "black")),
        axis.text.y = (element_text(size = 14, color = "black")),
        axis.title.x = (element_text(size = 16, color = "black")),
        axis.title.y = (element_text(size = 16, color = "black"))) -> br_nano_pca

#### Plot Hierarchial clustering ####
vsd <- vst(dds_br_nano_tiss, blind = F)
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$ID)
colnames(sampleDistMatrix) <- paste(vsd$Tissue)
colors <- colorRampPalette( rev(brewer.pal(9, "YlGnBu")) )(100)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col=colors,
         border_color = "gray",
         show_rownames = TRUE, 
         show_colnames = TRUE,
         angle_col = 45) -> br_nano_distance_matrix