##### Gene Type Assignemnt #####

###### Load packages ####
library(tidyverse)

#This should be the directory in which the gffcompare output files are stored
setwd("/mnt/Kanta/br_tissue/br_tissue_atlas/lncRNA/st_gtf/gffcompare_br_tiss")

#Read in your gffcompare tmap file
br_gffcomp <- read_delim("/mnt/Chikorita/Brassica_rapa/lncRNA/st_gtf/gffcompare_br_tiss/gffcompare_br_tiss.merge_br_tissue_atlas.gtf.tmap", 
                         delim = "\t", escape_double = FALSE, 
                         trim_ws = TRUE)

#### Intergenic lncRNAs ####
#First lets filter for intergenic lncRNAs - this is class code u
br_intergenic <- br_gffcomp %>% 
  dplyr::filter(class_code == "u" & len >= 200) # filter for intergenic lncRNAs

br_intergenic %>% 
  distinct(qry_id) %>% 
  write.table(., "/mnt/Chikorita/Brassica_rapa/lncRNA/st_gtf//gffcompare_br_tiss/br_intergenic_gffcompare_txIDs.txt", 
              row.names = F, col.names = F, quote = F) # get the transcript IDs to use so you can grep later
# these txIDs may have carriage returns and grepping won't work if they do! 
# make sure to check for them

#### Antisense lncRNAs ####
#Now lets filter for antisense RNAs and save those transcripts
br_antisense <- br_gffcomp %>% 
  dplyr::filter(class_code == "x" & len >= 200) # filter for antisense lncRNAs

br_antisense %>% 
  distinct(qry_id) %>% 
  write.table(., "/mnt/Chikorita/Brassica_rapa/lncRNA/st_gtf/gffcompare_br_tiss/br_antisense_gffcompare_txIDs.txt", 
              row.names = F, col.names = F, quote = F) # get the transcript IDs to use so you can grep later
# these txIDs may have carriage returns and grepping won't work if they do! 
# make sure to check for them


###### Intronic ####
#Now lets filter for intronic RNAs
br_intronic <- br_gffcomp %>%
  dplyr::filter(class_code == "i" & num_exons > 1 & len >= 200) # get intronic RNAs
# these intronic RNAs will be on both strands relative to the reference
#To further filter the intronic lncRNAs that are antisense intronic lncRNAs we need to pull in some more information i.e. our up-to-date genome annotation

br_merged_anno <- read.table("/mnt/Chikorita/Brassica_rapa/lncRNA/st_gtf/merge_br_tissue_atlas.gtf", sep = "\t", header = F)
# load in the stringtie merge annotation so that you can get the strand of the intronic RNA

#Clean up the annotation to make it easier to isolate transcript IDs
br_merged_anno_clean <- br_merged_anno %>%
  separate(col = V9, into = c("gene", "tx", "exon"),
           sep = ";") %>%
  dplyr::select(V7, tx) %>%
  distinct() # separate the last column of the stringtie merge annotation 
# to isolate the txIDs

# clean up the txID
br_merged_anno_clean$tx <- gsub(' transcript_id ', '', br_merged_anno_clean$tx)


#Now add a column that denotes the strand of the intronic lncRNA
intronic_with_strand <- br_intronic %>%
  left_join(., br_merged_anno_clean,
            by = c("qry_id" = "tx")) %>%
  distinct() %>%
  dplyr::rename(intronic_strand = V7) # join the strand to the intronic txIDs

#Now we need to figure out what strand the pcg is on
# load in the reference GTF 
br_ref_anno_mod <- read.table("/mnt/Chikorita/Brassica_rapa/coexpression/core_files/Brassica_rapa_ro18.SCU_BraROA_2.3.59.gtf", sep = "\t", header = F)


#lets clean up the reference annotation and extract the information that we really need
br_ref_anno_mod_clean <- br_ref_anno_mod %>%
  dplyr::filter(V3 == "transcript") %>%
  separate(col = V9, into = c("gene", "tx", "exon", "other", "other2", "other3"),
           sep = ";", extra = "merge") %>%
  dplyr::select(V7, tx) # separate

#clean up those transcript IDs
br_ref_anno_mod_clean$tx <- gsub('gene_id ', '', br_ref_anno_mod_clean$tx )
br_ref_anno_mod_clean$tx <- gsub(' ', '', br_ref_anno_mod_clean$tx )
#br_ref_anno_mod_clean$gene <- gsub('ID=', '', br_ref_anno_mod_clean$gene )

# Now, join together reference transcript strand with assembled transcript
intronic_with_strand %>%
  select(1,3,4,5,13) %>%
  as.data.frame()-> intronic_with_strand_alt

intronic_with_strand_alt %>%
  left_join(.,
            br_ref_anno_mod_clean, 
            by = c("ref_gene_id" = "tx"),
            copy = TRUE)%>%
  dplyr::rename(reference_strand = V7) %>%
  distinct() -> intron2

# filter for antisense intronic transcripts
intron3 <- intron2 %>%
  dplyr::filter(intronic_strand != reference_strand)

#Export those IDs, again watching out for carriage returns
setwd("/mnt/Chikorita/Brassica_rapa/lncRNA/st_gtf/gffcompare_br_tiss")
intron3 %>%
  distinct(qry_id) %>%
  write.table("br_antisense_intron_gffcompare_txIDs.txt",
              quote = F, row.names = F, col.names = F)

# Now we have already run Pfam and Rfam on the intergenic and antisense transcripts...
# And we have already run EDTA to predict TEs genome-wide and
# We ran bedtools intersect to see if any of our transcripts overlap with these other transcript types


#### Integrate  accessory information to assign transcript types ####
#Read in all the necessary output files
#read in intergenic gffcompare txids
br_intergenic_gffcompare_txIDs <- read_csv("/mnt/Chikorita/Brassica_rapa/lncRNA/st_gtf/gffcompare_br_tiss/br_intergenic_gffcompare_txIDs.txt", 
                                           col_names = FALSE) %>% setNames(c("tx_id"))

br_antisense_gffcompare_txIDs <- read_csv("/mnt/Chikorita/Brassica_rapa/lncRNA/st_gtf/gffcompare_br_tiss/br_all_antisense_gffcompare_txIDs.txt", 
                                          col_names = FALSE) %>% setNames(c("tx_id"))

#read in the the pfam hits
pfam_br_intergenic_hits <- read_delim("/mnt/Chikorita/Brassica_rapa/lncRNA/st_gtf/gffcompare_br_tiss/br_intergenic_gffcompare_txIDs.fa.transdecoder_dir/pfam_br_intergenic_hits.txt", 
                                      delim = "\t", escape_double = FALSE, 
                                      trim_ws = TRUE) %>% select(1,3) #lets just keep the transcript column

pfam_br_antisense_hits <- read_delim("/mnt/Chikorita/Brassica_rapa/lncRNA/st_gtf/gffcompare_br_tiss/br_all_antisense_gffcompare_txIDs.fa.transdecoder_dir/pfam_br_antisense_hits.txt", 
                                     delim = "\t", escape_double = FALSE, 
                                     trim_ws = TRUE) %>% select(1,3) #lets just keep the transcript column

#read in rfam hits
rfam_br_antisense <- read_table("//mnt/Chikorita/Brassica_rapa/lncRNA/st_gtf/gffcompare_br_tiss/kyle_test/br_all_antisense_rfam_intersect_all_output.txt", 
                                col_names = FALSE)

rfam_br_intergenic <- read_table("/mnt/Chikorita/Brassica_rapa/lncRNA/st_gtf/gffcompare_br_tiss/kyle_test/br_intergenic_rfam_all_overlap.txt", 
                                 col_names = FALSE)

#cpc2 output
cpc2_br_intergenic <- read_delim("/mnt/Chikorita/Brassica_rapa/lncRNA/st_gtf/gffcompare_br_tiss/CPC2_output/br_intergenic_gffcompare_txIDs_CPC2output.txt.txt", 
                                 delim = "\t", escape_double = FALSE, 
                                 trim_ws = TRUE)

cpc2_br_antisense <- br_all_antisense_gffcompare_txIDs_CPC2output_txt <- read_delim("/mnt/Chikorita/Brassica_rapa/lncRNA/st_gtf/gffcompare_br_tiss/CPC2_output/br_all_antisense_gffcompare_txIDs_CPC2output.txt.txt", 
                                                                                    delim = "\t", escape_double = FALSE, 
                                                                                    trim_ws = TRUE)

#TE overlap
br_intergenic_helitron_all_overlap <- read_delim("/mnt/Chikorita/Brassica_rapa/TE/Brassica_rapa_ro18.SCU_BraROA_2.3.dna.toplevel.fa.mod.EDTA.raw/br_intergenic_helitron_all_overlap.txt", 
                                                 delim = "\t", escape_double = FALSE, 
                                                 col_names = FALSE, trim_ws = TRUE)

br_antisense_helitron_all_overlap <- read_delim("/mnt/Chikorita/Brassica_rapa/TE/Brassica_rapa_ro18.SCU_BraROA_2.3.dna.toplevel.fa.mod.EDTA.raw/br_antisense_helitron_all_overlap.txt", 
                                                delim = "\t", escape_double = FALSE, 
                                                col_names = FALSE, trim_ws = TRUE)

br_intergenic_ltr_all_overlap <- read_delim("/mnt/Chikorita/Brassica_rapa/TE/Brassica_rapa_ro18.SCU_BraROA_2.3.dna.toplevel.fa.mod.EDTA.raw/br_intergenic_ltr_all_overlap.txt", 
                                            delim = "\t", escape_double = FALSE, 
                                            col_names = FALSE, trim_ws = TRUE)

br_antisense_ltr_all_overlap <- read_delim("/mnt/Chikorita/Brassica_rapa/TE/Brassica_rapa_ro18.SCU_BraROA_2.3.dna.toplevel.fa.mod.EDTA.raw/br_antisense_ltr_all_overlap.txt", 
                                           delim = "\t", escape_double = FALSE, 
                                           col_names = FALSE, trim_ws = TRUE)


###### Process rfam files ####
rfam_br_antisense %>%
  filter(X10 == ".") %>%
  select(X4) %>%
  unique() %>%
  mutate(Rfam_status = "no_rfam_hit") -> rfam_br_antisese_no_hits

rfam_br_antisense %>%
  filter(X10 != ".") %>%
  select(X4,X10)%>%
  group_by(X4) %>%
  summarise(Rfam_status = paste(X10, collapse = ",")) -> rfam_br_antisense_hits

rfam_br_intergenic %>%
  filter(X10 == ".") %>%
  select(X4) %>%
  unique()%>%
  mutate(Rfam_status = "no_rfam_hit") -> rfam_br_intergenic_no_hits

rfam_br_intergenic %>%
  filter(X10 != ".") %>%
  select(X4,X10)%>%
  group_by(X4) %>%
  summarise(Rfam_status = paste(X10, collapse = ",")) -> rfam_br_intergenic_hits  

rbind(rfam_br_intergenic_hits,rfam_br_intergenic_no_hits) -> rfam_br_intergenic_all
rbind(rfam_br_antisense_hits,rfam_br_antisese_no_hits) -> rfam_br_antisense_all


###### Process TE files ####
#Let's add in the final layer - TE associated transcripts
br_intergenic_helitron_all_overlap %>%
  filter(X10 != ".") %>% #183 transcripts overlap
  select(X4,X10) -> br_intergenic_helitron_tx
br_intergenic_ltr_all_overlap %>%
  filter(X10 != ".") %>% #150 transcripts overlap
  select(X4,X10) -> br_intergenic_ltr_tx

rbind(br_intergenic_helitron_tx, br_intergenic_ltr_tx) -> br_int_te_tx

br_antisense_helitron_all_overlap %>%
  filter(X10 != ".") %>% #51 transcripts overlap
  select(X4,X10) -> br_antisense_helitron_tx
br_antisense_ltr_all_overlap%>%
  filter(X10 != ".") %>% #36 transcripts overlap
  select(X4,X10) -> br_antisense_ltr_tx

rbind(br_antisense_helitron_tx,br_antisense_ltr_tx) -> br_anti_te_tx


##### Integrate CPC2 files ####
br_intergenic_gffcompare_txIDs %>%
  left_join(.,
            cpc2_br_intergenic, 
            by = c("tx_id" = "#ID")) %>%
  select(tx_id, label) %>%
  left_join(.,
            pfam_br_intergenic_hits,
            by = c("tx_id" = "tx_id_pfam_hit"))%>% 
  replace_na(list(hmm_name_type = 'no_pfam_hit')) %>%
  left_join(.,
            rfam_br_intergenic_all,
            by = c("tx_id" = "X4")) %>%
  setNames(c("Transcript_ID", "CPC2_label","Pfam_status","Rfam_status")) %>%
  mutate(lncRNA_class = "Intergenic") -> br_intergenic_meta


br_antisense_gffcompare_txIDs %>%
  left_join(.,
            cpc2_br_antisense, 
            by = c("tx_id" = "#ID")) %>%
  select(tx_id, label) %>%
  left_join(.,
            pfam_br_antisense_hits,
            by = c("tx_id" = "tx_id_pfam_hit"))%>% 
  replace_na(list(hmm_name_type = 'no_pfam_hit')) %>%
  left_join(.,
            rfam_br_antisense_all,
            by = c("tx_id" = "X4")) %>% 
  replace_na(list(Rfam = 'rfam_hit')) %>%
  setNames(c("Transcript_ID", "CPC2_label","Pfam_status","Rfam_status")) %>%
  mutate(lncRNA_class = "Antisense") -> br_antisense_meta


##### Integrate TE information #####
#Update the meta information with transposable element overlap info
br_intergenic_meta %>%
  left_join(.,
            br_int_te_tx,
            by = c("Transcript_ID" = "X4")) %>%
  replace_na(list(X10 = 'no_te_overlap')) %>%
  rename("TE_status" = "X10") -> br_intergenic_meta

br_antisense_meta %>%
  left_join(.,
            br_anti_te_tx,
            by = c("Transcript_ID" = "X4")) %>%
  replace_na(list(X10 = 'no_te_overlap')) %>%
  rename("TE_status" = "X10") -> br_antisense_meta

#save these files
write.table(br_intergenic_meta,"/mnt/Chikorita/Brassica_rapa/lncRNA/st_gtf/gffcompare_br_tiss/br_intergenic_meta.txt")

write.table(br_antisense_meta,"/mnt/Chikorita/Brassica_rapa/lncRNA/st_gtf/gffcompare_br_tiss/br_antisense_meta.txt")


#Do the antisense and intergenic independently 
#Let's see how many we have

#CPC2 non-coding, no_pfam_hit, no_rfam_hit, no_te_overlap
#CPC2_label == "noncoding" & Pfam_status == "no_pfam_hit" & Rfam_status == "no_rfam_hit" & TE_status == "no_te_overlap"  ~ "True lncRNA",

br_intergenic_meta %>%
  filter(CPC2_label == "noncoding" & Pfam_status == "no_pfam_hit" & Rfam_status == "no_rfam_hit" & TE_status == "no_te_overlap") %>% unique() #6078

br_antisense_meta %>%
  filter(CPC2_label == "noncoding" & Pfam_status == "no_pfam_hit" & Rfam_status == "no_rfam_hit" & TE_status == "no_te_overlap") %>% unique() #3452



#CPC2 coding, no_pfam_hit, no_rfam_hit
#CPC2_label == "coding" & Pfam_status == "no_pfam_hit" & Rfam_status == "no_rfam_hit" & TE_status == "no_te_overlap" ~ "Putative novel PCG"
br_intergenic_meta %>%
  filter(CPC2_label == "coding" & Pfam_status == "no_pfam_hit" & Rfam_status == "no_rfam_hit" & TE_status == "no_te_overlap") %>% unique() #869

br_antisense_meta %>%
  filter(CPC2_label == "coding" & Pfam_status == "no_pfam_hit" & Rfam_status == "no_rfam_hit" & TE_status == "no_te_overlap") %>% unique() #302



#CPC2_label == "coding" & Pfam_status == "no_pfam_hit" & Rfam_status == "no_rfam_hit" & TE_status != "no_te_overlap" ~ "Polycistronic transcript"
br_intergenic_meta %>%
  filter(CPC2_label == "coding" & Pfam_status == "no_pfam_hit" & Rfam_status == "no_rfam_hit" & TE_status != "no_te_overlap") %>% unique() #25

br_antisense_meta %>%
  filter(CPC2_label == "coding" & Pfam_status == "no_pfam_hit" & Rfam_status == "no_rfam_hit" & TE_status != "no_te_overlap") %>% unique() #7



#CPC2 noncoding, pfam_hit, no_rfam_hit, no te overlap
#CPC2_label == "noncoding" & Pfam_status != "no_pfam_hit" & Rfam_status == "no_rfam_hit" & TE_status == "no_te_overlap" ~ "Polycistronic transcript"
br_intergenic_meta %>%
  filter(CPC2_label == "noncoding" & Pfam_status != "no_pfam_hit" & Rfam_status == "no_rfam_hit" & TE_status == "no_te_overlap") %>% unique() #594

br_antisense_meta %>%
  filter(CPC2_label == "noncoding" & Pfam_status != "no_pfam_hit" & Rfam_status == "no_rfam_hit" & TE_status == "no_te_overlap") %>% unique() #43



#CPC2 noncoding, pfam_hit, no_rfam_hit, TE overlap
#CPC2_label == "noncoding" & Pfam_status != "no_pfam_hit" & Rfam_status == "no_rfam_hit" & TE_status == "no_te_overlap" ~ "Polycistronic transcript"
br_intergenic_meta %>%
  filter(CPC2_label == "noncoding" & Pfam_status != "no_pfam_hit" & Rfam_status == "no_rfam_hit" & TE_status != "no_te_overlap") %>% unique() #10

br_antisense_meta %>%
  filter(CPC2_label == "noncoding" & Pfam_status != "no_pfam_hit" & Rfam_status == "no_rfam_hit" & TE_status != "no_te_overlap") %>% unique() #3


#CPC2 coding, pfam_hit, no_rfam_hit, no te overlap
#CPC2_label == "coding" & Pfam_status != "no_pfam_hit" & Rfam_status == "no_rfam_hit" & TE_status == "no_te_overlap" ~ "Polycistronic transcript"
br_intergenic_meta %>%
  filter(CPC2_label == "coding" & Pfam_status != "no_pfam_hit" & Rfam_status == "no_rfam_hit" & TE_status == "no_te_overlap") %>% unique() #3182

br_antisense_meta %>%
  filter(CPC2_label == "coding" & Pfam_status != "no_pfam_hit" & Rfam_status == "no_rfam_hit" & TE_status == "no_te_overlap") %>% unique() #396



#CPC2 coding, pfam_hit, no_rfam_hit,te overlap
#CPC2_label == "coding" & Pfam_status != "no_pfam_hit" & Rfam_status == "no_rfam_hit" & TE_status == "no_te_overlap" ~ "Polycistronic transcript"
br_intergenic_meta %>%
  filter(CPC2_label == "coding" & Pfam_status != "no_pfam_hit" & Rfam_status == "no_rfam_hit" & TE_status != "no_te_overlap") %>% unique() #78

br_antisense_meta %>%
  filter(CPC2_label == "coding" & Pfam_status != "no_pfam_hit" & Rfam_status == "no_rfam_hit" & TE_status != "no_te_overlap") %>% unique() #28



#CPC2 non-coding, no_pfam_hit, rfam_hit
#CPC2_label == "noncoding" & Pfam_status == "no_pfam_hit" & Rfam_status != "no_rfam_hit" & TE_status == "no_te_overlap" ~ "Housekeeping ncRNA precursor"
br_intergenic_meta %>%
  filter(CPC2_label == "noncoding" & Pfam_status == "no_pfam_hit" & Rfam_status != "no_rfam_hit" & TE_status == "no_te_overlap") %>% unique() #1077

br_antisense_meta %>%
  filter(CPC2_label == "noncoding" & Pfam_status == "no_pfam_hit" & Rfam_status != "no_rfam_hit" & TE_status == "no_te_overlap") %>% unique()#77



#CPC2 non-coding, no_pfam_hit, rfam_hit, te_overlap
#CPC2_label == "noncoding" & Pfam_status == "no_pfam_hit" & Rfam_status != "no_rfam_hit" & TE_status == "no_te_overlap" ~ "Polycistronic transcript"
br_intergenic_meta %>%
  filter(CPC2_label == "noncoding" & Pfam_status == "no_pfam_hit" & Rfam_status != "no_rfam_hit" & TE_status != "no_te_overlap") %>% unique() #25

br_antisense_meta %>%
  filter(CPC2_label == "noncoding" & Pfam_status == "no_pfam_hit" & Rfam_status != "no_rfam_hit" & TE_status != "no_te_overlap") %>% unique()#7



#CPC2 coding, no_pfam_hit, rfam_hit, no_te_overlap 
#CPC2_label == "coding" & Pfam_status == "no_pfam_hit" & Rfam_status != "no_rfam_hit" & TE_status == "no_te_overlap" ~ Polycistronic transcript
br_intergenic_meta %>%
  filter(CPC2_label == "coding" & Pfam_status == "no_pfam_hit" & Rfam_status != "no_rfam_hit"& TE_status == "no_te_overlap") %>% unique() #1

br_antisense_meta %>%
  filter(CPC2_label == "coding" & Pfam_status == "no_pfam_hit" & Rfam_status != "no_rfam_hit"& TE_status == "no_te_overlap") %>% unique() #0



#CPC2 coding, no_pfam_hit, rfam_hit, te_overlap 
#CPC2_label == "coding" & Pfam_status == "no_pfam_hit" & Rfam_status != "no_rfam_hit" & TE_status == "no_te_overlap" ~ Polycistronic transcript
br_intergenic_meta %>%
  filter(CPC2_label == "coding" & Pfam_status == "no_pfam_hit" & Rfam_status != "no_rfam_hit"& TE_status != "no_te_overlap") %>% unique() #0

br_antisense_meta %>%
  filter(CPC2_label == "coding" & Pfam_status == "no_pfam_hit" & Rfam_status != "no_rfam_hit"& TE_status != "no_te_overlap") %>% unique() #0



#CPC2 coding, pfam_hit, rfam_hit, no_te_overlap 
#CPC2_label == "coding" & Pfam_status != "no_pfam_hit" & Rfam_status != "no_rfam_hit" & TE_status == "no_te_overlap" ~ Polycistronic transcript
br_intergenic_meta %>%
  filter(CPC2_label == "coding" & Pfam_status != "no_pfam_hit" & Rfam_status != "no_rfam_hit"& TE_status == "no_te_overlap") %>% unique() #13

br_antisense_meta %>%
  filter(CPC2_label == "coding" & Pfam_status != "no_pfam_hit" & Rfam_status != "no_rfam_hit"& TE_status == "no_te_overlap") %>% unique() #0



#CPC2 coding, pfam_hit, rfam_hit, te_overlap 
#CPC2_label == "coding" & Pfam_status != "no_pfam_hit" & Rfam_status != "no_rfam_hit" & TE_status == "no_te_overlap" ~ Polycistronic transcript
br_intergenic_meta %>%
  filter(CPC2_label == "coding" & Pfam_status != "no_pfam_hit" & Rfam_status != "no_rfam_hit"& TE_status != "no_te_overlap") %>% unique() #7

br_antisense_meta %>%
  filter(CPC2_label == "coding" & Pfam_status != "no_pfam_hit" & Rfam_status != "no_rfam_hit"& TE_status != "no_te_overlap") %>% unique() #0


#CPC2 noncoding, no pfam hit, no rfam hit, te_overlap
#CPC2_label == "noncoding" & Pfam_status == "no_pfam_hit" & Rfam_status == "no_rfam_hit" & TE_status != "no_te_overlap" ~ "TE-associated lncRNA"
br_intergenic_meta %>%
  filter(CPC2_label == "noncoding" & Pfam_status == "no_pfam_hit" & Rfam_status == "no_rfam_hit" & TE_status != "no_te_overlap")%>% unique() #187

br_antisense_meta %>%
  filter(CPC2_label == "noncoding" & Pfam_status == "no_pfam_hit" & Rfam_status == "no_rfam_hit" & TE_status != "no_te_overlap")%>% unique() #61




#CPC2 noncoding, pfam_hit, rfam_hit, no_te_overlap
#CPC2_label == "noncoding" & Pfam_status != "no_pfam_hit" & Rfam_status != "no_rfam_hit" & TE_status == "no_te_overlap" ~ Polycistronic transcript

br_intergenic_meta %>%
  filter(CPC2_label == "noncoding" & Pfam_status != "no_pfam_hit" & Rfam_status != "no_rfam_hit" & TE_status == "no_te_overlap") %>% unique() #2

br_antisense_meta %>%
  filter(CPC2_label == "noncoding" & Pfam_status != "no_pfam_hit" & Rfam_status != "no_rfam_hit" & TE_status == "no_te_overlap") %>% unique() #9

br_antisense_meta %>% select(1) %>% unique()



#CPC2 noncoding, pfam_hit, rfam_hit, te_overlap ~ Polycistronic transcript
#CPC2_label == "noncoding" & Pfam_status != "no_pfam_hit" & Rfam_status != "no_rfam_hit" & TE_status == "no_te_overlap"
br_intergenic_meta %>%
  filter(CPC2_label == "noncoding" & Pfam_status != "no_pfam_hit" & Rfam_status != "no_rfam_hit" & TE_status != "no_te_overlap") %>% unique() #1

br_antisense_meta %>%
  filter(CPC2_label == "noncoding" & Pfam_status != "no_pfam_hit" & Rfam_status != "no_rfam_hit" & TE_status != "no_te_overlap") %>% unique() #0



#Add lncRNA transcript category column
#Update on how exactly we are categorizing the identified transcripts

##### Assign gene type classifications ####

#USE THIS ONE
br_intergenic_meta %>%
  mutate(
    Type = case_when(
      CPC2_label == "noncoding" & Pfam_status == "no_pfam_hit" & Rfam_status == "no_rfam_hit" & TE_status == "no_te_overlap"  ~ "True lncRNA",
      CPC2_label == "noncoding" & Pfam_status == "no_pfam_hit" & Rfam_status == "no_rfam_hit" & TE_status != "no_te_overlap" ~ "TE-associated lncRNA",
      CPC2_label == "coding" & Pfam_status == "no_pfam_hit" & Rfam_status == "no_rfam_hit" & TE_status == "no_te_overlap" ~ "Putative novel PCG",
      CPC2_label == "noncoding" & Pfam_status != "no_pfam_hit" & Rfam_status == "no_rfam_hit" & TE_status == "no_te_overlap" ~ "Polycistronic transcript",
      CPC2_label == "coding" & Pfam_status != "no_pfam_hit" & Rfam_status == "no_rfam_hit" & TE_status == "no_te_overlap" ~ "Putative novel PCG",
      CPC2_label == "coding" & Pfam_status != "no_pfam_hit" & Rfam_status != "no_rfam_hit" & TE_status == "no_te_overlap" ~ "Polycistronic transcript",
      CPC2_label == "coding" & Pfam_status == "no_pfam_hit" & Rfam_status != "no_rfam_hit" & TE_status == "no_te_overlap" ~ "Polycistronic transcript",
      CPC2_label == "noncoding" & Pfam_status == "no_pfam_hit" & Rfam_status != "no_rfam_hit" & TE_status == "no_te_overlap" ~ "Housekeeping ncRNA precursor",
      CPC2_label == "noncoding" & Pfam_status != "no_pfam_hit" & Rfam_status != "no_rfam_hit" & TE_status == "no_te_overlap" ~ "Polycistronic transcript",
      CPC2_label == "coding" & Pfam_status == "no_pfam_hit" & Rfam_status == "no_rfam_hit" & TE_status != "no_te_overlap" ~ "Polycistronic transcript",
      CPC2_label == "noncoding" & Pfam_status != "no_pfam_hit" & Rfam_status != "no_rfam_hit" & TE_status == "no_te_overlap" ~ "Polycistronic transcript",
      CPC2_label == "coding" & Pfam_status != "no_pfam_hit" & Rfam_status != "no_rfam_hit"& TE_status != "no_te_overlap" ~ "Polycistronic transcript",
      CPC2_label == "noncoding" & Pfam_status == "no_pfam_hit" & Rfam_status != "no_rfam_hit" & TE_status != "no_te_overlap" ~ "Polycistronic transcript",
      CPC2_label == "coding" & Pfam_status != "no_pfam_hit" & Rfam_status == "no_rfam_hit" & TE_status != "no_te_overlap" ~ "Polycistronic transcript",
      CPC2_label == "noncoding" & Pfam_status != "no_pfam_hit" & Rfam_status == "no_rfam_hit" & TE_status != "no_te_overlap" ~ "Polycistronic transcript",
      CPC2_label == "noncoding" & Pfam_status != "no_pfam_hit" & Rfam_status != "no_rfam_hit" & TE_status != "no_te_overlap" ~ "Polycistronic transcript"
    )
  ) -> br_intergenic_meta_labels

#rows_with_na <- br_intergenic_meta_labels[apply(
#br_intergenic_meta_labels, 
#1, 
#function(x) any(is.na(x))
#), ]

br_antisense_meta %>%
  mutate(
    Type = case_when(
      CPC2_label == "noncoding" & Pfam_status == "no_pfam_hit" & Rfam_status == "no_rfam_hit" & TE_status == "no_te_overlap"  ~ "True lncRNA",
      CPC2_label == "noncoding" & Pfam_status == "no_pfam_hit" & Rfam_status == "no_rfam_hit" & TE_status != "no_te_overlap" ~ "TE-associated lncRNA",
      CPC2_label == "coding" & Pfam_status == "no_pfam_hit" & Rfam_status == "no_rfam_hit" & TE_status == "no_te_overlap" ~ "Putative novel PCG",
      CPC2_label == "noncoding" & Pfam_status != "no_pfam_hit" & Rfam_status == "no_rfam_hit" & TE_status == "no_te_overlap" ~ "Polycistronic transcript",
      CPC2_label == "coding" & Pfam_status != "no_pfam_hit" & Rfam_status == "no_rfam_hit" & TE_status == "no_te_overlap" ~ "Putative novel PCG",
      CPC2_label == "coding" & Pfam_status != "no_pfam_hit" & Rfam_status != "no_rfam_hit" & TE_status == "no_te_overlap" ~ "Polycistronic transcript",
      CPC2_label == "coding" & Pfam_status == "no_pfam_hit" & Rfam_status != "no_rfam_hit" & TE_status == "no_te_overlap" ~ "Polycistronic transcript",
      CPC2_label == "noncoding" & Pfam_status == "no_pfam_hit" & Rfam_status != "no_rfam_hit" & TE_status == "no_te_overlap" ~ "Housekeeping ncRNA precursor",
      CPC2_label == "noncoding" & Pfam_status != "no_pfam_hit" & Rfam_status != "no_rfam_hit" & TE_status == "no_te_overlap" ~ "Polycistronic transcript",
      CPC2_label == "coding" & Pfam_status == "no_pfam_hit" & Rfam_status == "no_rfam_hit" & TE_status != "no_te_overlap" ~ "Polycistronic transcript",
      CPC2_label == "noncoding" & Pfam_status != "no_pfam_hit" & Rfam_status != "no_rfam_hit" & TE_status == "no_te_overlap" ~ "Polycistronic transcript",
      CPC2_label == "coding" & Pfam_status != "no_pfam_hit" & Rfam_status != "no_rfam_hit"& TE_status != "no_te_overlap" ~ "Polycistronic transcript",
      CPC2_label == "noncoding" & Pfam_status == "no_pfam_hit" & Rfam_status != "no_rfam_hit" & TE_status != "no_te_overlap" ~ "Polycistronic transcript",
      CPC2_label == "coding" & Pfam_status != "no_pfam_hit" & Rfam_status == "no_rfam_hit" & TE_status != "no_te_overlap" ~ "Polycistronic transcript",
      CPC2_label == "noncoding" & Pfam_status != "no_pfam_hit" & Rfam_status == "no_rfam_hit" & TE_status != "no_te_overlap" ~ "Polycistronic transcript",
      CPC2_label == "noncoding" & Pfam_status != "no_pfam_hit" & Rfam_status != "no_rfam_hit" & TE_status != "no_te_overlap" ~ "Polycistronic transcript"
    )
  ) -> br_antisense_meta_labels


# if it overlaps with more than one transcipt feature it is polycistronic
# overlaps with one then it is that *_associated

rbind(br_intergenic_meta_labels, br_antisense_meta_labels) -> br_lncRNA_meta

#save the non-coding transcript meta information
write_tsv(br_lncRNA_meta, "/mnt/Chikorita/Brassica_rapa/coexpression/core_files/br_lncRNA_metaUpd.csv")

br_lncRNA_meta %>%
  select(Transcript_ID,lncRNA_class,Type)%>%
  group_by(lncRNA_class,Type)%>%
  summarize(n = n()) -> br_lncRNA_meta_tab 

with(br_lncRNA_meta_tab, sum(n[lncRNA_class == 'Antisense'])) #4859 Upd - 
with(br_lncRNA_meta_tab, sum(n[lncRNA_class == 'Intergenic']))#12149 Upd - 

br_lncRNA_meta_tab %>%
  pivot_wider(names_from = lncRNA_class, values_from = n)


ggplot(br_lncRNA_meta_tab, aes(x="", y=n, fill=Type))+
  geom_bar(width = 1, stat = "identity") +
  scale_fill_manual(values = c("#5f0f40","#9a031e","#e36414","#0f4c5c","#143601"))+
  facet_wrap(~lncRNA_class) +
  theme_minimal()+
  theme(axis.title = element_blank(),
        text = element_text(size = 18,color = "black"),
        axis.text = element_text(color = "black"))

#How did the number of Putative novel PCG and polycistronic transcripts change?
br_lncRNA_meta %>%
  ungroup()%>%
  group_by(Type)%>%
  summarize(n = n())

#How many of the true lncRNAs are putative housekeeping ncRNA precursors
br_lncRNA_meta %>%filter(Type =="Housekeeping ncRNA precursor") %>%
  mutate(
    `Housekeeping ncRNA Type` = case_when(
      grepl("rRNA",   Rfam_status) ~ "rRNA associated",
      grepl("tRNA",   Rfam_status) ~ "tRNA associated",
      str_detect(Rfam_status, "sno|SNORD|sca_ncR26") ~ "snoRNA associated",
      str_detect(Rfam_status, "mir|MIR") ~ "microRNA associated",
      str_detect(Rfam_status, "U1|U2|U3|U4|U6|U5") ~ "Spliceosome-associated",
      str_detect(Rfam_status, "SRP") ~ "SRP associated")) -> br_lncRNA_meta_housekeeping

br_lncRNA_meta_housekeeping %>%
  filter(lncRNA_class == "Intergenic") %>%
  group_by(`Housekeeping ncRNA Type`) %>%
  summarize(n = n())-> br_lncRNA_meta_housekeeping_int

br_lncRNA_meta_housekeeping %>%
  filter(lncRNA_class == "Antisense") %>%
  group_by(`Housekeeping ncRNA Type`) %>%
  summarize(n = n())-> br_lncRNA_meta_housekeeping_anti

br_lncRNA_meta_housekeeping_int

#rows_with_na <- br_lncRNA_meta_housekeeping[apply(
#br_lncRNA_meta_housekeeping, 
#1, 
#function(x) any(is.na(x))
#), ]

bp_int<- ggplot(br_lncRNA_meta_housekeeping_int, aes(x="", y=n, fill=`Housekeeping ncRNA Type`))+
  geom_bar(width = 1, stat = "identity")
bp_int

bp_anti<- ggplot(br_lncRNA_meta_housekeeping_anti, aes(x="", y=n, fill=`Housekeeping ncRNA Type`))+
  geom_bar(width = 1, stat = "identity")
bp_anti

pie_int <- bp_int + coord_polar("y", start=0)
pie_int + scale_fill_manual(values=c("#34091e","#55062d","#802754","#cc9ab5","#ecc5dd","#f3dbea")) + theme_minimal()

pie_anti <- bp_anti + coord_polar("y", start=0)
pie_anti + scale_fill_manual(values=c("#34091e","#55062d","#802754","#cc9ab5","#ecc5dd","#f3dbea")) + theme_minimal()



#### Pull it all together ####
#what are the gene types available
br_gene_type <- read_delim("/mnt/Chikorita/Brassica_rapa/coexpression/core_files/br_gene_type.txt", 
                           delim = "\t", escape_double = FALSE, 
                           trim_ws = TRUE) #these are all PCGs

br_lncRNA_meta <- read_delim("/mnt/Chikorita/Brassica_rapa/coexpression/core_files/br_lncRNA_metaUpd2.txt", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)

#wait our lncRNA meta is at the transcript level and we want it at the gene level so let's read in our little tx2gene to make this change
#combine these
br_lncRNA_meta %>%
  filter(True_type == "True lncRNA") %>%
  filter(lncRNA_class == "Intergenic") %>%
  dplyr::select(Transcript_ID) %>%
  mutate(`Gene type` = "intergenic lncRNA")-> int_lncRNA

br_lncRNA_meta %>%
  filter(True_type == "True lncRNA") %>%
  filter(lncRNA_class == "Antisense")%>%
  dplyr::select(Transcript_ID) %>%
  mutate(`Gene type` = "antisense lncRNA") -> anti_lncRNA

rbind(int_lncRNA, anti_lncRNA) -> br_lncRNA
br_lncRNA %>%
  left_join(.,
            tx2gene_br,
            by = c("Transcript_ID" = "TXNAME")) %>%
  dplyr::select(3,2) %>%
  unique() -> br_lncRNA_type

br_lncRNA_type %>%
  group_by(`Gene type`)%>%
  summarize(n = n())
tx2gene_br

br_gene_type %>%
  setNames(c("GENEID","Gene type")) -> br_gene_type

br_lncRNA_meta %>%
  dplyr::select(1,8)%>%
  filter(True_type != "True lncRNA")-> br_novel_tx

br_novel_tx %>%
  left_join(.,
            tx2gene_br,
            by = c("Transcript_ID" = "TXNAME"))%>% #so we have 5 genes that are na we need to figure out why
  dplyr::select(3,2) %>%
  unique() %>%
  setNames(c("GENEID","Gene type"))-> br_novel_tx_type


#Combine all of these into your gene name and type
rbind(br_gene_type, br_lncRNA_type,br_novel_tx_type) %>% na.omit() -> br_all_genes_of_int

br_all_genes_of_int %>%
  group_by(`Gene type`)%>%
  summarise(n = n())

##### Save gene type mapping file ####
write.csv(br_all_genes_of_int, "/mnt/Chikorita/Brassica_rapa/coexpression/core_files/br_all_genes_of_int_type.csv")