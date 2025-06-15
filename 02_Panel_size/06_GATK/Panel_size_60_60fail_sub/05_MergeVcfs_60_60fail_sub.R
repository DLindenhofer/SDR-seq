

library(dplyr)
library(stringr)
library(data.table)
library(readr)


#__________________
#Merge vcf files here

#Load in all files without any prior qualitiy filtering and combine
file_path <- "/g/steinmetz/projects/SDR_Shared_John_Dominik/SDR002/20240918_GATK_Subset/SDR002_60_60fail_sub/04_GATK_Out"
files <- list.files(file_path, pattern = "\\.vcf$", full.names = TRUE)

#Modify skip number depending on vcf input
read_and_preprocess <- function(file_path) {
  read_tsv(file_path, skip = 324, col_names = TRUE, show_col_types = FALSE) %>%
    dplyr::select(-QUAL, -FILTER, -INFO, -FORMAT)
}

# Initialize the combined data frame with the first file's data
combined_df <- read_tsv(skip = 324, files[1], show_col_types = FALSE) %>% mutate(QUAL = ".", INFO = ".")


# Iteratively left join the rest of the files
for(i in 2:length(files)) {
  next_df <- read_and_preprocess(files[i])
  combined_df <- full_join(combined_df, next_df, by = c("#CHROM", "POS", "ID", "REF", "ALT"))
}


#Load annotation file for gene annotation and combine
annot <- read.delim("/g/steinmetz/projects/SDR_Shared_John_Dominik/SDR002/20240918_GATK_Subset/SDR002_60_60fail_sub/05_Merged_VCFs/SDR002_genomic_annotation_file.txt") %>% dplyr::select(Chromosome, Start, End, ID = ID_full)
to_add <- combined_df %>% dplyr::select("#CHROM", POS, ID, REF, ALT)
to_add <- to_add %>% mutate(`#CHROM` = sub("_chrom", "\\1", `#CHROM`))
to_add <- to_add %>% left_join(annot, by = c("#CHROM" = "ID")) %>% dplyr::mutate(POS_Chrom = POS + Start -1)

#Replace

combined_df <- combined_df %>% mutate(`#CHROM` = sub("_chrom", "\\1", `#CHROM`))
combined_df <- combined_df %>% mutate(ID = `#CHROM`)
combined_df <- combined_df %>% mutate(`#CHROM` = to_add$Chromosome, 
                                      POS = to_add$POS_Chrom)
                                      

#Get only header file - adjust according to vcf files

header <- read_tsv(files[1], col_names = FALSE) 
header <- header[1:324,]


#Write out header and merged file 
write.table(combined_df, file="/g/steinmetz/projects/SDR_Shared_John_Dominik/SDR002/20240918_GATK_Subset/SDR002_60_60fail_sub/05_Merged_VCFs/SDR002_60_60fail_sub_merged.txt", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
write.table(header, file="/g/steinmetz/projects/SDR_Shared_John_Dominik/SDR002/20240918_GATK_Subset/SDR002_60_60fail_sub/05_Merged_VCFs/header.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")



#Read in again 

#combined_df <- read.delim("/g/steinmetz/projects/SDR_Shared_John_Dominik/SDR002/20240918_GATK_Subset/SDR002_60_60fail_sub/05_Merged_VCFs/SDR002_60_60fail_sub_merged.txt", header = TRUE)
#combined_df <- combined_df %>% rename("#CHROM" = X.CHROM) #Fix chromosme naming column 

#Merge files with cat in bash


#Do some filtering for variants

nCells = length(files)
cutoff = 0.001

all_variants <- as.data.table(combined_df[,c(10:length(combined_df))])
to_add <- combined_df[,c(1:9)]

#Get genotype
extract_vcf_info <- function(x) {
  sapply(strsplit(x, ":"), function(y) y[1])
}

all_variants_GT <- all_variants[, lapply(.SD, extract_vcf_info)]
all_variants_GT <- cbind(to_add, all_variants_GT)

#Count number of cells found per edit
count_GT <- all_variants_GT[,c(10:length(all_variants_GT))] %>% 
  apply(1, function(x) table(factor(x, levels = c("0/0", "0/1", "1/1")))) %>%
  as.data.frame()

count_GT <- as.data.frame(t(count_GT))

Variant_meta <- to_add %>% mutate(count_GT) 
Variant_meta <- Variant_meta %>% mutate(sum_variant = Variant_meta$`0/1` + Variant_meta$`1/1`)


#Filter out combined_df
combined_df_filtered <- cbind(combined_df, Variant_meta$sum_variant) 
combined_df_filtered <- combined_df_filtered %>% filter(Variant_meta$sum_variant > round(nCells*cutoff))
combined_df_filtered <- combined_df_filtered %>% dplyr::select(-"Variant_meta$sum_variant")

combined_df_filtered <- combined_df_filtered %>% mutate(QUAL = ".", INFO = ".")


write.table(combined_df_filtered, file="/g/steinmetz/projects/SDR_Shared_John_Dominik/SDR002/20240918_GATK_Subset/SDR002_60_60fail_sub/05_Merged_VCFs/SDR002_60_60fail_sub_merged_cut.txt", quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")


#Merge file with header using cat in bash

#Then uplaod to VEP ensemble website

#https://www.ensembl.org/info/docs/tools/vep/index.html


