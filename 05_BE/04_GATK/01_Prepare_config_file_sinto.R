


library(dplyr)
library(stringr)
library(here)


here()
#Make input array with all cells to submit them indiviudally

#Get file names
#From SDRranger script - define cells which are high qualitty to do variant calling on
#File input is split into multiple batches for better processing

file_names <- list.files(path = "01_Input_cells_txt", pattern = "\\.txt$", full.names = FALSE)


file_path_add <- "01_Input_cells_txt/"

chunks_to_run <- data.frame(Chunk_ID = file_names)


# Assuming chunks_to_run, file_path_add, and Chunk_ID are already defined
chunks_to_run <- chunks_to_run %>%
  mutate(path = str_c(file_path_add, Chunk_ID, sep = "")) %>% # Adjust sep as needed
  mutate(Chunk_ID = sub("\\.txt$", "", Chunk_ID))

#Write out config file for submitting sinto jobs as an array.
write.table(chunks_to_run, file="01_config_file_sinto.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
