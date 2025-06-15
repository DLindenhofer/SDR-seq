

library(dplyr)
library(tidyverse)
library(stringr)
library(here)



#Make input array with all cells to submit them indiviudally



#Get file names
file_names <- list.files(path = "/Volumes/steinmetz/projects/SDR_Shared_John_Dominik/SDR002/20240918_GATK_Subset/SDR002_240_240_sub/02_Input_cells_bam", pattern = "\\.bam$", full.names = FALSE)

file_path_add <- "/g/steinmetz/projects/SDR_Shared_John_Dominik/SDR002/20240918_GATK_Subset/SDR002_240_240_sub/02_Input_cells_bam/"



chunks_to_run <- data.frame(Chunk_ID = file_names)

#Remove the currently processed intermediate bam from preparation bash script
chunks_to_run <- chunks_to_run %>% filter(!str_detect(Chunk_ID, "withRG"))

# Assuming chunks_to_run, file_path_add, and Chunk_ID are already defined
chunks_to_run <- chunks_to_run %>%
  mutate(path = str_c(file_path_add, Chunk_ID, sep = "")) %>% # Adjust sep as needed
  mutate(Chunk_ID = sub("\\.bam$", "", Chunk_ID))

#write.table(chunks_to_run, file="/Volumes/steinmetz/projects/SDR_Shared_John_Dominik/SDR005/20240209_GATK/SDR001_Run1/02_input_array_GATK.txt", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

#Make smaller chunks for submission and submit individual batches as individual scriptes (array within script)


#Need to also make smaller chunks as otherwise too many files are handled

out_dir <- "/Volumes/steinmetz/projects/SDR_Shared_John_Dominik/SDR002/20240918_GATK_Subset/SDR002_240_240_sub/02_Input_GATK/"

# Set the batch size - samtools mileup gives an issue with high input reads somehow -try to stay below one million 
batch_size <- 1000

# Calculate the number of files needed
num_files <- ceiling(nrow(chunks_to_run) / batch_size)

# Loop through each batch and write to a file
for (file_num in 1:num_files) {
  # Calculate row indices for the current batch
  start_row <- (file_num - 1) * batch_size + 1
  end_row <- min(file_num * batch_size, nrow(chunks_to_run))
  
  # Create a file name for the current batch
  file_name <- paste0(out_dir, "Batch_", file_num, ".txt")
  
  # Write the current batch of rows to a file
  write.table(chunks_to_run[start_row:end_row, ], file = file_name, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
}
