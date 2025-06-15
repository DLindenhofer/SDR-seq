#' RepeatAnnot
#'
#' This function annotates a given dataset with repeat sequence information based on chromosomal coordinates.
#' It overlaps query genomic positions with repetitive sequence regions (such as those from RepeatMasker), adding repeat metadata to the original dataset.
#'
#' The function uses genomic ranges (GRanges) to efficiently find overlaps between the query data and the repeat annotations. The final output includes a new column with repeat IDs for each query.
#'
#' @param data A data frame containing genomic information with columns `Chromosome`, `POS`, and `ID_Geno` representing the genomic coordinates and unique identifiers for each row.
#' @param repeats A data frame containing repeat annotations, such as those from RepeatMasker, with columns `Chromosome`, `Start`, `End`, `Strand`, `Name`, `Family`, `Class`, `SW_Score`, and `Divergence`. 
#' @return The function returns the input `data` data frame with an additional `Repeat` column. This column contains repeat information for regions that overlap with the provided repeat annotations.
#' @import dplyr
#' @import tidyr
#' @import stringr
#' @import GenomicRanges
#' @import S4Vectors
#' @export



RepeatAnnot <- function(data, repeats){
  
  #Add colnames
  colnames(repeats) <- c("V1",
                         "SW_Score",
                         "Divergence", 
                         "Deletions",
                         "Insertions",
                         "Chromosome",
                         "Start",
                         "End", 
                         "V9",
                         "Strand",
                         "Name",
                         "Family",
                         "Class",
                         "Begin_in_repeat", 
                         "End_in_repeat", 
                         "Left_in_repeat",
                         "V17")
  
  #Calculate width
  repeats <- repeats %>% filter(!is.na(Start) | !is.na(End)) %>%
    mutate(Width = End - Start +1)
  
  #Select relevant columns, make ID  and GRanges
  repeats <- repeats %>% dplyr::select(Chromosome, Start, End, Strand, Name, Family, Class, SW_Score, Divergence, Width_Rep = Width)
  repeats <- repeats %>% mutate(Start_Rep = Start, End_Rep = End)
  repeats <- repeats  %>% mutate(ID_Rep = str_c(Chromosome, ":", Start, "-", End, "-", Name))
  
  
  repeats_gr <- makeGRangesFromDataFrame(repeats, keep.extra.columns = TRUE)
  
  #Prepare data file  
  data_mod <- data %>% dplyr::select(Chromosome, Start = POS, End = POS, ID_Geno)
  data_mod_gr <- makeGRangesFromDataFrame(data_mod, keep.extra.columns = TRUE)
  
  
  # Find overlaps between the two GRanges objects
  Overlap <- findOverlaps(data_mod_gr, repeats_gr)
  
  
  #_________________________
  #For ID_Rep
  
  # Initialize a CharacterList with empty lists for all rows in Geno_meta_i_gr
  Repeat <- CharacterList(vector("list", length(data_mod_gr)))
  
  
  # Replace the empty values in Repeats with the actual Repeats from overlaps
  # Split Repeatss by the queryHits (i.e., the index of Geno_meta_i_gr)
  overlap_repeats <- split(repeats_gr$ID_Rep[subjectHits(Overlap)], queryHits(Overlap))
  
  
  # Assign the Repeats to the correct positions in the Repeat list
  Repeat[as.numeric(names(overlap_repeats))] <- overlap_repeats
  
  
  # Add the Repeat column as metadata to Geno_meta_i_gr
  mcols(data_mod_gr ) <- DataFrame(mcols(data_mod_gr), Repeat)
  
  
  #_________________________
  #For Start_Rep
  
  # Initialize a CharacterList with empty lists for all rows in Geno_meta_i_gr
  Start_Rep <- CharacterList(vector("list", length(data_mod_gr)))
  
  
  # Replace the empty values in Repeats with the actual Repeats from overlaps
  # Split Repeatss by the queryHits (i.e., the index of Geno_meta_i_gr)
  overlap_start_rep <- split(repeats_gr$Start_Rep[subjectHits(Overlap)], queryHits(Overlap))
  
  
  # Assign the Repeats to the correct positions in the Repeat list
  Start_Rep[as.numeric(names(overlap_start_rep))] <- overlap_start_rep
  
  
  # Add the Repeat column as metadata to Geno_meta_i_gr
  mcols(data_mod_gr ) <- DataFrame(mcols(data_mod_gr), Start_Rep) 
  
  #_________________________
  #For End_Rep
  
  # Initialize a CharacterList with empty lists for all rows in Geno_meta_i_gr
  End_Rep <- CharacterList(vector("list", length(data_mod_gr)))
  
  
  # Replace the empty values in Repeats with the actual Repeats from overlaps
  # Split Repeatss by the queryHits (i.e., the index of Geno_meta_i_gr)
  overlap_end_rep <- split(repeats_gr$End_Rep[subjectHits(Overlap)], queryHits(Overlap))
  
  
  # Assign the Repeats to the correct positions in the Repeat list
  End_Rep[as.numeric(names(overlap_end_rep))] <- overlap_end_rep
  
  
  # Add the Repeat column as metadata to Geno_meta_i_gr
  mcols(data_mod_gr ) <- DataFrame(mcols(data_mod_gr), End_Rep) 
  
  
  #_________________________
  #For Width_Rep
  
  # Initialize a CharacterList with empty lists for all rows in Geno_meta_i_gr
  Width_Rep <- CharacterList(vector("list", length(data_mod_gr)))
  
  
  # Replace the empty values in Repeats with the actual Repeats from overlaps
  # Split Repeatss by the queryHits (i.e., the index of Geno_meta_i_gr)
  overlap_width_rep <- split(repeats_gr$Width_Rep[subjectHits(Overlap)], queryHits(Overlap))
  
  
  # Assign the Repeats to the correct positions in the Repeat list
  Width_Rep[as.numeric(names(overlap_width_rep))] <- overlap_width_rep
  
  
  # Add the Repeat column as metadata to Geno_meta_i_gr
  mcols(data_mod_gr ) <- DataFrame(mcols(data_mod_gr), Width_Rep) 
  
  
  # View the result
  data_mod_add <- data_mod_gr %>% as.data.frame() %>% tibble::as_tibble() %>% dplyr::select(ID_Geno, Repeat, Start_Rep, End_Rep, Width_Rep)
  
  
  #Convert the list columns to normal columns
  
  data_mod_add$Repeat <- sapply(data_mod_add$Repeat, function(x) paste(unlist(x), collapse = "|"))
  data_mod_add$Width_Rep <- sapply(data_mod_add$Width_Rep, function(x) paste(unlist(x), collapse = "|"))
  
  #Keep only the biggest Start 
  
  # Assuming your dataframe is called df and the column with lists is named 'list_column'
  data_mod_add$Start_Rep <- lapply(data_mod_add$Start_Rep, function(x) {
    if (length(x) > 0) {
      max(unlist(x))  # Unlist to convert list to numeric vector and then take the max
    } else {
      NA  # Return NA for empty lists
    }
  })
  
  # If you want to convert the list column back into a numeric column
  data_mod_add$Start_Rep <- as.numeric( data_mod_add$Start_Rep)
  
  
  
  #Keep only the smallest End
  
  # Assuming your dataframe is called df and the column with lists is named 'list_column'
  data_mod_add$End_Rep <- lapply(data_mod_add$End_Rep, function(x) {
    if (length(x) > 0) {
      min(unlist(x))  # Unlist to convert list to numeric vector and then take the min
    } else {
      NA  # Return NA for empty lists
    }
  })
  
  # If you want to convert the list column back into a numeric column
  data_mod_add$End_Rep <- as.numeric(data_mod_add$End_Rep)
  
  
  
  
  #Join in
  data <- data %>% left_join(data_mod_add, by = "ID_Geno")
  
  return(data)
  
  
}
