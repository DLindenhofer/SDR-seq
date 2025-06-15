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
                         "Left_in_repeat" )
  
  #Select relevant columns, make ID  and GRanges
  repeats <- repeats %>% dplyr::select(Chromosome, Start, End, Strand, Name, Family, Class, SW_Score, Divergence)
  repeats <- repeats  %>% mutate(ID_Rep = str_c(Chromosome, ":", Start, "-", End, "-", Name))
  
  
  repeats_gr <- makeGRangesFromDataFrame(repeats, keep.extra.columns = TRUE)
  
  #Prepare data file  
  data_mod <- data %>% dplyr::select(Chromosome, Start = POS, End = POS, ID_Geno)
  data_mod_gr <- makeGRangesFromDataFrame(data_mod, keep.extra.columns = TRUE)
  
  
  # Find overlaps between the two GRanges objects
  Overlap <- findOverlaps(data_mod_gr, repeats_gr)
  
  # Initialize a CharacterList with empty lists for all rows in Geno_meta_i_gr
  Repeat <- CharacterList(vector("list", length(data_mod_gr)))
  
  
  # Replace the empty values in Repeats with the actual Repeatss from overlaps
  # Split Repeatss by the queryHits (i.e., the index of Geno_meta_i_gr)
  overlap_repeats <- split(repeats_gr$ID_Rep[subjectHits(Overlap)], queryHits(Overlap))
  
  
  # Assign the Repeats to the correct positions in the Repeat list
  Repeat[as.numeric(names(overlap_repeats))] <- overlap_repeats
  
  
  # Add the Repeat column as metadata to Geno_meta_i_gr
  mcols(data_mod_gr ) <- DataFrame(mcols(data_mod_gr), Repeat)
  
  
  # View the result
  data_mod_add <- data_mod_gr %>% as.data.frame() %>% tibble::as_tibble() %>% dplyr::select(ID_Geno, Repeat)
  
  
  
  data <- data %>% left_join(data_mod_add, by = "ID_Geno")
  
  return(data)
  
  
}
