#' NoiseCallAmp
#'
#' This function processes amplicon sequencing data to identify and summarize base content, variant miscalls, and insertion/deletion frequencies. It integrates positional information and normalizes variant frequencies based on the base composition.
#' 
#' The function generates data on conversion frequencies (base-to-base changes) and insertion/deletion events, providing output files with detailed information on each type of event. The results can be used to plot or further analyze noise in amplicon sequencing data.
#'
#' @param Geno_meta A data frame containing SNP information for each amplicon with columns `Amp_ID`, `Start`, `End`, `Start_seq`, `End_seq`, `X.CHROM`, `Prim_leng_for`, `Prim_leng_rev`, `POS`, `REF`, `ALT`, `NA_count`, `REF_count`, `HET_count`, and `ALT_count`.
#' @param ref_seqs A data frame containing reference sequences for the amplicons with columns `ID` and `Sequence`.
#' @param outdir A string specifying the output directory to save the results.
#' @return Three CSV files saved in the specified `outdir`: `base_cont_df.csv`, `con_freq_df.csv`, and `inde_freq_df.csv`. Each file contains data about base content, conversion frequencies, and insertion/deletion frequencies, respectively.
#' @return A list containing three data frames: `base_cont_df`, `con_freq_df`, and `inde_freq_df`. Also writes these data frames as CSV files to the specified `outdir`.
#' @import dplyr
#' @import tidyr
#' @import stringr
#' @export


NoiseCallAmp <- function(Geno_meta, ref_seqs, nCells ){
  
  #Amplicons to loop for 
  
  to_plot <- unique(Geno_meta$Amp_ID)
  
  
  #Initiate list to store results in 
  
  con_freq_list <- list()
  inde_freq_list <- list()
  base_cont_list <- list()
  
  for(i in to_plot){
    
    Geno_meta_i <- Geno_meta %>% filter(Amp_ID == i)
    start_amp <- unique(Geno_meta_i$Start)
    end_amp <- unique(Geno_meta_i$End)
    start_seq <- unique(Geno_meta_i$Start_seq)
    end_seq <- unique(Geno_meta_i$End_seq)
    chromosome <- unique(Geno_meta_i$X.CHROM)
    prim_len_for <- unique(Geno_meta_i$Prim_leng_for)
    prim_len_rev <- unique(Geno_meta_i$Prim_leng_rev)
    
    
    #Select columns for plotting to reshape
    
    Geno_meta_i_sel <- Geno_meta_i %>% 
      dplyr::select(ID_Geno, NA_count, REF_count, HET_count, ALT_count)
    
    #Select columns for merging later to add positional information 
    Geno_meta_POS_i <- Geno_meta_i %>% dplyr::select(ID_Geno, POS, REF, ALT)
    
    
    #Add to both the positions with not called variants 
    
    #Determine most frequent NA count to add to ref seqs
    
    NA_count_REF_2_i <- as.numeric(names(which.max(table(Geno_meta_i_sel$NA_count))))
    REF_count_REF_2_i <- 1-NA_count_REF_2_i
    
    #Get all positions and reference sequence
    all_pos_i <- seq(start_amp, end_amp)
    to_fill_i <- setdiff(all_pos_i, Geno_meta_POS_i$POS)
    
    ref_seq_i <- ref_seqs %>% filter(ID == i)
    ref_seq_i <- as.character(ref_seq_i$Sequence)
    ref_seq_i <- strsplit(ref_seq_i, "")[[1]]
    
    
    #Make dataframe with all positions and counts to add
    to_add_refs_i <- data.frame(POS = all_pos_i, REF = ref_seq_i)
    
    #Filter for the ones that are missing 
    to_add_refs_i <- to_add_refs_i %>% filter(POS %in% to_fill_i)
    
    Geno_meta_i_sel_add <- to_add_refs_i %>% mutate(ID_Geno = str_c(chromosome, POS, REF, sep = "-"), 
                                                    NA_count = NA_count_REF_2_i, 
                                                    REF_count = REF_count_REF_2_i, 
                                                    HET_count = 0,
                                                    ALT_count = 0) %>% 
      dplyr::select(-POS, -REF)
    
    
    Geno_meta_POS_i_add <- to_add_refs_i %>% mutate(ID_Geno = str_c(chromosome, POS, REF, sep = "-"), 
                                                    ALT = REF) %>%
      dplyr::select(ID_Geno, POS, REF, ALT)
    
    
    #Add to both 
    Geno_meta_i_sel <- rbind(Geno_meta_i_sel, Geno_meta_i_sel_add)
    Geno_meta_POS_i <- rbind(Geno_meta_POS_i, Geno_meta_POS_i_add )
    
    
    #Make long for plotting 
    Geno_meta_i_long <- Geno_meta_i_sel %>% 
      pivot_longer(cols = -ID_Geno, names_to = "Geno", values_to = "frac")
    
    Geno_meta_i_long <- Geno_meta_i_long %>% left_join(Geno_meta_POS_i, by = "ID_Geno")
    
    #Check length again and if all POS are there
    #length(unique(Geno_meta_i_long$POS))
    
    #Clean up for variant miscalling plots
    #1. Remove HET variants
    #HET_var_i <- Geno_meta_i_sel %>% filter(HET_count+ALT_count > cut_off_High_Var) %>% pull(ID_Geno)
    #Geno_meta_i_long <- Geno_meta_i_long %>% filter(!(ID_Geno %in% HET_var_i))
    
    
    #2. Exclude primer regions
    #POS_for <- seq(start_amp, start_seq-1)
    #POS_rev <- seq(end_seq+1, end_amp)
    #Geno_meta_i_long <- Geno_meta_i_long %>% filter(!(POS %in% POS_for))
    #Geno_meta_i_long <- Geno_meta_i_long %>% filter(!(POS %in% POS_rev))
    
    #Check length again and if all POS are there
    length(unique(Geno_meta_i_long$POS)) 
    
    #3.Only select low frequency ones
    #LOW_var_i <- Geno_meta_i_sel %>% filter(HET_count+ALT_count < cut_off_Low_Var) %>% pull(ID_Geno)
    #Geno_meta_i_long <- Geno_meta_i_long %>% filter(ID_Geno %in% LOW_var_i)
    
    
    #3.1 - Now calculate the base content for normalization later
    ref_seq_i <- ref_seqs %>% filter(ID == i)
    ref_seq_i <- as.character(ref_seq_i$Sequence)
    ref_seq_i <- strsplit(ref_seq_i, "")[[1]]
    
    #Make dataframe with all positions and counts to add
    to_add_refs_i <- data.frame(POS = all_pos_i, REF = ref_seq_i)
    
    #Remove primer regions
    POS_for <- seq(start_amp, start_seq-1)
    POS_rev <- seq(end_seq+1, end_amp) 
    to_add_refs_i <- to_add_refs_i %>% filter(!(POS %in% POS_for))
    to_add_refs_i <- to_add_refs_i %>% filter(!(POS %in% POS_rev))
    
    base_cont <- to_add_refs_i %>% group_by(REF) %>%
      summarize(base_count = n()) %>% 
      mutate(base_sum = sum(base_count))
    
    base_cont_add <- base_cont %>% mutate(Amp_ID = i)
    
    
    #Add to list 
    
    base_cont_list[[i]] <- base_cont_add
    
    
    #4.Remove alleles with non called variants
    Geno_meta_i_long <- Geno_meta_i_long %>% filter(REF != ALT) 
    
    
    
    
    #Split into conversion and insertions/deletions
    Geno_meta_i_long_con <- Geno_meta_i_long %>% filter(str_count(REF) <2 &
                                                          str_count(ALT) < 2)
    
    
    Geno_meta_i_long_inde <- Geno_meta_i_long %>% filter(str_count(REF) >=2 |
                                                           str_count(ALT) >= 2)
    
    #A. For con
    
    #Make frequencies of miscalled variants
    Geno_meta_i_long_con <- Geno_meta_i_long_con %>% filter(Geno %in% c("HET_count", "ALT_count"))
    con_freq <- Geno_meta_i_long_con %>% group_by(REF, ALT) %>% summarize(frac = sum(frac), count_POS = n()/2)
    
    con_freq <- con_freq %>% mutate(nCells_Con = frac*nCells,
                                    Conversion = str_c(REF, ALT, sep = "_"), 
                                    Amp_ID = i) 
    
    con_freq <- con_freq %>% left_join(base_cont, by = c("REF" = "REF"))
    
    #Add to list 
    con_freq_list[[i]] <- con_freq
    
    
    #B. For inde
    
    Geno_meta_i_long_inde <- Geno_meta_i_long_inde %>% filter(Geno %in% c("HET_count", "ALT_count"))
    inde_freq <- Geno_meta_i_long_inde %>% group_by(REF, ALT) %>% summarize(frac = sum(frac), count_POS = n()/2)
    inde_freq <- inde_freq %>% mutate(nCells_Con = frac*nCells,
                                      Conversion = str_c(REF, ALT, sep = "_"), 
                                      Amp_ID = i) 
    
    inde_freq <- inde_freq %>% mutate(Size = str_count(ALT) - str_count(REF)) %>%
      mutate(Type = if_else(Size > 0, "Ins", "Del"))
    
    inde_freq <- inde_freq %>% mutate(REF2 = substr(REF, 1, 1))
    
    inde_freq <- inde_freq %>% left_join(base_cont, by = c("REF2" = "REF"))
    
    
    #Add to list 
    
    inde_freq_list[[i]] <- inde_freq
    
    
    
  }
  
  
  #Make df of lists
  
  base_cont_df <- bind_rows(base_cont_list)
  con_freq_df <- bind_rows(con_freq_list)
  inde_freq_df <- bind_rows(inde_freq_list)
  
  
  # Return the three data frames to the R environment
  return(list(base_cont_df = base_cont_df, con_freq_df = con_freq_df, inde_freq_df = inde_freq_df))
}

