#' Hiearchical removal of MAF entries
#' 
#' Remove MAF entries with the same gene/sample in an ordered fashion such that the most deleterious are retained
#' @name hiearchial_remove_trv_type
#' @param data_frame a data frame in long format with columns sample, gene, trv_type
#' @param file_type The type of file to act on one of 'MAF" or 'TGI'
#' @return a data frame with multiple mutations in the same sample/gene collapsed on the most deleterious 

hiearchial_remove_trv_type <- function(data_frame, file_type)
{
  #############################################################################################################################
  ####### Function to hiearchial collapse the data frame on sample/gene leaving only the most deleterious trv_type, ###########
  ####### the function expects a dataframe with column names sample, gene, trv_type                                 ###########
  #############################################################################################################################
  
  # reorder the trv_type in terms of deleterious effect and refactor the data frame
  if(toupper(file_type) == toupper('MGI'))
  {
    mutation_order <- c("nonsense", "frame_shift_del", "frame_shift_ins", "splice_site_del", "splice_site_ins", "splice_site", "nonstop", "in_frame_del", "in_frame_ins", "missense", "splice_region", "5_prime_flanking_region", "3_prime_flanking_region", "3_prime_untranslated_region", "5_prime_untranslated_region", "rna", "intronic", "silent")
    mutation_order <- c("nonsense", "frame_shift_del", "frame_shift_ins", "splice_site_del", "splice_site_ins", "splice_site", "nonstop", "in_frame_del", "in_frame_ins", "missense", "splice_region", "5_prime_flanking_region", "3_prime_flanking_region", "3_prime_untranslated_region", "5_prime_untranslated_region", "rna", "intronic", "silent", "negative", "positive", "unknown")    
  } else if (toupper(file_type) == toupper('MAF'))
  {
    mutation_order <- c("Nonsense_Mutation", "Frame_Shift_Ins", "Frame_Shift_Del", "In_Frame_Ins", "In_Frame_Del", "Nonstop_Mutation", "Splice_Site", "Missense_Mutation", "5\'Flank", "3\'Flank", "5\'UTR", "3\'UTR", "RNA", "Intron", "IGR", "Silent", "Targeted_Region")
  }
  
  data_frame$trv_type <- factor(data_frame$trv_type, levels=mutation_order)
  
  # sort the data frame so that the duplicated call will remove the proper trv_type
  data_frame <- data_frame[order(data_frame$sample, data_frame$gene, data_frame$trv_type),]
  
  # collapse the data on sample/gene
  data_frame <- data_frame[!duplicated(data_frame[, c("sample", "gene")]), ]
  
  return(data_frame)
}