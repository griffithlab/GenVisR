#' Hiearchical removal of MAF entries
#' 
#' Remove MAF entries with the same gene/sample in an ordered fashion such that
#' the most deleterious are retained
#' @name waterfall_hierarchyTRV
#' @param x a data frame in long format with columns sample, gene,
#' trv_type
#' @param file_type The type of file to act on one of 'MAF" or 'TGI'
#' @return a data frame with multiple mutations in the same sample/gene
#' collapsed on the most deleterious 

waterfall_hierarchyTRV <- function(x, file_type)
{
    # reorder the trv_type in terms of deleterious effect and refactor
    # the data frame
    if(toupper(file_type) == toupper('MGI'))
    {
        mutation_order <- c("nonsense", "frame_shift_del", "frame_shift_ins",
                            "splice_site_del", "splice_site_ins", "splice_site",
                            "nonstop", "in_frame_del", "in_frame_ins",
                            "missense", "splice_region",
                            "5_prime_flanking_region",
                            "3_prime_flanking_region",
                            "3_prime_untranslated_region",
                            "5_prime_untranslated_region", "rna", "intronic",
                            "silent")
    } else if (toupper(file_type) == toupper('MAF')) {
        mutation_order <- c("Nonsense_Mutation", "Frame_Shift_Ins",
                            "Frame_Shift_Del", "In_Frame_Ins", "In_Frame_Del",
                            "Nonstop_Mutation", "Splice_Site",
                            "Missense_Mutation", "5\'Flank", "3\'Flank",
                            "5\'UTR", "3\'UTR", "RNA", "Intron", "IGR",
                            "Silent", "Targeted_Region")
    }
    
    # Check that elements in trv_type are in the mutation order
    if(any(!x$trv_type %in% mutation_order))
    {
        stop("Detected an invalid mutation type, valid values for ", file_type,
             " are ", mutation_order)
    }
    x$trv_type <- factor(x$trv_type, levels=mutation_order)
    
    # sort the data frame so that the duplicated call will remove the
    # proper trv_type
    x <- x[order(x$sample, x$gene, x$trv_type),]
    
    # collapse the data on sample/gene
    x <- x[!duplicated(x[, c("sample", "gene")]), ]
  
  return(x)
}