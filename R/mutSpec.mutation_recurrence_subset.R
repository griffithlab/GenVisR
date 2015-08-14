#' Mutation Recurrence Cutoff
#' 
#' Subset a MAF file keeping only samples that meet a mutation recurrence cutoff
#' @name mutSpec.mutation_recurrence_subset
#' @param data_frame a data frame in long format with columns 'gene', 'trv_type'
#' @param recurrence_cutoff an integer, remove genes whose total mutations are
#' less than this amount
#' @return a subset data frame

mutSpec.mutation_recurrence_subset <- function(data_frame, recurrence_cutoff)
{
    # convert the dataframe to a table of mutation counts at the gene level
    mutation_recurrence <- table(data_frame[,c('gene', 'trv_type')],
                                 useNA="ifany")
    
    # add totals to the data frame (total mutations for each gene and total
    # mutations for each mutation type), remove the last row
    # (total mutations for each mutation)
    mutation_recurrence_total <- as.data.frame.matrix(addmargins(mutation_recurrence,
                                                                 FUN = list(Total = sum),
                                                                 quiet = TRUE))
    mutation_recurrence_total <- mutation_recurrence_total[-nrow(mutation_recurrence_total),]
    
    # If recurrence cutoff specified exceeds upper limit such that no plot
    # useful would be generated, reset recurrence cutoff
    if(max(mutation_recurrence_total$Total) < recurrence_cutoff)
    {
        message <- paste0('The recurrence cutoff specified exceeds the
                          recurrence seen in the data, resetting this value
                          to equal max recurrence:',
                          max(mutation_recurrence_total$Total))
        warning(message)
        recurrence_cutoff <- max(mutation_recurrence_total$Total)
    }
    
    # obtain a vector of gene names where the genes have mutation recurrence
    # greater than or equal to the recurrence cutoff
    gene_above_recurrence <- sort(row.names(subset(mutation_recurrence_total,
                                                   mutation_recurrence_total$Total >= recurrence_cutoff)))
    
    # add NA to the end of 'gene_above_recurrence' vector, allowing for all
    # samples having NA as a gene name to be retained in the subset below
    gene_above_recurrence <- c(gene_above_recurrence, NA)
    
    # subset the original data frame based on the following: keep gene if it is
    # in the gene vector in "mutation_recurrence_subset"
    subset_data_frame <- data_frame[(data_frame$gene %in% gene_above_recurrence), ]
  
  return(subset_data_frame)
}