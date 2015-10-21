#' fetch protein length
#' 
#' Retrieve protein length from ensembl database given enseml transcript id
#' @name lolliplot_transcriptID2length
#' @param transcriptID character string giving ensembl transcript id
#' @return length in residues of ensembl transcript id

lolliplot_transcriptID2length <- function(transcriptID)
{
    # Load in database and select dataset
    ensembl_mart <- biomaRt::useMart("ensembl")
    ensembl_mart <- biomaRt::useDataset("hsapiens_gene_ensembl", mart=ensembl_mart)
    
    # Apply various filters using vector of values
    filters <- c("ensembl_transcript_id")
    values <- as.list(c(transcriptID))
    
    # Select attributes to retrieve (protein domain, start, stop)
    attributes <- c("cds_length", "coding", "peptide")
    
    # Retrieve data
    result <- biomaRt::getBM(attributes=attributes, filters=filters,
                              values=values, mart=ensembl_mart)
    
    return(result)
}