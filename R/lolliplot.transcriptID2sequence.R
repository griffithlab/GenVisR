#' Convert trascript ID to protien length
#' 
#' Retrieve the length of an ensembl transcript ID in amino acids
#' @name lolliplot.transcriptID2sequence
#' @param transcriptID String specifying ensembl transcript ID
#' @param up uniprot data object from UniProt.ws
#' @return character vector giving protien sequence

lolliplot.transcriptID2sequence <- function(transcriptID, up)
{
    # select the database to grab the key from
    kt <- "ENSEMBL_TRANSCRIPT"
    
    # query the transcript ID
    keys <- c(transcriptID)
    
    # select data to retrieve back (AA length)
    columns <- c("SEQUENCE")
    
    # submit the query to return the length in amino acid residues
    result <- UniProt.ws::select(up, keys, columns, kt)
    AAsequence <- as.character(result[2])
    
    # Format into a data frame
    AAsequence <- as.data.frame(strsplit(AAsequence, "", fixed=TRUE))
    colnames(AAsequence) <- c("AA")
    AAsequence$coord <- rownames(AAsequence)
    
    return(AAsequence)
}