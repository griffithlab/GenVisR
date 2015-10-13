#' Convert trascript ID to uniprot ID
#'
#' Convert an ensembl transcript ID to the corresponding uniprot ID
#' @name lolliplot.transcriptID2uniprotID
#' @param transcriptID String specifying ensembl transcript ID
#' @param up uniprot data object from UniProt.ws
#' @return String specifying uniprot ID

lolliplot.transcriptID2uniprotID <- function(transcriptID, up)
{
    # select the database to grab the key from
    kt <- "ENSEMBL_TRANSCRIPT"

    # Query the transcript ID
    keys <- c(transcriptID)

    # Select data to retrieve back (uniprot ID)
    columns <- c("UNIPROTKB")

    # submit the query to return the uniprot ID
    result <- UniProt.ws::select(up, keys, columns, kt)
    uniprotID <- as.character(result[2])

    return(uniprotID)
}
