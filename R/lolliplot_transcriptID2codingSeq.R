#' fetch protein length
#' 
#' Retrieve protein length from ensembl database given enseml transcript id
#' @name lolliplot_transcriptID2codingSeq
#' @param transcriptID character string giving ensembl transcript id
#' @param species character string to use when searching for ensemblMart dataset
#' @return length in residues of ensembl transcript id

lolliplot_transcriptID2codingSeq <- function(transcriptID, species="hsapiens")
{
    # Load in mart
    ensembl_mart <- biomaRt::useMart("ensembl")
    
    # select proper data set given regexp print warnings if unexpected out occur
    index <- which(grepl(species, biomaRt::listDatasets(ensembl_mart)$dataset))
    if(length(index)>1)
    {
        memo <- paste0(species, " Matches more than one dataset for the",
                       " ensembl mart, please specify a species in the, ",
                       "following format: hsapiens")
        stop(memo)
    } else if(length(index)==0) {
        memo <- paste0(species, " does not appear to be supported by biomaRt",
                       "if you beleive this to be in error please modify", 
                       "you're input to to conform to this format: hsapiens")
        stop(memo)
    }
    ensembl_mart <- biomaRt::useDataset(as.character(biomaRt::listDatasets(ensembl_mart)$dataset[index]),
                                        mart=ensembl_mart)
    
    # Apply various filters using vector of values
    filters <- c("ensembl_transcript_id")
    values <- as.character(transcriptID)
    
    # Select attributes to retrieve coding dna sequence
    attributes <- as.list(c("coding"))
    
    # Retrieve data
    result <- biomaRt::getBM(attributes=attributes, filters=filters,
                              values=values, mart=ensembl_mart)
    
    return(result)
}