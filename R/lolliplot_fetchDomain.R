#' fetch protein domains
#' 
#' Retrieve protein domains given ensembl transcript ID
#' @name lolliplot_fetchDomain
#' @param transcriptID String specifying ensembl transcript id
#' @param species character string to use when searching for ensemblMart dataset
#' @param host Host to connect to.
#' @return data frame of protien domains and start/stop coordinates
#' @importFrom biomaRt useMart
#' @importFrom biomaRt listDatasets
#' @importFrom biomaRt useDataset
#' @importFrom biomaRt getBM
#' 
lolliplot_fetchDomain <- function(transcriptID,
                                  species="hsapiens",
                                  host="www.ensembl.org")
{
    # display message
    message("Querying biomaRt for protein domains")
    
    # Load in mart
    ensembl_mart <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL",
                                     host="www.ensembl.org")
    
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
    values <- as.list(c(transcriptID))
    
    # Select attributes to retrieve (protein domain, start, stop)
    attributes <- c("interpro_description",
                    "interpro_start",
                    "interpro_end")
    
    # Retrieve data
    result <- biomaRt::getBM(attributes=attributes, filters=filters,
                             values=values, mart=ensembl_mart)
    
    return(result)
}