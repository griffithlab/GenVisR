#' fetch protein length
#' 
#' Retrieve protein length from ensembl database given enseml transcript id
#' @name lolliplot_transcriptID2codingSeq
#' @param transcriptID character string giving ensembl transcript id
#' @param species character string to use when searching for ensemblMart dataset
#' @param host Host to connect to.
#' @return length in residues of ensembl transcript id
#' @importFrom biomaRt useMart
#' @importFrom biomaRt listDatasets
#' @importFrom biomaRt useDataset
#' @importFrom biomaRt getBM

lolliplot_transcriptID2codingSeq <- function(transcriptID,
                                             species="hsapiens",
                                             host="www.ensembl.org")
{
    # display mesage
    memo <- paste0("Using the following host: ", host, " for biomaRt queries",
                   " to change the ensembl annotation version alter this",
                   " parameter!")
    message(memo)
    message("Querying biomaRt for transcript sequence")

    # Load in mart
    ensembl_mart <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL",
                                     host=host)
    
    # select proper data set given regexp print warnings if unexpected out occur
    index <- which(grepl(species, biomaRt::listDatasets(ensembl_mart)$dataset))
    if(length(index)>1)
    {
        memo <- paste0(species, " Matches more than one dataset for the",
                       " ensembl mart, please specify a species in the, ",
                       "following format: hsapiens")
        stop(memo)
    } else if(length(index)==0) {
        valid_species <- toString(gsub("_gene_ensembl",
                                       "",
                                       biomaRt::listDatasets(ensembl_mart)$dataset))
        
        memo <- paste0(species, " does not appear to be supported by biomaRt",
                       " please specify one of the following species:",
                       valid_species)
        stop(memo)
    }
    ensembl_mart <- biomaRt::useDataset(as.character(biomaRt::listDatasets(ensembl_mart)$dataset[index]),
                                        mart=ensembl_mart)
    
    # Apply various filters using vector of values
    filters <- c("ensembl_transcript_id")
    ensg_id <- as.character(transcriptID)
    
    # Select attributes to retrieve coding dna sequence
    attributes <- as.list(c("coding","cds_length"))
    
    # Retrieve data
    result <- biomaRt::getBM(attributes=attributes, filters=filters,
                              values=ensg_id, mart=ensembl_mart)
    
    return(as.list(result))
}