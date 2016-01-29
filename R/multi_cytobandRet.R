#' Retrieve cytogenetic bands
#'
#' given a genome query UCSC for cytogenetic band locations
#' @name multi_cytobandRet
#' @param genome character string giving a UCSC genome
#' @return object of class data frame
#' @importFrom DBI dbConnect
#' @importFrom RMySQL MySQL
#' @importFrom DBI dbGetQuery
#' @importFrom DBI dbDisconnect

multi_cytobandRet <- function(genome)
{
    # connect to UCSC mySQL genome database
    conn <- DBI::dbConnect(RMySQL::MySQL(), user="genome",
                           host="genome-mysql.cse.ucsc.edu",
                           dbname=genome)
    
    # query for chromosome positions and cytogenetic staining
    result <- DBI::dbGetQuery(conn=conn,
                              statement="SELECT chrom,
                              chromStart, chromEnd, name,
                              gieStain FROM cytoBand")
    DBI::dbDisconnect(conn)
    
    return(result)
}
