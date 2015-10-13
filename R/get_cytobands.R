#' obtain cytogenetic bands
#'
#' given a genome query UCSC for cytogenetic band locations
#' @name get_cytobands
#' @param genome character string giving a UCSC genome
#' @return object of class data frame
#' @import RMySQL
#' @import DBI

get_cytobands <- function(genome)
{
    conn <- dbConnect(MySQL(), user="genome", host="genome-mysql.cse.ucsc.edu", dbname=genome)
    result <- dbGetQuery(conn=conn, statement="SELECT chrom, chromStart, chromEnd, name, gieStain FROM cytoBand")
    dbDisconnect(conn)
    return(result)
}
