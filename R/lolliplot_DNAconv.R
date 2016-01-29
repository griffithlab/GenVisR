#' Convert DNA character string
#' 
#' Convert a character string of nucleotides to amino acids or side chain class
#' @name lolliplot_DNAconv
#' @param x Character string of nucleotides to convert
#' @param to Character string specifying conversion to do, one of "codon", 
#' "residue", "sidechain"
#' @return Converted string of nucleotides as character vector

lolliplot_DNAconv <- function(x, to="residue")
{
    # check if given character string is a multiple of 3
    if(nchar(x)%%3 != 0)
    {
        memo <- paste0("Coding sequence retrieved for given ensembl transctipt",
                       ", is not a multiple of three. output may not be,",
                       " accurate!")
        warning(memo)
    }
    
    # split the character string into codons
    codon <- substring(x, seq(1,nchar(x), 3), seq(3, nchar(x), 3))
    if(toupper(to)=="CODON")
    {
        return(as.character(codon))
    }
    
    # convert the codons into amino acid residues
    residue <- sapply(codon, lolliplot_Codon2AA)
    if(toupper(to)=="RESIDUE")
    {
        return(as.character(residue))
    }
    
    # convert the residues into sidechain classifications
    sidechain <- sapply(residue, lolliplot_AA2sidechain)
    if(toupper(to)=="SIDECHAIN")
    {
        return(as.character(sidechain))
    }
    
    # return a warning if code gets this far
    memo <- paste0("did not recognize input to variable \"to\",",
                    " returning residue data")
    warning(memo)
    return(as.character(residue))
}