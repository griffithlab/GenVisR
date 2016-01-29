#' Convert AA to side chain classification
#' 
#' Given the 1 letter code an amino acid, return the side chian classification
#' @name lolliplot_AA2sidechain
#' @param x Character of length 1 giving the 1 letter amino acid code
#' @return Object of class character

lolliplot_AA2sidechain <- function(x)
{
    # Coerce all AA changes to uppercase and then apply switch statement
    x <- toupper(x)
    x <- switch(EXPR=x, "F"="Nonpolar", "L"="Nonpolar", "S"="Polar",
                "Y"="Polar", "C"="Polar", "W"="Nonpolar", "L"="Nonpolar",
                "P"="Nonpolar", "H"="Basic", "Q"="Polar", "R"="Basic",
                "I"="Nonpolar", "M"="Nonpolar", "T"="Polar", "N"="Polar",
                "K"="Basic", "S"="Polar", "R"= "Basic", "V"="Nonpolar",
                "A"="Nonpolar", "D"= "Acidic", "E"="Acidic", "G"="Polar")
    
    return(x)
}