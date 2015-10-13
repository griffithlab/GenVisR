#' Convert Codon to AA
#' 
#' Convert a Codon to the appropriate amino acid
#' @name lolliplot_Codon2AA
#' @param x Character string of length 1 giving the DNA codon to convert
#' @return Character corresponding to the residue for the given codon

lolliplot_Codon2AA <- function(x)
{
    # Convert codons to single AA code
    x <- toupper(x)
    x <- switch(x, "TTT"="F", "TTC"="F", "TTA"="L", "TTG"="L", "CTT"="L",
                "CTC"="L", "CTA"="L", "CTG"="L", "ATT"="I", "ATC"="I",
                "ATA"="I", "ATG"="M", "GTT"="V", "GTC"="V", "GTA"="V",
                "GTG"="V", "TCT"="S", "TCC"="S", "TCA"="S", "TCG"="S",
                "CCT"="P", "CCC"="P", "CCA"="P", "CCG"="P", "ACT"="T",
                "ACC"="T", "ACA"="T", "ACG"="T", "GTC"="A", "GCC"="A",
                "GCA"="A", "GCG"="A", "TAT"="Y", "TAC"="Y", "TAA"="OCHRE",
                "TAG"="AMBER", "CAT"="H", "CAC"="H", "CAA"="Q", "CAG"="Q",
                "AAT"="N", "AAC"="N", "AAA"="K", "AAG"="K", "GAT"="D",
                "GAC"="D", "GAA"="E", "GAG"="E", "TGT"="C", "TGC"="C",
                "TGA"="OPAL", "TGG"="W", "CGT"="R", "CGC"="R", "CGA"="R",
                "CGG"="R", "AGT"="S", "AGC"="S", "AGA"="R", "AGG"="R",
                "GGT"="G", "GGC"="G", "GGA"="G", "GGG"="G")
    
    return(x)
}