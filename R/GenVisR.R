#' GenVisR
#'
#' @name GenVisR
#' @docType package
NULL

#' Truncated BRCA MAF file
#'
#' A data set containing 50 samples corresponding to "Breast invasive carcinoma" originating from the TCGA project in .maf format (version 2.4):
#' https://wiki.nci.nih.gov/display/TCGA/TCGA+MAF+Files#TCGAMAFFiles-BRCA:Breastinvasivecarcinoma, /dccfiles_prod/tcgafiles/distro_ftpusers/anonymous/tumor/brca/gsc/genome.wustl.edu/illuminaga_dnaseq/mutations/genome.wustl.edu_BRCA.IlluminaGA_DNASeq.Level_2.5.3.0/genome.wustl.edu_BRCA.IlluminaGA_DNASeq.Level_2.5.3.0.somatic.maf
#'
#' @docType data
#' @keywords datasets
#' @name brcaMAF
#' @usage data(brcaMAF)
#' @format a data frame with 2773 observations and 55 variables
NULL

#' cytogenetic band information for select genomes
#'
#' A data set containing cytogenetic band information for all chromosomes in the following genomes "hg38", "hg19", "mm10", "mm9", "rn5", obtained from the UCSC sql database at genome-mysql.cse.ucsc.edu
#'
#' @docType data
#' @keywords datasets
#' @name cytoGeno
#' @usage data(cytoGeno)
#' @format a data frame with 3207 observations and 6 variables
NULL

#' hg19 chromosome boundaries
#'
#' A data set containg chromosome boundaries corresponding to hg19 in base pairs
#'
#' @docType data
#' @keywords datasets
#' @name hg19chr
#' @usage data(hg19chr)
#' @format a data frame with 24 observations and 3 variables
NULL


#' CN data for Luc2
#'
#' Downsampled copy number calls for Luc2 chromosome 14 originating from cnv-hmm "Genomic Landscape of Non-Small Cell Lung Cancer in Smokers and Never-Smokers" PMID:22980976
#'
#' @docType data
#' @keywords datasets
#' @name Luc2CNraw
#' @usage data(Luc2CNraw)
#' @format a data frame with 10000 observations and 5 variables
NULL

#' Truncated CN segments
#'
#' A data set in long format containing Copy Number segments for 4 samples corresponding to "lung cancer" from "Genomic Landscape of Non-Small Cell Lung Cancer in Smokers and Never-Smokers" PMID:22980976
#'
#' @docType data
#' @keywords datasets
#' @name LucCNseg
#' @usage data(LucCNseg)
#' @format a data frame with 3336 observations and 6 variables
NULL

#' Coverage for the pten gene
#'
#' A data set containing coverage data for pten obtained using "bedtools multicov" from "Convergent loss of PTEN leads to clinical resistance to a PI3K-alpha inhibitor" PMID:25409150
#'
#' @docType data
#' @keywords datasets
#' @name ptenCOV
#' @usage data(ptenCOV)
#' @format a data frame with 107338 observations and 4 variables
NULL
