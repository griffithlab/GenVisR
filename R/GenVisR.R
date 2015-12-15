#' GenVisR
#'
#' @name GenVisR
#' @docType package
NULL

#' Truncated BRCA MAF file
#'
#' A data set containing 50 samples corresponding to "Breast invasive carcinoma"
#' originating from the TCGA project in .maf format (version 2.4):
#' https://wiki.nci.nih.gov/display/TCGA/TCGA+MAF+Files#TCGAMAFFiles-BRCA:Breastinvasivecarcinoma, /dccfiles_prod/tcgafiles/distro_ftpusers/anonymous/tumor/brca/gsc/genome.wustl.edu/illuminaga_dnaseq/mutations/genome.wustl.edu_BRCA.IlluminaGA_DNASeq.Level_2.5.3.0/genome.wustl.edu_BRCA.IlluminaGA_DNASeq.Level_2.5.3.0.somatic.maf
#'
#' @docType data
#' @keywords datasets
#' @name brcaMAF
#' @usage data(brcaMAF)
#' @format a data frame with 2773 observations and 55 variables
#' @return Object of class data drame
NULL

#' Cytogenetic banding dataset
#'
#' A data set containing cytogenetic band information for all chromosomes in the
#' following genomes "hg38", "hg19", "mm10", "mm9", "rn5", obtained from the
#' UCSC sql database at genome-mysql.cse.ucsc.edu.
#'
#' @docType data
#' @keywords datasets
#' @name cytoGeno
#' @usage data(cytoGeno)
#' @format a data frame with 3207 observations and 6 variables
#' @return Object of class data frame
NULL

#' hg19 chromosome boundaries
#'
#' A data set containg chromosome boundaries corresponding to hg19.
#'
#' @docType data
#' @keywords datasets
#' @name hg19chr
#' @usage data(hg19chr)
#' @format a data frame with 24 observations and 3 variables
#' @return Object of class data frame
NULL


#' CN data for Luc2
#'
#' Downsampled copy number calls for Luc2 chromosome 14 originating from
#' cnv-hmm. Govindan et al. Cell. 2012, PMID:22980976
#'
#' @docType data
#' @keywords datasets
#' @name Luc2CNraw
#' @usage data(Luc2CNraw)
#' @format a data frame with 10000 observations and 5 variables
#' @return Object of class data frame
NULL

#' Truncated CN segments
#'
#' A data set in long format containing Copy Number segments for 4 samples
#' corresponding to "lung cancer" from Govindan et al. Cell. 2012, PMID:22980976
#'
#' @docType data
#' @keywords datasets
#' @name LucCNseg
#' @usage data(LucCNseg)
#' @format a data frame with 3336 observations and 6 variables
#' @return Object of class data frame
NULL

#' Coverage for PTEN
#'
#' A data set containing coverage data for PTEN from Juric et al. Nature 2015,
#' PMID:25409150
#'
#' @docType data
#' @keywords datasets
#' @name ptenCov
#' @usage data(ptenCov)
#' @format a data frame with 10734 observations and 4 variables
#' @return Object of class data frame
NULL

#' Identity snps
#' 
#' A data set containing locations of 24 identity snps originating from:
#' Pengelly et al. Genome Med. 2013, PMID 24070238
#'  
#' @docType data
#' @keywords datasets
#' @name SNPloci
#' @usage data(SNPloci)
#' @format a data frame with 24 observations and 3 variables
#' @return Object of class data frame
NULL

#' Tumor BAM
#' 
#' A data set containing aligned reads intersecting 24 identity snp locations
#' from GenVisR::SNPloci. Aligned reads are downsampled and originate from
#' tumor tissue corresponding to the HCC1395 breast cancer cell line.
#' @docType data
#' @keywords datasets
#' @name HCC1395_T
#' @usage data(HCC1395_T)
#' @format a named list with 24 elements and 6 variables
#' @return Object of class list
NULL

#' Normal BAM
#' 
#' A data set containing aligned reads intersecting 24 identity snp locations
#' from GenVisR::SNPloci. Aligned reads are downsamples and originate from
#' normal tissue corresponding to the HCC1395 breast cancer cell line.
#' @docType data
#' @keywords datasets
#' @name HCC1395_N
#' @usage data(HCC1395_N)
#' @format a named list with 24 elements and 6 variables
#' @return Object of class list
NULL
