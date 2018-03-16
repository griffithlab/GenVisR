# packges needed
library(ggplot2)

# define the objects for testing
# note Lolliplot can run in two modes, we need to check both,
# we use the internal dataset PIK3CA for this


# mode 1, amino acid changes are present
keep <- c("Chromosome", "Start_Position", "End_Position", "Reference_Allele",
          "Tumor_Seq_Allele2", "Tumor_Sample_Barcode", "Gene", "Variant_Classification")
dfObject <- PIK3CA[,keep]
colnames(dfObject) <- c("chromosome", "start", "stop", "reference", "variant",
                        "sample", "gene", "consequence")
toLolliplot.mode1.out <- toLolliplot(dfObject, verbose=FALSE)

# mode 2, amino acid changes are not present
keep <- c("Chromosome", "Start_Position", "End_Position", "Reference_Allele",
          "Tumor_Seq_Allele2", "Tumor_Sample_Barcode", "Gene", "Variant_Classification",
          "Transcript_ID", "HGVSp")
dfObject <- PIK3CA[,keep]
colnames(dfObject) <- c("chromosome", "start", "stop", "reference", "variant",
                        "sample", "gene", "consequence", "transcript", "proteinCoord")
toLolliplot.mode2.out <- toLolliplot(dfObject, verbose=FALSE)

################################################################################
###### Test LolliplotData class and associated functions in constructor ########
################################################################################

context("LolliplotData Constructor")

