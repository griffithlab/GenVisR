# packges needed
library(ggplot2)

# define the objects for testing
# note Lolliplot can run in two modes, we need to check both,
# we use the internal dataset PIK3CA for this

# mode 1, amino acid changes are present
keep <- c("Chromosome", "Start_Position", "End_Position", "Reference_Allele",
          "Tumor_Seq_Allele2", "Tumor_Sample_Barcode", "Gene", "Variant_Classification")
dfObject.mode1 <- PIK3CA[,keep]
colnames(dfObject.mode1) <- c("chromosome", "start", "stop", "reference", "variant",
                        "sample", "gene", "consequence")


# mode 2, amino acid changes are not present
keep <- c("Chromosome", "Start_Position", "End_Position", "Reference_Allele",
          "Tumor_Seq_Allele2", "Tumor_Sample_Barcode", "Gene", "Variant_Classification",
          "Transcript_ID", "HGVSp")
dfObject.mode2 <- PIK3CA[,keep]
colnames(dfObject.mode2) <- c("chromosome", "start", "stop", "reference", "variant",
                        "sample", "gene", "consequence", "transcript", "proteinCoord")

################################################################################
###### Test LolliplotData class and associated functions in constructor ########
################################################################################

context("LolliplotData Constructor")

###################### test toLolliplot ########################################

# Lolliplot can run in two modes depending on columns detected, test both
toLolliplot.mode1.out <- toLolliplot(dfObject.mode1, verbose=FALSE)
toLolliplot.mode2.out <- toLolliplot(dfObject.mode2, verbose=FALSE)

test_that("toLolliplot stops if a required column is not detected", {
    
    dfObject.mode1$sample <- NULL
    expect_error(toLolliplot(dfObject.mode1, verbose=FALSE))
})

test_that("toLolliplot works in verbose mode",{
    
    expect_message(toLolliplot(dfObject.mode1, verbose=TRUE))
})

##################### test retriveMart #########################################

# make sure there's an internet connection, if not skip tests
biomart_success <- try(biomaRt::useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org"))
biomart_success <- !class(biomart_success) == "try-error"
if(biomart_success){
    retrieveMart.out <- retrieveMart("hsapiens", host="www.ensembl.org", verbose=FALSE)
}

test_that("retrieveMart retrieves a valid mart for a given species", {
    
    skip_if_not(biomart_success, "mart recieved try-error, biomart is possibly down or there is no internet connection.")
    expect_true(class(retrieveMart.out) == "Mart")
})

test_that("retrieveMart errors if the species matches more than one dataset in the mart", {
    
    skip_if_not(biomart_success, "mart recieved try-error")
    expect_error(retrieveMart("m", host="www.ensembl.org", verbose=FALSE))
})

test_that("retreiveMart errors if the supplied species is not found", {
    
    skip_if_not(biomart_success, "mart recieved try-error")
    expect_error(retrieveMart("titan", host="www.ensembl.org", verbose=FALSE))
})

test_that("retreiveMart works in verbose mode", {
    
    expect_message(retrieveMart("hsapiens", host="www.ensembl.org", verbose=TRUE))
})

############### test annotateGene ##############################################
if(biomart_success){
    annotateGene.mode1.out <- annotateGene(toLolliplot.mode1.out,
                                           transcript="ENST00000263967",
                                           ensembl_mart=retrieveMart.out,
                                           verbose=FALSE) 
    annotateGene.mode2.out <- annotateGene(toLolliplot.mode2.out,
                                           transcript="ENST00000263967",
                                           ensembl_mart=retrieveMart.out,
                                           verbose=FALSE) 
}

test_that("annotateGene adds the ensembl, hgnc, and entrez gene identifiers for the supplied transcript", {
    
    skip_if_not(biomart_success, "mart recieved try-error")
    expect_equivalent(as.character(unique(annotateGene.mode1.out$biomaRt_ensembl_gene_id)), "ENSG00000121879")
    expect_equivalent(as.character(unique(annotateGene.mode1.out$biomaRt_hgnc_symbol)), "PIK3CA")
    expect_equivalent(as.character(unique(annotateGene.mode1.out$biomaRt_entrez_gene_id)), "5290")
})

test_that("annotateGene errors if no transcript is supplied", {
    
    expect_error(annotateGene(toLolliplot.mode1.out,
                              transcript=NULL,
                              ensembl_mart=retrieveMart.out,
                              verbose=FALSE))
})

test_that("annotateGene warns if transcript is not a character", {
    
    expect_warning(annotateGene(toLolliplot.mode1.out,
                                transcript=factor("ENST00000263967"),
                                ensembl_mart=retrieveMart.out,
                                verbose=FALSE))
})

test_that("annotateGene warns if more than 1 element is supplied to transcript", {
    
    expect_warning(annotateGene(toLolliplot.mode1.out,
                                transcript=rep("ENST00000263967", 2),
                                ensembl_mart=retrieveMart.out,
                                verbose=FALSE))
})

test_that("annotateGene errors if a valid ensembl transcript is not supplied", {
    
    expect_error(annotateGene(toLolliplot.mode1.out,
                              transcript="unicorn",
                              ensembl_mart=retrieveMart.out,
                              verbose=FALSE))
})

test_that("annotateGene works in verbose mode", {
    
    expect_message(annotateGene(toLolliplot.mode1.out,
                                transcript="ENST00000263967",
                                ensembl_mart=retrieveMart.out,
                                verbose=TRUE))
})

################ test filterByGene #############################################
if(biomart_success){
    filterByGene.mode1.out <- filterByGene(annotateGene.mode1.out, verbose=FALSE) 
    filterByGene.mode2.out <- filterByGene(annotateGene.mode2.out, verbose=FALSE) 
}

test_that("filterByGene correctly filters data by gene", {
    
    skip_if_not(biomart_success, "mart recieved try-error")
    annotateGene.mode1.out$gene <- as.character(annotateGene.mode1.out$gene)
    annotateGene.mode1.out$gene[1:length(annotateGene.mode1.out$gene)-1] <- c("newTranscript")
    filterByGene.mode1.out <- filterByGene(annotateGene.mode1.out, verbose=FALSE) 
    expect_equal(nrow(filterByGene.mode1.out), 1)
})

test_that("filterByGene can determine which value to filter by", {

    skip_if_not(biomart_success, "mart recieved try-error")
    annotateGene.mode1.out$gene <- as.character(annotateGene.mode1.out$gene)
    
    annotateGene.mode1.out$gene[1] <- as.character("PIK3CA")
    annotateGene.mode1.out$gene[2:nrow(annotateGene.mode1.out)] <- as.character("unicorn")
    filterByGene.mode1.out <- filterByGene(annotateGene.mode1.out, verbose=FALSE) 
    expect_equal(nrow(filterByGene.mode1.out), 1)
    
    annotateGene.mode1.out$gene[1] <- as.character("5290")
    annotateGene.mode1.out$gene[2:nrow(annotateGene.mode1.out)] <- as.character("unicorn")
    filterByGene.mode1.out <- filterByGene(annotateGene.mode1.out, verbose=FALSE) 
    expect_equal(nrow(filterByGene.mode1.out), 1)
    
    annotateGene.mode1.out$gene[1] <- as.character("ENSG00000121879")
    annotateGene.mode1.out$gene[2:nrow(annotateGene.mode1.out)] <- as.character("unicorn")
    filterByGene.mode1.out <- filterByGene(annotateGene.mode1.out, verbose=FALSE) 
    expect_equal(nrow(filterByGene.mode1.out), 1)
})

test_that("filterByGene errors if no gene is found to filter on", {
    
    skip_if_not(biomart_success, "mart recieved try-error")
    annotateGene.mode1.out$gene <- as.character(annotateGene.mode1.out$gene)
    annotateGene.mode1.out$gene[1:length(annotateGene.mode1.out$gene)] <- c("unicorn")
    expect_error(filterByGene(annotateGene.mode1.out, verbose=FALSE))
})

test_that("filterByGene prints status messages", {
    
    skip_if_not(biomart_success, "mart recieved try-error")
    expect_message(filterByGene(annotateGene.mode1.out, verbose=TRUE))
})

############### test annotateTranscript ########################################

if(biomart_success){
    annotateTranscript.mode1.out <- annotateTranscript(filterByGene.mode1.out, ensembl=retrieveMart.out, verbose=FALSE) 
    annotateTranscript.mode2.out <- annotateTranscript(filterByGene.mode2.out, ensembl=retrieveMart.out, verbose=FALSE) 
}

test_that("annotateTranscript annotates an ensembl transcript if none is present", {
    
    expect_true("biomaRt_ensembl_transcript_id" %in% colnames(annotateTranscript.mode1.out))
})

test_that("annotateTranscript does not attempt to annotate a transcript if one is already present", {
    
    expect_true("ensembl_transcript_id" %in% colnames(annotateTranscript.mode2.out))
})

test_that("annotateTranscript works in verbose mode", {
    
    expect_message(annotateTranscript(filterByGene.mode1.out, ensembl=retrieveMart.out, verbose=TRUE))
    expect_message(annotateTranscript(filterByGene.mode2.out, ensembl=retrieveMart.out, verbose=TRUE))
})

