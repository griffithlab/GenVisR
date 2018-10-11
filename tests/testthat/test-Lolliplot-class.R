# packges needed
library(ggplot2)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
BSgenome <- BSgenome.Hsapiens.UCSC.hg38

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

    skip_if_not(biomart_success, "mart recieved try-error")
    expect_true("biomaRt_ensembl_transcript_id" %in% colnames(annotateTranscript.mode1.out))
})

test_that("annotateTranscript does not attempt to annotate a transcript if one is already present", {

    skip_if_not(biomart_success, "mart recieved try-error")
    expect_true("ensembl_transcript_id" %in% colnames(annotateTranscript.mode2.out))
})

test_that("annotateTranscript works in verbose mode", {

    skip_if_not(biomart_success, "mart recieved try-error")
    expect_message(annotateTranscript(filterByGene.mode1.out, ensembl=retrieveMart.out, verbose=TRUE))
    expect_message(annotateTranscript(filterByGene.mode2.out, ensembl=retrieveMart.out, verbose=TRUE))
})

######################### test annotateProteinCoord ############################

if(biomart_success){

    annotateProteinCoord.mode1.out <- suppressWarnings(annotateProteinCoord(annotateTranscript.mode1.out,
                                                           ensembl=retrieveMart.out,
                                                           txdb=txdb,
                                                           BSgenome=BSgenome,
                                                           verbose=FALSE))

    annotateProteinCoord.mode2.out <- suppressWarnings(annotateProteinCoord(annotateTranscript.mode2.out,
                                                           ensembl=retrieveMart.out,
                                                           txdb=NULL,
                                                           BSgenome=NULL,
                                                           verbose=FALSE))
}

test_that("annotateProteinCoord adds protein coordinates if none are present", {

    skip_if_not(biomart_success, "mart recieved try-error")
    expect_true(class(annotateProteinCoord.mode1.out$proteinCoord) == "integer")

    # check 1 position of what the annotation file said to what genvisr said
    expect <- unique(annotateProteinCoord.mode2.out[annotateProteinCoord.mode2.out$ensembl_transcript_id == "ENST00000263967" &
                                                    annotateProteinCoord.mode2.out$chromosome == "chr3" &
                                                    annotateProteinCoord.mode2.out$start == "179218294" &
                                                    annotateProteinCoord.mode2.out$stop == "179218294",]$proteinCoord)
    actual <- unique(annotateProteinCoord.mode1.out[annotateProteinCoord.mode1.out$ensembl_transcript_id == "ENST00000263967" &
                                                    annotateProteinCoord.mode1.out$chromosome == "chr3" &
                                                    annotateProteinCoord.mode1.out$start == "179218294" &
                                                    annotateProteinCoord.mode1.out$stop == "179218294",]$proteinCoord)
    expect_equivalent(expect, actual)
})

test_that("annotateProteinCoord removes p. notation from protein coords if they are present", {

    skip_if_not(biomart_success, "mart recieved try-error")
    annotateTranscript.mode2.out <- annotateTranscript.mode2.out[1:5,]
    annotateProteinCoord.mode2.out <- suppressWarnings(annotateProteinCoord(annotateTranscript.mode2.out,
                                                           ensembl=retrieveMart.out,
                                                           txdb=NULL,
                                                           BSgenome=NULL,
                                                           verbose=FALSE))
    expect_true(is.numeric(annotateProteinCoord.mode1.out$proteinCoord))
})

test_that("annotateProteinCoord works in verbose mode", {
    skip_if_not(biomart_success, "mart recieved try-error")
    expect_message(suppressWarnings(annotateProteinCoord(annotateTranscript.mode1.out,
                                                         ensembl=retrieveMart.out,
                                                         txdb=txdb,
                                                         BSgenome=BSgenome,
                                                         verbose=TRUE)))
    expect_message(suppressWarnings(annotateProteinCoord(annotateTranscript.mode2.out,
                                                         ensembl=retrieveMart.out,
                                                         txdb=txdb,
                                                         BSgenome=BSgenome,
                                                         verbose=TRUE)))
})

##################### test filterByTranscript ##################################

if(biomart_success){
    filterByTranscript.mode1.out <- filterByTranscript(annotateProteinCoord.mode1.out,
                                                       transcript="ENST00000263967",
                                                       verbose=FALSE)
    filterByTranscript.mode2.out <- filterByTranscript(annotateProteinCoord.mode2.out,
                                                       transcript="ENST00000263967",
                                                       verbose=FALSE)
}

test_that("filterByTranscript leaves only one transcript", {
    skip_if_not(biomart_success, "mart recieved try-error")
    expect_equal(length(unique(filterByTranscript.mode1.out$ensembl_transcript_id)), 1)
    expect_equal(length(unique(filterByTranscript.mode2.out$ensembl_transcript_id)), 1)
})

test_that("filterByTranscript will choose a transcript if none is provided", {
    skip_if_not(biomart_success, "mart recieved try-error")
    filterByTranscript.mode1.out <- filterByTranscript(annotateProteinCoord.mode1.out,
                                                       transcript=NULL,
                                                       verbose=FALSE)
    filterByTranscript.mode2.out <- filterByTranscript(annotateProteinCoord.mode2.out,
                                                       transcript=NULL,
                                                       verbose=FALSE)
    expect_equal(length(unique(filterByTranscript.mode1.out$ensembl_transcript_id)), 1)
    expect_equal(length(unique(filterByTranscript.mode2.out$ensembl_transcript_id)), 1)
})

test_that("filterByTranscript warns of transcript input is not a character vector", {
    skip_if_not(biomart_success, "mart recieved try-error")
    expect_warning(filterByTranscript(annotateProteinCoord.mode1.out,
                                      transcript=factor("ENST00000263967"),
                                      verbose=FALSE))
})

test_that("filterByTranscript warns if transcript input has multiple elements", {
    skip_if_not(biomart_success, "mart recieved try-error")
    expect_warning(filterByTranscript(annotateProteinCoord.mode1.out,
                                      transcript=rep("ENST00000263967", 2),
                                      verbose=FALSE))
})

test_that("filterByTranscript warns if a given transcript is not found", {
    skip_if_not(biomart_success, "mart recieved try-error")
    expect_warning(filterByTranscript(annotateProteinCoord.mode1.out,
                                      transcript="unicorn",
                                      verbose=FALSE), "Did not find")
})

test_that("filterByTranscript works in verbose mode", {
    skip_if_not(biomart_success, "mart recieved try-error")
    expect_message(filterByTranscript(annotateProteinCoord.mode1.out,
                                      transcript="ENST00000263967",
                                      verbose=TRUE))
    expect_message(filterByTranscript(annotateProteinCoord.mode1.out,
                                      transcript=NULL,
                                      verbose=TRUE))
})

################ test calcMutFreq ##############################################


if(biomart_success){
    calcMutFreq.mode1.out <- calcMutFreq(filterByTranscript.mode1.out,
                                         verbose=FALSE)
    calcMutFreq.mode2.out <- calcMutFreq(filterByTranscript.mode2.out,
                                         verbose=FALSE)
}

test_that("calcMutFreq correctly calculates the mutation frequency", {
    skip_if_not(biomart_success, "mart recieved try-error")
    expected <- 13
    actual <- unique(calcMutFreq.mode1.out[calcMutFreq.mode1.out$chromosome == "chr3" &
                                           calcMutFreq.mode1.out$start == 179234297 &
                                           calcMutFreq.mode1.out$stop == 179234297 &
                                           calcMutFreq.mode1.out$reference == "A" &
                                           calcMutFreq.mode1.out$variant == "T",]$mutationFreq)
    expect_equal(expected, actual)

    expected <- 13
    actual <- unique(calcMutFreq.mode2.out[calcMutFreq.mode2.out$chromosome == "chr3" &
                                           calcMutFreq.mode2.out$start == 179234297 &
                                           calcMutFreq.mode2.out$stop == 179234297 &
                                           calcMutFreq.mode2.out$reference == "A" &
                                           calcMutFreq.mode2.out$variant == "T",]$mutationFreq)
    expect_equal(expected, actual)
})

test_that("calcMutFreq works in verbose mode", {
    skip_if_not(biomart_success, "mart recieved try-error")
    expect_message(calcMutFreq(filterByTranscript.mode1.out,
                               verbose=TRUE))
})

################# test constructTranscriptData fetches domains #################

if(biomart_success){
    constructTranscriptData.mode1.out <- constructTranscriptData(calcMutFreq.mode1.out,
                                                                 ensembl_mart=retrieveMart.out,
                                                                 verbose=FALSE)
    constructTranscriptData.mode2.out <- constructTranscriptData(calcMutFreq.mode2.out,
                                                                 ensembl_mart=retrieveMart.out,
                                                                 verbose=FALSE)
}

test_that("constructTranscriptData is able to retrieve protein domains", {

    skip_if_not(biomart_success, "mart recieved try-error")
    constructTranscriptData.mode1.out <- constructTranscriptData.mode1.out[constructTranscriptData.mode1.out$source == "domain",]
    expect_true(length(constructTranscriptData.mode1.out) >= 1)

    constructTranscriptData.mode2.out <- constructTranscriptData.mode2.out[constructTranscriptData.mode2.out$source == "domain",]
    expect_true(length(constructTranscriptData.mode2.out) >= 1)
})

test_that("constructTranscriptData is able to retrieve a protein length", {
    skip_if_not(biomart_success, "mart recieved try-error")
    constructTranscriptData.mode1.out <- constructTranscriptData.mode1.out[constructTranscriptData.mode1.out$source == "protein",]
    expect_true(constructTranscriptData.mode1.out$proteinStop >= 1)

    constructTranscriptData.mode2.out <- constructTranscriptData.mode2.out[constructTranscriptData.mode2.out$source == "protein",]
    expect_true(constructTranscriptData.mode2.out$proteinStop >= 1)
})

test_that("constructTranscriptData stops if more than one transcript is detected", {
    skip_if_not(biomart_success, "mart recieved try-error")
    calcMutFreq.mode1.out[1, "ensembl_transcript_id"] <- "ENST00000263968"
    expect_error(constructTranscriptData(calcMutFreq.mode1.out,
                                         ensembl_mart=retrieveMart.out,
                                         verbose=FALSE))
})

test_that("contructTranscriptData works in verbose mode", {
    skip_if_not(biomart_success, "mart recieved try-error")
    expect_message(constructTranscriptData(calcMutFreq.mode1.out,
                                           ensembl_mart=retrieveMart.out,
                                           verbose=TRUE))
})

############### test setTierTwo ################################################

if(biomart_success){
    setTierTwo.mode1.out <- setTierTwo(calcMutFreq.mode1.out, proteinData=constructTranscriptData.mode1.out,
                                       emphasize=NULL, verbose=FALSE)
    setTierTwo.mode2.out <- setTierTwo(calcMutFreq.mode2.out, proteinData=constructTranscriptData.mode2.out,
                                       emphasize=NULL, verbose=FALSE)
}

test_that("setTierTwo does not extend coordinates outside the protein", {
    skip_if_not(biomart_success, "mart recieved try-error")
    expect_true(max(setTierTwo.mode1.out$proteinCoordAdj) <= max(constructTranscriptData.mode1.out$proteinStop))
    expect_true(min(setTierTwo.mode1.out$proteinCoordAdj) <= min(constructTranscriptData.mode1.out$proteinStop))
})

test_that("setTierTwo increases the height of each different mutation that resides on the same amino acid", {
    skip_if_not(biomart_success, "mart recieved try-error")
    setTierTwo.mode1.out <- unique(setTierTwo.mode1.out[setTierTwo.mode1.out$proteinCoord=="345",
                                                        c("chromosome", "start", "stop", "reference",
                                                          "variant", "height", "mutationFreq")])
    setTierTwo.mode1.out <- setTierTwo.mode1.out[order(-setTierTwo.mode1.out$mutationFreq),]
    expect_true(head(setTierTwo.mode1.out$height, 1) < tail(setTierTwo.mode1.out$height, 1))
})

test_that("setTierTwo emphasize parameter works", {
    skip_if_not(biomart_success, "mart recieved try-error")

    setTierTwo.mode1.out <- setTierTwo(calcMutFreq.mode1.out, proteinData=constructTranscriptData.mode1.out,
                                       emphasize=c(70, 81), verbose=FALSE)
    setTierTwo.mode1.out <- setTierTwo.mode1.out[setTierTwo.mode1.out$tier == "second",]
    expect_true(nrow(setTierTwo.mode1.out) == 3)
})

test_that("setTierTwo works in verbose mode", {
    skip_if_not(biomart_success, "mart recieved try-error")
    expect_message(setTierTwo(calcMutFreq.mode1.out, proteinData=constructTranscriptData.mode1.out,
                              emphasize=NULL, verbose=TRUE))
})

################## test setTierOne #############################################

if(biomart_success){
    setTierOne.mode1.out <- setTierOne(setTierTwo.mode1.out, verbose=FALSE)
    setTierOne.mode2.out <- setTierOne(setTierTwo.mode2.out, verbose=FALSE)
}

test_that("first tier mutation heights increase as the mutation frequency increases", {
    skip_if_not(biomart_success, "mart recieved try-error")

    setTierOne.mode1.out <- setTierOne.mode1.out[setTierOne.mode1.out$tier == "first",]
    setTierOne.mode1.out <- setTierOne.mode1.out[order(setTierOne.mode1.out$mutationFreq),]
    expect_true(all(diff(setTierOne.mode1.out$height) >= 0))

    setTierOne.mode2.out <- setTierOne.mode2.out[setTierOne.mode2.out$tier == "first",]
    setTierOne.mode2.out <- setTierOne.mode2.out[order(setTierOne.mode2.out$mutationFreq),]
    expect_true(all(diff(setTierOne.mode2.out$height) >= 0))
})

test_that("setTierOne works in verbose mode", {
    skip_if_not(biomart_success, "mart recieved try-error")

    expect_message(setTierOne(setTierTwo.mode1.out, verbose=TRUE))
})

###################### test addLabel ###########################################

if(biomart_success){
    addLabel.mode1.out <- addLabel(setTierOne.mode1.out, verbose=FALSE)
    addLabel.mode2.out <- addLabel(setTierOne.mode2.out, verbose=FALSE)
}

test_that("addLabel adds only 1 label for each x/y coordinate", {
    skip_if_not(biomart_success, "mart recieved try-error")
    addLabel.mode1.out <- addLabel.mode1.out[addLabel.mode1.out$proteinCoordAdj == 118,]
    actual <- sum(is.na(addLabel.mode1.out$label))
    expected <- 4
    expect_equivalent(expected, actual)
})

test_that("addLabel is able to correctly add labels in the expected format", {
    skip_if_not(biomart_success, "mart recieved try-error")
    expected <- "p.Phe83Ser"
    actual <- as.character(addLabel.mode2.out[addLabel.mode2.out$proteinCoordAdj == 83]$label)
    expect_equivalent(expected, actual)

    expected <- "p.F83S"
    actual <- as.character(addLabel.mode1.out[addLabel.mode1.out$proteinCoordAdj == 83]$label)
    expect_equivalent(expected, actual)
})

test_that("addLabel works in verbose mode", {
    skip_if_not(biomart_success, "mart recieved try-error")
    expect_message(addLabel(setTierOne.mode1.out, verbose=TRUE))
})

###################### test setDomainHeights ###################################

if(biomart_success){
    setDomainHeights.mode1.out <- setDomainHeights(constructTranscriptData.mode1.out, verbose=FALSE)
    setDomainHeights.mode2.out <- setDomainHeights(constructTranscriptData.mode2.out, verbose=FALSE)
}

test_that("setDomainHeights domains are nested correctly", {
    
    skip_if_not(biomart_success, "mart recieved try-error")
    expected <- 8
    
    actual <- setDomainHeights.mode1.out[setDomainHeights.mode1.out$interproShortDesc == "PI3/4_kinase_cat_sf",]$nest
    expect_equivalent(expected, actual)

    expected <- 3
    actual <- setDomainHeights.mode1.out[setDomainHeights.mode1.out$interproShortDesc == "Ubiquitin-like_domsf",]$nest
    expect_equivalent(expected, actual)
})

test_that("setDomainHeights adjusts heights correctly based on the nest value", {

    skip_if_not(biomart_success, "mart recieved try-error")

    setDomainHeights.mode1.out <- setDomainHeights.mode1.out[order(setDomainHeights.mode1.out$nest),]
    expect_true(all(diff(setDomainHeights.mode1.out$maxHeight) >= 0))

    setDomainHeights.mode1.out <- setDomainHeights.mode1.out[order(setDomainHeights.mode1.out$nest),]
    expect_true(all(diff(setDomainHeights.mode1.out$minHeight) <= 0))
})

test_that("setDomainHeights works in verbose mode", {

    skip_if_not(biomart_success, "mart recieved try-error")

    expect_message(setDomainHeights(constructTranscriptData.mode1.out, verbose=TRUE))
})

############ test Lolliplotdata constructor ####################################

if(biomart_success){

    LolliplotData.mode1.out <- suppressWarnings(LolliplotData(dfObject.mode1, transcript="ENST00000263967",
                                                              species="hsapiens", host="www.ensembl.org",
                                                              txdb=txdb, BSgenome=BSgenome, emphasize=NULL,
                                                              verbose=FALSE))
    LolliplotData.mode2.out <- suppressWarnings(LolliplotData(dfObject.mode2, transcript="ENST00000263967",
                                                              species="hsapiens", host="www.ensembl.org",
                                                              txdb=NULL, BSgenome=NULL, emphasize=NULL,
                                                              verbose=FALSE))
}



test_that("LolliplotData outputs a S4 class object", {

    skip_if_not(biomart_success, "mart recieved try-error")

    expect_s4_class(LolliplotData.mode1.out, "LolliplotData")
    expect_s4_class(LolliplotData.mode2.out, "LolliplotData")
})

################################################################################
###### test the LolliplotPlots class and associated constructor functions ######
################################################################################

context("LolliplotPlots Constructor")

##################### test buildDensityPlot ####################################

test_that("buildDensityPlot constructs the expected plot", {

    skip_on_bioc()
    skip_if_not(biomart_success, "mart recieved try-error")

    buildDensityPlot.out <- buildDensityPlot(LolliplotData.mode1.out, plotALayers=NULL, verbose=FALSE)
    vdiffr::expect_doppelganger("Lolliplot Density Plot", grid::grid.draw(buildDensityPlot.out))
})

test_that("buildDensityPlot is able to add layers to the plot", {

     skip_on_bioc()
     skip_if_not(biomart_success, "mart recieved try-error")

     test_layer <- list(ggplot2::geom_text(aes(label="TEST LAYER", x=500, y=.025), colour="red", size=10))
     buildDensityPlot.out <- buildDensityPlot(LolliplotData.mode1.out, plotALayers=test_layer, verbose=FALSE)
     vdiffr::expect_doppelganger("Lolliplot Density Plot Layer", grid::grid.draw(buildDensityPlot.out))
})

test_that("buildDensityPlot works in verbose mode", {

    skip_if_not(biomart_success, "mart recieved try-error")

    expect_message(buildDensityPlot(LolliplotData.mode1.out, plotALayers=NULL, verbose=TRUE))
})

test_that("buildDensityPlot errors if plotALayers is not a list", {

    skip_if_not(biomart_success, "mart recieved try-error")

    test_layer <- ggplot2::geom_text(aes(label="TEST LAYER", x=500, y=.025), colour="red", size=10)
    expect_error(buildDensityPlot(LolliplotData.mode1.out, plotALayers=test_layer, verbose=FALSE))
})

####################### test buildLolliplot ####################################

test_that("buildLolliplot constructs a plot", {

    skip_on_bioc()
    skip_if_not(biomart_success, "mart recieved try-error")

    buildLolliplot.out <- buildLolliplot(LolliplotData.mode1.out, labelAA=FALSE,
                                         DomainPalette=NULL, MutationPalette=NULL,
                                         plotBLayers=NULL, verbose=FALSE)
    vdiffr::expect_doppelganger("Lolliplot Plot Base", grid::grid.draw(buildLolliplot.out))
})

test_that("buildLolliplot successfullly adds layers", {

    skip_on_bioc()
    skip_if_not(biomart_success, "mart recieved try-error")

    test_layer <- list(ggplot2::geom_text(aes(label="TEST LAYER", x=500, y=1), colour="red", size=10))

    buildLolliplot.out <- buildLolliplot(LolliplotData.mode1.out, labelAA=FALSE,
                                         DomainPalette=NULL, MutationPalette=NULL,
                                         plotBLayers=test_layer, verbose=FALSE)
    vdiffr::expect_doppelganger("Lolliplot Plot Base add layer", grid::grid.draw(buildLolliplot.out))
})

test_that("buildLolliplot successfullly adds labels", {

    skip_on_bioc()
    skip_if_not(biomart_success, "mart recieved try-error")

    buildLolliplot.out <- buildLolliplot(LolliplotData.mode1.out, labelAA=TRUE,
                                         DomainPalette=NULL, MutationPalette=NULL,
                                         plotBLayers=NULL, verbose=FALSE)
    vdiffr::expect_doppelganger("Lolliplot Plot Base add labels", grid::grid.draw(buildLolliplot.out))
})

test_that("buildLolliplot successfullly changes colors for domains", {

    skip_on_bioc()
    skip_if_not(biomart_success, "mart recieved try-error")

    domain_palette <- c("blue", "red", "green", "yellow", "orange", "purple",
                        "indianred", "salmon", "plum", "thistle", "seagreen1",
                        "slateblue1", "slategrey")

    buildLolliplot.out <- buildLolliplot(LolliplotData.mode1.out, labelAA=FALSE,
                                         DomainPalette=domain_palette, MutationPalette=NULL,
                                         plotBLayers=NULL, verbose=FALSE)
    vdiffr::expect_doppelganger("Lolliplot Plot Base add domain palette", grid::grid.draw(buildLolliplot.out))
})

test_that("buildLolliplot successfullly changes colors for mutations", {

    skip_on_bioc()
    skip_if_not(biomart_success, "mart recieved try-error")

    mutation_palette <- c("blue", "red", "green", "yellow")

    buildLolliplot.out <- buildLolliplot(LolliplotData.mode1.out, labelAA=FALSE,
                                         DomainPalette=NULL, MutationPalette=mutation_palette,
                                         plotBLayers=NULL, verbose=FALSE)
    vdiffr::expect_doppelganger("Lolliplot Plot Base add mutation palette", grid::grid.draw(buildLolliplot.out))
})

test_that("buildLolliplot works in verbose mode", {

    skip_if_not(biomart_success, "mart recieved try-error")

    expect_message(buildLolliplot(LolliplotData.mode1.out, labelAA=FALSE,
                                  DomainPalette=NULL, MutationPalette=NULL,
                                  plotBLayers=NULL, verbose=TRUE))
})

test_that("buildLolliplot warns if not enough colors are supplied to DomainPalette", {

    skip_if_not(biomart_success, "mart recieved try-error")

    domain_palette <- c("blue", "red", "green")

    expect_warning(buildLolliplot(LolliplotData.mode1.out, labelAA=FALSE,
                                  DomainPalette=domain_palette, MutationPalette=NULL,
                                  plotBLayers=NULL, verbose=FALSE))

})

test_that("buildLolliplot warns if not enough colors are supplied to MutationPalette", {

    skip_if_not(biomart_success, "mart recieved try-error")

    mutation_palette <- c("blue", "red")

    expect_warning(buildLolliplot(LolliplotData.mode1.out, labelAA=FALSE,
                                  DomainPalette=NULL, MutationPalette=mutation_palette,
                                  plotBLayers=NULL, verbose=FALSE))

})

test_that("buildLolliplot errors if plotBLayers is not a list", {

    skip_if_not(biomart_success, "mart recieved try-error")

    test_layer <- ggplot2::geom_text(aes(label="TEST LAYER", x=500, y=1), colour="red", size=10)

    expect_error(buildLolliplot(LolliplotData.mode1.out, labelAA=FALSE,
                                DomainPalette=NULL, MutationPalette=NULL,
                                plotBLayers=test_layer, verbose=FALSE))
})

###################### test LolliplotPlots constructor #########################

if(biomart_success){
    LolliplotPlots.out <- LolliplotPlots(LolliplotData.mode1.out, labelAA=FALSE,
                                         DomainPalette=NULL, MutationPalette=NULL,
                                         plotALayers=NULL, plotBLayers=NULL, verbose=FALSE)
}

test_that("LolliplotPlots constructor outputs a s4 class object", {

    skip_if_not(biomart_success, "mart recieved try-error")

    expect_s4_class(LolliplotPlots.out, "LolliplotPlots")

})

################################################################################
################## test Lolliplot constructor ##################################

context("Lolliplot Final Plot")

if(biomart_success){
    arrangeLolliplotPlot.out <- arrangeLolliplotPlot(LolliplotPlots.out, sectionHeights=NULL, verbose=FALSE)
}

test_that("arrangeLolliplotPlot constructs a plot", {

    skip_on_bioc()
    skip_if_not(biomart_success, "mart recieved try-error")

    vdiffr::expect_doppelganger("Final Lolliplot Base", grid::grid.draw(arrangeLolliplotPlot.out))

})

test_that("arrangeLolliplotPlot alters section heights", {

    skip_on_bioc()
    skip_if_not(biomart_success, "mart recieved try-error")

    arrangeLolliplotPlot.out <- arrangeLolliplotPlot(LolliplotPlots.out, sectionHeights=c(.1, .9), verbose=FALSE)

    vdiffr::expect_doppelganger("Final Lolliplot alter section height", grid::grid.draw(arrangeLolliplotPlot.out))
})

test_that("arrangeLolliplotPlot correctly warns if section heights does not match the number of plots", {

    skip_if_not(biomart_success, "mart recieved try-error")
    expect_warning(arrangeLolliplotPlot(LolliplotPlots.out, sectionHeights=c(1), verbose=FALSE))
})

test_that("arrangeLolliplotPlot correctly warns if section heights is not numeric", {

    skip_if_not(biomart_success, "mart recieved try-error")
    expect_warning(arrangeLolliplotPlot(LolliplotPlots.out, sectionHeights=c("A"), verbose=FALSE))
})

############################ test Lolliplot and accessors#######################

if(biomart_success){
    Lolliplot.out <- suppressWarnings(Lolliplot(dfObject.mode1, transcript="ENST00000263967",
                                                species="hsapiens", host="www.ensembl.org",
                                                txdb=txdb, BSgenome=BSgenome, emphasize=NULL,
                                                DomainPalette=NULL, MutationPalette=NULL,
                                                labelAA=TRUE, plotALayers=NULL, plotBLayers=NULL,
                                                sectionHeights=NULL, verbose=FALSE))
}

test_that("Lolliplot constructor outputs a S4 class object", {
    expect_s4_class(Lolliplot.out, "Lolliplot")
})

################################ draw plot #####################################

context("Lolliplot accessors")

test_that("drawPlot constructs a Lolliplot", {

    skip_on_bioc()
    skip_if_not(biomart_success, "mart recieved try-error")
    vdiffr::expect_doppelganger("drawPlot Lolliplot", drawPlot(Lolliplot.out))
})

################################# getGrob ######################################

test_that("getGrob outputs error if index is out of bounds", {

    expect_error(getGrob(Lolliplot.out, index=10))
})

test_that("getGrob successfully retrieves grob objects", {

    expect_s3_class(getGrob(Lolliplot.out, index=1), "gtable")
    expect_s3_class(getGrob(Lolliplot.out, index=2), "gtable")
    expect_s3_class(getGrob(Lolliplot.out, index=3), "gtable")
})

################################ getData #######################################

test_that("getData outputs error if no name or index is given", {

    expect_error(getData(Lolliplot.out))

})

test_that("getData outputs error if index exceeds the number of slots", {

    expect_error(getData(Lolliplot.out, index=10))

})

test_that("getData outputs error if supplied name is not a valid slot name", {

    expect_error(getData(Lolliplot.out, name="shouldNotexist"))

})

test_that("getData retrieves specified slot data correctly", {

    expect_s3_class(getData(Lolliplot.out, index=1), "data.table")
    expect_equivalent(getData(Lolliplot.out, name="primaryData"), getData(Lolliplot.out, index=1))

})
