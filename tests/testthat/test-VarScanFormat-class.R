# Get the disk location for test files
testFileDir <- system.file("extdata", package="GenVisR")
testFile <- Sys.glob(paste0(testFileDir, "/HCC1395.varscan.tsv"))

# Define the object for testing
varscanObject <- VarScanFormat(testFile)

################################################################################
##################### test VarScanFormat class construction ####################
################################################################################

context("VarScanFormat")

test_that("VarScanFormat can construct object from a file path", {
    expect_s4_class(varscanObject, "VarScanFormat")
})

################################################################################
############################# test accessor methods ############################
################################################################################

test_that("accessor method getLohData extracts the proper columns", {
    expectedCol <- c("chromosome", "position", "tumor_var_freq", 
                     "normal_var_freq", "sample")
    extractedCol <- colnames(getLohData(varscanObject))
    expect_true(all(extractedCol %in% expectedCol))
})





