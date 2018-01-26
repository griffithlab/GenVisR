# Get the disk location for test files
testFileDir <- system.file("extdata", package="GenVisR")
testFile <- Sys.glob(paste0(testFileDir, "/HCC1395.varscan.tsv"))

# Define the object for testing
varscanObject <- VarScanFormat(testFile, varscanType = "LOH")

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



################################################################################
########### test the getLohData method in lohSpec/combinedCnLohPlot ############
################################################################################

getLohData.out <- getLohData(varscanObject, verbose=FALSE, lohSpec=TRUE, germline=FALSE)
test_that("accessor method getLohData extracts the proper columns with heterozygous calls", {
    # test that it is a data.table
    expect_is(getLohData.out, "data.table")
    
    # test that it has the proper columns
    expectedCol <- c("chromosome", "position", "tumor_var_freq", 
                     "normal_var_freq", "sample")
    extractedCol <- colnames(getLohData.out)
    expect_true(all(extractedCol %in% expectedCol))
    
    # test that there are no coordinates with normal VAF less than 0.4 or greater than 0.6
     
})





