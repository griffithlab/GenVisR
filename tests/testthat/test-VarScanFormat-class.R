# Get the disk location for test files
testFileDir <- system.file("extdata", package="GenVisR")
lohTestFile <- Sys.glob(paste0(testFileDir, "/HCC1395.varscan.tsv"))
cnvTestFile <- Sys.glob(paste0(testFileDir, "/"))

# Define the object for testing
lohVarscanObject <- VarScanFormat(lohTestFile, varscanType = "LOH")
cnvVarscanObject <- VarScanFormat()
################################################################################
##################### test VarScanFormat class construction ####################
################################################################################

context("VarScanFormat")

test_that("VarScanFormat can construct object from a file path", {
    expect_s4_class(lohVarscanObject, "VarScanFormat")
})

test_that("VarScanFormat errors if both path and varscanData are NULL", {
    expect_error(VarScanFormat(path=NULL, varscanData=NULL))
})

test_that("VariantCallFormat warns if conversion to a data.table is required", {
    dataset <- data.frame(data.table::fread(lohTestFile))
    expect_message(VarScanFormat(varscanData=dataset, varscanType="LOH"))
})

test_that("VarScanFormat errors if specified varscanType is not supported", {
    dataset <- data.table::fread(lohTestFile)
    expect_error(VarScanFormat(path=NULL, varscanData=dataset, varscanType="CNA"))
})

test_that("VarScanFormat prints a message if the LOH VAF data is a percentage", {
    expect_warning(VarScanFormat(path=lohTestFile, varscanType="LOH"))
})

################################################################################
############################# test accessor methods ############################
################################################################################


################################################################################
########### test the getLohData method in lohSpec/combinedCnLohPlot ############
################################################################################

getLohData.out <- getLohData(lohVarscanObject, verbose=FALSE, lohSpec=TRUE, germline=FALSE)
test_that("method getLohData extracts the proper columns with heterozygous calls", {
    # test that it is a data.table
    expect_is(getLohData.out, "data.table")
    
    # test that it has the proper columns
    expectedCol <- c("chromosome", "position", "tumor_var_freq", 
                     "normal_var_freq", "sample")
    extractedCol <- colnames(getLohData.out)
    expect_true(all(extractedCol %in% expectedCol))
    
    # test that there are no coordinates with normal VAF less than 0.4 or greater than 0.6
    expect_true(length(which(getLohData.out$normal_var_freq < 0.4 | getLohData.out$normal_var_freq > 0.6))==0)
})




