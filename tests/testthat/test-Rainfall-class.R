# packges needed
library(ggplot2)

# get the disk location for test files
testFileDir <- system.file("extdata", package="GenVisR")
testFile <- Sys.glob(paste0(testFileDir, "/*.vep"))

# define the object for testing
vepObject <- VEP(testFile)
#toRainFall.out <- toRainfall(vepObject, BSgenome=NULL, verbose=F)

################################################################################
###### Test Rainfall class and associated functions in constructor #############
################################################################################