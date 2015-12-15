test_that("cnView_qual warns if input to x is not a data frame", {
    x <- data.frame(coordinate=c(100), chromosome=c(1), cn=0)
    x <- as.matrix(x)
    y <- NULL
    z <- NULL
    genome <- "hg19"
    expect_warning(cnView_qual(x, y, z, genome), "not appear to be")
})

test_that("cnView_qual checks for proper column names in x", {
    x <- data.frame(incorrect=c(100), chromosome=c(1), cn=0)
    y <- NULL
    z <- NULL
    genome <- "hg19"
    expect_error(cnView_qual(x, y, z, genome), "Did not detect") 
    
    x <- data.frame(coordinate=c(100), incorrect=c(1), cn=0)
    expect_error(cnView_qual(x, y, z, genome), "Did not detect") 
    
    x <- data.frame(coordinate=c(100), chromosome=c(1), incorrect=0)
    expect_error(cnView_qual(x, y, z, genome), "Did not detect") 
})

test_that("cnView_qual adds chr prefix if it is not detected", {
    x <- data.frame(coordinate=c(100), chromosome=c(1), cn=0)
    y <- NULL
    z <- NULL
    genome <- "hg19"
    expect_message(cnView_qual(x, y, z, genome),"the prefix chr")
    
    out <- cnView_qual(x, y, z, genome)
    expect_equivalent(as.character(out[[1]]$chromosome), "chr1")
})

test_that("cnView_qual outputs status message if chr prefix is detected", {
    x <- data.frame(coordinate=c(100), chromosome=c("chr1"), cn=0)
    y <- NULL
    z <- NULL
    genome <- "hg19"
    expect_message(cnView_qual(x, y, z, genome), "Detected chr")
})

test_that("cnView_qual outputs error message if an unrecognized prefix occurs in the chromosome column of x", {
    x <- data.frame(coordinate=c(100, 200), chromosome=c("chr1", 1), cn=c(0, 5))
    y <- NULL
    z <- NULL
    genome <- "hg19"
    expect_error(cnView_qual(x, y, z, genome), "mixed prefixes")    
})

test_that("cnView_qual checks that input to y is a data frame", {
    x <- data.frame(coordinate=c(100), chromosome=c(1), cn=0)
    y <- data.frame(chrom=c("chr1"), chromStart=c(1), chromEnd=c(10), name=c("test"), gieStain=c("gneg"))
    y <- as.matrix(y)
    z <- NULL
    genome <- "hg19"
    expect_message(cnView_qual(x, y, z, genome), "not appear to be")
})

test_that("cnView_qual checks the input to y contains correct column names", {
    x <- data.frame(coordinate=c(100), chromosome=c(1), cn=0)
    y <- data.frame(incorrect=c("chr1"), chromStart=c(1), chromEnd=c(10), name=c("test"), gieStain=c("gneg"))
    z <- NULL    
    expect_error(cnView_qual(x, y, z, genome), "not detect correct")
    
    y <- data.frame(chrom=c("chr1"), incorrect=c(1), chromEnd=c(10), name=c("test"), gieStain=c("gneg"))
    expect_error(cnView_qual(x, y, z, genome), "not detect correct")
    
    y <- data.frame(chrom=c("chr1"), chromStart=c(1), incorrect=c(10), name=c("test"), gieStain=c("gneg"))
    expect_error(cnView_qual(x, y, z, genome), "not detect correct")
    
    y <- data.frame(chrom=c("chr1"), chromStart=c(1), chromEnd=c(10), incorrect=c("test"), gieStain=c("gneg"))
    expect_error(cnView_qual(x, y, z, genome), "not detect correct")
    
    y <- data.frame(chrom=c("chr1"), chromStart=c(1), chromEnd=c(10), name=c("test"), incorrect=c("gneg"))
    expect_error(cnView_qual(x, y, z, genome), "not detect correct")
})

test_that("cnView_qual attempts to check if a genome supplied is not in UCSC terms", {
    x <- data.frame(coordinate=c(100), chromosome=c(1), cn=0)
    y <- NULL
    z <- NULL
    genome <- "grch38"
    expect_warning(cnView_qual(x, y, z, genome), "please specify a genome")
})

test_that("cnView_qual checks if input to z is not a data frame", {
    x <- data.frame(coordinate=c(100), chromosome=c(1), cn=0)
    y <- NULL
    z <- data.frame(chromosome=c(1), start=c(500), end=c(1000), segmean=c("4"))
    z <- as.matrix(z)
    genome <- "hg19"
    expect_message(cnView_qual(x, y, z, genome), "does not appear")    
})

test_that("cnView_qual checks if input to z does not contain correct column names", {
    x <- data.frame(coordinate=c(100), chromosome=c(1), cn=0)
    y <- NULL
    z <- data.frame(incorrect=c(1), start=c(500), end=c(1000), segmean=c("4"))
    genome <- "hg19"
    expect_error(cnView_qual(x, y, z, genome), "not detect correct")
    
    z <- data.frame(chromosome=c(1), incorrect=c(500), end=c(1000), segmean=c("4"))
    expect_error(cnView_qual(x, y, z, genome), "not detect correct")
    
    z <- data.frame(chromosome=c(1), start=c(500), incorrect=c(1000), segmean=c("4"))
    expect_error(cnView_qual(x, y, z, genome), "not detect correct")
    
    z <- data.frame(chromosome=c(1), start=c(500), end=c(1000), incorrect=c("4"))
    expect_error(cnView_qual(x, y, z, genome), "not detect correct")
})

test_that("cnView_qual adds chr prefix if it is not detected in z", {
    x <- data.frame(coordinate=c(100), chromosome=c("chr1"), cn=0)
    y <- NULL
    z <- data.frame(chromosome=c(1), start=c(500), end=c(1000), segmean=c("4"))
    genome <- "hg19"
    expect_message(cnView_qual(x, y, z, genome),"the prefix chr")
    
    out <- cnView_qual(x, y, z, genome)
    expect_equivalent(as.character(out[[1]]$chromosome), "chr1")
})

test_that("cnView_qual outputs status message if chr prefix is detected in z", {
    x <- data.frame(coordinate=c(100), chromosome=c("chr1"), cn=0)
    y <- NULL
    z <- data.frame(chromosome=c("chr1"), start=c(500), end=c(1000), segmean=c("4"))
    genome <- "hg19"
    expect_message(cnView_qual(x, y, z, genome), "chromosome column of z")
})

test_that("cnView_qual outputs error message if an unrecognized prefix occurs in the chromosome column of z", {
    x <- data.frame(coordinate=c(100), chromosome=c("chr1"), cn=c(0))
    y <- NULL
    z <- data.frame(chromosome=c("chr1", 1), start=c(500, 2000), end=c(1000, 3000), segmean=c("4", "4"))
    genome <- "hg19"
    expect_error(cnView_qual(x, y, z, genome), "mixed prefixes in the chromosome column of z")    
})