test_that("cnSpec_qual checks if genome is not of UCSC naming convention", {
    x <- data.frame(chromosome=c('chr1'), start=c(100), end=c(200), segmean=c(4), sample=c("samp1"))
    y <- NULL
    genome <- 'GRCH38'
    
    expect_warning(cnSpec_qual(x, y, genome), "UCSC terms")
})

test_that("cnSpec_qual checks if input to y is not a data frame", {
    x <- data.frame(chromosome=c('chr1'), start=c(100), end=c(200), segmean=c(4), sample=c("samp1"))
    y <- data.frame(chromosome=c('chr1'), start=c(0), end=c(1000000))
    y <- as.matrix(y)
    genome <- 'hg19'
    
    expect_message(cnSpec_qual(x, y, genome), "y is not a data frame")
})

test_that("cnSpec_qual checks for correct columns in input to y", {
    x <- data.frame(chromosome=c('chr1'), start=c(100), end=c(200), segmean=c(4), sample=c("samp1"))
    y <- data.frame(incorrect=c('chr1'), start=c(0), end=c(1000000))
    genome <- 'hg19'
    expect_error(cnSpec_qual(x, y, genome), "not detect correct")
    
    y <- data.frame(chromosome=c('chr1'), incorrect=c(0), end=c(1000000))
    expect_error(cnSpec_qual(x, y, genome), "not detect correct")
    
    y <- data.frame(chromosome=c('chr1'), start=c(0), incorrect=c(1000000))
    expect_error(cnSpec_qual(x, y, genome), "not detect correct")   
})

test_that("cnSpec_qual checks if input to x is not a data frame", {
    x <- data.frame(chromosome=c('chr1'), start=c(100), end=c(200), segmean=c(4), sample=c("samp1"))
    x <- as.matrix(x)
    y <- NULL
    genome <- 'hg19'
    
    expect_message(cnSpec_qual(x, y, genome), "not a data frame")
})

test_that("cnSpec_qual checks if input to x contains correct column names", {
    x <- data.frame(incorrect=c('chr1'), start=c(100), end=c(200), segmean=c(4), sample=c("samp1"))
    x <- as.matrix(x)
    y <- NULL
    genome <- 'hg19'
    expect_error(cnSpec_qual(x, y, genome), "not detect correct")
    
    x <- data.frame(chromosome=c('chr1'), incorrect=c(100), end=c(200), segmean=c(4), sample=c("samp1"))
    expect_error(cnSpec_qual(x, y, genome), "not detect correct")
    
    x <- data.frame(chromosome=c('chr1'), start=c(100), incorrect=c(200), segmean=c(4), sample=c("samp1"))
    expect_error(cnSpec_qual(x, y, genome), "not detect correct")
    
    x <- data.frame(chromosome=c('chr1'), start=c(100), end=c(200), incorrect=c(4), sample=c("samp1"))
    expect_error(cnSpec_qual(x, y, genome), "not detect correct")
    
    x <- data.frame(chromosome=c('chr1'), start=c(100), end=c(200), segmean=c(4), incorrect=c("samp1"))
    expect_error(cnSpec_qual(x, y, genome), "not detect correct")
})