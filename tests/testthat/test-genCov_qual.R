suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19))
genome <- BSgenome.Hsapiens.UCSC.hg19
suppressPackageStartupMessages(library(TxDb.Hsapiens.UCSC.hg19.knownGene))
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
gr <- GRanges(seqnames=c("chr1"), ranges=IRanges(start=c(1), end=c(2)), strand=strand(c("+")))

test_that("genCov_qual correctly identifies if input to x is not of proper class", {
    x <- data.frame(end=c(1, 2), cov=c(50, 50))
    expect_warning(genCov_qual(x=x, txdb=txdb, gr=gr, genome=genome), "appear to be a list")
    
    x <- as.matrix(x)
    expect_warning(genCov_qual(x=x, txdb=txdb, gr=gr, genome=genome), "appear to be a list")    
})

test_that("genCov_qual correctly identifies if elements in the list x are not of propper class", {
    x <- data.frame(end=c(1, 2), cov=c(50, 50))
    x <- list('sample1'=as.matrix(x))
    
    expect_warning(genCov_qual(x=x, txdb=txdb, gr=gr, genome=genome), "not of type data frame")        
})

test_that("genCov_qual correctly identifies if a column name is missing in x", {
    x <- data.frame(incorrect=c(1, 2), cov=c(50, 50))
    x <- list("test"=x)
    
    expect_error(genCov_qual(x=x, txdb=txdb, gr=gr, genome=genome), "x are missing column names")        
})

test_that("genCov_qual correctly identifies if input to txdb in not a TxDb object", {
    x <- data.frame(end=c(1, 2), cov=c(50, 50))
    x <- list("test"=x)
    incorrect <- 'incorrect'
    
    expect_error(genCov_qual(x=x, txdb=incorrect, gr=gr, genome=genome), "txdb does not appear")        
})

test_that("genCov_qual correctly identifies if input to gr in not a Granges object", {
    x <- data.frame(end=c(1, 2), cov=c(50, 50))
    x <- list("test"=x)
    incorrect <- 'incorrect'
    
    expect_error(genCov_qual(x=x, txdb=txdb, gr=incorrect, genome=genome), "gr does not appear")        
})

test_that("genCov_qual correctly identifies if input to genome in not a BSgenome object", {
    x <- data.frame(end=c(1, 2), cov=c(50, 50))
    x <- list("test"=x)
    incorrect <- 'incorrect'
    
    expect_error(genCov_qual(x=x, txdb=txdb, gr=gr, genome=incorrect), "genome does not appear")        
})