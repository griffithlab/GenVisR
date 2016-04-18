test_that("cnView_subsetChr correctly subset on a desired chromosome", {
    x <- data.frame(chromosome=c("chr1", "chr1", "chr2"), test=c("a", "b", "c"))
    out <- multi_subsetChr(x, "chr1")
    expect_equal(nrow(out), 2)
})

test_that("cnView_subsetChr returns all chromosomes if specified", {
    x <- data.frame(chromosome=c("chr1", "chr1", "chr2"), test=c("a", "b", "c"))
    out <- multi_subsetChr(x, "all")
    expect_equal(nrow(out), 3)    
})

test_that("cnView_subsetChr warns if chromosome column is not of class factor", {
    x <- data.frame(chromosome=c("chr1", "chr1", "chr2"), test=c("a", "b", "c"))
    x$chromosome <- as.character(x$chromosome)
    expect_warning(multi_subsetChr(x, "chr1"), "not a factor")
})

test_that("cnView_subsetChr produces an error if requested subset is not possible", {
    x <- data.frame(chromosome=c("chr1", "chr1", "chr2"), test=c("a", "b", "c"))
    expect_error(multi_subsetChr(x, "chr4000"), "does not match levels")
})