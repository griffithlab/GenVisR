test_that("cnSpec plot data is formated as expected", {
    x <- data.frame(chromosome=c('chr1', 'chr1'), start=c(1000, 5000), end=c(5000, 10000), segmean=c(0, 4), sample=c('samp1', 'samp2'))
    y <- data.frame(chromosome=c('chr1', 'chr1'), start=c(0), end=c(10000))
    out <- cnSpec(x, y, out="data")
    expect <- data.frame(chromosome=as.factor(rep(1, 12)),
                         start=rep(c(1000, 5000, 0, 0, 10000, 10000), 2),
                         end=rep(c(5000, 10000, 0, 0, 10000, 10000), 2),
                         sample=c(rep('samp1', 6), rep('samp2', 6)),
                         cn=c(0, rep(NA, 6), 4, rep(NA, 4)))
    expect_equal(out[[1]], expect)
})