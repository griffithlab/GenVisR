test_that("cnSpec plot data is formated as expected", {
    x <- data.frame(chromosome=c('chr1', 'chr1'), start=c(1000, 5000), end=c(5000, 10000), segmean=c(0, 4), sample=c('samp1', 'samp2'))
    y <- data.frame(chromosome=c('chr1', 'chr1'), start=c(0), end=c(10000))
    out <- cnSpec(x, y, out="data")
    expect_a <- data.frame(chromosome=as.factor(rep(1, 2)),
                         start=c(1000, 5000),
                         end=c(5000, 10000),
                         sample=as.factor(c('samp1', 'samp2')),
                         cn=c(0, 4))
    expect_b <- data.frame(chromosome=as.factor(rep(1, 4)),
                           start=rep(c(0, 10000), 2),
                           end=rep(c(0, 10000), 2),
                           sample=c('samp1', 'samp1', 'samp2', 'samp2'))

    expect_equivalent(out[[1]], expect_a)
    expect_equivalent(out[[2]], expect_b)
    
})