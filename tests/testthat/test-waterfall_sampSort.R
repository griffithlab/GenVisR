test_that("waterfall_sampSort properly sorts samples based on a hierarchy", {
    x <- data.frame(sample=c('a', 'a', 'a', 'b', 'b', 'c'), gene=c('x', 'y', 'z', 'x', 'y', 'z'), trv_type=rep('missense', 6))
    out <- waterfall_sampSort(x)
    expec <- c('a', 'c', 'b')
    expect_equal(out, expec)
    
    # Test that convert to boolean properly works on data frame with one gene
    x <- data.frame(sample=c('a', 'b', 'c'), gene=c('x', 'x', 'x'), trv_type=rep('missense', 3))
    out <- waterfall_sampSort(x)
    expec <- c('a', 'b', 'c')
    expect_equal(out, expec)
    
    # Test that samples added in via the plotSamples parameter in waterfall are
    # put last
    x <- data.frame(sample=c('TESTA', 'A'), gene=c(NA, 'x'), trv_type=c(NA, 'missense'))
    out <- waterfall_sampSort(x)
    expec <- c('A', 'TESTA')
    #expect_equal(out, expec)
})

test_that("waterfall_sampSort removes samples in sampOrder not found in primary input", {
    x <- data.frame(sample=c('a', 'a', 'a', 'b', 'b', 'c'), gene=c('x', 'y', 'z', 'x', 'y', 'z'), trv_type=rep('missense', 6))
    sample_order <- c("a", "b", "c", "d")
    out <- suppressWarnings(waterfall_sampSort(x, sampOrder = sample_order))
    expec <- c('a', 'b', 'c')
    expect_equal(out, expec)
})