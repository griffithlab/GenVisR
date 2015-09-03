test_that("waterfall_sampSort properly sorts samples based on a hierarchy", {
    x <- data.frame(sample=c('a', 'a', 'a', 'b', 'b', 'c'), gene=c('x', 'y', 'z', 'x', 'y', 'z'), trv_type=rep('missense', 6))
    out <- waterfall_sampSort(x)
    expec <- c('a', 'c', 'b')
    expect_equal(out, expec)
})