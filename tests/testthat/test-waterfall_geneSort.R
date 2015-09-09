test_that("waterfall_geneSort properly orders genes", {
    x <- data.frame(sample=rep('a', 10), gene=c(rep(c('b', 'c', 'd'), 3), 'e'), trv_type=rep('missense', 10))
    out <- waterfall_geneSort(x) 
    expec <- c('e', 'b', 'c', 'd')
    expect_equal(out, expec)
})