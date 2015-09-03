test_that("waterfall_rmvSilent removes silent mutations", {
    x <- data.frame(sample=rep('a', 3), gene=rep('b', 3), trv_type=c('Silent', 'siLenT', 'missense'))
    out <- droplevels(waterfall_rmvSilent(x))    
    expec <- data.frame(sample=rep('a', 3), gene=c(NA, NA, 'b'), trv_type=c(NA, NA, 'missense'))
    expect_equal(out, expec)
})