test_that("waterfall_calcMutFreq correctly counts mutations", {
    x <- data.frame(sample=c('a', 'a', 'a', 'b', 'b'), gene=rep('b', 5), trv_type=rep('silent', 5))
    out <- waterfall_calcMutFreq(x)
    expect_equal(out[out$sample == 'a' & out$trv_type == 'Synonymous', 'mutation_total'], 3)
    expect_equal(out[out$sample == 'b' & out$trv_type == 'Synonymous', 'mutation_total'], 2)
    expect_equal(out[out$sample == 'a' & out$trv_type == 'Non Synonymous', 'mutation_total'], 0)
    expect_equal(out[out$sample == 'b' & out$trv_type == 'Non Synonymous', 'mutation_total'], 0)
})

test_that("waterfall_calcMutFreq correctly identifies Synonymous, Non Synonymous mutations", {
    x <- data.frame(sample=c('a', 'a', 'a'), gene=rep('b', 3), trv_type=c('not silent', 'Silent', 'siLenT'))
    out <- waterfall_calcMutFreq(x)
    expect_equal(out[out$sample == 'a' & out$trv_type == 'Synonymous', 'mutation_total'], 2)
    expect_equal(out[out$sample == 'a' & out$trv_type == 'Non Synonymous', 'mutation_total'], 1)   
})