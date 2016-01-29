test_that("waterfall_geneRecurCutoff properly detects if recurrence cutoff exceeds max recurrence", {
    x <- data.frame(sample=factor(seq(1:11)), gene=c(rep(c('b', 'c'), 5), 'd'), trv_type=rep('missense', 11))
    recurrence_cutoff <- .51
    expect_warning(waterfall_geneRecurCutoff(x, recurrence_cutoff), "cutoff specified exceeds")
})

test_that("waterfall_geneRecurCutoff properly removes entries not meeting a recurrence cutoff", {
    x <- data.frame(sample=factor(seq(1:11)), gene=c(rep(c('b', 'c'), 5), 'd'), trv_type=rep('missense', 11))
    recurrence_cutoff <- .51
    out <- suppressWarnings(waterfall_geneRecurCutoff(x, recurrence_cutoff))
    expect_equal(nrow(out), 10)    
})