test_that("waterfall_hierarchyTRV leaves only one entry for a gene/mutation", {
    # Check MGI file
    x <- data.frame(sample=rep(c('a', 'b'), 4), gene=rep(c('c', 'd'), 4), trv_type=rep(c('nonsense', 'missense', 'rna', 'silent'), 2))
    file_type <- 'MGI'
    out <- waterfall_hierarchyTRV(x, file_type)
    expect_equal(as.character(out[out$sample == 'a' & out$gene == 'c', 'trv_type']), 'nonsense')
    expect_equal(as.character(out[out$sample == 'b' & out$gene == 'd', 'trv_type']), 'missense')
    
    # Check MAF file
    x <- data.frame(sample=rep(c('a', 'b'), 4), gene=rep(c('c', 'd'), 4), trv_type=rep(c('Splice_Site', 'IGR', 'Silent', 'Nonsense_Mutation'), 2))
    file_type <- 'MAF'
    out <- waterfall_hierarchyTRV(x, file_type)
    expect_equal(as.character(out[out$sample == 'a' & out$gene == 'c', 'trv_type']), 'Splice_Site')
    expect_equal(as.character(out[out$sample == 'b' & out$gene == 'd', 'trv_type']), 'Nonsense_Mutation')    
})

test_that("waterfall_hierarchyTRV recognizes an invalid mutation type", {
    # Check for MGI
    x <- data.frame(sample=rep(c('a', 'b'), 4), gene=rep(c('c', 'd'), 4), trv_type=rep(c('missense', 'rna', 'silent', 'invalid'), 2))
    file_type <- 'MGI'
    expect_error(waterfall_hierarchyTRV(x, file_type), "invalid mutation type")
    
    # Check for MAF
    x <- data.frame(sample=rep(c('a', 'b'), 4), gene=rep(c('c', 'd'), 4), trv_type=rep(c('Splice_Site', 'IGR', 'Silent', 'invalid'), 2))
    file_type <- 'MAF'
    expect_error(waterfall_hierarchyTRV(x, file_type), "invalid mutation type")
})