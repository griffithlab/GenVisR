test_that("waterfall_MGI2anno checks for proper column names",{
    # Check without label column
    x <- data.frame(incorrect='samp1', gene_name='egfr', trv_type='missense')
    label_col <- NULL
    expect_error(waterfall_MGI2anno(x, label_col), "Did not detect")
    
    # Check with label column
    x <- data.frame(sample='samp1', gene_name='egfr', trv_type='missense', incorrect='R504K')
    label_col <- 'label'
    expect_error(waterfall_MGI2anno(x, label_col), "Did not detect")   
})

test_that("waterfall_MGI2anno output format is correct", {
    # check without label column
    x <- data.frame(sample='samp1', gene_name='egfr', trv_type='missense')
    label_col <- NULL
    out <- waterfall_MGI2anno(x, label_col)
    expect_equal(sort(colnames(out)), sort(c('sample', 'gene', 'trv_type')))
    
    # check with label column
    x <- data.frame(sample='samp1', gene_name='egfr', trv_type='missense', label='R504K')
    label_col <- 'label'
    out <- waterfall_MGI2anno(x, label_col)
    expect_equal(sort(colnames(out)), sort(c('sample', 'gene', 'trv_type', label_col)))
    
    # Check class
    expect_is(out, "data.frame")
})