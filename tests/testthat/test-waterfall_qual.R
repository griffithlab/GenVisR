test_that("waterfall_qual returns a list with three elements", {
    x <- data.frame(sample='samp1', gene_name='egfr', trv_type='missense')
    y <- NULL
    z <- NULL
    file_type <- 'MGI'
    label_col <- NULL
    out <- waterfall_qual(x, y, z, file_type, label_col)
    
    expect_equal(length(out), 3)
    expect_is(out, 'list')
})

test_that("waterfall_qual checks input to file_type is valid", {
    x <- data.frame(sample='samp1', gene_name='egfr', trv_type='missense')
    y <- NULL
    z <- NULL
    file_type <- 'TGI'
    label_col <- NULL
    
    expect_error(waterfall_qual(x, y, z, file_type, label_col), "Unrecognized")
})

test_that("waterfall_qual verifies input is of proper class", {
    
    # Check x
    x <- c(sample='samp1', gene_name='egfr', trv_type='missense')
    y <- NULL
    z <- NULL
    file_type <- 'MGI'
    label_col <- NULL
    
    expect_error(waterfall_qual(x, y, z, file_type, label_col),
                 'Did not detect')
    
    # Check y
    x <- data.frame(sample='samp1', gene_name='egfr', trv_type='missense')
    y <- c(sample='samp1', variable='gender', value='male') 
    z <- NULL
    
    expect_error(waterfall_qual(x, y, z, file_type, label_col),
                 'Did not detect')
    
    # Check z
    x <- data.frame(sample='samp1', gene_name='egfr', trv_type='missense')
    y <- NULL
    z <- c(sample='samp1', mut_burden='5')
    expect_error(waterfall_qual(x, y, z, file_type, label_col),
                 'Did not detect')
})

test_that("waterfall_qual verifies correct columns are in input", {
    
    # Check x
    x <- data.frame(incorrect='samp1', gene_name='egfr', trv_type='missense')
    y <- NULL
    z <- NULL
    file_type <- 'MGI'
    label_col <- NULL
    
    expect_error(waterfall_qual(x, y, z, file_type, label_col),
                 'Did not detect')
    
    x <- data.frame(sample='samp1', incorrect='egfr', trv_type='missense')
    y <- NULL
    z <- NULL
    file_type <- 'MGI'
    label_col <- NULL
    
    expect_error(waterfall_qual(x, y, z, file_type, label_col),
                 'Did not detect')
    
    x <- data.frame(sample='samp1', gene_name='egfr', incorrect='missense')
    y <- NULL
    z <- NULL
    file_type <- 'MGI'
    label_col <- NULL
    
    expect_error(waterfall_qual(x, y, z, file_type, label_col),
                 'Did not detect')
    
    # Check y
    x <- data.frame(sample='samp1', gene_name='egfr', trv_type='missense')
    y <- data.frame(incorrect='samp1', variable='gender', value='male')
    z <- NULL
    file_type <- 'MGI'
    label_col <- NULL
    
    expect_error(waterfall_qual(x, y, z, file_type, label_col),
                 'Did not detect')
    
    x <- data.frame(sample='samp1', gene_name='egfr', trv_type='missense')
    y <- data.frame(sample='samp1', incorrect='gender', value='male')
    z <- NULL
    file_type <- 'MGI'
    label_col <- NULL
    
    expect_error(waterfall_qual(x, y, z, file_type, label_col),
                 'Did not detect')
    
    x <- data.frame(sample='samp1', gene_name='egfr', trv_type='missense')
    y <- data.frame(sample='samp1', variable='gender', incorrect='male')
    z <- NULL
    file_type <- 'MGI'
    label_col <- NULL
    
    expect_error(waterfall_qual(x, y, z, file_type, label_col),
                 'Did not detect')
    
    # Check z
    x <- data.frame(sample='samp1', gene_name='egfr', trv_type='missense')
    y <- NULL
    z <- c(incorrect='samp1', mut_burden='5')
    file_type <- 'MGI'
    label_col <- NULL
    
    expect_error(waterfall_qual(x, y, z, file_type, label_col),
                 'Did not detect')
    
    x <- data.frame(sample='samp1', gene_name='egfr', trv_type='missense')
    y <- NULL
    z <- c(sample='samp1', incorrect='5')
    file_type <- 'MGI'
    label_col <- NULL
    
    expect_error(waterfall_qual(x, y, z, file_type, label_col),
                 'Did not detect')    
})




