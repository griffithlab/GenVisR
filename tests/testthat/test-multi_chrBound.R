test_that("multi_chrBound checks if input is empty", {
    x <- data.frame()
    expect_error(multi_chrBound(x), "input has 0 rows")
})

test_that("multi_chrBound correctly formats data", {
    x <- GenVisR::cytoGeno[GenVisR::cytoGeno$genome == 'hg19',]
    out <- multi_chrBound(x)
    expect_equal(nrow(out[out$chromosome == 'chr1',]), 2)
    
    expect_equal(out$start, out$end)
})