test_that("lolliplot_reduceLolli successfully reduces the number of lollis to the max limit",{
    x <- data.frame(mutation_coord=c(5, 5, 5, 2, 2), label=c("a", "b", "c", "d", "e"))
    max <- 2
    out <- lolliplot_reduceLolli(x, max=max)
    expect_equal(nrow(out[out$mutation_coord == 5,]), 2)
})

test_that("lolliplot_reduceLolli reduces nothing if max is set to NULL", {
    x <- data.frame(mutation_coord=c(5, 5, 5, 2, 2), label=c("a", "b", "c", "d", "e"))
    max <- NULL
    out <- lolliplot_reduceLolli(x, max=max)
    expect_equal(nrow(out), 5)
})