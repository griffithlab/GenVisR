test_that("lolliplot_dodgeCoordX succesfully dodges coordinates as expected", {
    x <- data.frame(x=c(1, 2), y=c(1, 1))
    out <- lolliplot_dodgeCoordX(x)
    expect_true(out[1] < 1)
    expect_true(out[2] > 2)
    
    x <- data.frame(x=c(1, 1), y=c(1, 1))
    out <- lolliplot_dodgeCoordX(x)
    expect_equal(out[1], 1)
})

test_that("lolliplot_dodgeCoordX does not apply a force field if length of vector is 1", {
    x <- data.frame(x=c(1), y=c(1))
    out <- lolliplot_dodgeCoordX(x)
    expect_equal(out, 1)
})