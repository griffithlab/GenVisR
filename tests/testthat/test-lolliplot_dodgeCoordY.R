test_that("lolliplot_dodgeCoordY succesfully dodges y coordinates as expected", {
    x <- data.frame(coord_x_dodge=c(1, 2))
    out <- lolliplot_dodgeCoordY(x, track="top")
    expect_equal(out[1], out[2])
    
    x <- data.frame(coord_x_dodge=c(1, 1))
    out <- lolliplot_dodgeCoordY(x, track="top")
    expect_true(out[1] != out[2])
    expect_true(out[1] > 0)
    
    out <- lolliplot_dodgeCoordY(x, track="bottom")
    expect_true(out[1] != out[2])
    expect_true(out[1] < 0)
})