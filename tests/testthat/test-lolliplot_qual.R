test_that("lolliplot_qual checks that input to x is a data frame", {
    x <- data.frame(transcript_name=c("ENST00000269305"), gene=c("TP53"), amino_acid_change=c("p.K440L"))
    x <- as.matrix(x)
    y <- NULL
    z <- NULL
    expect_message(lolliplot_qual(x, y, z), "not a data frame")
})

test_that("lolliplot_qual checks for correct column names in x", {
    x <- data.frame(incorrect=c("ENST00000269305"), gene=c("TP53"), amino_acid_change=c("p.K440L"))
    y <- NULL
    z <- NULL 
    expect_error(lolliplot_qual(x, y, z), "not detect correct columns")
    
    x <- data.frame(transcript_name=c("ENST00000269305"), incorrect=c("TP53"), amino_acid_change=c("p.K440L"))
    expect_error(lolliplot_qual(x, y, z), "not detect correct columns")
    
    x <- data.frame(transcript_name=c("ENST00000269305"), gene=c("TP53"), incorrect=c("p.K440L"))
    expect_error(lolliplot_qual(x, y, z), "not detect correct columns")
})

test_that("lolliplot_qual checks that only one transcript is present", {
    x <- data.frame(transcript_name=c("ENST00000269305", "ENST00000222222"), gene=c("TP53", "TP53"), amino_acid_change=c("p.K440L", "p.T318G"))
    y <- NULL
    z <- NULL 
    expect_error(lolliplot_qual(x, y, z), "more than 1 transcript")
})

test_that("lolliplot_qual checks that input to y is a data frame", {
    x <- data.frame(transcript_name=c("ENST00000269305"), gene=c("TP53"), amino_acid_change=c("p.K440L"))
    y <- data.frame(transcript_name=c("ENST00000269305"), amino_acid_change=c("p.K440L"))
    y <- as.matrix(y)
    z <- NULL  
    expect_message(lolliplot_qual(x, y, z), "not a data frame")
})

test_that("lolliplot_qual checks for correct column names in input to y", {
    x <- data.frame(transcript_name=c("ENST00000269305"), gene=c("TP53"), amino_acid_change=c("p.K440L"))
    y <- data.frame(incorrect=c("ENST00000269305"), amino_acid_change=c("p.K440L"))
    z <- NULL
    expect_error(lolliplot_qual(x, y, z), "not detect correct")
    
    x <- data.frame(transcript_name=c("ENST00000269305"), gene=c("TP53"), amino_acid_change=c("p.K440L"))
    y <- data.frame(transcript_name=c("ENST00000269305"), incorrect=c("p.K440L"))
    z <- NULL
    expect_error(lolliplot_qual(x, y, z), "not detect correct")
})

test_that("lolliplot_qual checks that input to z is a data frame", {
    x <- data.frame(transcript_name=c("ENST00000269305"), gene=c("TP53"), amino_acid_change=c("p.K440L"))
    y <- NULL
    z <- data.frame(description=c("test"), start=c(5), stop=c(3))
    z <- as.matrix(z)
    expect_warning(lolliplot_qual(x, y, z), "not a data frame")
})

test_that("lolliplot_qual checks for correct column names in input to y", {
    x <- data.frame(transcript_name=c("ENST00000269305"), gene=c("TP53"), amino_acid_change=c("p.K440L"))
    y <- NULL
    z <- data.frame(incorrect=c("test"), start=c(5), stop=c(3))
    expect_error(lolliplot_qual(x, y, z), "not detect correct")
    
    z <- data.frame(description=c("test"), incorrect=c(5), stop=c(3))
    expect_error(lolliplot_qual(x, y, z), "not detect correct")
    
    z <- data.frame(description=c("test"), start=c(5), incorrect=c(3))
    expect_error(lolliplot_qual(x, y, z), "not detect correct")
})

