test_that("TvTi_qual correctly identifies if input to x is not of proper class", {
    x <- data.frame(sample=c("a", "b"), reference=c("a", "t"), variant=c("g", "c"))
    x <- as.matrix(x)
    y <- NULL
    file_type="MGI"
    expect_warning(TvTi_qual(x, y=y, file_type=file_type), "not an object of class data frame")
})

test_that("TvTi_qual correctly identifies if duplicate rows are present in input to x", {
    x <- data.frame(sample=c("a", "a"), reference=c("a", "a"), variant=c("c", "c"))
    y <- NULL
    file_type="MGI"
    expect_warning(TvTi_qual(x, y=y, file_type=file_type), "Detected duplicate rows")
})

test_that("TvTi_qual correctly identifies if input to y is not of proper class", {
    x <- data.frame(sample=c("a", "b"), reference=c("a", "t"), variant=c("g", "c"))
    y <- c("A->C or T->G"=1/6, "A->G or T->C"=1/6, "A->T or T->A"=1/6, "G->A or C->T"=1/6, "G->C or C->G"=1/6, "G->T or C->A"=1/6)
    y <- as.matrix(y)
    file_type <- "MGI"
    expect_error(TvTi_qual(x, y=y, file_type=file_type), "not an object of class data frame")
})

test_that("TvTi_qual Checks for proper column names if input to y is a data frame", {
    x <- data.frame(sample=c("a", "b"), reference=c("a", "t"), variant=c("g", "c"))
    y <- data.frame(Prop=rep(1/6, 6), incorrect=rep("A->C or T->G", 6))
    file_type <- "MGI"
    expect_error(TvTi_qual(x, y=y, file_type=file_type), "Did not detect correct column names")
})

test_that("TvTi_qual checks that proportion values are of numeric type", {
    x <- data.frame(sample=c("a", "b"), reference=c("a", "t"), variant=c("g", "c"))
    y <- c("A->C or T->G"="1/6", "A->G or T->C"="1/6", "A->T or T->A"="1/6", "G->A or C->T"="1/6", "G->C or C->G"="1/6", "G->T or C->A"="1/6")
    file_type <- "MGI"
    expect_error(TvTi_qual(x, y=y, file_type=file_type), "not of type double or numeric")
})

test_that("TvTi_qual recognizes if input to y is not of proper class", {
    x <- data.frame(sample=c("a", "b"), reference=c("a", "t"), variant=c("g", "c"))
    y <- c("A->C or T->G"="1/6", "A->G or T->C"="1/6", "A->T or T->A"="1/6", "G->A or C->T"="1/6", "G->C or C->G"="1/6", "G->T or C->A"="1/6")
    y <- as.matrix(y)
    file_type <- "MGI"
    expect_error(TvTi_qual(x, y=y, file_type=file_type), "input to y is not")
})

test_that("TvTi_qual detects if proper column names are not supplied to x", {
    x <- data.frame(incorrect=c("a", "b"), reference=c("a", "t"), variant=c("g", "c"))
    y <- data.frame(Prop=rep(1/6, 6), trans_tranv=rep("A->C or T->G", 6))
    file_type <- "MGI"
    expect_error(TvTi_qual(x, y=y, file_type=file_type), "not find all columns")
    
    x <- data.frame(Tumor_Sample_Barcode=c("a", "b"), incorrect=c("a", "a"), Tumor_Seq_Allele1=c("a", "a"), Tumor_Seq_Allele2=c("a", "a"))
    file_type <- "MAF"
    expect_error(TvTi_qual(x, y=y, file_type=file_type), "not find all columns")
})

test_that("TvTi_qual detects unexpected nucleotide codes in input supplied to x", {
    x <- data.frame(sample=c("a", "b"), reference=c("a", "I"), variant=c("g", "c"))
    y <- data.frame(Prop=rep(1/6, 6), trans_tranv=rep("A->C or T->G", 6))
    file_type <- "MGI"   
    
    expect_error(TvTi_qual(x, y=y, file_type=file_type), "Unrecognized Base")
    
    x <- data.frame(sample=c("a", "b"), reference=c("a", "I"), variant=c("I", "c"))
    expect_error(TvTi_qual(x, y=y, file_type=file_type), "Unrecognized Base")
})

test_that("TvTi_qual checks for proper transition/transversion levels in input to y", {
    x <- data.frame(sample=c("a", "b"), reference=c("a", "t"), variant=c("g", "c"))
    y <- data.frame(Prop=rep(1/6, 6), trans_tranv=rep("A->C or T->G", 6))
    file_type <- "MGI"
    expect_error(TvTi_qual(x, y=y, file_type=file_type), "all combinations of transitions/transversions")
})

test_that("TvTi_qual checks that proportions supplied in y sum to one", {
    x <- data.frame(sample=c("a", "b"), reference=c("a", "t"), variant=c("g", "c"))
    y <- c("A->C or T->G"=1/5, "A->G or T->C"=1/6, "A->T or T->A"=1/6, "G->A or C->T"=1/6, "G->C or C->G"=1/6, "G->T or C->A"=1/6)
    file_type <- "MGI"
    expect_error(TvTi_qual(x, y=y, file_type=file_type), "should equal 1")
})
