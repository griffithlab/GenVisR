test_that("lolliplot_constructGene properly nest's overlappping genomic features", {
    domain_data <- data.frame(desc=c("domainA", "domainB"), start=c(5, 10), end=c(20, 15))
    gene <- "test"
    length <- 50
    out <- lolliplot_constructGene(gene, domain_data, length)  
    
    expect_equivalent(out[out$Domain == 'domainA',]$nest, 1)
    expect_equivalent(out[out$Domain == 'domainB',]$nest, 2)
    expect_equivalent(out[out$Domain == 'test',]$nest, 1)
})

test_that("lolliplot_constructGene properly identifies if a domain exceeds the length of the protein", {
    domain_data <- data.frame(desc=c("domainA"), start=c(5), end=c(20))
    gene <- "test"
    length <- 15
    expect_warning(lolliplot_constructGene(gene, domain_data, length), "is exceeding the length")
})

test_that("lolliplot_constructGene properly identifies if a domain is before the start of the protein", {
    domain_data <- data.frame(desc=c("domainA"), start=c(0), end=c(10))
    gene <- "test"
    length <- 15
    expect_warning(lolliplot_constructGene(gene, domain_data, length), "less than the start of the protein")
})

test_that("lolliplot_constructGene identifies if a start position is greater than an end position", {
    domain_data <- data.frame(desc=c("domainA"), start=c(10), end=c(5))
    gene <- "test"
    length <- 20
    expect_warning(lolliplot_constructGene(gene, domain_data, length), "Found a start position greater")
})