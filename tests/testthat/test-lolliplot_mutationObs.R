test_that("lolliplot_mutationObs correctly removes input variants annotated as within an intron", {
    x <- data.frame(gene=c("TP53", "TP53"), amino_acid_change=c('e0+2', 'p.K101N'), transcript_name=c("ENST00000269305", "ENST00000269305"))
    expect_message(lolliplot_mutationObs(x,
                                         track="top",
                                         fill_value=NULL,
                                         label_column=NULL,
                                         rep.fact=5000,
                                         rep.dist.lmt=500,
                                         attr.fact=.1,
                                         adj.max=.1,
                                         adj.lmt=.5,
                                         iter.max=10000), "not within a residue")
})

test_that("lolliiplot_mutationObs correctly identifies if there are no mutations to plot", {
    x <- data.frame(gene=c("TP53"), amino_acid_change=c('e0+2'), transcript_name=c("ENST00000269305"))
    expect_error(lolliplot_mutationObs(x,
                                       track="top",
                                       fill_value=NULL,
                                       label_column=NULL,
                                       rep.fact=5000,
                                       rep.dist.lmt=500,
                                       attr.fact=.1,
                                       adj.max=.1,
                                       adj.lmt=.5,
                                       iter.max=10000), "Did not detect any residues")
})

test_that("lolliplot_mutationObs correctly extract coordinates from p. notation", {
    x <- data.frame(gene=c("TP53"), amino_acid_change=c('p.F426A'), transcript_name=c("ENST00000269305"))
    out <- lolliplot_mutationObs(x, track="top", fill_value=NULL, label_column=NULL,
                                 rep.fact=5000, rep.dist.lmt=500, attr.fact=.1, 
                                 adj.max=.1, adj.lmt=.5, iter.max=10000)
    expect_equal(out$mutation_coord, 426)
})

test_that("lolliplot_mutationObs correctly identifies amino acid changes in c. notation", {
    x <- data.frame(gene=c("TP53"), amino_acid_change=c('c.F4260A'), transcript_name=c("ENST00000269305"))
    expect_error(lolliplot_mutationObs(x, track="top", fill_value=NULL, label_column=NULL,
                                       rep.fact=5000, rep.dist.lmt=500, attr.fact=.1,
                                       adj.max=.1, adj.lmt=.5, iter.max=10000), "c. notation is not currently supported")
})

test_that("lolliplot_mutationObs correctly notifies if amino_acid_change notation is not recognized", {
    x <- data.frame(gene=c("TP53"), amino_acid_change=c('4260'), transcript_name=c("ENST00000269305"))
    expect_error(lolliplot_mutationObs(x, track="top", fill_value=NULL, label_column=NULL,
                                       rep.fact=5000, rep.dist.lmt=500, attr.fact=.1,
                                       adj.max=.1, adj.lmt=.5, iter.max=10000), "Could not determine notation type")
})