#' Construct a lolliplot
#'
#' Given a data frame construct a plot displaying mutations on a transcript
#' framework.
#' @name lolliplot
#' @param x Object of class data frame with rows representing mutations. The
#' data frame must contain columns with the following names "transcript_name",
#' "gene", and "amino_acid_change". Values in the "transcript_name" column must
#' represent an ensembl transcript id and values in the "amino_acid_change"
#' column must be in p.notation (see details).
#' @param y Object of class data frame with rows representing mutations. The
#' data frame must contain columns with the following names "transcript_name"
#' and "amino_acid_change". Values in the "transcript_name" column must
#' represent an ensembl transcript id and values in the "amino_acid_change"
#' column must be in p. notation (optional, see details).
#' @param z Object of class data frame with rows representing regions of
#' interest. The data frame must contain columns with the following names
#' "description", "start", "stop" (optional see details).
#' @param fillCol Character string specifying the column name of the argument
#' supplied to parameter x on which to colour the lollis representing mutations
#' (see details).
#' @param labelCol Character string specifying the column name of the argument
#' supplied to parameter x from which to extract and display text corresponding
#' to mutations (see details).
#' @param txtAngle Integer specifying the angle of label text to be plotted if
#' an argument is supplied to the labelCol parameter.
#' @param txtSize Integer specifying the size of label text to be plotted if an
#' argument is supplied to the labelCol parameter.
#' @param pntSize Integer specifying the size of lolli points representing
#' mutations.
#' @param proteinColour Character string specifying the background colour of the
#' protein.
#' @param obsA.rep.fact Numeric value representing the repulsive factor for the
#' lollis plotted, which were derived from the argument supplied to parameter x
#' (see details and vignette).
#' @param obsA.rep.dist.lmt Numberic value representing the repulsive distance
#' limit for the lollis plotted, which were derived from the argument supplied
#' to parameter x (see details and vignette).
#' @param obsA.attr.fact Numeric value representing the attraction factor for
#' the lollis plotted, which were derived from the argument supplied to
#' parameter x (see details and vignette).
#' @param obsA.adj.max Numeric value representing the max position adjustment
#' for the lollis plotted, which were derived from the argument supplied to
#' parameter x (see details and vignette).
#' @param obsA.adj.lmt Numeric value representing the adjustment limit for the
#' lollis plotted, which were derived from the argument supplied to parameter x
#' (see details and vignette).
#' @param obsA.iter.max Integer representing the number of iterations of
#' position adjustments for the lollis plotted, which were derived from the 
#' argument supplied to parameter x (see details and vignette).
#' @param obsB.rep.fact Numeric value representing the repulsive factor for the
#' lollis plotted, which were derived from the argument supplied to parameter y
#' (see details and vignette).
#' @param obsB.rep.dist.lmt Numberic value representing the repulsive distance
#' limit for the lollis plotted, which were derived from the argument supplied
#' to parameter y (see details and vignette).
#' @param obsB.attr.fact Numeric value representing the attraction factor for
#' the lollis plotted, which were derived from the argument supplied to
#' parameter y (see details and vignette).
#' @param obsB.adj.max Numeric value representing the max position adjustment
#' for the lollis plotted, which were derived from the argument supplied to
#' parameter y (see details and vignette).
#' @param obsB.adj.lmt Numeric value representing the adjustment limit for the
#' lollis plotted, which were derived from the argument supplied to parameter y
#' (see details and vignette).
#' @param obsB.iter.max Integer representing the number of iterations of
#' position adjustments for the lollis plotted, which were derived from the 
#' argument supplied to parameter y (see details and vignette).
#' @param sideChain Boolean specifying if amino acid sidechain data should be
#' plotted in lieu of protein domains (see details).
#' @param species A valid species from which to retrieve protein domain and
#' sequence data for a given transcript (see details).
#' @param maxLolliStack Integer specifying the cutoff for the maximum number of
#' lollis allowed to be stacked at a single position.
#' @param plotLayer Valid ggplot2 layer to be added to the plot.
#' @param paletteA Character vector specifying colours for protein domains,
#' valid only if sideChain==FALSE.
#' @param paletteB Character vector specifying colours for lollis representing
#' mutations, valid only if argument is supplied to fillCol.
#' @param host Host to connect to for biomaRt queries (see details).
#' @param out Character vector specifying the the object to output, one of
#' "data", "grob", or "plot", defaults to "plot" (see returns).
#' @details lolliplot is a function designed to display mutation information in
#' the context of a protien identified by an ensembl transcript id. The
#' lolliplot function will query ensembl via biomart to retrieve sequence and
#' domain information in order to construct a representation of a protein and
#' therefore requires an internet connection. A value must be supplied to the
#' species parameter (defaults to hsapiens) in order for a successful biomart
#' query. Valid arguments to this field are those species with datasts available
#' via ensembl. please specify species in lowercase without a period
#' (i.e. hsapiens instead of H.sapiens), lolliplot will inform the user of
#' available species if input to the species parameter is not recognized.
#' Further lolliplot will build a protein framework based on sequence data
#' obtained from biomaRt, by default this will default to the latest ensembl
#' version. In order for the most accurate representation the annotation version
#' of the mutations given to lolliplot should match the annotation version used 
#' by biomaRt. The annotation version used by biomaRt can be changed via the 
#' host paramter (see vignette for more details).
#' 
#' lolliplot is capable of plotting two seperate sets of data on the protein
#' representation specified by parameters `x` and `y`, the data supplied to
#' these parameters will be plotted on the top and bottom of the protein
#' respectively. Note that input to these parameters is expected to correspond
#' to a single ensembl transcript and that values in the "amino_acid_change"
#' columns are required to be in p. notation (i.e. p.V600E). Further lolliplot
#' is able to plot custom domain annotation if supplied via the parameter `z`,
#' this will override domain information obtained from biomart.
#' 
#' lolliplot uses a forcefield model from the package FField to attract and 
#' repulse lollis. The parameters for this force field model are set to
#' reasonable defaults however may be adjusted via the obsA... and obsB...
#' family of parameters. Please see the package FField available on cran for
#' a description of these parameters. Note that the time to construct the
#' lolliplot will in large part depend on the number of mutations and the values
#' supplied to the forcefield parameters.
#' 
#' @examples
#' # Create input data
#' data <- brcaMAF[brcaMAF$Hugo_Symbol == 'TP53',c('Hugo_Symbol', 'amino_acid_change_WU')]
#' data <- as.data.frame(cbind(data, 'ENST00000269305'))
#' colnames(data) <- c('gene', 'amino_acid_change', 'transcript_name')
#'
#' # Call lolliplot
#' lolliplot(data)
#' @return One of the following, a list of dataframes containing data to be
#' plotted, a grob object, or a plot.
#' @export

lolliplot <- function(x, y=NULL, z=NULL, fillCol=NULL, labelCol=NULL,
                      txtAngle=45, txtSize=5, pntSize=4,
                      proteinColour='#999999', obsA.rep.fact=5000,
                      obsA.rep.dist.lmt=500, obsA.attr.fact=.1, obsA.adj.max=.1,
                      obsA.adj.lmt=.5, obsA.iter.max=50000, obsB.rep.fact=5000,
                      obsB.rep.dist.lmt=500, obsB.attr.fact=.1, obsB.adj.max=.1,
                      obsB.adj.lmt=.5, obsB.iter.max=50000,
                      sideChain=FALSE, species="hsapiens",
                      maxLolliStack=NULL, plotLayer=NULL, paletteA=NULL,
                      paletteB=NULL, host="www.ensembl.org", out="plot")
{
    # Perform quality check
    input <- lolliplot_qual(x, y, z)
    x <- input[[1]]
    y <- input[[2]]
    z <- input[[3]]

    # extract transcript id and subset data y on that id if it exists
    transcriptID <- as.character(x$transcript_name[1])
    if(!is.null(y))
    {
        y <- y[y$transcript_name == transcriptID,]
    }

    # extract HUGO gene name
    gene <- as.character(x$gene[1])

    # Obtain length of protein
    result <- lolliplot_transcriptID2codingSeq(transcriptID,
                                               species=species,
                                               host=host)
    codingSeq <- result$coding
    cdsLen <- result$cds_length
        
    
    # Get the sequence length in AA, perform quality checks along the way
    residueSeq <- lolliplot_DNAconv(codingSeq, to="residue")    
    # If it is requested grab the sidechain information and bind to residues
    if(sideChain==TRUE)
    {
        sidechain <- lolliplot_DNAconv(codingSeq, to="sidechain")
        AAsequence <- cbind(sidechain, residueSeq)
        AAsequence <- as.data.frame(AAsequence)
        AAsequence$coord <- seq(from=1, to=nrow(AAsequence))
    } else {
        AAsequence <- NULL
    }
    
    # if there are any stop codons remove them as they are not considered part
    # of the protein
    if(any(residueSeq %in% c("OPAL", "OCHRE", "AMBER")))
    {
        stopRes <- c("OPAL", "OCHRE", "AMBER")
        residueSeq <- residueSeq[-which(residueSeq %in% stopRes)]
        if(!is.null(AAsequence))
        {
            AAsequence <- AAsequence[-which(AAsequence$residueSeq %in% stopRes),]
        }
    }
    
    # grab the length of the protein in Amino Acids
    proteinLength <- length(residueSeq)    
    
    # if z is specified plot that instead of fetching the domain information
    if(!is.null(z))
    {
        geneData <- lolliplot_constructGene(gene, z, proteinLength)
    } else {
        # extract protien domain data
        protein_domain <- lolliplot_fetchDomain(transcriptID,
                                                species=species,
                                                host=host)
        
        # construct gene from data collected
        geneData <- lolliplot_constructGene(gene, protein_domain, proteinLength)
    }
 
    # construct data frame of observed mutations for top track
    observed_mutation <- lolliplot_mutationObs(x, 'top', fillCol, labelCol,
                                               obsA.rep.fact, obsA.rep.dist.lmt,
                                               obsA.attr.fact, obsA.adj.max,
                                               obsA.adj.lmt, obsA.iter.max)
    observed_mutation <- lolliplot_reduceLolli(observed_mutation,
                                               maxLolliStack)

    # construct data frame of observed mutations for bottom track
    if(!is.null(y))
    {
        observed_mutation2 <- lolliplot_mutationObs(y, 'bottom', fillCol,
                                                    labelCol, obsB.rep.fact,
                                                    obsB.rep.dist.lmt,
                                                    obsB.attr.fact,
                                                    obsB.adj.max, obsB.adj.lmt,
                                                    obsB.iter.max)
        observed_mutation2 <- lolliplot_reduceLolli(observed_mutation2,
                                                    maxLolliStack)
    } else {
        observed_mutation2 <- NULL
    }
    
    # construct the lolliplot
    plot <- lolliplot_buildMain(geneData, length, observed_mutation,
                                observed_mutation2,fillCol, labelCol,
                                txtAngle, txtSize, pntSize,
                                proteinColour, AAsequence,
                                plot_sidechain=sideChain, layers=plotLayer,
                                paletteA=paletteA, paletteB=paletteB)

    # Decide what to output
    dataOut <- list("gene"=geneData,
                    "observed_mutation"=observed_mutation,
                    "observed_mutation2"=observed_mutation2)
    output <- multi_selectOut(data=dataOut, plot=plot, out=out)
    return(output)
}
