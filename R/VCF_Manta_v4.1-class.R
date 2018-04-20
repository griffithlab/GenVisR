################################################################################
##################### Public/Private Class Definitions #########################

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Public Class !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#

#' Class VCF_Manta_v4.1
#' 
#' An S4 class to represent data in vcf version 4.1 format, inherits from the 
#' VCF_Virtual class
#' @name VCF_Manta_v4.1
#' @rdname VCF_Manta_v4.1-class
#' @slot header data.table object containing header information
#' @slot meta data.table object containing meta information lines
#' @slot vcfHeader data.table object containing header for vcf data
#' @slot vcfData data.table object containing vcf data lines
#' @slot sample data.table object containing sample information
#' @include VCF_Virtual-class.R
#' @import methods
setClass("VCF_Manta_v4.1",
         contains="VCF_Virtual",
         validity=function(object){
             cnames <- c("chromosome", "position", "chromosome2", "position2", "direction",
                         "REF", "ALT", "svtype", "total_read_support", "FILTER", "sample", 
                         "ID", "INFO", "FORMAT", "tumorSample", "paired")
             
             ## Check the columns
             sampleCol <- which(!colnames(object@vcfData) %in% cnames)
             if (length(sampleCol) > 0) {
                 memo <- paste0("Columns in the input data.table are missing. Required ",
                                "columns are: chromosome, position, chromosome2, position2, direction,", 
                                "REF, ALT, svtype, total_read_support, FILTER, sample ", 
                                "ID, INFO, FORMAT, tumorSample, paired")
                 message(memo)
             }
             return(TRUE)
             
         })

#' Constuctor for the VCF_Manta_v4.1 sub-class
#' 
#' @name VCF_Manta_v4.1
#' @rdname VCF_Manta_v4.1-class
#' @param vcfData data.table object containing a VCF file conforming to the 
#' version 4.1 specifications
#' @param vcfHeader Object of class list containing character vectors for vcf
#' header information
#' @param paired Boolean object specifying if the svCaller was ran in paired mode
#' @param tumorColumn String specifying the name of the sample column with read support information
#' @importFrom data.table data.table
VCF_Manta_v4.1 <- function(vcfData, vcfHeader, paired, tumorColumn) {

    ## Set the data descriptions for the object
    if (length(vcfHeader)==0) {
        finalDescription <- data.table::data.table()
    } else {
        description <- lapply(vcfHeader, function(x){
            descriptionFieldIndex <- which(grepl("Description", x))
            x <- x[descriptionFieldIndex]
            
            split1 <- unlist(strsplit(unlist(strsplit(x, ",")), "="))
            id <- split1[grep("<ID", split1)+1]
            description <- gsub(">", "", split1[grep("Description", split1)+1])
            description <- gsub("\"", "", description, fixed=TRUE)
            
            x <- paste("ID=", id, "|", "Description=", description, sep="")
            return(x)
        })
        description <- unique(unlist(description))
        
        # convert these results to a data.table after splitting into two columns
        d <- unlist(strsplit(description, split="|", fixed=TRUE))
        id <- d[grep("ID=", d)]
        id <- unlist(strsplit(id, split="ID="))[
            grep("[A-Z]", unlist(strsplit(id, split="ID=")))]
        description <- d[grep("Description=", d)]
        description <- unlist(strsplit(description, split="Description="))[
            grep("[A-Z]", unlist(strsplit(description, split="Description=")))]
        
        finalDescription <- data.table(name=id, description=description)
    }
    
    ## Get the samples
    sample <- data.table(sample=unique(vcfData$sample))
    
    ## Assign the column names
    cnames <- c("chromosome", "position", "ID", "REF", "ALT", "QUAL",
                "FILTER", "INFO", "FORMAT")
    colnames(vcfData)[1:9] <- cnames

    ## Check if the sample is paired
    if (paired == TRUE) {
        ## Check if the sv data is paired based on the input files/dataset
        cnames <- c("chromosome", "position", "ID", "REF", "ALT", "QUAL",
                    "FILTER", "INFO", "FORMAT", "sample")
        cols <- length(colnames(vcfData[,-which(colnames(vcfData) %in% cnames), with=FALSE]))
        if (cols !=2) {
            memo <- paste("Are you sure the data is paired. There are", cols, 
                          "columns with sample read support data",
                          "in the input sv data when there should be 2.")
            stop(memo)
        }
        
        ## Check if the tumorColumn variable actually specifies a sample column to use
        num <- which(colnames(vcfData[,tumorColumn,with=FALSE]) %in% 
                         colnames(vcfData[,-which(colnames(vcfData) %in% cnames), with=FALSE]))
        if (length(num) != 1) {
            memo <- paste("The column designated as the tumor sample does not", 
                          "correspond to sample read support. The valid values to use",
                          "for the tumorColumn variable are:", 
                          paste(which(!colnames(vcfData) %in% cnames), collapse=" or "))
            stop(memo)
        }
        
        ## Check if the tumorColumn variable is NULL
        if (is.null(tumorColumn)) {
            memo <- paste0("Input was designated as paired but the tumor/diseased sample ",
                           "was not designated. If the samples are paired, please ",
                           "use the tumorColumn variable to identify which column in the vcf datasets has the ",
                           "read support for calls in the tumor sample.")
        }
        if (is.null(tumorColumn) == FALSE) {
            cnames <- c("chromosome", "position", "ID", "REF", "ALT", "QUAL",
                        "FILTER", "INFO", "FORMAT", colnames(vcfData)[tumorColumn], "sample")
            vcfData <- vcfData[,which(colnames(vcfData) %in% cnames), with=FALSE]
            colnames(vcfData)[10] <- "tumorSample"
        }
    }
    if (paired == FALSE) {
        ## Check if the sv data is not paired based on the input files/dataset
        cnames <- c("chromosome", "position", "ID", "REF", "ALT", "QUAL",
                    "FILTER", "INFO", "FORMAT", "sample")
        cols <- length(colnames(vcfData[,-which(colnames(vcfData) %in% cnames), with=FALSE]))
        if (cols !=1) {
            memo <- paste("Are you sure the data is NOT paired. There are", cols, 
                          "columns with sample read support data",
                          "in the input sv data when there should be 1.")
            stop(memo)
        }
        
        if (!is.null(tumorColumn)) {
            memo <- paste("The sv caller output was designated as not paired. ", 
                          "But a value was assigned to the tumorColumn variable.", 
                          "This value will be ignored.")
        }
        tumorColumn <- which(!colnames(vcfData) %in% cnames)
        colnames(vcfData)[tumorColumn] <- "tumorSample"
    }
    vcfData$paired <- paired
    
    ## Get the structural variant type from the INFO column
    ## Get the chr and position of the second breakpoint
    ## Get the read support
    ## Get the direction of the translocation
    temp <- suppressWarnings(data.table::rbindlist(apply(vcfData, 1, function(x, paired) {
        ## SV type
        svtype <- unlist(strsplit(as.character(x["INFO"]), split=";"))[
            grep("SVTYPE", unlist(strsplit(as.character(x["INFO"]), split=";")))]
        svtype <- strsplit(svtype, "SVTYPE=")[[1]][2]
        x$svtype <- svtype[1]
        
        ## 2nd breakpoint and get the BND direction
        if (svtype == "BND" | svtype == "TRA") {
            alt <- strsplit(as.character(x$ALT), "[[:punct:]]")[[1]]
            chromosome2 <- alt[2]
            position2 <- alt[3]
            dir <- strsplit(as.character(x$ALT), split="")[[1]]
            if (dir[length(dir)] == "]"){
                final <- "N]P]"
            }
            if (dir[length(dir)] == "[") { 
                final <- "N[P["
            }
            if (dir[1] == "]") {
                final <- "]P]N"
            }
            if (dir[1] == "[") { 
                final <- "[P[N"
            }
            dir <- final      
        } 
        else if (svtype == "INV" | svtype == "DEL" | svtype == "DUP" | svtype == "INS") {
            chromosome2 <- x$chromosome
            tmp <- unlist(strsplit(as.character(x$INFO), split=";")[[1]])
            tmp <- tmp[grep("END=", tmp, fixed=TRUE)][1]
            position2 <- strsplit(tmp, split="END=")[[1]][2]
            dir <- svtype
        }
        x$chromosome2 <- as.character(chromosome2)
        x$position2 <- as.numeric(position2)
        x$direction <- dir
        
        ## Get the read support 
        format <- as.character(x$FORMAT)
        ## Read support for the sample
        a <- as.character(x$tumorSample)
        if (format == "PR") {
            ## Get read support for the sample
            sample1_temp <- strsplit(a, split=":")[[1]]
            sample1_chr1_PR <- as.numeric(strsplit(sample1_temp[1], split=",")[[1]][1])
            sample1_chr1_SR <- 0
            sample1_chr2_PR <- as.numeric(strsplit(sample1_temp[1], split=",")[[1]][2])
            sample1_chr2_SR <- 0
        }
        if (format == "PR:SR") {
            ## Get read support for the first sample
            sample1_temp <- strsplit(a, split=":")[[1]]
            sample1_chr1_PR <- as.numeric(strsplit(sample1_temp[1], split=",")[[1]][1])
            sample1_chr1_SR <- as.numeric(strsplit(sample1_temp[2], split=",")[[1]][1])
            sample1_chr2_PR <- as.numeric(strsplit(sample1_temp[1], split=",")[[1]][2])
            sample1_chr2_SR <- as.numeric(strsplit(sample1_temp[2], split=",")[[1]][2])
        }
        ## Get the total read support
        x$total_read_support <- sum(sample1_chr1_PR, sample1_chr1_SR,
                                          sample1_chr2_PR, sample1_chr2_SR)
        x <- as.data.table(t(cbind(x)))
        return(x)
    }, paired=paired)))
    vcfData <- temp[,c("chromosome", "position", "chromosome2", "position2", "direction",
            "REF", "ALT", "svtype", "total_read_support", "FILTER", "sample", 
            "ID", "INFO", "FORMAT", "tumorSample", "paired")]
    
    ## Get the ID for each SV call
    vcfData$ID <- data.table::rbindlist(lapply(strsplit(as.character(vcfData$ID), split=":"), function(x) {
        y <- data.table(paste(x[1], x[2], x[3], x[4], sep = "_"))
        return(y)
    }))
    
    ## Get the svtype
    svType <- data.table(unique(vcfData$svtype))
    colnames(svType) <- "svtype"
    
    ## Initialize the object
    new("VCF_Manta_v4.1", description=finalDescription, sample=sample, 
        vcfData=vcfData, svType=svType)
}