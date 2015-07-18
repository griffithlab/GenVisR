# GenVisR
Intuitively visualizing and interpreting data from high-throughput genomic technologies continues to be challenging. GenVisr (Genomic Visualizations in R) attempts to alleviate this burden by providing highly customizable publication-quality graphics focused primarily on a cohort level (i.e., multiple samples/patients). GGgenome attempts to maintain a high degree of flexibility while leveraging the abilities of ggplot2 and bioconductor to achieve this goal.

![alt text](https://github.com/griffithlab/GenVisR/blob/master/vignettes/GenVisR_intro_files/figure-latex/unnamed-chunk-3-1.pdf)
![alt text](https://github.com/griffithlab/GenVisR/blob/master/vignettes/GenVisR_intro_files/figure-latex/unnamed-chunk-7-1.pdf)
![alt text](https://github.com/griffithlab/GenVisR/blob/master/vignettes/GenVisR_intro_files/figure-latex/unnamed-chunk-14-1.pdf)
![alt text](https://github.com/griffithlab/GenVisR/blob/master/vignettes/GenVisR_intro_files/figure-latex/unnamed-chunk-16-1.pdf)


## Installation Instructions

Note requires R >= 3.2.1, see Description file for a full list of dependencies

* Install devtools
```R
install.packages("devtools")
```

* Install Bioconductor dependencies
```R
source("http://bioconductor.org/biocLite.R")
biocLite(c("GenomicRanges", "biomaRt", "GenomicFeatures", "UniProt.ws"))
```

* It is suggested but not required to install knitr and rmarkdown for vignette creation
```R
install.packages(c("rmarkdown", "knitr"))
```

* Install GenVisR
```R
devtools::install_github("griffithlab/GenVisR")
```
