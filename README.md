# GenVisR
Intuitively visualizing and interpreting data from high-throughput genomic technologies continues to be challenging. GenVisr (Genomic Visualizations in R) attempts to alleviate this burden by providing highly customizable publication-quality graphics focused primarily on a cohort level (i.e., multiple samples/patients). GenVisr attempts to maintain a high degree of flexibility while leveraging the abilities of ggplot2 and bioconductor to achieve this goal.

## Installation Instructions

Note requires R >= 3.2.1, see Description file for a full list of dependencies

* Install devtools
```R
install.packages("devtools")
```

* Install Bioconductor dependencies
```R
source("http://bioconductor.org/biocLite.R")
biocLite(c("AnnotationDbi", "biomaRt", "Biostrings", "GenomicFeatures", "GenomicRanges", "Rsamtools"))
```

* It is suggested but not required to install knitr and rmarkdown for vignette creation and testthat for unit tests
```R
install.packages(c("rmarkdown", "knitr", "testthat"))
```

* Install GenVisR
```R
devtools::install_github("griffithlab/GenVisR")
```

## Tips for developers
* It is recommended to use Rstudio for development, Clone the repo at (https://github.com/griffithlab/GenVisR.git) and open the "GenVisR.Rproj" folder in Rstudio.
* It is necessary to manual install all packages (since you are cloning and building, not installing), all packages within the imports, Depends, and Suggests flags in the "DESCRIPTION" file (GenVisR/DESCRIPTION) should be installed! The r-package devtools should be installed as well!
* Build and Reload the package at this step (Rstudio shorcut: Shift+Cmd+B)
* GenVisR uses roxygen2 for documentation, to update the Reference Manuscript edit the roxygen2 flags (i.e. @param, @details, etc.) in the code. These are within the R function files in GenVisR/R/*.R, Then run devtools::document() or build and reload the package.
* To update the vignette edit the R markdown file in the vignettes subdirectory (GenVisR/vignettes/GenVisR_intro.Rmd) and press the knit button in Rstudio
* Code for GenVisR is within the R subdirectory (~/GenVisR/R/), make changes there. Be sure to build and reload the package from Rstudio to test changes (Rstudio -> Build -> Build & Reload || Shift+CMD+B)
* test cases are in the test subdirectory (~/GenVisR/tests/testthat/) and use the testthat package, tests must be named with the character string "test" in the beginning of the name (i.e. for the waterfall function the test file would look like test-waterfall.R). To run tests run devtools::test()
