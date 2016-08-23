SCRAT: Single-Cell Regulome Analysis Tool
====

## SCRAT Installation

SCRAT software can be installed via Github.
Users should have R installed on their computer before installing SCRAT. R can be downloaded here: http://www.r-project.org/  (64-bit R is recommended).
Users should first install the SCRAT data packages by running following commands in R:
```{r }
if (!require("devtools"))
  install.packages("devtools")
devtools::install_github("SCRATdatahg19","zji90")
devtools::install_github("SCRATdatahg38","zji90")
devtools::install_github("SCRATdatamm10","zji90")
devtools::install_github("SCRATdatamm9","zji90")
```

To install the latest version of SCRAT package via Github, run following commands in R:
```{r }
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("GenomicAlignments")
for (pkg in c("shiny","ggplot2","reshape2","pheatmap","scatterD3","DT","gplots","mclust","tsne","devtools"))
      if (!require(pkg)) install.packages(pkg)
devtools::install_github("zji90/SCRAT")
```
To launch user interface after installation, run following commands in R:
```{r }
library(SCRAT)
SCRATui()
```

## SCRAT Online User Interface

SCRAT is also available online at https://zhiji.shinyapps.io/scrat/

## SCRAT User Manual

SCRAT user manual is available on Github: https://github.com/zji90/SCRATdata/blob/master/manual.pdf

## SCRAT Example Data

SCRAT example data from GM12878 and HEK293T single-cell ATAC-seq are available on Github: https://github.com/zji90/SCRATdata

## Q an A

Why SCRAT installation fails on my computer?

If both 32-bit and 64-bit R are installed on your computer, it is possible that SCRAT fails to be installed. In this case please uninstall the 32-bit R and only keep the 64-bit R and then try installing SCRAT again.

## Contact the Author
Author: Zhicheng Ji, Weiqiang Zhou, Hongkai Ji

Report bugs and provide suggestions by sending email to:

Maintainer: Zhicheng Ji (zji4@jhu.edu)

Or open a new issue on this Github page
