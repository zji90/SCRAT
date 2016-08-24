SCRAT: Single-Cell Regulome Analysis Tool
====

## SCRAT Online User Interface

SCRAT is also available online at https://zhiji.shinyapps.io/scrat/

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

To install the latest version of SCRAT package via Github, run following commands in R (R version 3.2.5 and up is recommended. See Q & A if you want to install SCRAT on an older version of R):
```{r }
source("https://raw.githubusercontent.com/zji90/SCRATdata/master/installcode.R")
```
To launch user interface after installation, run following commands in R:
```{r }
library(SCRAT)
SCRATui()
```

## SCRAT User Manual

SCRAT user manual is available on Github: https://github.com/zji90/SCRATdata/blob/master/manual.pdf

## SCRAT Example Data

SCRAT example data from GM12878 and HEK293T single-cell ATAC-seq are available on Github: https://github.com/zji90/SCRATdata

## Q & A

How to install SCRAT with older versions of R (e.g., R 3.1.3)

Since older version of R does not support https connection, to install the SCRAT package, instead of using:
```{r }
source("https://raw.githubusercontent.com/zji90/SCRATdata/master/installcode.R")
```

Please use the following commands:
```{r }
source("http://bioconductor.org/biocLite.R")
biocLite("GenomicAlignments")
if(!require("shiny")) install.packages("shiny",dependencies=T)
if(!require("ggplot2")) install.packages("ggplot2",dependencies=T)
if(!require("reshape2")) install.packages("reshape2",dependencies=T)
if(!require("gplots")) install.packages("gplots",dependencies=T)
if(!require("pheatmap")) install.packages("pheatmap",dependencies=T)
if(!require("scatterD3")) install.packages("scatterD3",dependencies=T)
if(!require("DT")) install.packages("DT",dependencies=T)
if(!require("mclust")) install.packages("mclust",dependencies=T)
if(!require("tsne")) install.packages("tsne",dependencies=T)
if(!require("devtools")) install.packages("devtools",dependencies=T)
if(!require("png")) install.packages("png",dependencies=T)
devtools::install_github("zji90/SCRAT")
```

## Contact the Author
Author: Zhicheng Ji, Weiqiang Zhou, Hongkai Ji

Report bugs and provide suggestions by sending email to:

Maintainer: Zhicheng Ji (zji4@jhu.edu)

Or open a new issue on this Github page
