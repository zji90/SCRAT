SCRAT: Tools For Single-Cell ANalysis
====

## Overview
Given single-cell RNA-seq data and true experiment time of cells or pseudo-time cell ordering, SCRAT provides convenient functions for users to assign genes into different gene expression patterns such as constant, monotone increasing and increasing then decreasing. SCRAT then performs GO enrichment analysis to analysis the functional roles of genes with same or similar patterns.

## SCRAT Online User Interface
SCRAT user interface can be directly launched online without installing any software package: https://zhiji.shinyapps.io/SCRAT. PLEASE NOTE: Currently the online version only allows one concurrent user. If the online user interface shows "please wait" for a long time, probably another user is using the online interface and please come back at another time. Users are recommended to install SCRAT on their own computers with following procedures.

## SCRAT Installation

SCRAT software can be installed via Github.
Users should have R installed on their computer before installing GSCA. R can be downloaded here: http://www.r-project.org/.
To install the latest version of GSCA package via Github, run following commands in R:
```{r }
if (!require("devtools"))
  install.packages("devtools")
devtools::install_github("SCRAT","zji90")
```
To launch user interface after installation, run following commands in R:
```{r }
library(SCRAT)
SCRATui()
```
For users with R programming experience, command line tools are also available in SCRAT R package. Please check the manual package included in the package for details.

## SCRAT Datasets
The example datasets for SCRAT paper can be downloaded on Github: https://github.com/zji90/SCRATdata

## Contact the Author
Author: Zhicheng Ji, Hongkai Ji

Report bugs and provide suggestions by sending email to:

Maintainer: Zhicheng Ji (zji4@jhu.edu)

Or open a new issue on this Github page
