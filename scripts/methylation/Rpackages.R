if (!any(rownames(installed.packages()) == "BiocManager")){
  install.packages("BiocManager")
}
library(BiocManager)


if (!any(rownames(installed.packages()) == "bsseq")){
  BiocManager::install("bsseq")
}
library(bsseq)

if (!any(rownames(installed.packages()) == "DSS")){
  BiocManager::install("DSS")
}
library(DSS)

if (!any(rownames(installed.packages()) == "dplyr")){
  install.packages("dplyr")
}
library(dplyr)

if (!any(rownames(installed.packages()) == "ggplot2")){
  install.packages("ggplot2")
}
library(ggplot2)

if (!any(rownames(installed.packages()) == "kableExtra")){
  install.packages("kableExtra")
}
library(kableExtra)

if (!any(rownames(installed.packages()) == "knitr")){
  install.packages("knitr")
}
library(knitr)

