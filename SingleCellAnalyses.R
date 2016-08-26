## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(out.width='750px', out.height='750px', dpi=300,
                      fig.height=7, fig.width=7)
knitr::opts_knit$set(root.dir="~/FestivalWorkshopSC/BrainAtlas")

## ----Check for data files, eval=TRUE, echo=TRUE--------------------------
setwd("~/FestivalWorkshopSC/BrainAtlas")
file.exists("cell_metadata.csv")
file.exists("genes_counts.csv")
file.exists("README.txt")

## ----Read in Counts, eval=TRUE, echo=TRUE--------------------------------
counts <- read.csv("genes_counts.csv", stringsAsFactors = FALSE, header=TRUE, row.names = 1)
str(counts[,1:20]) # restrict to the first 20 columns (cells) 

## ----Read in Cell Metadata, eval=TRUE, echo=TRUE-------------------------
cells <- read.csv("cell_metadata.csv", stringsAsFactors = FALSE, header = TRUE)
str(cells)

## ----Peek at README.txt file, eval=TRUE, echo=TRUE, engine="bash"--------
# This is a bash command, to be executed at the command line (not within R);
# Alternatively, simply open the README.txt in your favorite text editor to view its contents
head README.txt

## ----Check for desired packages, eval=TRUE, echo=TRUE, results="hide", message=FALSE, warning=FALSE----
require(scde)     #bioconductor
require(monocle)  #bioconductor
require(sincell)  #bioconductor
require(scDD)     #github
require(ggplot2)  #cran
require(devtools) #cran

## ----install bioconductor package, echo=TRUE, eval=FALSE, results="hide", message=FALSE----
## source("http://bioconductor.org/biocLite.R")
## biocLite("monocle")

## ----install cran packages, echo=TRUE, eval=FALSE------------------------
## install.packages(devtools)

## ----install github packages, echo=TRUE, eval = FALSE--------------------
## install.packages("devtools")
## devtools::install_github("kdkorthauer/scDD")

## ----PCA, eval = TRUE, echo = TRUE---------------------------------------
# extract top 500 variable genes
gene.var <- apply(counts, 1, function(x) var(log(x[x>0])))
counts.top500 <- counts[which(rank(-gene.var)<=500),]

counts.pca <- prcomp(log(counts.top500+1),
                   center = TRUE,
                   scale. = TRUE) 
summary(counts.pca)$importance[,1:5]
plot(counts.pca, type="l", main="Top 10 PCs")

color_class <- rainbow(length(unique(cells$major_class)))
plot(counts.pca$rotation[,1], counts.pca$rotation[,2], 
      xlab="PC 1", ylab="PC 2", col=color_class[as.numeric(factor(cells$major_class))], pch=20,
      main="PCA plot of cells colored by derived major class")

color_class <- rainbow(length(unique(cells$layer_dissectoin)))
plot(counts.pca$rotation[,1], counts.pca$rotation[,2], 
      xlab="PC 1", ylab="PC 2", col=color_class[as.numeric(factor(cells$layer_dissectoin))], pch=20,
      main="PCA plot of cells colored by Dissection Layer")

## ----Detection Rate, eval = TRUE, echo = TRUE----------------------------
detectionRate <- apply(counts, 2, function(x) sum(x > 0) / length(x))
hist(detectionRate)

## ----Highly Variable Genes, eval = TRUE, echo = TRUE---------------------
library(sincell)

## ----DE, eval = TRUE, echo = TRUE----------------------------------------
library(scde)
library(scDD)

## ----Pseudotime, eval = TRUE, echo = TRUE--------------------------------
library(monocle)

## ----Oscillating genes, eval = TRUE, echo = TRUE-------------------------
library(Oscope)

## ----Extract Code Snippets, eval = TRUE, echo = FALSE--------------------
# This snippet just generates a .R file that only contains the code snippets within this document
library(knitr)
purl("~/Desktop/scRNAseq/FestivalWorkshop2016/FestivalWorkshopVignettes/SingleCellAnalyses.Rmd",  output="~/Desktop/scRNAseq/FestivalWorkshop2016/FestivalWorkshopVignettes/SingleCellAnalyses.R")

