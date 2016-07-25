Dismantling the bulk: examining neuronal heterogeneity using single-cell techniques
================
Sara Linker, Apua Paquola, Roger Lasken, and Keegan Korthauer
7/25/2016

Welcome to the Festival of Genomics workshop on single-cell analyses
--------------------------------------------------------------------

This is an R Markdown document that contains instructions and code for the examples used in todays workshop. The first few steps will check that you have all the packages and data files necessary to carry out all of the analyses.

Getting Started
---------------

### Check that the Brain Atlas data files are present

The following code chunk assumes that Brain Atlas files have been downloaded and placed in a folder in your home directory entitled "FestivalWorkshopSC". If you have downloaded these files to another location, either create a new folder and move the three .csv data files and one documentation text file there, or substitute the file path to where they are currently located for "~/FestivalWorkshopSC".

``` check
setwd("~/FestivalWorkshopSC/")
file.exists("columns-cells.csv")
file.exists("fpkm_table.csv")
file.exists("rows-genes.csv")
file.exists("README.txt")
```

If any of the preceding lines return `FALSE`, double check that you have set the correct working directory and that all four download files have been placed in that folder. If they are missing, you can download the files [here](http://celltypes.brain-map.org/api/v2/well_known_file_download/502999251).

### Read the Brain Atlas data files into R

The `read.csv` function in R is useful for reading in .csv (comma-separated value) files. First, we'll read in the main data file \`fpkm\_table.csv'.

``` read
fpkm <- read.csv("fpkm_table.csv", stringsAsFactors = FALSE, header=TRUE)
```

read in gene and cell metadata files ...

### Check that the desired R packages have been installed

Once a list of desired R packages is finalized, can check that they are installed with

``` check
require(edgeR)
require(scde)
require(monocle)
```

If any of these commands return a message that includes "there is no package called...", then the package is missing and needs to be installed. For Bioconductor packages, for example `edgeR`, this can be done with the following code:

``` install
source("http://bioconductor.org/biocLite.R")
biocLite("monocle")
```
