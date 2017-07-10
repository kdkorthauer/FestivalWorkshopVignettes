## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(out.width='750px', out.height='750px', dpi=300,
                      fig.height=7, fig.width=7)
knitr::opts_knit$set(root.dir="~/FestivalWorkshopSC/BrainAtlas")

## ----Check for data files, eval=TRUE, echo=TRUE--------------------------
setwd("~/FestivalWorkshopSC/BrainAtlas")
file.exists("cell_metadata.csv")
file.exists("genes_counts.csv")
file.exists("genes_rpkm.csv")
file.exists("ercc_counts.csv")
file.exists("README.txt")

## ----Read in Counts and RPKM, eval=TRUE, echo=TRUE-----------------------
counts <- read.csv("genes_counts.csv", stringsAsFactors = FALSE, header=TRUE, row.names = 1)
str(counts[,1:20]) # restrict to the first 20 columns (cells)

rpkm <- read.csv("genes_rpkm.csv", stringsAsFactors = FALSE, header=TRUE, row.names = 1)
str(rpkm[,1:20]) # restrict to the first 20 columns (cells) 
rpkm <- log(rpkm + 1, 2) # put on the log scale

## ----Read in Cell Metadata, eval=TRUE, echo=TRUE-------------------------
cells <- read.csv("cell_metadata.csv", stringsAsFactors = FALSE, header = TRUE)
str(cells)

## ----Read in ERCC, eval=TRUE, echo=TRUE----------------------------------
ercc <- read.csv("ercc_counts.csv", stringsAsFactors = FALSE, header = TRUE, row.names=1)
str(ercc[,1:20]) # restrict to the first 20 columns (cells)

# remove the tdTomato row
whichTomato <- grep("tdTomato", rownames(ercc))
ercc <- ercc[-whichTomato,]

# combine the counts and erccs into the same data frame
all.equal(colnames(counts), colnames(ercc)) # check that the cells are in the same order across the two datasets
counts <- rbind(counts, ercc)  

## # This is a bash command, to be executed at the command line (not within R);

## ----Check for desired packages, eval=TRUE, echo=TRUE, results="hide", message=FALSE, warning=FALSE----
require(scde)         #bioconductor
require(monocle)      #bioconductor
require(scran)        #bioconductor
require(scater)       #bioconductor
require(Biobase)      #bioconductor
require(EBSeq)        #bioconductor
require(scDD)         #bioconductor
require(SummarizedExperiment) # bioconductor
require(ggplot2)      #cran
require(RColorBrewer) #cran

## ----install bioconductor package, echo=TRUE, eval=FALSE, results="hide", message=FALSE----
## source("http://bioconductor.org/biocLite.R")
## biocLite("scDD")

## ----install cran packages, echo=TRUE, eval=FALSE------------------------
## install.packages(ggplot2)

## ----examine data, echo=TRUE, eval = TRUE--------------------------------
# Take the dimensions of the count and normalized objects
dim(rpkm)
dim(counts)

#both objects contain the same genes in the same order
head(rownames(rpkm))
head(rownames(counts))

#as well as the same samples in the same order
head(colnames(rpkm))
head(colnames(counts))


## ----hist counts, eval=TRUE, echo=TRUE-----------------------------------
hist(as.vector(as.matrix(counts)))
hist(as.vector(as.matrix(rpkm)))

## ----table metadata, eval=TRUE, echo=TRUE--------------------------------
#There is a row for each column in the gene expression matrices
dim(cells)
#Each column represents a different variable associated with the dataset
colnames(cells)
#For example, the columns cre contains the cre-drivers used. We can see how many cells of each line are present in the dataset
table(cells$cre)
barplot(table(cells$cre))

## ----PCA, eval = TRUE, echo = TRUE---------------------------------------
# extract top 1000 variable genes
gene.var <- apply(rpkm, 1, function(x) var(x[x>0]))
rpkm.top1000 <- rpkm[which(rank(-gene.var)<=1000),]

rpkm.pca <- prcomp(rpkm.top1000,
                   center = TRUE,
                   scale. = TRUE) 
summary(rpkm.pca)$importance[,1:5]
plot(rpkm.pca, type="l", main="Top 10 PCs")

#Plot the first PCs
cells.inform <- cells[,c(2,3,4,7,11,12,14)]
df <- data.frame(PC1 = rpkm.pca$rotation[,1],
                 PC2 = rpkm.pca$rotation[,2],
                 cells.inform
                )

ggplot(df, aes(PC1, PC2, colour = major_class))+
  geom_point()+
  labs(title = "PCA plot of cells colored by Major Class")

#it's often advantageous to run a quick association of the top components of variation with your known variables
#Association with PC1
model.pc1 <- anova(lm(PC1 ~. , df[,-c(2)]))
#Association with PC2
model.pc2 <- anova(lm(PC2 ~. , df[,-c(1)]))
summary(model.pc2)

#many of these look significant
ggplot(df, aes(PC1, PC2, colour = mRNA_percent))+
  geom_point()+
  labs(title = "PCA plot of cells colored by alignment to mRNA")


ggplot(df, aes(PC1, PC2, colour = cre))+
  geom_point()+
  labs(title = "PCA plot of cells colored by Cre reporter")
  
rm(gene.var, rpkm.pca, rpkm.top1000, rpkm); gc()  # clean up workspace

## ----Preprocess, eval = TRUE, echo = TRUE--------------------------------
library(scater)
library(scran)
rownames(cells) <- cells$long_name

# construct a SCESet that also contains the gene counts, ercc counts, and metadata
eset <- newSCESet(countData = counts, phenoData = AnnotatedDataFrame(cells))

rm(counts); gc() # remove counts matrix to free up memory (the counts are now stored in the eset object)

## ----QC, eval=TRUE, echo=TRUE--------------------------------------------
# QC to compare the level of dropout in endogeneous genes to ERCC spike ins in raw data
eset <- calculateQCMetrics(eset, feature_controls=list(ERCC=which(grepl("^ERCC", rownames(eset)))))

#designate which rows contain spikeins instead of genes (for HVG analysis)
setSpike(eset) <- "ERCC"

plotQC(eset, type = "exprs-freq-vs-mean")

## ----filter, eval=TRUE, echo=TRUE----------------------------------------
# first, filter out genes that are almost always zero (at least 50 out of 1679 cells must have nonzero expression)
keep <- rowSums(counts(eset) > 0) >= 50
eset <- eset[keep,] 
sum(keep)

# next, filter out genes with very low average nonzero expression across all cells (requres a mean of at least 5 counts across all cells)
keep <- rowMeans(counts(eset)) >= 5
eset <- eset[keep,]
sum(keep)

# normalize counts for library size using the pool & deconvolve method of Lun et al. (2016)
eset <- computeSumFactors(eset, sizes=c(20, 40, 60, 80))
summary(sizeFactors(eset))

# use the size factors calculated above to normalize the counts - these get placed in the 'exprs' slot
eset <- normalize.SCESet(eset)

## ----plot size factors, eval=TRUE, echo=TRUE-----------------------------
plot(sizeFactors(eset), colSums(counts(eset))/1e6, log="xy",
    ylab="Library Size (Total Counts in Millions)", xlab="Pooled Size Factor Estimate",
    main="Normalization factor versus library size")

## ----Detection Rate, eval = TRUE, echo = TRUE----------------------------
detectionRate <- apply(counts(eset), 2, function(x) sum(x > 0) / length(x))
hist(detectionRate)

plot(detectionRate, df$PC1, pch=20, col="grey", ylab="PC 1")

rm(detectionRate, df); gc()  #clean up workspace

## ----Indentify Variable Genes, eval=TRUE, echo=TRUE----------------------
var.fit <- trendVar(eset, trend="loess", use.spikes=FALSE, span=0.2)
var.out <- decomposeVar(eset, var.fit)

# plot the mean versus variance of log-expression, along with the technical variance fit
plot(var.out$mean, var.out$total, pch=16, cex=0.6, xlab="Mean log-expression", 
    ylab="Variance of log-expression")
o <- order(var.out$mean)
lines(var.out$mean[o], var.out$tech[o], col="dodgerblue", lwd=2)

## ----Color by dropout, eval=TRUE, echo=TRUE------------------------------
# function to create a gradient of colors
color.gradient <- function(x, colors=c("orange2", "blue"), colsteps=500) {
  return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )
}
pzero <- apply(counts(eset), 1, function(x) sum(x > 0) / length(x))

# replot with color by proportion of zero
plot(var.out$mean, var.out$total, pch=16, cex=0.6, xlab="Mean log-expression", 
    ylab="Variance of log-expression", col=color.gradient(pzero))
colscale <- c(0,0.25,0.5,0.75,1)
legend(11.5,40, title="Proportion cells zero", legend=paste0(1-round(quantile(pzero, colscale),2)),
        col=color.gradient(quantile(pzero, colscale)), pch=16)

rm(pzero)    # clean up workspace

## ----Indentify Variable Genes with spikes, eval=TRUE, echo=TRUE----------
var.fit.spike <- trendVar(eset, trend="loess", span=0.3, use.spikes=TRUE)
var.out.spike <- decomposeVar(eset, var.fit.spike)

# plot the mean versus variance of log-expression, along with the technical variance fit
plot(var.out.spike$mean, var.out.spike$total, pch=16, cex=0.6, xlab="Mean log-expression", 
    ylab="Variance of log-expression")
o <- order(var.out$mean)
lines(var.out$mean[o], var.out$tech[o], col="dodgerblue", lwd=2)
o <- order(var.out.spike$mean)
lines(var.out.spike$mean[o], var.out.spike$tech[o], col="red", lwd=2)
points(var.fit.spike$mean, var.fit.spike$var, col="red", pch=16)

## ----Extract Variable Genes, eval=TRUE, echo=TRUE------------------------
# extract and examine the top 1000 genes by biological variance
top.hvg <- order(var.out.spike$bio, decreasing=TRUE)[1:1000]
head(var.out.spike[top.hvg,])

# construct a new eset object that only contains the highly variable genes for downstream analysis
eset.hvg <- eset[top.hvg,]

# plot distribution of the top 25 highly variable genes
top25 <- top.hvg[1:25]
boxplot(t(exprs(eset)[top25,]), las=2, ylab="Normalized log-expression", col="dodgerblue", main="Top 25 Highly Variable Genes")

## ----Replot HVG, eval=TRUE, echo=TRUE------------------------------------
# plot the mean versus variance of log-expression, along with the technical variance fit
plot(var.out$mean, var.out$total, pch=16, cex=0.6, xlab="Mean log-expression", 
    ylab="Variance of log-expression")
o <- order(var.out.spike$mean)
lines(var.out.spike$mean[o], var.out.spike$tech[o], col="red", lwd=2)
points(var.out$mean[top.hvg], var.out$total[top.hvg], col="darkorange2", cex=0.6, pch=16)
points(var.out$mean[top25], var.out$total[top25], col="red4", cex=0.7, pch=15)
legend(11.5, 40, legend=c("Top 25 HVGs", "Top 1000 HVGs"), pch = c(15, 16), 
       col= c("red4", "darkorange2"), cex=1.1)

rm(var.out, var.out.spike, var.fit, var.fit.spike) # clean up workspace

## ----heatmap HVG, eval=TRUE, echo=TRUE-----------------------------------
# extract matrix of hvg expression for plotting
m <- exprs(eset)[top25,]

# plot heatmap 
library(RColorBrewer)
heatmap(m/apply(m,1,max),zlim=c(0,1),labCol=NA, col=brewer.pal(9,"YlOrRd"),
        scale="none",ColSideColors=rainbow(5)[as.numeric(factor(phenoData(eset)$major_class))])
par(lend = 1)           # square line ends for the color legend
legend("topleft", inset=-0.04,
    legend = unique(phenoData(eset)$major_class),
    col = rainbow(5)[as.numeric(as.factor(unique(phenoData(eset)$major_class)))],
    lty= 1, lwd = 5, cex=0.6
)

rm(m)  # clean up workspace

## ----scde prep, eval = TRUE, echo = TRUE---------------------------------
library(scde)

# find DE genes between excitatory and inhibitory neuronal subtypes
# first, find which cells are neuron (inhibitory or excitatory)
which.neur <- which(pData(eset.hvg)$major_class %in% c("Inhibitory", "Excitatory"))

# next, construct the group factor labeling for neuron subtype using the index of neuron cells `which.neur' to subset the SCESet object
group <- factor(pData(eset.hvg[,which.neur])$major_class)

# re-create the group factor for the subsampled cell indices
group <- factor(pData(eset.hvg[,which.neur])$major_class)
names(group) <- sampleNames(eset.hvg[,which.neur])

# save counts matrix as integer class; note SCDE requires raw (unnormalized) counts as input. 
cts <- apply(counts(eset.hvg[,which.neur]),2,function(x) {storage.mode(x) <- 'integer'; x}) 

## ----scde fit, eval=TRUE, echo=TRUE--------------------------------------
# fit error models using the scde.error.models function (the following steps have already been run in the interest of time - it is computationally intensive)
########## err.mod <- scde.error.models(counts = cts, groups = group, n.cores = 10, min.nonfailed = 30,
########## 						             verbose=1, save.crossfit.plots=FALSE, threshold.segmentation = TRUE,
##########                         min.size.entries = 500, save.model.plots = FALSE, linear.fit=FALSE,
##########                         min.pairs.per.cell=20, max.pairs=10000)

# estimate the prior for gene expression 
########## prior.mod <- scde.expression.prior(models = err.mod, counts = cts, length.out = 400, show.plot = FALSE, max.quantile=0.9999)

# test for differential expression for all genes using the scde.expression.difference function
########## exp.diff <- scde.expression.difference(models = err.mod, counts = cts, prior = prior.mod, groups = group, n.cores = 1)


## ----scde DE, eval=TRUE, echo=TRUE---------------------------------------
# load the pre-computed results objects err.mod, prior.mod, and exp.diff (all contained
# in the file "scde.results.RData").  Assumes the file is saved to the current working
# directory
load("scde.results.RData")

# view the top of the main results table
head(exp.diff)

# add a column with p-values and multiplicity-corrected pvalues
exp.diff$Pval <- 2*pnorm(-abs(exp.diff$Z))
exp.diff$cPval <- p.adjust(exp.diff$Pval, method="fdr")

## ----scde results, eval=TRUE, echo=TRUE----------------------------------
# count how many genes have an fdr-adjusted p-value less than 0.05
sum(exp.diff$cPval < 0.05)

# reorder the genes by pvalue
exp.diff <- exp.diff[order(exp.diff$cPval),]
head(exp.diff)

# create a .csv file that lists all the DE genes with fdr-adjusted p-value less than 0.05
write.csv(exp.diff[exp.diff$cPval < 0.05, ], file = "scdeResults_DEgenes.csv", row.names = TRUE, quote = FALSE)

## ----scde plot de, eval=TRUE, echo =TRUE---------------------------------
# visualize the results for a particular DE gene
exp.diff[rownames(exp.diff)=="Gad1",]
scde.test.gene.expression.difference("Gad1", models = err.mod, counts = cts, prior = prior.mod)

## ----scde plot ee, eval=TRUE, echo=TRUE----------------------------------
exp.diff[rownames(exp.diff)=="Rgs17",]
scde.test.gene.expression.difference("Rgs17", models = err.mod, counts = cts, prior = prior.mod)

## ----scde plot Vip, eval=TRUE, echo=TRUE---------------------------------
exp.diff[rownames(exp.diff)=="Vip",]
scde.test.gene.expression.difference("Vip", models = err.mod, counts = cts, prior = prior.mod)

# reset plot window
dev.off()

## ----heatmap scde, eval=TRUE, echo=TRUE----------------------------------
# extract matrix of hvg expression for plotting
m <- exprs(eset.hvg)[rownames(eset.hvg) %in% rownames(exp.diff[exp.diff$cPval < 0.05,]), colnames(eset.hvg) %in% names(group)]

# plot heatmap 
library(RColorBrewer)
heatmap(m/apply(m,1,max),zlim=c(0,1),labCol=NA, col=brewer.pal(9,"YlOrRd"),
        scale="none",ColSideColors=rainbow(2)[as.numeric(as.factor(group))])
par(lend = 1)           # square line ends for the color legend
legend("topleft", inset=-0.01,
    legend = unique(group),
    col = rainbow(2)[as.numeric(as.factor(unique(group)))],
    lty= 1, lwd = 5, cex=0.6
)

rm(m, exp.diff, err.mod, prior.mod, cts) # clean up workspace

## ----scdd setup, eval = TRUE, echo = TRUE--------------------------------
library(scDD)
library(Biobase)
library(SummarizedExperiment)

# construct object to send to scDD; we'll reuse group object created for SCDE analysis and create a new (raw-scale) normalized counts matrix
condition <- as.numeric(as.factor(group))
names(condition) <- names(group)
cts <- exp(exprs(eset.hvg[,which.neur])) - 1 #scdd expects raw counts (not log-scaled)
se.scdd <- SummarizedExperiment(assays=list("NormCounts"=cts),
                                colData=data.frame(condition))
rm(cts) # clean up workspace

# Create list of prior parameters for model fitting 
prior_param=list(alpha=0.01, mu0=0, s0=0.01, a0=0.01, b0=0.01)

## ----scdd fit, eval=TRUE, echo=TRUE--------------------------------------
# find DE genes between excitatory and inhibitory neuronal subtypes
set.seed(6767) # set random seed for reprodicbility
se.scdd <- scDD(se.scdd, prior_param=prior_param, permutations=0, testZeroes=FALSE)

head(results(se.scdd))

# how many genes were significantly differentially distributed?
sum(results(se.scdd)$nonzero.pvalue.adj < 0.05)

# what categories to the DD genes belong to?
table(results(se.scdd)$DDcategory)

## ----scdd gad1, eval=TRUE, echo=TRUE-------------------------------------
results(se.scdd)[results(se.scdd)$gene == "Gad1",]
sideViolin(assays(se.scdd)$NormCounts[rownames(se.scdd) ==  "Gad1",], colData(se.scdd)$condition,
            title.gene="Gad1")

## ----scdd Rgs17, eval=TRUE, echo=TRUE------------------------------------
results(se.scdd)[results(se.scdd)$gene == "Rgs17",]
sideViolin(assays(se.scdd)$NormCounts[rownames(se.scdd) ==  "Rgs17",], colData(se.scdd)$condition,
            title.gene="Rgs17")

## ----scdd Vip, eval=TRUE, echo=TRUE--------------------------------------
results(se.scdd)[results(se.scdd)$gene == "Vip",]
sideViolin(assays(se.scdd)$NormCounts[rownames(se.scdd) ==  "Vip",], colData(se.scdd)$condition,
            title.gene="Vip")

## ----scdd Pla2g7, eval=TRUE, echo=TRUE-----------------------------------
results(se.scdd)[results(se.scdd)$gene == "Pla2g7",]
sideViolin(assays(se.scdd)$NormCounts[rownames(se.scdd) ==  "Pla2g7",], colData(se.scdd)$condition,
            title.gene="Pla2g7")

## ----scdd Scg2, eval=TRUE, echo=TRUE-------------------------------------
results(se.scdd)[results(se.scdd)$gene == "Scg2",]
sideViolin(assays(se.scdd)$NormCounts[rownames(se.scdd) ==  "Scg2",], colData(se.scdd)$condition,
            title.gene="Scg2")

## ----heatmap DM, eval=TRUE, echo=TRUE------------------------------------
# extract matrix of hvg expression for plotting
DMgenes <- results(se.scdd)$gene[results(se.scdd)$DDcategory == "DM"] 
m <- exprs(eset.hvg)[rownames(eset.hvg) %in% DMgenes, colnames(eset.hvg) %in% names(group)]

# plot heatmap 
library(RColorBrewer)
heatmap(m/apply(m,1,max),zlim=c(0,1),labCol=NA, col=brewer.pal(9,"YlOrRd"),
        scale="none",ColSideColors=rainbow(2)[as.numeric(as.factor(group))])
par(lend = 1)           # square line ends for the color legend
legend("topleft", inset=-0.01,
    legend = unique(group),
    col = rainbow(2)[as.numeric(as.factor(unique(group)))],
    lty= 1, lwd = 5, cex=0.6
)

rm(m, se.scdd)  # clean up workspace

## ----Monocle Object, eval = TRUE, echo = TRUE----------------------------
library(monocle)
# construct a CellDataSet object with our SCESet object that contains only the top 1000 highly variable genes
cts <- exp(exprs(eset)) - 1 #scdd expects raw counts (not log-scaled)
genes <- data.frame(gene_short_name=rownames(eset))
rownames(genes) <- genes$gene_short_name
pheno <-  phenoData(eset)@data
rownames(pheno) <- pheno$long_name
cset <- newCellDataSet(as(as.matrix(cts), "sparseMatrix"), 
                       phenoData = new("AnnotatedDataFrame", 
                                       data=pheno),
                       featureData = new("AnnotatedDataFrame", data = genes),
                       expressionFamily=negbinomial.size())
class(cset)
cset <- estimateSizeFactors(cset)
cset <- estimateDispersions(cset)

## ----Pseudotime algorithm, eval = TRUE, echo = TRUE----------------------
# Run Monocle on subset of Inhibitory neurons
cset.inhibitory <- cset[,phenoData(cset)$major_class=="Inhibitory" & 
                         phenoData(cset)$layer_dissectoin %in% c("lower", "upper")]
cset.inhibitory <- setOrderingFilter(cset.inhibitory, ordering_genes=rownames(cset.inhibitory)[1:100])
cset.inhibitory <- reduceDimension(cset.inhibitory) # Reduce dimensionality
cset.inhibitory <- orderCells(cset.inhibitory) # Order cells

# Run Monocle on subset of Excitatory neurons
cset.excitatory <- cset[,phenoData(cset)$major_class=="Excitatory" & 
                         phenoData(cset)$layer_dissectoin %in% c("L1", "L2/3", "L4", "L5", "L6", "L6a", "L6b")]
cset.excitatory <- setOrderingFilter(cset.excitatory, ordering_genes=rownames(cset.excitatory)[1:100])
cset.excitatory <- reduceDimension(cset.excitatory) # Reduce dimensionality
cset.excitatory <- orderCells(cset.excitatory) # Order cells

rm(cset) # clean up workspace

## ----Pseudotime plotting, eval = TRUE, echo = TRUE-----------------------
# plotting by various factors
plot_cell_trajectory(cset.inhibitory, color_by="layer_dissectoin") # plot spanning tree
plot_cell_trajectory(cset.inhibitory, color_by = "cre")

# Excitatory
plot_cell_trajectory(cset.excitatory, color_by="layer_dissectoin") # plot spanning tree
plot_cell_trajectory(cset.excitatory, color_by="cre") # plot spanning tree

## ----sesh, eval=TRUE, echo=TRUE------------------------------------------
sessionInfo()

## ----Extract Code Snippets, eval = FALSE, echo = FALSE-------------------
## # This snippet just generates a .R file that only contains the code snippets within this document
## # Not executed; need to run separately to update the .R file after this file is modified ...
## library(knitr)
## knitr:::purl("~/Desktop/scRNAseq/FestivalWorkshop2016/FestivalWorkshopVignettes/SingleCellAnalyses.Rmd",
##      output="~/Desktop/scRNAseq/FestivalWorkshop2016/FestivalWorkshopVignettes/SingleCellAnalyses.R")

