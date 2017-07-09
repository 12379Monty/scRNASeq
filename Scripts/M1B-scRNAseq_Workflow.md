### A step-by-step workflow for low-level analysis of single-cell RNA-seq data 

<!-- Part II - Analysis of cell types in the brain -->

<style type="text/css"> body, td { font-size: 14px; } code.r{ font-size: 12px; } pre { font-size: 12px } </style>

[bioconductor workflows simpleSingleCell](https://www.bioconductor.org/help/workflows/simpleSingleCell/)



<!-- ***************************************************** -->


<!-- ***************************************************** -->


<!-- ***************************************************** -->



<!-- ***************************************************** -->


## Analysis of cell types in the brain


### Overview

We proceed to a more heterogeneous dataset from a study of cell types in the mouse brain (Zeisel et al. 2015). This contains approximately 3000 cells of varying types such as oligodendrocytes, microglia and neurons. Individual cells were isolated using the Fluidigm C1 microfluidics system and library preparation was performed on each cell using a UMI-based protocol. After sequencing, expression was quantified by counting the number of UMIs mapped to each gene. Count data for all endogenous genes, mitochondrial genes and spike-in transcripts were obtained from [http://linnarssonlab.org/cortex](http://linnarssonlab.org/cortex).

### Count loading

The count data are distributed across several files, so some work is necessary to consolidate them into a single matrix. We define a simple utility function for loading data in from each file. (We stress that this function is only relevant to the current dataset, and should not be used for other datasets. This kind of effort is generally not required if all of the counts are in a single file and separated from the metadata.)



```r
readFormat <- function(infile) { 
    # First column is empty.
    metadata <- read.delim(infile, stringsAsFactors=FALSE, header=FALSE, nrow=10)[,-1] 
    rownames(metadata) <- metadata[,1]
    metadata <- metadata[,-1]
    metadata <- as.data.frame(t(metadata))
    # First column after row names is some useless filler.
    counts <- read.delim(infile, stringsAsFactors=FALSE, header=FALSE, row.names=1, skip=11)[,-1] 
    counts <- as.matrix(counts)
    return(list(metadata=metadata, counts=counts))
}
```
<!-- ******************************************************** -->

Using this function, we read in the counts for the endogenous genes, ERCC spike-ins and mitochondrial genes.


```r
 endo.data <- readFormat(file.path(EXT_DATA, "expression_mRNA_17-Aug-2014.txt"))
 spike.data <- readFormat(file.path(EXT_DATA, "expression_spikes_17-Aug-2014.txt"))
 mito.data <- readFormat(file.path(EXT_DATA, "expression_mito_17-Aug-2014.txt"))
```
<!-- ******************************************************** -->

We also need to rearrange the columns for the mitochondrial data, as the order is not consistent with the other files.


```r
 m <- match(endo.data$metadata$cell_id, mito.data$metadata$cell_id)
 mito.data$metadata <- mito.data$metadata[m,]
 mito.data$counts <- mito.data$counts[,m]
```
<!-- ******************************************************** -->
 
In this particular data set, some genes are represented by multiple rows corresponding to alternative genomic locations. We sum the counts for all rows corresponding to a single gene for ease of interpretation.


```r
 stopifnot(identical(endo.data$metadata$cell_id, spike.data$metadata$cell_id)) # should be the same.
 stopifnot(all(endo.data$metadata$cell_id==mito.data$metadata$cell_id)) # should now be the same.
 
 raw.names <- sub("_loc[0-9]+$", "", rownames(endo.data$counts))
 new.counts <- rowsum(endo.data$counts, group=raw.names, reorder=FALSE)
 endo.data$counts <- new.counts
```
<!-- ******************************************************** -->

The counts are then combined into a single matrix for constructing a **SCESet** object. For convenience, metadata for all cells are stored in the same object for later access.
 

```r
 Cortex.Counts.frm <- rbind(endo.data$counts, mito.data$counts, spike.data$counts)
 
 metadata <- AnnotatedDataFrame(endo.data$metadata)
 Cortex.sce <- newSCESet(countData=Cortex.Counts.frm, phenoData=metadata)
 dim(Cortex.sce)
```

```
## Features  Samples 
##    19896     3005
```
<!-- ******************************************************** -->

We also add annotation identifying rows that correspond to each class of features.


```r
 ## ------------------------------------------------------------------------
 nrows <- c(nrow(endo.data$counts), nrow(mito.data$counts), nrow(spike.data$counts))
 is.spike <- rep(c(FALSE, FALSE, TRUE), nrows)
 is.mito <- rep(c(FALSE, TRUE, FALSE), nrows)
```

<!-- ******************************************************** -->

### Quality control on the cells

The original authors of the study have already removed low-quality cells prior to data publication. Nonetheless, we compute some quality control metrics to check whether the remaining cells are satisfactory.


```r
 Cortex.sce <- calculateQCMetrics(Cortex.sce, feature_controls=list(Spike=is.spike, Mt=is.mito)) 
 setSpike(Cortex.sce) <- "Spike"

 # save this
 save(Cortex.sce, file=file.path(WRKDIR, 'Data', 'Cortex.sce'))
```
<!-- ******************************************************** -->

We examine the distribution of library sizes and numbers of expressed genes across cells (Figure 15).

 
#### Figure 15:
* Histograms of library sizes (left) and number of expressed genes (right) for all cells in the brain dataset.

```r
 load(file.path(WRKDIR, 'Data', 'Cortex.sce'))

 par(mfrow=c(1,2))
 hist(Cortex.sce$total_counts/1e3, xlab="Library sizes (thousands)", main="", 
     breaks=20, col="grey80", ylab="Number of cells")
 hist(Cortex.sce$total_features, xlab="Number of expressed genes", main="", 
     breaks=20, col="grey80", ylab="Number of cells")
```

![plot of chunk libplotbrain](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M1B/libplotbrain-1.png)
<!-- ******************************************************** -->

We also examine the distribution of the proportions of UMIs assigned to mitochondrial genes or spike-in transcripts (Figure 16). The spike-in proportions here are more variable than in the HSC dataset. This may reflect a greater variability in the total amount of endogenous RNA per cell when many cell types are present.
 

#### Figure 16:
* Histogram of the proportion of UMIs assigned to mitochondrial genes (left) or spike-in transcripts (right) across all cells in the brain dataset.


```r
 par(mfrow=c(1,2))
 hist(Cortex.sce$pct_counts_feature_controls_Mt, xlab="Mitochondrial proportion (%)", 
     ylab="Number of cells", breaks=20, main="", col="grey80")
 hist(Cortex.sce$pct_counts_feature_controls_Spike, xlab="ERCC proportion (%)",
     ylab="Number of cells", breaks=20, main="", col="grey80")
```

![plot of chunk UMIplot](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M1B/UMIplot-1.png)
<!-- ******************************************************** -->

We remove small outliers in Figure 15 and large outliers in Figure 16, using a MAD-based threshold as previously described.


```r
 libsize.drop <- isOutlier(Cortex.sce$total_counts, nmads=3, type="lower", log=TRUE)
 feature.drop <- isOutlier(Cortex.sce$total_features, nmads=3, type="lower", log=TRUE)
 mito.drop <- isOutlier(Cortex.sce$pct_counts_feature_controls_Mt, nmads=3, type="higher")
 spike.drop <- isOutlier(Cortex.sce$pct_counts_feature_controls_Spike, nmads=3, type="higher")
```
<!-- ******************************************************** -->

Removal of low-quality cells is then performed by combining the filters for all of the metrics. The vast majority of cells are retained, which suggests that the original quality control procedures were generally adequate.



```r
 Cortex.sce <- Cortex.sce[,!(libsize.drop | feature.drop | spike.drop | mito.drop)]
 data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop), 
     ByMito=sum(mito.drop), BySpike=sum(spike.drop), Remaining=ncol(Cortex.sce))
```

```
##         ByLibSize ByFeature ByMito BySpike Remaining
## Samples         8         3     87       8      2902
```
<!-- ******************************************************** -->

### Cell cycle classification

Application of **cyclone** to the brain dataset suggests that most of the cells are in G1 phase (Figure 17). However, the intepretation of this result requires some caution due to the differences between the test and training datasets. The classifier was trained on C1 SMARTer data (Scialdone et al. 2015) and accounts for the biases in that protocol. The brain dataset uses UMI counts, which has an entirely different set of biases, e.g., 3’-end coverage only, no length bias, no amplification noise. These new biases (and the absence of expected biases) may interfere with accurate classification of some cells.

 

 
#### Figure 17:
* Cell cycle phase scores from applying the pair-based classifier on the brain dataset, where each point represents a cell.


```r
 mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", 
           package="scran"))

 suppressMessages(require(org.Mm.eg.db))
 anno <- select(org.Mm.eg.db, keys=rownames(Cortex.sce), 
                keytype="SYMBOL", column="ENSEMBL")
 ensembl <- anno$ENSEMBL[match(rownames(Cortex.sce), anno$SYMBOL)]
 assignments <- cyclone(Cortex.sce, mm.pairs, gene.names=ensembl)
 plot(assignments$score$G1, assignments$score$G2M, xlab="G1 score", 
      ylab="G2/M score", pch=16)
 abline(0,1, col='red')
 abline(h=.5, v=0.5, col='blue', lty=2) 
```

![plot of chunk phasePlot](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M1B/phasePlot-1.png)
<!-- ******************************************************** -->


An additional complication is that many neuronal cell types are expected to lie in the G0 resting phase, which is distinct from the other phases of the cell **cycle** (Coller, Sang, and Roberts 2006). Application of cyclone to these cells may be suboptimal if each cell must be assigned into one of the G1, S or G2/M phases. To avoid problems from misclassification, we will not perform any processing of this dataset by cell cycle phase. This is unlikely to be problematic for this analysis, as the cell cycle effect will be relatively subtle compared to the obvious differences between cell types in a diverse population. Thus, the former is unlikely to distort the conclusions regarding the latter.
 

### Removing uninteresting genes

Low-abundance genes are removed by applying a simple mean-based filter. We use a lower threshold for UMI counts compared to that used for read counts. This is because the number of transcript molecules will always be lower than the number of reads generated from such molecules. While some information and power will be lost due to the decrease in the size of the counts, this is mitigated by a concomitant reduction in the variability of the counts. Specifically, the use of UMIs eliminates technical noise due to amplification biases (Islam et al. 2014).


#### Figure 18:
* Histogram of log-average counts for all genes in the brain dataset. The filter threshold is represented by the blue line.

Figure 18 suggests that our choice of threshold is appropriate. The filter removes the bulk of lowly expressed genes while preserving the peak of moderately expressed genes.


```r
 ave.counts <- calcAverage(Cortex.sce)
 keep <- ave.counts >= 0.1
 
 hist(log10(ave.counts), breaks=100, main="", col="grey",
     xlab=expression(Log[10]~"average count"))
 abline(v=log10(0.1), col="blue", lwd=2, lty=2)
```

![plot of chunk abhist](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M1B/abhist-1.png)
<!-- ******************************************************** -->

The mean-based filter is applied to the dataset by subsetting sce as previously described. Despite the reduced threshold, the number of retained genes is lower than that in the HSC dataset, simply because the library sizes are much smaller with UMI counts.



```r
 Cortex.sce <- Cortex.sce[keep,]
 nrow(Cortex.sce)
```

```
## Features 
##    10765
```
<!-- ******************************************************** -->

 
Some datasets also contain strong heterogeneity in mitochondrial RNA content, possibly due to differences in mitochondrial copy number or activity between cell types. This heterogeneity will cause mitochondrial genes to dominate the top set of results, e.g., for identification of correlated HVGs. However, these genes are largely uninteresting given that most studies focus on nuclear regulation. As such, we filter them out prior to further analysis. Other candidates for removal include pseudogenes or ribosome-associated genes, which might not be relevant for characterising cell types but can still interfere with the interpretation of the results.



```r
 Cortex.sce <- Cortex.sce[!fData(Cortex.sce)$is_feature_control_Mt,]
```
<!-- ******************************************************** -->

 
### Normalization of cell-specific biases

Normalization of cell-specific biases is performed using the deconvolution method in the **computeSumFactors** function. Here, we cluster similar cells together and normalize the cells in each cluster using the deconvolution method. This improves normalization accuracy by reducing the number of DE genes between cells in the same cluster. Scaling is then performed to ensure that size factors of cells in different clusters are comparable.

 

```r
 clusters <- quickCluster(Cortex.sce)
 Cortex.sce <- computeSumFactors(Cortex.sce, cluster=clusters)
 summary(sizeFactors(Cortex.sce))
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.1288  0.4509  0.8223  1.0000  1.3495  4.9446
```
<!-- ******************************************************** -->

 
Compared to the HSC analysis, more scatter is observed around the trend between the total count and size factor for each cell (Figure 19). This is consistent with an increased amount of DE between cells of different types, which compromises the accuracy of library size normalization (Robinson and Oshlack 2010). In contrast, the size factors are estimated based on median ratios and are more robust to the presence of DE between cells.


#### Figure 19:
* Size factors from deconvolution, plotted against library sizes for all cells in the brain dataset. Axes are shown on a log-scale.


```r
 plot(sizeFactors(Cortex.sce), Cortex.sce$total_counts/1e3, log="xy",
     ylab="Library size (thousands)", xlab="Size factor")
```

![plot of chunk plotNormFactors](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M1B/plotNormFactors-1.png)
<!-- ******************************************************** -->

We also compute size factors specific to the spike-in set, as previously described.
 

```r
 Cortex.sce <- computeSpikeFactors(Cortex.sce, type="Spike", general.use=FALSE)
```
<!-- ******************************************************** -->

Finally, normalized log-expression values are computed for each endogenous gene or spike-in transcript using the appropriate size factors.

 

```r
 Cortex.sce <- normalize(Cortex.sce)
```
<!-- ******************************************************** -->
 

### Checking for important technical factors

Larger experiments contain more technical factors that need to be investigated. In this dataset, factors include the sex of the animal from which the cells were extracted, the age of the animal, the tissue of origin for each cell, and the total spike-in count in each cell. Figure 20 shows that the tissue of origin explains a substantial proportion of the variance for a subset of genes. This is probably because each tissue contains a different composition of cell types, leading to systematic differences in gene expression between tissues. The other factors explain only a small proportion of the variance for most genes and do not need to be incorporated into our downstream analyses.


#### Figure 20:
* Density plot of the percentage of variance explained by each factor across all genes in the brain dataset. For each gene, the percentage of the variance of the normalized log-expression values that is explained by the (log-transformed) total spike-in counts, the sex or age of the mouse, or the tissue of origin is calculated. Each curve corresponds to one factor and represents the distribution of percentages across all genes.


```r
 fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))
 plotExplanatoryVariables(Cortex.sce, variables=c("counts_feature_controls_Spike", 
     "log10_counts_feature_controls_Spike", "sex", "tissue", "age")) + fontsize
```

![plot of chunk plotExplanatoryVariables](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M1B/plotExplanatoryVariables-1.png)
<!-- **************************************************************** -->

#### Replot as pairs

```r
 fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))
 plotExplanatoryVariables(Cortex.sce, method='pairs',
    variables=c("counts_feature_controls_Spike", "log10_counts_feature_controls_Spike", 
                "sex", "tissue", "age")) + fontsize
```

![plot of chunk pairsExplanatoryVariables](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M1B/pairsExplanatoryVariables-1.png)
<!-- **************************************************************** -->

Nonetheless, we demonstrate how to account for uninteresting technical factors by using sex as an example. We set up a design matrix with the sex of the animal as the explanatory factor for each cell. This ensures that any sex-specific changes in expression will be modelled in our downstream analyses. We do not block on the tissue of origin, despite the fact that it explains more of the variance than sex in Figure 20. This is because the tissue factor is likely to be associated with genuine differences between cell types, so including it in the model might regress out interesting biological effects.


```r
 design <- model.matrix(~Cortex.sce$sex)
```
<!-- **************************************************************** -->

Other relevant factors include the chip or plate on which the cells were processed and the batch in which the libraries were sequenced. Blocking on these factors may be necessary to account for batch effects that are often observed in scRNA-seq data (Hicks, Teng, and Irizarry 2015; Tung et al. 2016).

### Identifying correlated HVGs

We identify HVGs that may be involved in driving population heterogeneity. This is done by fitting a trend to the technical variances for the spike-in transcripts. We compute the biological component of the variance for each endogenous gene by subtracting the fitted value of the trend from the total variance. To account for uninteresting factors, we supply **design** to **trendVar** to regress out any technical differences due to sex. HVGs are then detected in the same manner as described in the HSC data analysis.


```r
 var.fit.des <- trendVar(Cortex.sce, trend="loess", span=0.4, design=design)
 var.out.des <- decomposeVar(Cortex.sce, var.fit.des)
 hvg.out.des <- var.out.des[which(var.out.des$FDR <= 0.05 & var.out.des$bio >= 0.5),]
 nrow(hvg.out.des)
```

```
## [1] 1647
```
<!-- **************************************************************** -->
 
Alternatively, for data sets containing multiple batches, a more robust strategy is to perform trend fitting and variance decomposition separately for each batch. This accommodates differences in the mean-variance trends between batches, especially if a different amount of spike-in RNA was added to the cells in each batch. We demonstrate the second approach below by treating each sex as a different “batch”.


```r
 collected <- list()
 for (block in levels(as.factor(Cortex.sce$sex))) {
     cur.Cortex.sce <- Cortex.sce[,Cortex.sce$sex==block]
     cur.Cortex.sce <- normalize(cur.Cortex.sce) 
     var.fit <- trendVar(cur.Cortex.sce, trend="loess", span=0.4)
     collected[[block]] <- decomposeVar(cur.Cortex.sce, var.fit)
 }
```
<!-- **************************************************************** -->

Figure 21 suggests that the trend is fitted accurately to the technical variances for each sex. Errors in fitting are negligible due to the precision of the variance estimates in a large dataset containing thousands of cells. The technical and total variances are also much smaller than those in the HSC dataset. This is due to the use of UMIs which reduces the noise caused by variable PCR amplification. Furthermore, the spike-in trend is consistently lower than the variances of the endogenous genes. This reflects the heterogeneity in gene expression across cells of different types. It also means the previous strategy of fitting a trend to the endogenous variances would not be appropriate here (or necessary, given the quality of the spike-in trend).

 

<!-- **************************************************************** -->
 
### Figure 21:
* Variance of normalized log-expression values against the mean for each gene, calculated across all cells from male (left) or female mice (right). In each plot, the red line represents the mean-dependent trend in the technical variance of the spike-in transcripts (also highlighted as red points).


```r
 par(mfrow=c(1,2))
 for (block in c("1", "-1")) { 
     var.out <- collected[[block]]
     plot(var.out$mean, var.out$total, pch=16, cex=0.6, xlab="Mean log-expression", 
         ylab="Variance of log-expression", main=ifelse(block=="1", "Male", "Female"))
     points(var.out$mean[isSpike(Cortex.sce)], var.out$total[isSpike(Cortex.sce)], col="red", pch=16)
     o <- order(var.out$mean)
     lines(var.out$mean[o], var.out$tech[o], col="red", lwd=2)
 }
```

![plot of chunk hvgplot](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M1B/hvgplot-1.png)
<!-- **************************************************************** -->

Statistics are combined across the two sexes using the **combineVar** function. HVGs are identified as genes with large positive biological components, and are saved to file for future reference. Note that some of the p-values are reported as zero due to numerical imprecision.


```r
 var.out <- do.call(combineVar, collected)
 hvg.out <- var.out[which(var.out$FDR <= 0.05 & var.out$bio >= 0.5),]
 hvg.out <- hvg.out[order(hvg.out$bio, decreasing=TRUE),] 
 nrow(hvg.out)
```

```
## [1] 1712
```

```r
 write.table(file=file.path(TABLE_DIR,"brain_hvg.tsv"), 
      hvg.out, sep="\t", quote=FALSE, col.names=NA)
 head(hvg.out)
```

```
##          mean     total       bio      tech p.value FDR
## Plp1 3.942784 16.268240 16.015984 0.2522561       0   0
## Mal  2.287297 10.898808 10.309165 0.5896426       0   0
## Trf  2.115009  9.958051  9.341630 0.6164211       0   0
## Mbp  2.276203  7.977436  7.388900 0.5885360       0   0
## Mog  1.813190  8.046103  7.369380 0.6767230       0   0
## Apod 1.749638  7.844978  7.181341 0.6636372       0   0
```
<!-- **************************************************************** -->

#### Figure 22:
* Violin plots of normalized log-expression values for the top 10 HVGs in the brain dataset. For each gene, each point represents the log-expression value for an individual cell.


```r
 plotExpression(Cortex.sce, rownames(hvg.out)[1:10], alpha=0.05, jitter="jitter") + fontsize
```

![plot of chunk hvgvioplotbrain](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M1B/hvgvioplotbrain-1.png)
<!-- **************************************************************** -->

To identify genes involved in defining subpopulations, correlations in the expression profiles are calculated between pairs of HVGs. Correlated HVGs are defined as those with significant correlations to one or more other HVGs at a FDR of 5%. Here, the proportion of HVGs with significant correlations is much higher than in the HSC dataset, indicating that strong substructure is present.


```r
 set.seed(100)
 var.cor <- correlatePairs(Cortex.sce, design=design, subset.row=rownames(hvg.out), per.gene=TRUE)
 write.table(file="brain_cor.tsv", var.cor, sep="\t", quote=FALSE, row.names=FALSE)
 head(var.cor)
```

```
##   gene       rho      p.value          FDR limited
## 1 Plp1 0.7489356 2.565215e-06 5.961472e-06    TRUE
## 2  Mal 0.7489356 2.604259e-06 5.961472e-06    TRUE
## 3  Trf 0.7161424 2.592422e-06 5.961472e-06    TRUE
## 4  Mbp 0.7376700 2.618207e-06 5.961472e-06    TRUE
## 5  Mog 0.7034789 2.594387e-06 5.961472e-06    TRUE
## 6 Apod 0.6885070 2.654768e-06 5.961472e-06    TRUE
```

```r
 sig.cor <- var.cor$FDR <= 0.05
 summary(sig.cor)
```

```
##    Mode   FALSE    TRUE 
## logical       1    1711
```
<!-- **************************************************************** -->
 
### Further data exploration with the correlated HVGs

We first remove the sex effect using the **removeBatchEffect** function from the **limma** package (Ritchie et al. 2015). This ensures that any sex-specific differences will not dominate the visualization of the expression profiles. In this manner, we maintain consistency with the use of design in the previous steps. (However, if an analysis method can accept a design matrix, blocking on nuisance factors in the design matrix is preferable to manipulating the expression values with **removeBatchEffect.** This is because the latter does not account for the loss of residual degrees of freedom, nor the uncertainty of estimation of the blocking factor terms.) We store these sex-corrected expression values in the **norm_exprs** field of the **SCESet** object for later use.


```r
 suppressMessages(require(limma))
 adj.exprs <- exprs(Cortex.sce)
 adj.exprs <- removeBatchEffect(adj.exprs, batch=Cortex.sce$sex)
 norm_exprs(Cortex.sce) <- adj.exprs 
```
<!-- **************************************************************** -->

We perform dimensionality reduction on the correlated HVGs to check if there is any substructure. Cells separate into clear clusters in the t-SNE plot (Figure 23), corresponding to distinct subpopulations. This is consistent with the presence of multiple cell types in the diverse brain population.

 
#### Figure 23:
* _t_-SNE plots constructed from the normalized and corrected log-expression values of correlated HVGs for cells in the brain dataset. Each point represents a cell and is coloured according to its expression of the top HVG (left) or _Mog_ (right).


```r
 DUMMY <- 0
 chosen <- var.cor$gene[sig.cor]
 top.hvg <- rownames(hvg.out)[1]
 tsne1 <- plotTSNE(Cortex.sce, exprs_values="norm_exprs", colour_by=top.hvg,
     perplexity=10, rand_seed=100, feature_set=chosen) + fontsize
 tsne2 <- plotTSNE(Cortex.sce, exprs_values="norm_exprs", colour_by="Mog",
     perplexity=10, rand_seed=100, feature_set=chosen) + fontsize
 multiplot(tsne1, tsne2, cols=2)
```

![plot of chunk tsneplot](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M1B/tsneplot-1.png)
<!-- **************************************************************** -->
 


The PCA plot is less effective at separating cells into many different clusters (Figure 24). This is because the first two principal components (PCs) are driven by strong differences between specific subpopulations, which reduces the resolution of more subtle differences between some of the other subpopulations. Nonetheless, some substructure is still visible.

 
#### Figure 24:
* PCA plots constructed from the normalized and corrected log-expression values of correlated HVGs for cells in the brain dataset. Each point represents a cell and is coloured according to its expression of the top HVG (left) or _Mog_ (right).


```r
 pca1 <- plotPCA(Cortex.sce, exprs_values="norm_exprs", colour_by=top.hvg) + fontsize
 pca2 <- plotPCA(Cortex.sce, exprs_values="norm_exprs", colour_by="Mog") + fontsize
 multiplot(pca1, pca2, cols=2)
```

![plot of chunk PCAplot](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M1B/PCAplot-1.png)
<!-- **************************************************************** -->

For both methods, we colour each cell based on the expression of a particular gene. This is a useful strategy for visualizing changes in expression across the lower-dimensional space. It can also be used to characterise each cluster if the selected genes are known markers for particular cell types. For example, Mog can be used to identify clusters corresponding to oligodendrocytes.

### Clustering cells into putative subpopulations

The normalized and sex-adjusted log-expression values for correlated HVGs are used to cluster cells into putative subpopulations. Specifically, we perform hierarchical clustering on the Euclidean distances between cells, using Ward’s criterion to minimize the total variance within each cluster. This yields a dendrogram that groups together cells with similar expression patterns across the chosen genes. An alternative approach is to cluster on a matrix of distances derived from correlations (e.g., as in **quickCluster**). This is more robust to noise and normalization errors, but is also less sensitive to subtle changes in the expression profiles.

Clusters are explicitly defined by applying a dynamic tree cut (Langfelder, Zhang, and Horvath 2008) to the dendrogram. This exploits the shape of the branches in the dendrogram to refine the cluster definitions, and is more appropriate than **cutree** for complex dendrograms. Greater control of the empirical clusters can be obtained by manually specifying **cutHeight** in **cutreeDynamic.**


```r
 suppressMessages(require(dynamicTreeCut))

 chosen.exprs <- norm_exprs(Cortex.sce)[chosen,]
 my.dist <- dist(t(chosen.exprs))
 my.tree <- hclust(my.dist, method="ward.D2")
 
 my.clusters <- unname(cutreeDynamic(my.tree, distM=as.matrix(my.dist), verbose=0))
```
<!-- **************************************************************** -->

Figure 25 contains a clear block-like pattern, representing systematic differences between clusters of cells with distinct expression profiles. This is consistent with the presence of well-defined subpopulations that were previously observed in the dimensionality reduction plots.

 
### Figure 25:
* Heatmap of mean-centred normalized and corrected log-expression values for correlated HVGs in the brain dataset. Dendrograms are formed by hierarchical clustering on the Euclidean distances between genes (row) or cells (column). Column colours represent the cluster to which each cell is assigned after a dynamic tree cut. Heatmap colours are capped at a maximum absolute log-fold change of 5."----
`

```r
 suppressMessages(require(gplots))

 heat.vals <- chosen.exprs - rowMeans(chosen.exprs)
 clust.col <- rainbow(max(my.clusters))
 heatmap.2(heat.vals, col=bluered, symbreak=TRUE, trace='none', cexRow=0.3,
     ColSideColors=clust.col[my.clusters], Colv=as.dendrogram(my.tree), 
     breaks=seq(-5, 5, length.out=21))
```

![plot of chunk heatMap](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M1B/heatMap-1.png)
<!-- **************************************************************** -->

 
This heatmap can be stored at a greater resolution for detailed inspection later. It is advisable to verify the separatedness of the clusters using metrics such as the silhouette width or gap statistic (see the **cluster** package for details). The same statistics can also be used to gauge the optimal parameter values (e.g., cut height, number of clusters) that maximize the separation between clusters.


```r
 pdf(file.path(TABLE_DIR,"brain_heat.pdf"), width=10, height=100) # Large 'height' to show all gene names.
 heatmap.2(heat.vals, col=bluered, symbreak=TRUE, trace='none', cexRow=0.3,
     ColSideColors=clust.col[my.clusters], Colv=as.dendrogram(my.tree), 
     lwid=c(1, 4), lhei=c(1, 50), # Avoid having a very large colour key.
     breaks=seq(-5, 5, length.out=21))
 dev.off()
```

```
## png 
##   2
```
<!-- **************************************************************** -->
 
A complementary approach is to perform PCA on the log-expression values for correlated HVGs and cluster on the first few PCs. This assumes that random technical noise in each gene will be represented by later PCs, while biological substructure involving coregulated groups of genes will contribute more variance and be represented by earlier PCs. The **denoisePCA** function will remove later PCs until the total discard variance is equal to the sum of technical components for all genes used in the PCA. This eliminates technical noise and enriches for biological signal in the remaining PCs, allowing for improved resolution during clustering.



```r
 # Using the technical components computed with a design matrix.
 pcs <- denoisePCA(Cortex.sce, design=design, technical=var.fit.des$trend, subset.row=chosen) 
 dim(reducedDimension(pcs)) # Cluster on this instead of t(chosen.exprs)
```

```
## [1] 2902  126
```
<!-- **************************************************************** -->

 
When examining very heterogeneous datasets, it can be useful to repeat the HVG detection and clustering using only the cells within a particular cluster. This can be achieved by subsetting **sce** according to a particular level of **my.clusters**, and applying the same functions that were previously described. Doing so may identify a different set of correlated HVGs that define heterogeneity within the cluster, as opposed to those that define differences between clusters. This would allow fine-scale structure within each cluster to be explored at greater resolution. For simplicity, though, we will only use the broad clusters corresponding to clear subpopulations in this workflow.


### Detecting marker genes between subpopulations

Once putative subpopulations are identified by clustering, we can identify marker genes for each cluster using the **findMarkers** function. This fits a linear model to the log-expression values for each gene, using methods in the **limma** package (Law et al. 2014; Ritchie et al. 2015). The aim is to test for DE in each cluster compared to the others while blocking on uninteresting factors in **design**. The top DE genes are likely to be good candidate markers as they can effectively distinguish between cells in different clusters.


```r
 markers <- findMarkers(Cortex.sce, my.clusters, design=design)
```
 
For each cluster, the DE results of the relevant comparisons are consolidated into a single output table. This allows a set of marker genes to be easily defined by taking the top DE genes from each pairwise comparison between clusters. For example, to construct a marker set for cluster 1 from the top 10 genes of each comparison, one would filter **marker.set** to retain rows with **Top** less than or equal to 10. Other statistics are also reported for each gene, including the adjusted p-values (see below) and the log-fold changes relative to every other cluster.


```r
 old.digits <- options()$digits
 options(digits=3)
 
 marker.set <- markers[["1"]]
 print(kable(head(marker.set, 10)))
```

```
## 
## 
## | Top|Gene    | FDR|      2|      3|      4|      5|     6|      7|      8|
## |---:|:-------|---:|------:|------:|------:|------:|-----:|------:|------:|
## |   1|Htr3a   |   0| -0.029| -0.112|  0.006| -0.057| -3.16| -0.143| -0.106|
## |   1|Prkar1b |   0| -2.251| -2.673| -2.175| -0.369| -2.99| -1.629| -0.147|
## |   1|Mllt11  |   0| -2.015| -2.564| -2.763|  0.039| -2.90| -1.387|  0.327|
## |   1|Clu     |   0| -0.422| -1.391| -0.416| -0.276| -1.72| -0.653| -5.577|
## |   1|Snap25  |   0| -3.263| -4.599| -4.072| -0.716| -4.39| -3.328| -0.237|
## |   1|Syt11   |   0|  1.255|  1.529|  1.968|  3.677|  1.42|  0.933|  2.061|
## |   2|Kcnip1  |   0| -0.015| -0.043|  0.010| -0.004| -1.45| -0.154| -0.039|
## |   2|Syt1    |   0| -3.668| -3.691| -3.244| -0.967| -4.30| -2.766| -0.483|
## |   2|Ndrg4   |   0| -3.007| -3.441| -3.180| -0.383| -4.02| -2.298| -0.350|
## |   2|Snhg11  |   0| -4.269| -3.737| -3.964| -0.892| -3.80| -3.435| -0.399|
```

```r
 options(digits=old.digits)

 write.table(marker.set, file=file.path(TABLE_DIR,"brain_marker_1.tsv"), 
         sep="\t", quote=FALSE, col.names=NA)
```
We save the list of candidate marker genes for further examination. The **overlapExprss** function may also be useful here, to prioritize candidates where there is clear separation between the distributions of expression values of different clusters.

We visualize the expression profiles of the top candidates to verify that the DE signature is robust. Figure 26 indicates that most of the top markers have strong and consistent up- or downregulation in cells of cluster 1 compared to some or all of the other clusters. Thus, cells from the subpopulation of interest can be identified as those that express the upregulated markers and do not express the downregulated markers.
 
#### Figure 26:
* Heatmap of mean-centred normalized and corrected log-expression values for the top set of markers for cluster 1 in the brain dataset. Column colours represent the cluster to which each cell is assigned, as indicated by the legend.


```r
 top.markers <- marker.set$Gene[marker.set$Top <= 10]
 top.exprs <- norm_exprs(Cortex.sce)[top.markers,,drop=FALSE]
 heat.vals <- top.exprs - rowMeans(top.exprs)
 heatmap.2(heat.vals, col=bluered, symbreak=TRUE, trace='none', cexRow=0.6,
     ColSideColors=clust.col[my.clusters], Colv=as.dendrogram(my.tree), dendrogram='none')
 legend("bottomleft", col=clust.col, legend=sort(unique(my.clusters)), pch=16)
```

![plot of chunk heatmapmarker](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M1B/heatmapmarker-1.png)

Many of the markers in Figure 26 are not uniquely up- or downregulated in the chosen cluster. Testing for unique DE tends to be too stringent as it overlooks important genes that are expressed in two or more clusters. For example, in a mixed population of CD4+-only, CD8+-only, double-positive and double-negative T cells, neither Cd4 or Cd8 would be detected as subpopulation-specific markers because each gene is expressed in two subpopulations. With our approach, both of these genes will be picked up as candidate markers as they will be DE between at least one pair of subpopulations. A combination of markers can then be chosen to characterize a subpopulation, which is more flexible than trying to find uniquely DE genes.

It must be stressed that the (adjusted) p-values computed here cannot be properly interpreted as measures of significance. This is because the clusters have been empirically identified from the data. **limma** does not account for the uncertainty of clustering, which means that the p-values are much lower than they should be. However, this is not a concern in other analyses where the groups are pre-defined. For such analyses, the FDR-adjusted p-value can be directly used to define significant genes for each DE comparison, though some care may be required to deal with plate effects (Hicks, Teng, and Irizarry 2015; Tung et al. 2016). 

The *8SCESet** object can also be easily transformed for use in other DE analysis methods. For example, the **convertTo** function can be used to construct a **DGEList** for input into the **edgeR** pipeline (Robinson, McCarthy, and Smyth 2010). This allows users to construct their own marker detection pipeline, though we find that direct use of **findMarkers** is usually sufficient.


```r
 library(edgeR)
 y <- convertTo(Cortex.sce, type="edgeR")
```

### Additional comments

Having completed the basic analysis, we save the SCESet object with its associated data to file. This is especially important here as the brain dataset is quite large. If further analyses are to be performed, it would be inconvenient to have to repeat all of the pre-processing steps described above.


```r
 save(file=file.path(WRKDIR,'Data',"Cortex.sce"), Cortex.sce)
```

 ### Parameter settings:
   * WRKDIR = /mnt100/home/Dropbox/SingleCell/Jun2017/R
   * FN = M1B-scRNAseq_Workflow
   * Scripts = Scripts
   * EXT_DATA = /mnt100/home/Dropbox/SingleCell//Data
   * TABLE_DIR = /mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/tables/M1B 
   * RUN DATE = Wed Jun 21 22:38:55 2017
 
 

```
## R version 3.4.0 (2017-04-21)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 16.04.2 LTS
## 
## Matrix products: default
## BLAS: /usr/lib/openblas-base/libblas.so.3
## LAPACK: /usr/lib/libopenblasp-r0.2.18.so
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8       
##  [4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
## [10] LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
##  [1] grid      stats4    parallel  methods   stats     graphics  grDevices utils     datasets 
## [10] base     
## 
## other attached packages:
##  [1] gplots_3.0.1          dynamicTreeCut_1.63-1 org.Mm.eg.db_3.4.1    AnnotationDbi_1.38.1 
##  [5] scran_1.4.4           BiocParallel_1.10.1   destiny_2.4.0         mvoutlier_2.0.8      
##  [9] sgeostat_1.0-27       Rtsne_0.13            readxl_1.0.0          R.utils_2.5.0        
## [13] R.oo_1.21.0           R.methodsS3_1.7.1     BiocStyle_2.4.0       scater_1.4.0         
## [17] ggplot2_2.2.1         readr_1.1.1           irlba_2.2.1           Matrix_1.2-10        
## [21] plyr_1.8.4            tidyr_0.6.3           reshape2_1.4.2        stringr_1.2.0        
## [25] vcd_1.4-3             GenomicRanges_1.28.3  GenomeInfoDb_1.12.2   IRanges_2.10.2       
## [29] S4Vectors_0.14.3      Biobase_2.36.2        BiocGenerics_0.22.0   edgeR_3.18.1         
## [33] limma_3.32.2          RColorBrewer_1.1-2    data.table_1.10.4     magrittr_1.5         
## [37] dplyr_0.7.0           knitr_1.16            rmarkdown_1.6        
## 
## loaded via a namespace (and not attached):
##   [1] backports_1.1.0         Hmisc_4.0-3             RcppEigen_0.3.3.3.0    
##   [4] igraph_1.0.1            lazyeval_0.2.0          sp_1.2-4               
##   [7] shinydashboard_0.6.1    splines_3.4.0           digest_0.6.12          
##  [10] htmltools_0.3.6         viridis_0.4.0           gdata_2.18.0           
##  [13] checkmate_1.8.2         memoise_1.1.0           cluster_2.0.6          
##  [16] matrixStats_0.52.2      xts_0.9-7               colorspace_1.3-2       
##  [19] rrcov_1.4-3             RCurl_1.95-4.8          tximport_1.4.0         
##  [22] lme4_1.1-13             survival_2.41-3         zoo_1.8-0              
##  [25] glue_1.1.0              gtable_0.2.0            zlibbioc_1.22.0        
##  [28] XVector_0.16.0          MatrixModels_0.4-1      car_2.1-4              
##  [31] kernlab_0.9-25          prabclus_2.2-6          DEoptimR_1.0-8         
##  [34] SparseM_1.77            VIM_4.7.0               scales_0.4.1           
##  [37] mvtnorm_1.0-6           DBI_0.7                 GGally_1.3.1           
##  [40] Rcpp_0.12.11            sROC_0.1-2              htmlTable_1.9          
##  [43] viridisLite_0.2.0       xtable_1.8-2            laeken_0.4.6           
##  [46] proxy_0.4-17            foreign_0.8-68          mclust_5.3             
##  [49] Formula_1.2-1           DT_0.2                  htmlwidgets_0.8        
##  [52] FNN_1.1                 fpc_2.1-10              acepack_1.4.1          
##  [55] modeltools_0.2-21       reshape_0.8.6           XML_3.98-1.8           
##  [58] flexmix_2.3-14          nnet_7.3-12             locfit_1.5-9.1         
##  [61] labeling_0.3            rlang_0.1.1             munsell_0.4.3          
##  [64] cellranger_1.1.0        tools_3.4.0             RSQLite_1.1-2          
##  [67] pls_2.6-0               evaluate_0.10           cvTools_0.3.2          
##  [70] yaml_2.1.14             robustbase_0.92-7       caTools_1.17.1         
##  [73] nlme_3.1-131            mime_0.5                quantreg_5.33          
##  [76] biomaRt_2.32.1          compiler_3.4.0          pbkrtest_0.4-7         
##  [79] beeswarm_0.2.3          e1071_1.6-8             statmod_1.4.30         
##  [82] smoother_1.1            tibble_1.3.3            robCompositions_2.0.3  
##  [85] pcaPP_1.9-61            stringi_1.1.5           highr_0.6              
##  [88] lattice_0.20-35         trimcluster_0.1-2       nloptr_1.0.4           
##  [91] lmtest_0.9-35           cowplot_0.7.0           bitops_1.0-6           
##  [94] httpuv_1.3.3            latticeExtra_0.6-28     R6_2.2.2               
##  [97] KernSmooth_2.23-15      gridExtra_2.2.1         vipor_0.4.5            
## [100] gtools_3.5.0            boot_1.3-19             MASS_7.3-47            
## [103] assertthat_0.2.0        rhdf5_2.20.0            rprojroot_1.2          
## [106] rjson_0.2.15            GenomeInfoDbData_0.99.0 diptest_0.75-7         
## [109] mgcv_1.8-17             hms_0.3                 rpart_4.1-11           
## [112] class_7.3-14            minqa_1.2.4             TTR_0.23-1             
## [115] shiny_1.0.3             base64enc_0.1-3         ggbeeswarm_0.5.3
```

