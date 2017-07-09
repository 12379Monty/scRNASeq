### A step-by-step workflow for low-level analysis of single-cell RNA-seq data 

<!-- Part I - Analysis of haematopoietic stem cells -->


<style type="text/css"> body, td { font-size: 14px; } code.r{ font-size: 12px; } pre { font-size: 12px } </style>

[bioconductor workflows simpleSingleCell](https://www.bioconductor.org/help/workflows/simpleSingleCell/)



<!-- ***************************************************** -->


<!-- ***************************************************** -->


<!-- ***************************************************** -->



<!-- ***************************************************** -->


## Analysis of haematopoietic stem cells

### Overview

To introduce most of the concepts of scRNA-seq data analysis, we use a
relatively simple dataset from a study of haematopoietic stem cells (HSCs)
(N. K. Wilson et al. 2015). Single mouse HSCs were isolated into microtiter
plates and libraries were prepared for 96 cells using the Smart-seq2 protocol.
A constant amount of spike-in RNA from the External RNA Controls Consortium (ERCC)
was also added to each cell’s lysate prior to library preparation. High-throughput
sequencing was performed and the expression of each gene was quantified by
counting the total number of reads mapped to its exonic regions.
Similarly, the quantity of each spike-in transcript was measured by counting
the number of reads mapped to the spike-in reference sequences.
Counts for all genes/transcripts in each cell were obtained from 
the NCBI Gene Expression Omnibus (GEO) as a supplementary file
under the accession number
[GSE6153](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE61533).

For simplicity, we forego a description of the read processing steps
required to generate the count matrix, i.e., read alignment and counting into features.
These steps have been described in some detail elsewhere (Love et al. 2015;
Y. Chen, Lun, and Smyth 2016), and are largely the same for bulk and single-cell data.
The only additional consideration is that the spike-in information must be
included in the pipeline. Typically, spike-in sequences can be included as
additional FASTA files during genome index building prior to alignment, while
genomic intervals for both spike-in transcripts and endogenous genes can be
concatenated into a single GTF file prior to counting. For users favouring an
R-based approach to read alignment and counting, we suggest using the methods in the
**Rsubread** package (Liao, Smyth, and Shi 2013; Liao, Smyth, and Shi 2014).
Alternatively, rapid quantification of expression with alignment-free methods such
as kallisto (Bray et al. 2016) or Salmon (Patro, Duggal, and Kingsford 2015)
can be performed using the functions **runKallisto** and **runSalmon**
in the **scater** package.

### Count loading

The first task is to load the count matrix into memory. In this case,
some work is required to retrieve the data from the Gzip-compressed Excel format.
Each row of the matrix represents an endogenous gene or a spike-in transcript,
and each column represents a single HSC. For convenience, the counts
for spike-in transcripts and endogenous genes are stored in a
**SCESet** object from the **scater** package (McCarthy et al. 2016).


```r
 ## ---- eval=on.bioc, echo=FALSE, results='hide'---------------------------
 all.urls <- c("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE61533&format=file&file=GSE61533%5FHTSEQ%5Fcount%5Fresults%2Exls%2Egz", 
"https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE29087&format=file&file=GSE29087%5FL139%5Fexpression%5Ftab%2Etxt%2Egz",
"https://storage.googleapis.com/linnarsson-lab-www-blobs/blobs/cortex/expression_mRNA_17-Aug-2014.txt",
"https://storage.googleapis.com/linnarsson-lab-www-blobs/blobs/cortex/expression_mito_17-Aug-2014.txt",
"https://storage.googleapis.com/linnarsson-lab-www-blobs/blobs/cortex/expression_spikes_17-Aug-2014.txt",
"http://www.ebi.ac.uk/teichmann-srv/espresso/static/counttable_es.csv", 
"http://www.nature.com/nbt/journal/v33/n2/extref/nbt.3102-S7.xlsx")

 all.basenames <- basename(all.urls)
 all.basenames[1] <- "GSE61533_HTSEQ_count_results.xls.gz"
 all.basenames[2] <- "GSE29087_L139_expression_tab.txt.gz"
 all.modes <- rep("w", length(all.urls))
 all.modes[!grepl("(txt|csv)$", all.basenames)] <- "wb"

 for (x in seq_along(all.urls)) { 
     download.file(all.urls[x], file.path(EXT_DATA, all.basenames[x]), mode=all.modes[x])
 }
```
<!-- ***************************************************** -->


```r
 #### Was able to make this work...
 #### debug/undebug(calculateQCMetrics) fixed 'mvoutlier' problem!?#


 gunzip(file.path(EXT_DATA,"GSE61533_HTSEQ_count_results.xls.gz"), 
        remove=FALSE, overwrite=TRUE)

 GSE61533.Counts.frm <- as.data.frame(
   read_excel(file.path(EXT_DATA,'GSE61533_HTSEQ_count_results.xls'), sheet=1))
 rownames(GSE61533.Counts.frm) <- GSE61533.Counts.frm$ID
 GSE61533.Counts.frm <- GSE61533.Counts.frm[,-1]

 print(kable(GSE61533.Counts.frm[1:5, 1:5]))
 
 # Save this
 save(GSE61533.Counts.frm, file=file.path(WRKDIR,'Data', 'GSE61533.Counts.frm'))

 GSE61533.sce <- newSCESet(countData=GSE61533.Counts.frm)
 dim(GSE61533.sce)

 ## find spike-ins and mitochonrdial sequences
 is.spike <- grepl("^ERCC", rownames(GSE61533.sce))
 is.mito <- grepl("^mt-", rownames(GSE61533.sce))

 ## Get QC
 GSE61533.sce <- calculateQCMetrics(GSE61533.sce, 
                 feature_controls=list(ERCC=is.spike, Mt=is.mito))
 setSpike(GSE61533.sce) <- "ERCC"

 head(colnames(pData(GSE61533.sce)))

 # save
 save(GSE61533.sce, file.path(WRKDIR, 'Data', 'GSE61533.sce'))
```
<!-- **************************************************************** -->

### Quality control on the cells

Low-quality cells need to be removed to ensure that technical effects do not
distort downstream analysis results. Two common measures of cell quality are the 
library size and the number of expressed features in each library.
The library size is defined as the total sum of counts across *all* features, 
i.e., genes and spike-in transcripts. Cells with relatively small library sizes
are considered to be of low quality as the RNA has not been efficiently captured
(i.e., converted into cDNA and amplified) during library preparation.
The number of expressed features in each cell is defined as the number of
features with non-zero counts for that cell. Any cell with very few expressed genes
is likely to be of poor quality as the diverse transcript population has not
been successfully captured. The distributions of both of these metrics are shown 
in Figure 1.

 
<!-- 
 ## ----libplothsc, fig.height=6, fig.width=12, 
fig.cap="**Figure 1:** Histograms of library sizes (left) and number of expressed genes (right) for all cells in the HSC dataset."----
-->



<!-- ****************************************************** -->

<!--
#### Get results w/o having relying (mostly) on obscure and convoluted methods of the scater package
-->


```r
 load(file.path(WRKDIR, 'Data', 'GSE61533.Counts.frm'))

 Ctl.flg <- grepl("^ERCC", rownames(GSE61533.Counts.frm)) | 
            grepl("^mt-", rownames(GSE61533.Counts.frm))

 GSE61533.total_counts.vec <- colSums(GSE61533.Counts.frm)

 load(file.path(WRKDIR, 'Data', 'GSE61533.sce'))
 lowerDetectionLimit <- GSE61533.sce@lowerDetectionLimit
 GSE61533.expr.mtx <- exprs(GSE61533.sce)

 GSE61533.total_features.vec <- colSums(GSE61533.expr.mtx[!Ctl.flg,] > lowerDetectionLimit)
 
 par(mfrow=c(1,2),oma=c(0,0,2,0))

 hist(GSE61533.total_counts.vec/1e6, xlab="Library sizes (millions)", main="", 
    breaks=20, col="grey80", ylab="Number of cells")

 hist(GSE61533.total_features.vec, xlab="Number of expressed genes", main="", 
    breaks=20, col="grey80", ylab="Number of cells")

 mtext(side=3, outer=T, 
  "Fig. 1a: Library sizes (left) and number of expressed genes (right) for all cells in the HSC dataset.")
```

![plot of chunk libplothsc2](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M1A/libplothsc2-1.png)
<!-- ****************************************************** -->

Picking a threshold for these metrics is not straightforward as their absolute
values depend on the protocol and biological system. For example, sequencing
to greater depth will lead to more reads, regardless of the quality of the cells.
To obtain an adaptive threshold, we assume that most of the dataset consists of
high-quality cells. We remove cells with log-library sizes that are more than
3 median absolute deviations (MADs) below the median log-library size.
(A log-transformation improves resolution at small values, especially
when the MAD of the raw values is comparable to or greater than the median.)
We also remove cells where the log-transformed number of expressed
genes is 3 MADs below the median. 


Another measure of quality is the proportion of reads mapped to genes in the
mitochondrial genome. High proportions are indicative of poor-quality cells
(Islam et al. 2014; Ilicic et al. 2016), possibly because of increased apoptosis
and/or loss of cytoplasmic RNA from lysed cells. Similar reasoning applies to the
proportion of reads mapped to spike-in transcripts. The quantity of spike-in
RNA added to each cell should be constant, which means that the proportion should
increase upon loss of endogenous RNA in low-quality cells.
The distributions of mitochondrial and spike-in proportions across all 
cells are shown in Figure 2.

<!-- 
## ----controlplothsc, fig.width=12, fig.height=6, fig.cap="**Figure 2:** Histogram of the proportion of reads mapped to mitochondrial genes (left) or spike-in transcripts (right) across all cells in the HSC dataset."----
-->



<!--
#### Get results w/o relying (mostly) on obscure and convoluted methods of the scater package
-->


```r
 load(file.path(WRKDIR, 'Data', 'GSE61533.Counts.frm'))

 ERCC.flg <- grepl("^ERCC", rownames(GSE61533.Counts.frm)) 
 Mt.flg <- grepl("^mt-", rownames(GSE61533.Counts.frm))

 pctCountsFeatureCtl_Mt.vec <- 100*colSums(GSE61533.Counts.frm[Mt.flg,])/
                                 colSums(GSE61533.Counts.frm)
 pctCountsFeatureCtl_ERCC.vec <- 100*colSums(GSE61533.Counts.frm[ERCC.flg,])/
                                 colSums(GSE61533.Counts.frm)

 par(mfrow=c(1,2),oma=c(0,0,2,0))
 hist(pctCountsFeatureCtl_Mt.vec, xlab="Mitochondrial proportion (%)",
     ylab="Number of cells", breaks=20, main="", col="grey80")
 hist(pctCountsFeatureCtl_ERCC.vec, xlab="ERCC proportion (%)",
     ylab="Number of cells", breaks=20, main="", col="grey80")
 mtext(side=3, outer=T,
 "Fig. 2a: Proportion of reads mapped to mitochondrial genes (left) or spike-in transcripts (right) across all cells in the HSC dataset.")
```

![plot of chunk controlplothsc2](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M1A/controlplothsc2-1.png)
<!-- ****************************************************** -->

Again, the ideal threshold for these proportions depends on the cell
type and the experimental protocol. Cells with more mitochondria or more mitochondrial
activity may naturally have larger mitochondrial proportions. Similarly, cells with more
endogenous RNA or that are assayed with protocols using less spike-in RNA will have
lower spike-in proportions. If we assume that most cells in the dataset are of
high quality, then the threshold can be set to remove any large outliers from the
distribution of proportions. We use the MAD-based definition of outliers to remove
putative low-quality cells from the dataset.


<!-- ****************************************************** -->

<!-- ### Repeat w/o ...
-->

```r
 libsize.drop <- (log10(GSE61533.total_counts.vec) <
   median(log10(GSE61533.total_counts.vec) - 3*mad(log10(GSE61533.total_counts.vec) )))
 #colnames(GSE61533.sce)[libsize.drop]

 feature.drop <- (log10(GSE61533.total_features.vec) <
   median(log10(GSE61533.total_features.vec) - 3*mad(log10(GSE61533.total_features.vec) )))
 #colnames(GSE61533.sce)[feature.drop]

 mito.drop <- (pctCountsFeatureCtl_Mt.vec >
   median(pctCountsFeatureCtl_Mt.vec + 3*mad(pctCountsFeatureCtl_Mt.vec )))
 #colnames(GSE61533.sce)[mito.drop]

 spike.drop <- (pctCountsFeatureCtl_ERCC.vec >
   median(pctCountsFeatureCtl_ERCC.vec + 3*mad(pctCountsFeatureCtl_ERCC.vec )))
 #colnames(GSE61533.sce)[spike.drop]
```
<!-- ****************************************************** -->

Subsetting by column will retain only the high-quality cells that pass each
filter described above. We examine the number of cells removed by each filter as
well as the total number of retained cells. Removal of a substantial proportion of
cells (> 10%) may be indicative of an overall issue with data quality.
It may also reflect genuine biology in extreme cases (e.g., low numbers of
expressed genes in erythrocytes) for which the filters described here are
inappropriate.



```r
GSE61533.sce <- GSE61533.sce[,!(libsize.drop | feature.drop | mito.drop | spike.drop)]
cat('\nDropped Samples:\n')
```

```
## 
## Dropped Samples:
```

```r
data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop),
    ByMito=sum(mito.drop), BySpike=sum(spike.drop), Remaining=ncol(GSE61533.sce))
```

```
##         ByLibSize ByFeature ByMito BySpike Remaining
## Samples         2         2      6       3        86
```

An alternative approach to quality control is to perform a principal components
analysis (PCA) based on the quality metrics for each cell, e.g., 
the total number of reads, the total number of features and the proportion of
mitochondrial or spike-in reads. Outliers on a PCA plot may be indicative of
low-quality cells that have aberrant technical properties compared to the (presumed)
majority of high-quality cells. In Figure 3, no obvious outliers are present,
which is consistent with the removal of suspect cells in the preceding
quality control steps.

<!-- 
## ----pcaqualplothsc, fig.cap="**Figure 3:** PCA plot for cells in the HSC dataset, constructed using quality metrics. The first and second components are shown on each axis, along with the percentage of total variance explained by each component. Bars represent the coordinates of the cells on each axis."----
-->


```r
fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))
plotPCA(GSE61533.sce, pca_data_input="pdata") + fontsize +
ggtitle("Fig. 3: PCA plot for cells in the HSC dataset, constructed using quality metrics.")
```

![plot of chunk qcPCA](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M1A/qcPCA-1.png)

Methods like PCA-based outlier detection and support vector machines can
provide more power to distinguish low-quality cells from high-quality counterparts
(Ilicic et al. 2016). This is because they are able to detect subtle patterns across
many quality metrics simultaneously. However, this comes at some cost to
interpretability, as the reason for removing a given cell may not always be obvious.
Thus, for this workflow, we will use the simple approach whereby each quality metric
is considered separately. Users interested in the more sophisticated
approaches are referred to the **scater** and **cellity** packages.

### Classification of cell cycle phase

We use the prediction method described by Scialdone et al. (2015) to classify
cells into cell cycle phases based on the gene expression data. Using a training dataset,
the sign of the difference in expression between two genes was computed for each pair
of genes. Pairs with changes in the sign across cell cycle phases were chosen as markers.
Cells in a test dataset can then be classified into the appropriate phase,
based on whether the observed sign for each marker pair is consistent
with one phase or another. This approach is implemented in the **cyclone**
function using a pre-trained set of marker pairs for mouse data.
The result of phase assignment for each cell in the HSC dataset is shown in Figure 4.
(Some additional work is necessary to match the gene symbols in
the data to the Ensembl annotation in the pre-trained marker set.)


<!-- 
## ----phaseplothsc, message=FALSE, fig.cap="**Figure 4:** Cell cycle phase scores from applying the pair-based classifier on the HSC dataset, where each point represents a cell."----
-->


```r
 set.seed(100)

 mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))

 library(org.Mm.eg.db)
 anno <- select(org.Mm.eg.db, keys=rownames(GSE61533.sce), keytype="SYMBOL", column="ENSEMBL")
 ensembl <- anno$ENSEMBL[match(rownames(GSE61533.sce), anno$SYMBOL)]

 assignments <- cyclone(GSE61533.sce, mm.pairs, gene.names=ensembl)

 plot(assignments$score$G1, assignments$score$G2M, xlab="G1 score", ylab="G2/M score", pch=16)
 abline(0,1,col='red')
 abline(h=0.5, v=0.5, col='blue', lty=2)
 title("Fig. 4: Cell cycle phase scores where each point represents a cell.")
```

![plot of chunk cellCycleScores](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M1A/cellCycleScores-1.png)

Cells are classified as being in G1 phase if the G1 score is above 0.5 and 
greater than the G2/M score; in G2/M phase if the G2/M score is above 0.5 and
greater than the G1 score; and in S phase if neither score is above 0.5.
Here, the vast majority of cells are classified as being in G1 phase.
We will focus on these cells in the downstream analysis. Cells in other phases are
removed to avoid potential confounding effects from cell cycle-induced differences.
Alternatively, if a non-negligible number of cells are in other phases,
we can use the assigned phase as a blocking factor in downstream analyses.
This protects against cell cycle effects without discarding information.



```r
 GSE61533_G1.sce <- GSE61533.sce[,assignments$phases=="G1"]
```
Pre-trained classifiers are available in **scran** for human and mouse data.
While the mouse classifier used here was trained on data from embryonic stem cells,
it is still accurate for other cell types (Scialdone et al. 2015).
This may be due to the conservation of the transcriptional program associated
with the cell cycle (Bertoli, Skotheim, and Bruin 2013; Conboy et al. 2007).
The pair-based method is also a non-parametric procedure that is robust to most
technical differences between datasets. However, it will be less accurate for
data that are substantially different from those used in the training set, e.g.,
due to the use of a different protocol. In such cases, users can construct a custom
classifier from their own training data using the **sandbag** function.
This will also be necessary for other model organisms where pre-trained
classifiers are not available.

### Filtering out low-abundance genes

Low-abundance genes are problematic as zero or near-zero counts do not
contain enough information for reliable statistical inference
(Bourgon, Gentleman, and Huber 2010). In addition, the discreteness of the counts
may interfere with downstream statistical procedures, e.g., by compromising the
accuracy of continuous approximations. Here, low-abundance genes are defined
as those with an average count below a filter threshold of 1.
These genes are likely to be dominated by drop-out events (Brennecke et al. 2013),
which limits their usefulness in later analyses. Removal of these genes mitigates
discreteness and reduces the amount of computational work without major loss of
information.


```r
 ave.counts <- calcAverage(GSE61533_G1.sce)
 keep <- ave.counts >= 1
 sum(keep)
```

```
## [1] 13977
```

To check whether the chosen threshold is suitable, we examine the
distribution of log-means across all genes (Figure 5).
The peak represents the bulk of moderately expressed genes while the rectangular
component corresponds to lowly expressed genes. The filter threshold should cut the
distribution at some point along the rectangular component to remove the majority
of low-abundance genes.

<!--
## ----abhisthsc, fig.cap="**Figure 5:** Histogram of log-average counts for all genes in the HSC dataset. The filter threshold is represented by the blue line."----
-->


```r
 hist(log10(ave.counts), breaks=100, main="", col="grey80",
     xlab=expression(Log[10]~"average count"))
 abline(v=log10(1), col="blue", lwd=2, lty=2)
 title(paste("Fig. 5: of log-average counts for all genes in the HSC dataset.",
  "\nThe filter threshold is represented by the blue line."))
```

![plot of chunk histLogAveCount](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M1A/histLogAveCount-1.png)

We also look at the identities of the most highly expressed genes (Figure 6).
This should generally be dominated by constitutively expressed transcripts, such as
those for ribosomal or mitochondrial proteins. The presence of other classes of 
features may be cause for concern if they are not consistent with expected biology.
For example, a top set containing many spike-in transcripts suggests that too much
spike-in RNA was added during library preparation, while the absence of ribosomal
proteins and/or the presence of their pseudogenes are indicative of
suboptimal alignment.

<!--
## ----topgenehsc, fig.height=9, fig.width=6, fig.cap="**Figure 6:** Percentage of total counts assigned to the top 50 most highly-abundant features in the HSC dataset. For each feature, each bar represents the percentage assigned to that feature for a single cell, while the circle represents the average across all cells. Bars are coloured by the total number of expressed features in each cell, while circles are coloured according to whether the feature is labelled as a control feature."----
-->

#### Figure 6: Percentage of total counts assigned to the top 50 most highly-abundant features in the HSC dataset

For each feature, each bar represents the percentage assigned to that feature for 
a single cell, while the circle represents the average across all cells.
Bars are coloured by the total number of expressed features in each cell,
while circles are coloured according to whether the feature is labelled as
a control feature.


```r
 plotQC(GSE61533_G1.sce, type = "highest-expression", n=50) + fontsize
```

![plot of chunk top50Features](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M1A/top50Features-1.png)

An alternative approach to gene filtering is to select genes that have
non-zero counts in at least n cells. This provides some more protection against genes
with outlier expression patterns, i.e., strong expression in only one or two cells. 
Such outliers are typically uninteresting as they can arise from 
*amplification artifacts* that are not replicable across cells.
(The exception is for studies involving rare cells where the outliers
may be biologically relevant.) An example of this filtering approach is
shown below for n set to 10, though smaller values may be necessary
to retain genes expressed in rare cell types.


```r
 numcells <- nexprs(GSE61533_G1.sce, byrow=TRUE)
 alt.keep <- numcells >= 10
 sum(alt.keep)
```

```
## [1] 11988
```

The relationship between the number of expressing cells and the mean is shown in
Figure 7. The two statistics tend to be well-correlated so filtering on either
should give roughly similar results.
<!--
## ----geneplothsc, fig.cap="**Figure 7:** Number of expressing cells against the log-mean expression for each gene in the HSC dataset. Spike-in transcripts are highlighted in red."----
-->

```r
 smoothScatter(log10(ave.counts), numcells, xlab=expression(Log[10]~"average count"), 
     ylab="Number of expressing cells")
 is.ercc <- isSpike(GSE61533_G1.sce, type="ERCC")
 points(log10(ave.counts[is.ercc]), numcells[is.ercc], col="red", pch=16, cex=0.5)

 title(paste("Fig.  7: Number of expressing cells against the log-mean expression",
   "\nSpike-in transcripts are highlighted in red."))
```

![plot of chunk numExprPlot](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M1A/numExprPlot-1.png)

In general, we prefer the mean-based filter as it tends to be less aggressive.
A gene will be retained as long as it has sufficient expression in **any subset of cells**
Genes expressed in fewer cells require higher levels of expression in those cells to be
retained, but this is not undesirable as it avoids selecting uninformative genes 
(with low expression in few cells) that contribute little to downstream analyses,
e.g., HVG detection or clustering. In contrast, the **“at least n”** filter depends heavily
on the choice of n. With n = 10, a gene expressed in a subset of 9 cells would be
filtered out, regardless of the level of expression in those cells.
This may result in the failure to detect rare subpopulations that are present
at frequencies below n. While the mean-based filter will retain more outlier-driven genes,
this can be handled by choosing methods that are robust to outliers
in the downstream analyses.

Thus, we apply the mean-based filter to the data by subsetting the SCESet object
as shown below. This removes all rows corresponding to endogenous genes or
spike-in transcripts with abundances below the specified threshold.


```r
 GSE61533_G1.sce <- GSE61533_G1.sce[keep,] 
```

### Normalization of cell-specific biases

#### Using the deconvolution method to deal with zero counts

Read counts are subject to differences in capture efficiency and sequencing depth
between cells (Stegle, Teichmann, and Marioni 2015). Normalization is required
to eliminate these cell-specific biases prior to downstream quantitative analyses.
This is often done by assuming that most genes are not differentially expressed (DE)
between cells. Any systematic difference in count size across the non-DE majority
of genes between two cells is assumed to represent bias and is removed by scaling.
More specifically, “size factors” are calculated that represent the extent to
which counts should be scaled in each library.

Size factors can be computed with several different approaches, e.g., using
the **estimateSizeFactorsFromMatrix** function in the **DESeq2** package
(Anders and Huber 2010; Love, Huber, and Anders 2014), or with the 
**calcNormFactors** function (Robinson and Oshlack 2010) in the **edgeR** package.
However, single-cell data can be problematic for these bulk data-based methods
due to the dominance of low and zero counts. To overcome this, we pool counts from 
many cells to increase the count size for accurate size factor estimation
(Lun, Bach, and Marioni 2016). Pool-based size factors are then *“deconvolved”*
into cell-based factors for cell-specific normalization.


```r
 GSE61533_G1.sce <- computeSumFactors(GSE61533_G1.sce, sizes=seq(20, 80, 5))
 summary(sizeFactors(GSE61533_G1.sce))
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.4142  0.8159  0.9316  1.0000  1.1827  1.9895
```

In this case, the size factors are tightly correlated with the library sizes
for all cells (Figure 8). This suggests that the systematic differences between
cells are primarily driven by differences in capture efficiency or sequencing depth.
Any DE between cells would yield a non-linear trend between the total count and size factor,
and/or increased scatter around the trend. This does not occur here as strong DE is
unlikely to exist within a homogeneous population of cells.

<!--
## ----normplothsc, fig.cap="**Figure 8:** Size factors from deconvolution, plotted against library sizes for all cells in the HSC dataset. Axes are shown on a log-scale."----
-->


```r
plot(sizeFactors(GSE61533_G1.sce), GSE61533_G1.sce$total_counts/1e6, log="xy",
    ylab="Library size (millions)", xlab="Size factor")
title(paste("Fig. 8: Size factors from deconvolution, plotted against library sizes",
   "\n Axes are shown on a log-scale."))
```

![plot of chunk plotSizeFactors](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M1A/plotSizeFactors-1.png)

#### Computing separate size factors for spike-in transcripts


Size factors computed from the counts for endogenous genes are usually not appropriate
for normalizing the counts for spike-in transcripts. Consider an experiment without
library quantification, i.e., the amount of cDNA from each library is not equalized prior
to pooling and multiplexed sequencing. Here, cells containing more RNA have greater counts
for endogenous genes and thus larger size factors to scale down those counts. However,
the same amount of spike-in RNA is added to each cell during library preparation.
This means that the counts for spike-in transcripts are not subject to the effects of
RNA content. Attempting to normalize the spike-in counts with the gene-based size factors
will lead to over-normalization and incorrect quantification of expression. Similar
reasoning applies in cases where library quantification is performed. For a constant
total amount of cDNA, any increases in endogenous RNA content will suppress the coverage of
spike-in transcripts. As a result, the bias in the spike-in counts will be opposite to
that captured by the gene-based size factor.

To ensure normalization is performed correctly, we compute a separate set of size
factors for the spike-in set. For each cell, the spike-in-specific size factor is
defined as the total count across all transcripts in the spike-in set.
This assumes that none of the spike-in transcripts are differentially expressed,
which is reasonable given that the same amount and composition of spike-in RNA should
have been added to each cell. (See below for a more detailed discussion on spike-in
normalization.) These size factors are stored in a separate field of the **SCESet**
object by setting **general.use=FALSE** in **computeSpikeFactors**.
This ensures that they will only be used with the spike-in transcripts but
not the endogenous genes.


```r
 GSE61533_G1.sce <- computeSpikeFactors(GSE61533_G1.sce, type="ERCC", general.use=FALSE)
```
#### Applying the size factors to normalize gene expression

The count data are used to compute normalized log-expression values for use in
downstream analyses. Each value is defined as the log-ratio of each count to the
size factor for the corresponding cell, after adding a prior count of 1 to
avoid undefined values at zero counts. Division by the size factor ensures that any
cell-specific biases are removed. If spike-in-specific size factors are present
in **sce**, they will be automatically applied to normalize the spike-in
transcripts separately from the endogenous genes.



```r
 GSE61533_G1.sce <- normalize(GSE61533_G1.sce)
```

The log-transformation provides some measure of variance stabilization (Law et al. 2014),
so that high-abundance genes with large variances do not dominate downstream analyses.
The computed values are stored as an **exprs** matrix in addition to the other assay elements.


#### Checking for important technical factors

We check whether there are technical factors that contribute substantially to the heterogeneity of gene expression. If so, the factor may need to be regressed out to ensure that it does not inflate the variances or introduce spurious correlations. For this dataset, the simple experimental design means that there are no plate or batch effects to examine. Instead, we use the (log-transformed) total count for the spike-in transcripts as a proxy for the relative bias in each sample. This bias is purely technical in origin, given that the same amount of spike-in RNA should have been added to each cell. Thus, any association of gene expression with this factor is not biologically interesting and should be removed.

For each gene, we calculate the percentage of the variance of the expression values that is explained by the spike-in totals (Figure 9). The percentages are generally small (1-3%), indicating that the expression profiles of most genes are not strongly associated with this factor. This result is consistent with successful removal of cell-specific biases by scaling normalization. Thus, the spike-in total does not need to be explicitly modelled in our downstream analyses.



#### Figure 9: 
* Density plot of the % of variance explained by the (log-transformed) total spike-in counts across all genes For each gene, the percentage of the variance of the normalized log-expression values across cells that is explained by each factor is calculated. Each curve corresponds to one factor and represents the distribution of percentages across all genes.




```r
 plotExplanatoryVariables(GSE61533_G1.sce, variables=c("counts_feature_controls_ERCC", 
     "log10_counts_feature_controls_ERCC")) + fontsize
```

![plot of chunk plotExplanatoryVariables](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M1A/plotExplanatoryVariables-1.png)

Note that the use of the spike-in total as an accurate proxy for the relative technical bias assumes that no library quantification is performed. Otherwise, the coverage of the spike-in transcripts would be dependent on the total amount of endogenous RNA in each cell. (Specifically, if the same amount of cDNA is used for sequencing per cell, any increase in the amount of endogenous RNA will suppress the coverage of the spike-in transcripts.) This means that the spike-in totals could be confounded with genuine biological effects associated with changes in RNA content.

### Identifying HVGs from the normalized log-expression

We identify HVGs to focus on the genes that are driving heterogeneity across the population of cells. This requires estimation of the variance in expression for each gene, followed by decomposition of the variance into biological and technical components. HVGs are then identified as those genes with the largest biological components. This avoids prioritizing genes that are highly variable due to technical factors such as sampling noise during RNA capture and library preparation.

Ideally, the technical component would be estimated by fitting a mean-variance trend to the spike-in transcripts using the **trendVar** function. Recall that the same set of spike-ins was added in the same quantity to each cell. This means that the spike-in transcripts should exhibit no biological variability, i.e., any variance in their counts should be technical in origin. Given the mean abundance of a gene, the fitted value of the trend can be used as an estimate of the technical component for that gene. The biological component of the variance can then be calculated by subtracting the technical component from the total variance of each gene with the **decomposeVar** function.

In practice, this strategy is compromised by the small number of spike-in transcripts, the uneven distribution of their abundances and (for low numbers of cells) the imprecision of their variance estimates. This makes it difficult to accurately fit a complex mean-dependent trend to the spike-in variances. An alternative approach is to fit the trend to the variance estimates of the endogenous genes, using the **use.spikes=FALSE** setting as shown below. This assumes that the majority of genes are not variably expressed, such that the technical component dominates the total variance for those genes. The fitted value of the trend is then used as an estimate of the technical component. Obviously, this is the only approach that can be used if no spike-ins were added in the experiment.



```r
 var.fit <- trendVar(GSE61533_G1.sce, trend="loess", use.spikes=FALSE, span=0.2)
 var.out <- decomposeVar(GSE61533_G1.sce, var.fit)
```

We assess the suitability of the trend fitted to the endogenous variances by examining whether it is consistent with the spike-in variances (Figure 10). The trend passes through or close to most of the spike-in variances, indicating that our assumption (that most genes have low levels of biological variability) is valid. This strategy exploits the large number of endogenous genes to obtain a stable trend, with the spike-in transcripts used as diagnostic features rather than in the trend fitting itself. However, if our assumption did not hold, we would instead fit the trend directly to the spike-in variances with the default **use.spikes=TRUE**. This sacrifices stability to reduce systematic errors in the estimate of the biological component for each gene. (In such cases, tinkering with the trend fitting parameters may yield a more stable curve – see **?trendVar** for more details.)


#### Figure 10:
* Variance of normalized log-expression values for each gene in the HSC dataset, plotted against the mean log-expression. The blue line represents the mean-dependent trend fitted to the variances of the endogenous genes. Variance estimates for spike-in transcripts are highlighted in red.


```r
plot(var.out$mean, var.out$total, pch=16, cex=0.6, xlab="Mean log-expression", 
    ylab="Variance of log-expression")
o <- order(var.out$mean)
lines(var.out$mean[o], var.out$tech[o], col="dodgerblue", lwd=2)
cur.spike <- isSpike(GSE61533_G1.sce)
points(var.out$mean[cur.spike], var.out$total[cur.spike], col="red", pch=16)
```

![plot of chunk hvgplothsc](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M1A/hvgplothsc-1.png)

HVGs are defined as genes with biological components that are significantly greater than zero at a false discovery rate (FDR) of 5%. These genes are interesting as they drive differences in the expression profiles between cells, and should be prioritized for further investigation. In addition, we only consider a gene to be a HVG if it has a biological component greater than or equal to 0.5. For transformed expression values on the log2 scale, this means that the average difference in true expression between any two cells will be at least 2-fold. (This reasoning assumes that the true log-expression values are Normally distributed with variance of 0.5. The root-mean-square of the difference between two values is treated as the average log2-fold change between cells and is equal to unity.) We rank the results by the biological component to focus on genes with larger biological variability.


```r
 hvg.out <- var.out[which(var.out$FDR <= 0.05 & var.out$bio >= 0.5),]
 hvg.out <- hvg.out[order(hvg.out$bio, decreasing=TRUE),] 
 nrow(hvg.out)
```

```
## [1] 145
```

```r
 write.table(file=file.path(TABLE_DIR,"hsc_hvg.tsv"),
             hvg.out, sep="\t", quote=FALSE, col.names=NA)
 head(hvg.out)
```

```
##              mean     total       bio      tech      p.value          FDR
## Fos      6.411763 20.184481 12.080744  8.103738 2.060249e-12 1.304512e-09
## Rgs1     5.214931 20.269800  9.118208 11.151592 8.631784e-06 1.431438e-03
## Dusp1    6.694369 16.081746  8.860774  7.220972 1.252186e-09 4.714312e-07
## Ctla2a   8.653032  9.508748  7.344402  2.164347 1.553494e-36 3.327867e-33
## Ppp1r15a 6.544443 14.946762  7.229263  7.717498 7.781068e-07 1.642277e-04
## Sult1a1  5.612108 17.383183  7.015586 10.367597 1.211462e-04 1.339339e-02
```

We recommend checking the distribution of expression values for the top HVGs to ensure that the variance estimate is not being dominated by one or two outlier cells (Figure 11).


#### Figure 11:
*  Violin plots of normalized log-expression values for the top 10 HVGs in the HSC dataset. Each point represents the log-expression value in a single cell.


```r
plotExpression(GSE61533_G1.sce, rownames(hvg.out)[1:10]) + fontsize
```

![plot of chunk hvgvioplothsc](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M1A/hvgvioplothsc-1.png)

There are many other strategies for defining HVGs, e.g., by using the coefficient of variation (Brennecke et al. 2013; Kołodziejczyk et al. 2015; Kim et al. 2015), with the dispersion parameter in the negative binomial distribution (McCarthy, Chen, and Smyth 2012), or as a proportion of total variability (Vallejos, Marioni, and Richardson 2015). Some of these methods are available in **scran** – for example, see **DM** or **technicalCV2** for calculations based on the coefficient of variation. Here, we use the variance of the log-expression values because the log-transformation protects against genes with strong expression in only one or two cells. This ensures that the set of top HVGs is not dominated by genes with (mostly uninteresting) outlier expression patterns.

### Identifying correlated gene pairs with Spearman’s rho

Another useful procedure is to identify the HVGs that are highly correlated with one another. This distinguishes between HVGs caused by random noise and those involved in driving systematic differences between subpopulations. Correlations between genes are quantified by computing Spearman’s rho, which accommodates non-linear relationships in the expression values. Gene pairs with significantly large positive or negative values of rho are identified using the **correlatePairs** function. We only apply this function to the set of HVGs, because these genes have large biological components and are more likely to exhibit strong correlations driven by biology. In contrast, calculating correlations for all possible gene pairs would require too much computational time and increase the severity of the multiple testing correction. It may also prioritize uninteresting genes that have strong correlations but low variance, e.g., tightly co-regulated house-keeping genes.


```r
set.seed(100)
var.cor <- correlatePairs(GSE61533_G1.sce, subset.row=rownames(hvg.out))
write.table(file="hsc_cor.tsv", var.cor, sep="\t", quote=FALSE, row.names=FALSE)
head(var.cor)
```

```
##      gene1   gene2       rho      p.value         FDR limited
## 1   mt-Nd2 mt-Rnr1 0.6068805 1.999998e-06 0.002982854    TRUE
## 2     Egr1     Jun 0.5304143 1.999998e-06 0.002982854    TRUE
## 3    Pdia6   Hspa5 0.5177365 1.999998e-06 0.002982854    TRUE
## 4      Fos    Egr1 0.5107468 1.999998e-06 0.002982854    TRUE
## 5 Ppp1r15a   Zfp36 0.4986357 1.999998e-06 0.002982854    TRUE
## 6    Zfp36    Ier2 0.4929894 1.999998e-06 0.002982854    TRUE
```

The significance of each correlation is determined using a permutation test. For each pair of genes, the null hypothesis is that the expression profiles of two genes are independent. Shuffling the profiles and recalculating the correlation yields a null distribution that is used to obtain a p-value for each observed correlation value (Phipson and Smyth 2010). Correction for multiple testing across many gene pairs is performed by controlling the FDR at 5%. Correlated gene pairs can be directly used for experimental validation with orthogonal techniques (e.g., fluorescence-activated cell sorting, immunohistochemistry or RNA fluorescence in situ hybridization) to verify that these expression patterns are genuinely present across the cell population.




```r
sig.cor <- var.cor$FDR <= 0.05
summary(sig.cor)
```

```
##    Mode   FALSE    TRUE 
## logical   10390      50
```
Larger sets of correlated genes are assembled by treating genes as nodes in a graph and each pair of genes with significantly large correlations as an edge. In particular, an undirected graph is constructed using methods in the **RBGL** package. Highly connected subgraphs are then identified and defined as gene sets. This provides a convenient summary of the pairwise correlations between genes.



```r
 suppressMessages(require(RBGL))
 g <- ftM2graphNEL(cbind(var.cor$gene1, var.cor$gene2)[sig.cor,], 
     W=NULL, V=NULL, edgemode="undirected")
 cl <- highlyConnSG(g)$clusters
 cl <- cl[order(lengths(cl), decreasing=TRUE)]
 head(cl)
```

```
## [[1]]
## [1] "Pdia6"   "Hspd1"   "Phgdh"   "Pik3ip1" "Ncl"    
## 
## [[2]]
## [1] "Egr1"  "Fos"   "Zfp36" "Ier2" 
## 
## [[3]]
## [1] "mt-Nd2"  "Sh3bgrl" "mt-Rnr1"
## 
## [[4]]
## [1] "Vbp1"   "Cacybp" "Hnrnpu"
## 
## [[5]]
## [1] "Sqstm1" "Igfbp4" "Cct3"  
## 
## [[6]]
## [1] "Ppp1r15a" "Junb"
```
Significant correlations provide evidence for substructure in the dataset, i.e., subpopulations of cells with systematic differences in their expression profiles. The number of significantly correlated HVG pairs represents the strength of the substructure. If many pairs were significant, this would indicate that the subpopulations were clearly defined and distinct from one another. For this particular dataset, a relatively low number of HVGs exhibit significant correlations. This suggests that any substructure in the data will be modest, which is expected given that rigorous selection was performed to obtain a homogeneous population of HSCs (N. K. Wilson et al. 2015).

The **correlatePairs** function can also return gene-centric output by setting **per.gene=TRUE**. This calculates a combined p-value (Simes 1986) for each HVG that indicates whether it is significantly correlated to at least one other HVG. From a statistical perspective, this is a more natural approach to correcting for multiple testing when genes, rather than pairs of genes, are of interest. We define correlated HVGs as those that are detected here at a FDR of 5%.


```r
 var.cor <- correlatePairs(GSE61533_G1.sce, subset.row=rownames(hvg.out), per.gene=TRUE)
 head(var.cor)
```

```
##       gene        rho      p.value         FDR limited
## 1      Fos  0.5173376 0.0005759994 0.009279991   FALSE
## 2     Rgs1  0.4236178 0.0129599870 0.059657083   FALSE
## 3    Dusp1  0.4500231 0.0051839948 0.031319969   FALSE
## 4   Ctla2a -0.3899920 0.0478079522 0.116380211   FALSE
## 5 Ppp1r15a  0.5094035 0.0005759994 0.009279991   FALSE
## 6  Sult1a1  0.4107930 0.0204479796 0.078025185   FALSE
```

```r
 sig.cor <- var.cor$FDR <= 0.05
 summary(sig.cor)
```

```
##    Mode   FALSE    TRUE 
## logical     116      29
```

### Using correlated HVGs for further data exploration

We visualize the expression profiles of the correlated HVGs with a heatmap (Figure 12). All expression values are mean-centred for each gene to highlight the relative differences in expression between cells. If any subpopulations were present, they would manifest as rectangular “blocks” in the heatmap, corresponding to sets of genes that are systematically up- or down-regulated in specific groups of cells. This is not observed in Figure 12, consistent with the lack of strong substructure. There may be a subpopulation of Fos and Jun-negative cells, but it is poorly defined given the small numbers of cells and genes involved.


#### Figure 12:
* Heatmap of mean-centred normalized log-expression values for correlated HVGs in the HSC dataset. Dendrograms are formed by hierarchical clustering on the Euclidean distances between genes (row) or cells (column).

```r
 chosen <- var.cor$gene[sig.cor]
 norm.exprs <- exprs(GSE61533_G1.sce)[chosen,,drop=FALSE]
 heat.vals <- norm.exprs - rowMeans(norm.exprs)

 suppressMessages(require(gplots))
 heat.out <- heatmap.2(heat.vals, col=bluered, symbreak=TRUE, trace='none', cexRow=0.6)
```

![plot of chunk heatmap](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M1A/heatmap-1.png)

We also apply dimensionality reduction techniques to visualize the relationships between cells. This is done by constructing a PCA plot from the normalized log-expression values of the correlated HVGs (Figure 13). Cells with similar expression profiles should be located close together in the plot, while dissimilar cells should be far apart. We only use the correlated HVGs in **plotPCA** because any substructure should be most pronounced in the expression profiles of these genes. Even so, no clear separation of cells into distinct subpopulations is observed.


#### Figure 13:
* PCA plot constructed from normalized log-expression values of correlated HVGs, where each point represents a cell in the HSC dataset. First and second components are shown, along with the percentage of variance explained. Bars represent the coordinates of the cells on each axis. Each cell is coloured according to its total number of expressed features


```r
 plotPCA(GSE61533_G1.sce, exprs_values="exprs", colour_by="total_features",
    feature_set=chosen) + fontsize
```

![plot of chunk PCA_corrHVG](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M1A/PCA_corrHVG-1.png)

On a related note, we only show the first two components that contribute most to the variance in Figure 13. Additional components can be visualized by increasing the **ncomponents** argument in **plotPCA** to construct pairwise plots. The percentage of variance explained by each component can also be obtained by running **plotPCA** with **return_SCESet=TRUE**, and then calling **reducedDimension** on the returned object. This information may be useful for selecting high-variance components (possibly corresponding to interesting underlying factors) for further examination.

Another widely used approach is the t-stochastic neighbour embedding (t-SNE) method (Van der Maaten and Hinton 2008). t-SNE tends to work better than PCA for separating cells in more diverse populations. This is because the former can directly capture non-linear relationships in high-dimensional space, whereas the latter must represent them (suboptimally) as linear components. However, this improvement comes at the cost of more computational effort and complexity. In particular, t-SNE is a stochastic method, so users should run the algorithm several times to ensure that the results are representative, and then set a seed to ensure that the chosen results are reproducible. It is also advisable to test different settings of the “perplexity” parameter as this will affect the distribution of points in the low-dimensional space. This is demonstrated below in Figure 14, though no consistent substructure is observed in all plots.

#### Figure 14:
*  _t_-SNE plots constructed from normalized log-expression values of correlated HVGs, using a range of perplexity values. In each plot, each point represents a cell in the HSC dataset. Bars represent the coordinates of the cells on each axis. Each cell is coloured according to its total number of expressed features.", fig.width=12, fig.height=6----


```r
 set.seed(100)
 out5 <- plotTSNE(GSE61533_G1.sce, exprs_values="exprs", perplexity=5, colour_by="total_features", 
     feature_set=chosen) + fontsize + ggtitle("Perplexity = 5")
 out10 <- plotTSNE(GSE61533_G1.sce, exprs_values="exprs", perplexity=10, colour_by="total_features",
     feature_set=chosen) + fontsize + ggtitle("Perplexity = 10")
 out20 <- plotTSNE(GSE61533_G1.sce, exprs_values="exprs", perplexity=20, colour_by="total_features",
     feature_set=chosen) + fontsize + ggtitle("Perplexity = 20")
 multiplot(out5, out10, out20, cols=3)
```

![plot of chunk tsneplot](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M1A/tsneplot-1.png)

There are many other dimensionality reduction techniques that we do not consider here but could also be used, e.g., multidimensional scaling, diffusion maps. These have their own advantages and disadvantages – for example, diffusion maps (see **plotDiffusionMap**) place cells along a continuous trajectory and are suited for visualizing graduated processes like differentiation (Angerer et al. 2016). For each visualization method, additional cell-specific information can be incorporated into the colour, size or shape of each point. Here, cells are coloured by the total number of expressed features to demonstrate that this metric does not drive any systematic differences across the population. The **selectorPlot** function from **scran** can also be used to interactively select groups of cells in two-dimensional space. This facilitates data exploration as visually identified subpopulations can be directly selected for further examination.

Finally, putative subpopulations can be computationally defined by cutting the dendrogram in **heat.out$colDendrogram** with **cutree** to form clusters. We do not attempt this here as the substructure is too weak for reliable clustering. In fact, users should generally treat clustering results with some caution. If the differences between cells are subtle, the assignment of cells into clusters may not be robust. Moreover, different algorithms can yield substantially different clusters by focusing on different aspects of the data. Experimental validation of the clusters is critical to ensure that the putative subpopulations actually exist.

### Additional comments


Once the basic analysis is completed, it is often useful to save the **SCESet** object to file with the **saveRDS** function. The object can then be easily restored into new R sessions using the **readRDS** function. This allows further work to be conducted without having to repeat all of the processing steps described above.

##### DONT REALLY NEED read/saveRDS!!!


```r
save(file=file.path(WRKDIR, 'Data', "GSE61533_G1.sce"), GSE61533_G1.sce)
```
A variety of methods are available to perform more complex analyses on the processed expression data. For example, cells can be ordered in pseudotime (e.g., for progress along a differentiation pathway) with **monocle** (Trapnell et al. 2014) or **TSCAN** (Z. Ji and Ji 2016); cell-state hierarchies can be characterized with the **sincell** package (Julia, Telenti, and Rausell 2015); and oscillatory behaviour can be identified using **Oscope** (Leng et al. 2015). HVGs can be used in gene set enrichment analyses to identify biological pathways and processes with heterogeneous activity, using packages designed for bulk data like **topGO** or with dedicated single-cell methods like **scde** (J. Fan et al. 2016). Full descriptions of these analyses are outside the scope of this workflow, so interested users are advised to consult the relevant documentation.


### Parameter settings:
  * WRKDIR = /mnt100/home/Dropbox/SingleCell/Jun2017/R
  * FN = M1A-scRNAseq_Workflow
  * Scripts = Scripts
  * EXT_DATA = /mnt100/home/Dropbox/SingleCell//Data
  * TABLE_DIR = /mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/tables/M1A 
  * RUN DATE = Wed Jun 21 02:51:04 2017



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
##  [1] gplots_3.0.1         RBGL_1.52.0          graph_1.54.0         org.Mm.eg.db_3.4.1  
##  [5] AnnotationDbi_1.38.1 scran_1.4.4          BiocParallel_1.10.1  destiny_2.4.0       
##  [9] mvoutlier_2.0.8      sgeostat_1.0-27      Rtsne_0.13           readxl_1.0.0        
## [13] R.utils_2.5.0        R.oo_1.21.0          R.methodsS3_1.7.1    BiocStyle_2.4.0     
## [17] scater_1.4.0         ggplot2_2.2.1        readr_1.1.1          irlba_2.2.1         
## [21] Matrix_1.2-10        plyr_1.8.4           tidyr_0.6.3          reshape2_1.4.2      
## [25] stringr_1.2.0        vcd_1.4-3            GenomicRanges_1.28.3 GenomeInfoDb_1.12.2 
## [29] IRanges_2.10.2       S4Vectors_0.14.3     Biobase_2.36.2       BiocGenerics_0.22.0 
## [33] edgeR_3.18.1         limma_3.32.2         RColorBrewer_1.1-2   data.table_1.10.4   
## [37] magrittr_1.5         dplyr_0.7.0          knitr_1.16           rmarkdown_1.6       
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
##  [58] flexmix_2.3-14          nnet_7.3-12             dynamicTreeCut_1.63-1  
##  [61] locfit_1.5-9.1          labeling_0.3            rlang_0.1.1            
##  [64] munsell_0.4.3           cellranger_1.1.0        tools_3.4.0            
##  [67] RSQLite_1.1-2           pls_2.6-0               evaluate_0.10          
##  [70] cvTools_0.3.2           yaml_2.1.14             robustbase_0.92-7      
##  [73] caTools_1.17.1          nlme_3.1-131            mime_0.5               
##  [76] quantreg_5.33           biomaRt_2.32.1          compiler_3.4.0         
##  [79] pbkrtest_0.4-7          beeswarm_0.2.3          e1071_1.6-8            
##  [82] statmod_1.4.30          smoother_1.1            tibble_1.3.3           
##  [85] robCompositions_2.0.3   pcaPP_1.9-61            stringi_1.1.5          
##  [88] highr_0.6               lattice_0.20-35         trimcluster_0.1-2      
##  [91] nloptr_1.0.4            lmtest_0.9-35           cowplot_0.7.0          
##  [94] bitops_1.0-6            httpuv_1.3.3            latticeExtra_0.6-28    
##  [97] R6_2.2.2                KernSmooth_2.23-15      gridExtra_2.2.1        
## [100] vipor_0.4.5             gtools_3.5.0            boot_1.3-19            
## [103] MASS_7.3-47             assertthat_0.2.0        rhdf5_2.20.0           
## [106] rprojroot_1.2           rjson_0.2.15            GenomeInfoDbData_0.99.0
## [109] diptest_0.75-7          mgcv_1.8-17             hms_0.3                
## [112] rpart_4.1-11            class_7.3-14            minqa_1.2.4            
## [115] TTR_0.23-1              shiny_1.0.3             base64enc_0.1-3        
## [118] ggbeeswarm_0.5.3
```

