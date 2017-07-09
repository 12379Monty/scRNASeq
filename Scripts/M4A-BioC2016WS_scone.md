### BioC 2016 workshop - Analysis of single-cell RNA-seq data with R and Bioconductor

### I - Quality control (QC) and normalization with scone

<style type="text/css"> body, td { font-size: 14px; } code.r{ font-size: 12px; } pre { font-size: 12px } </style>

[BioC 2016 workshop](https://github.com/drisso/bioc2016singlecell)



<!-- ***************************************************** -->


<!-- ***************************************************** -->



<!-- ***************************************************** -->



<!-- ***************************************************** -->

 

<!--
%\VignetteEngine{knitr::rmarkdown}
%\VignetteIndexEntry{scone Vignette}
-->



### Introduction

This is the first part of the Bioc2016 workshop "Analysis of single-cell RNA-seq data with R and Bioconductor."

In this part we will cover single-cell RNA-Seq quality control (QC) and normalization with the *[scone](https://github.com/YosefLab/scone)* package.  The package is available on [Github](https://github.com/YosefLab/scone) and will be submitted to Bioconductor in the near future.

Single-cell RNA sequencing (scRNA-Seq) technologies are opening the way for transcriptome-wide profiling across diverse and complex mammalian tissues, facilitating unbiased identification of novel cell sub-populations and their functional roles. As in other high-throughput assays, a fraction of the heterogeneity observed in scRNA-Seq data results from batch effects and other technical artifacts. In particular, these protocols’ reliance on minuscule amounts of starting mRNA can lead to widespread  “drop-out effects,” in which expressed transcripts are missed. Due to the biases inherent to these assays, data normalization is an essential step prior to any downstream analyses. Furthermore, due to wide-range of scRNA-Seq study designs used in the field, we cannot expect to find a one-size-fits-all solution to these problems.

`scone` supports a rational, data-driven framework for assessing the efficacy of various normalization workflows, encouraging users to explore trade-offs inherent to their data set prior to finalizing a data normalization strategy. We provide an interface for running multiple normalization workflows in parallel. We also offer tools for ranking workflows and visualizing trade-offs. We import some common normalization modules used in traditional bulk sequencing, and provide support for integrating user-specified normalization modules.

### An Example Dataset

The first two parts of this workshop will analyze a small set of cells lying along a developmental trajectory. We will start from raw data objects obtained from a standard transcriptome alignment pipeline. Raw data and important reference data can be loaded directly from the workshop package. 


```r
 suppressMessages(require(bioc2016singlecell))

 ## Load Example Data
 FILES <- ls()
 data(ws_input)
 setdiff(ls(), c('FILES', FILES))
```

```
## [1] "batch"  "bio"    "counts" "de"     "hk"     "qc"
```

The `counts` data frame contains feature-level read counts from tophat alignments of 96 single-cell libraries to the mm10 reference genome [trapnell2009]. `qc` contains library alignment metrics obtained from [Picard](http://broadinstitute.github.io/picard/).


```r
 colnames(qc)
```

```
##  [1] "NREADS"               "NALIGNED"             "RALIGN"               "TOTAL_DUP"           
##  [5] "PRIMER"               "PCT_RIBOSOMAL_BASES"  "PCT_CODING_BASES"     "PCT_UTR_BASES"       
##  [9] "PCT_INTRONIC_BASES"   "PCT_INTERGENIC_BASES" "PCT_MRNA_BASES"       "MEDIAN_CV_COVERAGE"  
## [13] "MEDIAN_5PRIME_BIAS"   "MEDIAN_3PRIME_BIAS"
```

`de` and `hk` are positive and negative control gene sets derived from population studies. `batch` and `bio` are factors labeling batch and time point respectively. Consider the joint distribution of these factors:


```r
 ## Joint distribution of batches and biological conditions (time after induction)
 table(batch,bio)
```

```
##       bio
## batch  time0 time1 time2 time3 time4 time5
##   Y01     10     0     0     0     0     0
##   Y02      0     0     0     6     0     0
##   Y03      0     0     0     8     0     0
##   Y04      8     0     0     0     0     0
##   Y05A     0    10     0     0     0     0
##   Y05B     0     9     0     0     0     0
##   Y06      0     0     5     0     0     0
##   Y07A     0     0     9     0     0     0
##   Y07B     0     0    10     0     0     0
##   Y08      0     6     0     0     0     0
##   Y10      0     0     0     0     2     0
##   Y11      0     0     0     0     4     0
##   Y12A     0     0     0     0     0     2
##   Y12B     0     0     0     0     0     2
##   Y13      0     0     0     0     0     5
```

Notice that each time-point is composed of multiple batches. This nested design is common among scRNA-Seq studies due to limitations on the number of cells that can be harvested/sequenced concurrently.

### Technical Variability and Batch Effects

We will first visualize cell level quality readouts, starting with the ratio of reads aligned to the genome:


```r
 # Color scheme
 cc <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"),brewer.pal(9,"Set3"))

 # Barplot of read proportion mapping to the genome
 barplot(qc$RALIGN[order(batch)], 
         col=cc[batch][order(batch)], 
         border=cc[batch][order(batch)], 
         main="Percentage of mapped reads, colored by batch")

 legend("bottomleft", legend=levels(batch), fill=cc,cex=0.4)
```

![plot of chunk ralign](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M4A/ralign-1.png)

We see that there aren't important differences between batches, while there are a few cells of particularly low alignment efficiency. These cells can be removed via filtering procedures. We can alternatively consider the number of reads in each library:


```r
 # Barplot of total read number
 barplot(qc$NREADS[order(batch)], 
         col=cc[batch][order(batch)], 
         border=cc[batch][order(batch)], 
         main="Total number of reads, colored by batch")
 legend("topright", legend=levels(batch), fill=cc, cex=0.4)
```

![plot of chunk nreads](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M4A/nreads-1.png)

We see that coverage varies substantially between batches, and though some of this technical variability may be addressed with filtering, we will still need to carefully normalize the data to make batches comparable.

It can be very helpful to visualize distributions of single metrics, but it's important to note that these metrics are often correlated. Sometimes it may be more useful to consider principal components of the quality matrix, identifying principal factors of library variation:


```r
 qpc = prcomp(qc,center = T,scale. = TRUE)
 barplot(cumsum(qpc$sdev^2)/sum(qpc$sdev^2), 
         border="gray", 
         xlab="PC", 
         ylab="proportion of variance", 
         main="Quality PCA")
```

![plot of chunk qpc](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M4A/qpc-1.png)

Even though 14 quality metrics have been quantified, PCA shows us that only a small number of PCs are needed to described a majority of the variance (e.g. 3 to explain 74%). We will now visualize the distribution of the first PC in the context of batch:


```r
 barplot(qpc$x[,1][order(batch)], 
         col=cc[batch][order(batch)], 
         border=cc[batch][order(batch)], 
         main="Quality PC1, colored by batch")
 legend("bottomleft", legend=levels(batch), fill=cc, cex=0.8)
```

![plot of chunk qpc_view](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M4A/qpc_view-1.png)

This PC appears to correlate with batch (particularly Y04), but it also highlights 5 outlier libraries that are significantly different from the rest. While we won't discuss the application of these PCs to filtering in this workshop, we will show how they can be applied to normalization.

### Drop-out Characteristics

Next we consider a uniquely single-cell problem: drop-outs. `hk` contains a list of genes that are believed to be ubiquitously and uniformly expressed across the target tissue. Because we assume these genes are expressed in all cells, we can label all zeroes as drop-out events. Below we model detection failures as a logistic function of mean expression, in line with the standard logistic model for drop-outs employed by the field.


```r
 # Mean log10(x+1) expression
 mu_obs = rowMeans(log10(counts[hk,]+1))

 # Drop-outs
 drop_outs = counts[hk,] == 0

 # Logistic Regression Model of Failure
 ref.glms = list()
 for (si in 1:dim(drop_outs)[2]){
   fit = glm(cbind(drop_outs[,si],1 - drop_outs[,si]) ~ mu_obs,family=binomial(logit))
   ref.glms[[si]] = fit$coefficients
 }
```

The list `ref.glm` contains the intercept and slope of each fit. We can now visualize the fit curves and the corresponding Area Under the Curves (AUCs):


```r
 par(mfrow=c(1,1))

 # Plot Failure Curves and Calculate AUC
 plot(NULL, main = "False Negative Rate Curves",
      ylim = c(0,1),xlim = c(0,6), ylab = "Failure Probability", 
      xlab = "Mean log10 Expression")
 x = (0:60)/10
 AUC = NULL
 for(si in 1:ncol(counts)){
   y = 1/(exp(-ref.glms[[si]][1] - ref.glms[[si]][2] * x) + 1)
   AUC[si] = sum(y)/10
   lines(x, 1/(exp(-ref.glms[[si]][1] - ref.glms[[si]][2] * x) + 1), 
         type = 'l', lwd = 2, col = cc[batch][si])
 }
```

<img src="/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M4A/fnr_vis-1.png" title="plot of chunk fnr_vis" alt="plot of chunk fnr_vis" width="800px" height="400px" />

```r
 # Barplot of FNR AUC
 barplot(AUC[order(batch)], 
         col=cc[batch][order(batch)], 
         border=cc[batch][order(batch)], 
         main="FNR AUC, colored by batch")
 legend("topleft", legend=levels(batch), fill=cc, cex=0.4)
```

<img src="/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M4A/fnr_vis-2.png" title="plot of chunk fnr_vis" alt="plot of chunk fnr_vis" width="800px" height="400px" />

You may have noticed that the 5 outlier cells seen here are the same as the 5 outliers highlighted in the PCA of quality metrics above. Drop-out summaries such as the AUC can very useful for assessing single-cell library quality.

### The `scone` Workflow

So far we have only described potential problems with our data set. Now we will take steps to address them! The basic qc and normalization pipeline we will use today will allow us to: 

* Filter out poor libraries using the `metric_sample_filter` function.
*  Run and score many different normalization workflows (different combinations of normalization modules) using the main `scone` function.
* Browse top-ranked methods and visualize trade-offs with the `biplot_colored` function.

The SCONE's normalization workflow template is composed of 3 modules:

1) Data imputation: replacing zero-abundance values with expected values under a drop-out model. We will only consider identity imputation (no imputation) in this example.
2) Scaling or quantile normalization: either i) normalization that scales each sample's transcriptome abundances by a single factor or ii) more complex offsets that match quantiles across samples. Examples: TMM or DESeq scaling factors, upper quartile normalization, or full-quantile normalization.
3) Regression-based approaches for removing unwanted correlated variation from the data. Examples: RUVg [@risso2014] or regression on Quality Principal Components described above.

### Step 1: Sample Filtering with `metric_sample_filter`

The most basic sample filtering function in `scone` is the `metric_sample_filter`. It takes as input an expression matrix. The output depends on arguments provided, but generally consists of a list of 4 logicals designating each sample as having failed (T) or passed (FALSE) a threshold-based filter on 4 metrics

* Number of reads.
* Ratio of reads aligned to the genome. Requires the `ralign` argument.
* "Transcriptome breadth" - Defined here as the proportion of high-quality genes detected in the sample. Requires the `gene_filter` argument.
* FNR AUC. Requires the `pos_controls` argument.

If required arguments are missing for any of the 4, the function will simply return NA instead of the corresponding logical.


```r
 mfilt_report <- metric_sample_filter(expr = as.matrix(counts),
                                      nreads = qc$NREADS,ralign = qc$RALIGN,
                                      suff_nreads = 10^5,
                                      suff_ralign = 90,
                                      pos_controls = hk,
                                      zcut = 3,mixture = F, plot = T)
```

![plot of chunk metric_sample_filter](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M4A/metric_sample_filter-1.png)

In the call above, we have set the following parameters:

* suff_nreads = 10^5. Sets max (or "sufficient") threshold for nreads filter.
* suff_ralign = 90. Sets max threshold for ralign filter.
* zcut = 3. Filter leniency (see below).
* mixture = FALSE. Mixture modeling will not be used (see below).
* plot = T. Plot distributions of metrics before and after filtering.

### On Threshold Selection

Let's take a closer look at the computation behind the ralign filter. In selecting the threshold value 90, `metric_sample_filter` is taking 4 values into account:

1) `hard_ralign`, the default "hard" threshold at 15 - rather forgiving...
2) 3 (`zcut`) times the standard deviation below the mean `ralign` value.
3) 3 (`zcut`) times the MAD below the median `ralign` value.
4) `suff_ralign`, the sufficient threshold set to 90.


```r
 hist(qc$RALIGN, breaks = 0:100)

 # Hard threshold
 abline(v = 15, col = "yellow", lwd = 2) 

 # 3 (zcut) standard deviations below the mean ralign value
 abline(v = mean(qc$RALIGN) - 3*sd(qc$RALIGN), col = "green", lwd = 2) 

 # 3 (zcut) MADs below the median ralign value
 abline(v = median(qc$RALIGN) - 3*mad(qc$RALIGN), col = "red", lwd = 2)

 # Sufficient threshold
 abline(v = 90, col = "grey", lwd = 2)

 # Final threshold is the minimum of i
 # 1) the sufficient threshold and 2) the max of all others
 thresh = min(90,
              max(c(15,
                    mean(qc$RALIGN) - 3*sd(qc$RALIGN),
                    median(qc$RALIGN) - 3*mad(qc$RALIGN))
                  )
              )
 abline(v = thresh, col = "blue", lwd = 2, lty = 2)

 legend("topleft",
        legend = c("Hard","SD","MAD","Sufficient","Final"),
        lwd = 2, col = c("yellow","green","red","grey","blue"),
        lty = c(1,1,1,1,2), cex = .5)
```

<img src="/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M4A/thresh-1.png" title="plot of chunk thresh" alt="plot of chunk thresh" width="600px" height="400px" />

We see here that the 3rd threshold exceeds the sufficient threshold, and so `metric_sample_filter` settles for the sufficient threshold. Note that if `mixture=T` an additional criterion is included: bi-modal metric distributions are fit with to a two component mixture model, and a threshold is defined with respect to the mean and standard deviation of the "better" component. As `metric_sample_filter` relies on a maximum of candidate thresholds, we recommend users treat this function as a liberal filter.

### Applying the sample filter

With `metric_sample_filter` output in hand, filtering out the poor samples is fairly straightforward:


```r
 # Which thresholds are missing? (breadth)
 is_na_filt = unlist(lapply(is.na(mfilt_report),any)) 

 # Identify samples failing any threshold
 m_sampfilter = !apply(simplify2array(mfilt_report[!is_na_filt]),1,any)

 # Filter Samples
 fcounts = counts[,m_sampfilter]
 fqc = qc[m_sampfilter,]
 fbatch = batch[m_sampfilter]
 fbio = bio[m_sampfilter]

 # Simple gene filter
 filterCount <- function(counts, nRead=5, nCell=5){
  filter <- apply(counts, 1, function(x) length(x[x>=nRead])>=nCell)
  return(filter)
 }
 genefilter <- filterCount(fcounts)

 # Filter genes
 fcounts = fcounts[genefilter,]
 fhk = hk[hk %in% rownames(fcounts)]
 fde = de[de %in% rownames(fcounts)]
```

### Step 2: Run and Score Normalization Workflows with `scone`

As described earlier, not only does `scone` normalize expression data, but it also provides a framework for evaluation the performance of various normalization workflows. In order to run the `scone` function, we will need to decide which workflows we will want to compare. 

#### Selecting Workflow Parameters with `run=F`

The arguments to `scone` allow for a lot of flexibility, but sometimes you will need to run very specific combinations of modules. For this purpose, `scone` can be run in `run=F` mode, returning only a data frame of workflows to be performed.

#### scone now requires a SconeExperiment Object!!!! 
#### STOP EXECUTION HERE!!!


```r
 params <- scone(expr = as.matrix(fcounts),
                 scaling = c(none = identity,deseq = DESEQ_FN, 
                             tmm = TMM_FN, fq = FQT_FN),
             ######## uqp = UQ_FN_POS, --> DONT KNOW WHERE THAT IS!!!
                ruv_negcon = fhk, k_ruv = 3,
                qc = as.matrix(fqc), k_qc = 3,
                bio = fbio,adjust_bio = "yes",
                batch = fbatch,adjust_batch = "yes",
                run = F)
 head(params)
```

In the call above, we have set the following parameters:

* scaling = list(...). This argument contains a list of scaling normalization functions that will be applied, including the identity (no-op), DESeq scaling, TMM normalization, scaling by the upper quartile of positive counts, and full-quantile normalization. 
* ruv_negcon = fhk. A list of genes to be used as negative controls for RUVg normalization.
* k_ruv = 3. The maximum number of RUVg factors to consider.
* k_qc = 3. The maximum number of quality PCs (QPCs) to be included in a linear model, analogous to RUVg normalization. The qc argument must be provided.
* adjust_bio = "yes." Time points will be included in RUVg or QPC normalization and restored to the residual of the regression. The bio argument must be provided.
* adjust_batch = "yes." Batch will be included in RUVg or QPC normalization. The batch argument must be provided.

These arguments translate to the following set of options:


```r
 apply(params,2,unique)
```

While adjusting for biology prior to downstream analysis is generally problematic, it can be particularly dangerous in the single-cell context: there is significantly more sample-intrinsic biological variability in scRNA-Seq data, variability that cannot be captured by exposure. In our nested design scenario we would only consider adjustments for biology when batch effects are included. We can produce an updated `params` data frame reflecting this choice:


```r
 nis_screened = (params$adjust_biology == "bio") & 
                (params$adjust_batch != "batch")

 params = params[!is_screened,]
```

### Calling `scone` with `run=T`

Now that we have selected our workflows, we can run `scone` in `run=T` mode. This mode offers a few additional arguments, including an optional `params` argument to pass any results from the `run=F` mode. In order to understand these arguments, we must first understand the 8 metrics used to evaluate each normalization. The first 6 metrics rely on a reduction of the normalized data down to 3 dimensions via PCA (default). Each metric is taken to have a positive (higher is better) or negative (lower is better) signature.

* BIO_SIL. The average silhouette width of clusters defined by `bio`, defined with respect to a Euclidean distance metric over the first 3 expression PCs. Positive signature.
* BATCH_SIL. The average silhouette width of clusters defined by `batch`, defined with respect to a Euclidean distance metric over the first 3 expression PCs. Negative signature.
* PAM_SIL. The maximum average silhouette width of clusters defined by PAM clustering, defined with respect to a Euclidean distance metric over the first 3 expression PCs. Positive signature.
* EXP_QC_COR. Maximum squared Spearman correlation between first 3 expression PCs and first `k_qc` QPCs. Negative signature.
* EXP_UV_COR. Maximum squared Spearman correlation between first 3 expression PCs and first 3 PCs of the negative control (specified by `eval_negcon` or `ruv_negcon` by default) sub-matrix of the original (raw) data. Negative signature.
* EXP_WV_COR. Maximum squared Spearman correlation between first 3 expression PCs and first 3 PCs of the positive control (specified by `eval_poscon`) sub-matrix of the original (raw) data. Positive signature.
* RLE_MED. The mean squared median Relative Log Expression (RLE). Negative signature.
* RLE_IQR. The mean inter-quartile range (IQR) of the RLE. Negative signature.


```r
 res <- scone(expr = as.matrix(fcounts),
              scaling = c(none = identity, deseq = DESEQ_FN, tmm = TMM_FN, 
                             fq = FQT_FN),
              #######uqp = UQ_FN_POS, --> DONT KNOW WHERE THAT IS
              ruv_negcon = fhk, k_ruv = 3,
              qc = as.matrix(fqc), k_qc = 3,
              bio = fbio,adjust_bio = "yes",
              batch = fbatch,adjust_batch = "yes",
              run = T,params = params,
              eval_poscon = fde, eval_kclust = 2:3, conditional_pam = T)
```

In the call above, we have set the following parameters:

* eval_poscon = fde. A list of genes to be used as positive controls for evaluation.
* eval_kclust = 2:3. For PAM_SIL, range of k (# of clusters) to use when computing maximum average silhouette width of PAM clusterings.
* conditional_pam = T. For PAM_SIL, apply separate PAM clusterings to each time point rather than across all time points. Average is weighted by time point group size.

Running the command will take a couple minutes. The output will contain a list of four elements:


```r
 names(res)
```

`normalized_data` contains a list of normalized expression data (log-scale); each expression matrix is named according to the same convention as seen in the `params` argument.

`metrics` contains the 8 raw metrics for each normalization. `scores` contains metrics multiplied by their signature - or "scores" - as well as a 9th column that contains the mean score for that normalization. Normalization workflows in `normalized_data`,`metrics`, and `scores` are sorted in decreasing order by mean score. 


```r
 head(res$scores)
```

### Step 3: Selecting a normalization for downstream analysis

Based on our sorting criteria, it would appear that `none,fq,ruv_k=1,no_bio,no_batch` performs well compared to other normalization workflows. A useful way to visualize this method with respect to others is the `biplot_colored` function


```r
 pc_obj = prcomp(res$scores[,-ncol(res$scores)],center = T,scale = F)
 bp_obj = biplot_colored(pc_obj,y = -res$scores[,ncol(res$scores)],expand = .6)
```

We have colored each point above according the corresponding method's mean score. We can see that workflows cluster largely by three scores: BIO_SIL, EXP_WV_COR, and EXP_UV_COR. Since clustering by biology is playing a strong role, we can highlight the points corresponding the biology adjustment:


```r
 bp_obj = biplot_colored(pc_obj,y = -res$scores[,ncol(res$scores)],expand = .6)

 points(bp_obj[grepl(",bio,batch",rownames(bp_obj)),], pch = 1, col = "red", cex = 1)
 points(bp_obj[grepl(",bio,batch",rownames(bp_obj)),], pch = 1, col = "red", cex = 1.5)
```

We see that the intermediate-performance cluster to the bottom-left is made up only of these methods. This example highlights the danger of adjusting for biology, even in the nested design scenario. Although biology-adjusted data may reflect prior differences more strongly, other metrics would suggest that the data is poorly behaved. We can also consider cases where only batch is accounted for:


```r
 bp_obj = biplot_colored(pc_obj,y = -res$scores[,ncol(res$scores)],expand = .6)

 points(bp_obj[grepl(",no_bio,batch",rownames(bp_obj)),], pch = 1, col = "red", cex = 1)
 points(bp_obj[grepl(",no_bio,batch",rownames(bp_obj)),], pch = 1, col = "red", cex = 1.5)
```

Here we see that the poorly performing cluster on the right is made up solely of these methods. This follows from the simple fact that batch is nested within biology, and adjustments for the former can remove critical biological variation. Finally we will visualize the top-performing method and it's relation to un-normalized data: 


```r
 bp_obj = biplot_colored(pc_obj,y = -res$scores[,ncol(res$scores)],expand = .6)

 points(t(bp_obj[1,]), pch = 1, col = "red", cex = 1)
 points(t(bp_obj[1,]), pch = 1, col = "red", cex = 1.5)

 points(t(bp_obj[rownames(bp_obj) == rownames(params)[1],]), 
        pch = 1, col = "blue", cex = 1)
 points(t(bp_obj[rownames(bp_obj) == rownames(params)[1],]), 
        pch = 1, col = "blue", cex = 1.5)

 arrows(bp_obj[rownames(bp_obj) == rownames(params)[1],][1],
        bp_obj[rownames(bp_obj) == rownames(params)[1],][2],
        bp_obj[1,][1],
        bp_obj[1,][2],
        lty = 2, lwd = 2)
```

The arrow traces a line from the "no-op" normalization to the top-ranked normalization in SCONE. We see that SCONE has selected a method in-between the two extremes, reducing the signal of unwanted variation (as defined by `fhk`) while preserving biological signal (as defined by `fbio` and `fde`).

### Parameter settings:
  * WRKDIR = /mnt100/home/Dropbox/SingleCell/Jun2017/R
  * FN = M4A-BioC2016WS_scone
  * Scripts = Scripts
  * RUN DATE = Sun Jun 25 03:52:57 2017



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
## [1] parallel  stats4    methods   stats     graphics  grDevices utils     datasets  base     
## 
## other attached packages:
##  [1] bioc2016singlecell_0.0.2   devtools_1.13.2            BiocStyle_2.4.0           
##  [4] RColorBrewer_1.1-2         scone_1.0.0                SummarizedExperiment_1.6.3
##  [7] DelayedArray_0.2.7         matrixStats_0.52.2         Biobase_2.36.2            
## [10] GenomicRanges_1.28.3       GenomeInfoDb_1.12.2        IRanges_2.10.2            
## [13] S4Vectors_0.14.3           BiocGenerics_0.22.0        BiocParallel_1.10.1       
## [16] knitr_1.16                 rmarkdown_1.6             
## 
## loaded via a namespace (and not attached):
##   [1] backports_1.1.0          aroma.light_3.6.0        NMF_0.20.6              
##   [4] igraph_1.0.1             plyr_1.8.4               lazyeval_0.2.0          
##   [7] shinydashboard_0.6.1     splines_3.4.0            scater_1.4.0            
##  [10] ggplot2_2.2.1            gridBase_0.4-7           digest_0.6.12           
##  [13] foreach_1.4.3            htmltools_0.3.6          viridis_0.4.0           
##  [16] gdata_2.18.0             magrittr_1.5             memoise_1.1.0           
##  [19] cluster_2.0.6            doParallel_1.0.10        mixtools_1.1.0          
##  [22] limma_3.32.2             Biostrings_2.44.1        annotate_1.54.0         
##  [25] bayesm_3.0-2             R.utils_2.5.0            rARPACK_0.11-0          
##  [28] colorspace_1.3-2         dplyr_0.7.0              tximport_1.4.0          
##  [31] RCurl_1.95-4.8           jsonlite_1.5             hexbin_1.27.1           
##  [34] genefilter_1.58.1        survival_2.41-3          zoo_1.8-0               
##  [37] iterators_1.0.8          glue_1.1.0               registry_0.3            
##  [40] gtable_0.2.0             zlibbioc_1.22.0          XVector_0.16.0          
##  [43] compositions_1.40-1      kernlab_0.9-25           prabclus_2.2-6          
##  [46] DEoptimR_1.0-8           scales_0.4.1             DESeq_1.28.0            
##  [49] mvtnorm_1.0-6            DBI_0.7                  edgeR_3.18.1            
##  [52] rngtools_1.2.4           miniUI_0.1.1             Rcpp_0.12.11            
##  [55] viridisLite_0.2.0        xtable_1.8-2             mclust_5.3              
##  [58] DT_0.2                   htmlwidgets_0.8          httr_1.2.1              
##  [61] FNN_1.1                  gplots_3.0.1             fpc_2.1-10              
##  [64] modeltools_0.2-21        XML_3.98-1.8             R.methodsS3_1.7.1       
##  [67] flexmix_2.3-14           nnet_7.3-12              locfit_1.5-9.1          
##  [70] dynamicTreeCut_1.63-1    rlang_0.1.1              reshape2_1.4.2          
##  [73] AnnotationDbi_1.38.1     munsell_0.4.3            tools_3.4.0             
##  [76] visNetwork_1.0.3         RSQLite_1.1-2            evaluate_0.10           
##  [79] stringr_1.2.0            yaml_2.1.14              robustbase_0.92-7       
##  [82] caTools_1.17.1           purrr_0.2.2.2            EDASeq_2.10.0           
##  [85] mime_0.5                 scran_1.4.4              R.oo_1.21.0             
##  [88] biomaRt_2.32.1           compiler_3.4.0           beeswarm_0.2.3          
##  [91] plotly_4.7.0             tibble_1.3.3             statmod_1.4.30          
##  [94] geneplotter_1.54.0       stringi_1.1.5            highr_0.6               
##  [97] GenomicFeatures_1.28.3   RSpectra_0.12-0          lattice_0.20-35         
## [100] trimcluster_0.1-2        Matrix_1.2-10            tensorA_0.36            
## [103] data.table_1.10.4        bitops_1.0-6             httpuv_1.3.3            
## [106] rtracklayer_1.36.3       R6_2.2.2                 latticeExtra_0.6-28     
## [109] hwriter_1.3.2            ShortRead_1.34.0         gridExtra_2.2.1         
## [112] KernSmooth_2.23-15       vipor_0.4.5              codetools_0.2-15        
## [115] boot_1.3-19              energy_1.7-0             MASS_7.3-47             
## [118] gtools_3.5.0             assertthat_0.2.0         rhdf5_2.20.0            
## [121] rjson_0.2.15             pkgmaker_0.22            rprojroot_1.2           
## [124] withr_1.0.2              RUVSeq_1.10.0            GenomicAlignments_1.12.1
## [127] Rsamtools_1.28.0         GenomeInfoDbData_0.99.0  diptest_0.75-7          
## [130] grid_3.4.0               tidyr_0.6.3              class_7.3-14            
## [133] segmented_0.5-2.1        shiny_1.0.3              ggbeeswarm_0.5.3
```

