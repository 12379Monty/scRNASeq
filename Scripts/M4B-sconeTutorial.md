### scone Tutorial


<style type="text/css"> body, td { font-size: 14px; } code.r{ font-size: 12px; } pre { font-size: 12px } </style>

[scone tutorial](https://github.com/YosefLab/scone/tree/master/vignettes)



<!-- ***************************************************** -->


<!-- ***************************************************** -->



<!-- ***************************************************** -->



<!-- ***************************************************** -->

 

<!--
%\VignetteEngine{knitr::rmarkdown}
%\VignetteIndexEntry{scone Vignette}
-->


 
### Introduction
 
 Single-cell RNA sequencing (scRNA-Seq) technologies are opening the way for
 transcriptome-wide profiling across diverse and complex mammalian tissues,
 facilitating unbiased identification of novel cell sub-populations and
 discovery of novel cellular function. As in other high-throughput analyses, a
 large fraction of the variability observed in scRNA-Seq data results from batch
 effects and other technical artifacts [@hicks2015]. In particular, a unique
 reliance on minuscule amounts of starting mRNA can lead to widespread “drop-out
 effects,” in which expressed transcripts are missed during library preparation
 and sequencing. Due to the biases inherent to these assays, data normalization
 is an essential step prior to many downstream analyses. As we face a growing
 cohort of scRNA-Seq technologies, diverse biological contexts, and novel
 experimental designs, we cannot reasonably expect to find a one-size-fits-all
 solution to data normalization.
 
 `scone` supports a rational, data-driven framework for assessing the efficacy
 of various normalization workflows, encouraging users to explore trade-offs
 inherent to their data prior to finalizing their data normalization strategy.
 We provide an interface for running multiple normalization workflows in
 parallel, and we offer tools for ranking workflows and visualizing
 study-specific trade-offs.
 
 This package was originally developed to address normalization problems
 specific to scRNA-Seq expression data, but it should be emphasized that its use
 is not limited to scRNA-Seq data normalization. Analyses based on other
 high-dimensional data sets - including bulk RNA-Seq data sets - can utilize
 tools implemented in the `scone` package.
 
### Human Neurogenesis
 
 We will demonstrate the basic `scone` workflow by using an early scRNA-Seq data
 set [@pollen2014]. We focus on a set of 65 human cells sampled from four
 biological conditions: Cultured neural progenitor cells ("NPC") derived from
 pluripotent stem cells, primary cortical samples at gestation weeks 16 and 21
 ("GW16" and "GW21" respectively) and late cortical samples cultured for 3 weeks
 ("GW21+3"). Gene-level expression data for these cells can be loaded directly
 from the `scRNAseq` package on 
 [Bioconductor](http://bioconductor.org/packages/scRNAseq/).
 

```r
 suppressMessages(require(scRNAseq))
 
 ## ----- Load Example Data -----
 data(fluidigm)
 
 # Set assay to RSEM estimated counts
 assay(fluidigm) = assays(fluidigm)$rsem_counts
```
 
 The `rsem_counts` assay contains expected gene-level read counts via RSEM
 [@li2011] quantification of 130 single-cell libraries aligned to the hg38
 RefSeq transcriptome. The data object also contains library transcriptome
 alignment metrics obtained from
 [Picard](http://broadinstitute.github.io/picard/) and other basic tools.
 

```r
 ## ----- List all QC fields -----
 
 # List all qc fields (accessible via colData())
 metadata(fluidigm)$which_qc
```

```
##  [1] "NREADS"                       "NALIGNED"                     "RALIGN"                      
##  [4] "TOTAL_DUP"                    "PRIMER"                       "INSERT_SZ"                   
##  [7] "INSERT_SZ_STD"                "COMPLEXITY"                   "NDUPR"                       
## [10] "PCT_RIBOSOMAL_BASES"          "PCT_CODING_BASES"             "PCT_UTR_BASES"               
## [13] "PCT_INTRONIC_BASES"           "PCT_INTERGENIC_BASES"         "PCT_MRNA_BASES"              
## [16] "MEDIAN_CV_COVERAGE"           "MEDIAN_5PRIME_BIAS"           "MEDIAN_3PRIME_BIAS"          
## [19] "MEDIAN_5PRIME_TO_3PRIME_BIAS"
```
 
 All cell-level metadata, such as cell origin and sequence coverage ("low" vs
 "high" coverage) can be accessed using `colData()`:
 

```r
 # Joint distribution of "biological condition"" and "coverage type""
 table(colData(fluidigm)$Coverage_Type,
       colData(fluidigm)$Biological_Condition)
```

```
##       
##        GW16 GW21 GW21+3 NPC
##   High   26    8     16  15
##   Low    26    8     16  15
```
 
 Each cell had been sequenced twice, at different levels of coverage. In this
 vignette we will focus on the high-coverage data. Before we get started we will
 do some preliminary filtering to remove the low-coverage replicates and
 undetected gene features:
 

```r
 # Preliminary Sample Filtering: High-Coverage Only
 is_select = colData(fluidigm)$Coverage_Type == "High"
 fluidigm = fluidigm[,is_select]
 
 # Retain only detected transcripts
 fluidigm = fluidigm[which(apply(assay(fluidigm) > 0,1,any)),]
```
 
### Visualizing Technical Variability and Batch Effects
 
 One of our alignment quality readouts is the fraction of reads aligned to the
 transcriptome. We can use simple bar plots to visualize how this metric relates
 to the biological batch.
 

```r
 # Define a color scheme
 cc <- c(brewer.pal(9, "Set1"))
 
 # One batch per Biological Condition
 batch = factor(colData(fluidigm)$Biological_Condition)
 
 # Alignment Quality Metrics
 qc = colData(fluidigm)[,metadata(fluidigm)$which_qc]
 
 # Barplot of read proportion mapping to human transcriptome
 ralign = qc$RALIGN
 o = order(ralign)[order(batch[order(ralign)])] # Order by batch, then value
 
 barplot(ralign[o], col=cc[batch][o], 
         border=cc[batch][o], main="Percentage of reads mapped")
 legend("bottomleft", legend=levels(batch), fill=cc,cex=0.4)
```

![plot of chunk ralign](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M4B/ralign-1.png)
 
 We can see modest differences between batches, and we see that there is one
 GW21 cell with a particularly low alignment rate relative to the rest of the
 GW21 batch. These types of observations can inform us of "poor-quality"
 libraries or batches. We may alternatively consider the number of reads for
 each library:
 

```r
 # Barplot of total read number
 nreads = qc$NREADS
 o = order(nreads)[order(batch[order(nreads)])] # Order by batch, then value
 
 barplot(nreads[o], col=cc[batch][o], 
         border=cc[batch][o], main="Total number of reads")
 legend("topright", legend=levels(batch), fill=cc, cex=0.4)
```

![plot of chunk nreads](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M4B/nreads-1.png)
 
 We see that read coverage varies substantially between batches as well as
 within batches. These coverage differences and other technical features can
 induce non-intuitive biases upon expression estimates. Though some biases can
 be addressed with simple library-size normalization and cell-filtering, demand
 for greater cell numbers may require more sophisticated normalization methods
 in order to compare multiple batches of cells. Batch-specific biases are
 impossible to address directly in this study as biological origin and sample
 preparation are completely confounded.
 
 While it can be very helpful to visualize distributions of single quality
 metrics it should be noted that QC metrics are often correlated. In some cases
 it may be more useful to consider Principal Components (PCs) of the quality
 matrix, identifying latent factors of protocol variation:
 

```r
 ## ----- PCA of QC matrix -----
 qpc = prcomp(qc,center = TRUE,scale. = TRUE)
 barplot((qpc$sdev^2)/sum(qpc$sdev^2), border="gray", 
         xlab="PC", ylab="Proportion of Variance", main="Quality PCA")
```

![plot of chunk qpc](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M4B/qpc-1.png)
 
 Even though 19 different QC metrics have been quantified in this analysis, PCA
 shows us that only a small number of PCs are needed to described a majority of
 the QC variance (e.g. 3 to explain 76%). We will now visualize the distribution
 of the first PC in the context of batch:
 

```r
 # Barplot of PC1 of the QC matrix
 qc1 = as.vector(qpc$x[,1])
 o = order(qc1)[order(batch[order(qc1)])]
 
 barplot(qc1[o], col=cc[batch][o], 
         border=cc[batch][o], main="Quality PC1")
 legend("bottomright", legend=levels(batch), 
        fill=cc, cex=0.8)
```

![plot of chunk qpc_view](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M4B/qpc_view-1.png)
 
 This first PC appears to represent both inter-batch and intra-batch sample
 heterogeneity, similar the the total number of reads. If this latent factor
 reflects variation in sample preparation, we may expect expression artifacts to
 trace this factor as well: in other words, we should be very skeptical of genes
 for which expression correlates strongly with the first PC of quality metrics.
 In this vignette we will show how latent factors like this can be applied to
 the normalization problem.
 
### Drop-out Characteristics
 
 Before we move on to normalization, let's briefly consider a uniquely
 single-cell problem: "drop-outs." One of the greatest challenges in modeling
 drop-out effects is modeling both i) technical drop-outs and ii) biological
 expression heterogeneity. One way to simplify the problem is to focus on genes
 for which we have strong prior belief in true expression. The `scone` package
 contains lists of genes that are believed to be ubiquitously and even uniformly
 expressed across human tissues. If we assume these genes are truly expressed in
 all cells, we can label all zero abundance observations as drop-out events. We
 model detection failures as a logistic function of mean expression, in line
 with the standard logistic model for drop-outs employed by the field:
 

```r
 # Extract Housekeeping Genes
 data(housekeeping)
 hk = intersect(housekeeping$V1,rownames(assay(fluidigm)))
 
 # Mean log10(x+1) expression
 mu_obs = rowMeans(log10(assay(fluidigm)[hk,]+1))
 
 # Assumed False Negatives
 drop_outs = assay(fluidigm)[hk,] == 0
 
 # Logistic Regression Model of Failure
 ref.glms = list()
 for (si in 1:dim(drop_outs)[2]){
   fit = glm(cbind(drop_outs[,si],1 - drop_outs[,si]) ~ mu_obs,
             family=binomial(logit))
   ref.glms[[si]] = fit$coefficients
 }
```
 
 The list `ref.glm` contains the intercept and slope of each fit. We can now
 visualize the fit curves and the corresponding Area Under the Curves (AUCs):
 

```r
 par(mfrow=c(1,2))
 
 # Plot Failure Curves and Calculate AUC
 plot(NULL, main = "False Negative Rate Curves",
      ylim = c(0,1),xlim = c(0,6), 
      ylab = "Failure Probability", xlab = "Mean log10 Expression")
 x = (0:60)/10
 AUC = NULL
 for(si in 1:ncol(assay(fluidigm))){
   y = 1/(exp(-ref.glms[[si]][1] - ref.glms[[si]][2] * x) + 1)
   AUC[si] = sum(y)/10
   lines(x, 1/(exp(-ref.glms[[si]][1] - ref.glms[[si]][2] * x) + 1),
         type = 'l', lwd = 2, col = cc[batch][si])
 }
 
 # Barplot of FNR AUC
 o = order(AUC)[order(batch[order(AUC)])]
 
 barplot(AUC[o], col=cc[batch][o], border=cc[batch][o], main="FNR AUC")
 legend("topright", legend=levels(batch), fill=cc, cex=0.4)
```

<img src="/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M4B/fnr_vis-1.png" title="plot of chunk fnr_vis" alt="plot of chunk fnr_vis" width="800px" height="400px" />
 
 Model-based metrics such as these may be more interpretable with respect to
 upstream sample preparation, and can be very useful for assessing single-cell
 library quality.
 
### The `scone` Workflow
 
 So far we have only described potential problems with single-cell expression
 data. Now we will take steps to address problems with our example data set. The
 basic QC and normalization pipeline we will use in this vignette allows us to:
 
 * Filter out poor libraries using the `metric_sample_filter` function.
 * Run and score many different normalization workflows 
   (different combinations of normalization modules)
   using the main `scone` function.
 * Browse top-ranked methods and visualize trade-offs with the
   `biplot_color` and `sconeReport` function.
 
 In order to run many different workflows, SCONE relies on a normalization
 workflow template composed of 3 modules:
 
 1) Data imputation: replacing zero-abundance values with expected values under
 a prior drop-out model. As we will see below, this module may be used as a
 modifier for module 2, without passing imputed values forward to downstream
 analyses. 2) Scaling or quantile normalization: either i) normalization that
 scales each sample's transcriptome abundances by a single factor or ii) more
 complex offsets that match quantiles across samples. Examples: TMM or DESeq
 scaling factors, upper quartile normalization, or full-quantile normalization.
 3) Regression-based approaches for removing unwanted correlated variation from
 the data, including batch effects. Examples: RUVg [@risso2014] or regression on
 Quality Principal Components described above.
 
### Sample Filtering with `metric_sample_filter`
 
 The most basic sample filtering function in `scone` is the
 `metric_sample_filter`. The function takes a consensus approach, retaining
 samples that pass multiple data-driven criteria.
 
 `metric_sample_filter` takes as input an expression matrix. The output depends
 on arguments provided, but generally consists of a list of 4 logicals
 designating each sample as having failed (TRUE) or passed (FALSE)
 threshold-based filters on 4 sample metrics:
 
 * Number of reads.
 * Ratio of reads aligned to the genome. 
   Requires the `ralign` argument.
 * "Transcriptome breadth" - Defined here as the proportion of "high-quality"
   genes detected in the sample. Requires the `gene_filter` argument.
 * FNR AUC. Requires the `pos_controls` argument.
 
 If required arguments are missing for any of the 4, the function will simply
 return NA instead of the corresponding logical.
 

```r
 # Initial Gene Filtering: 
 # Select "common" transcripts based on proportional criteria.
 num_reads = quantile(assay(fluidigm)[assay(fluidigm) > 0])[4]
 num_cells = 0.25*ncol(fluidigm)
 is_common = rowSums(assay(fluidigm) >= num_reads ) >= num_cells
 
 # Metric-based Filtering
 mfilt = metric_sample_filter(assay(fluidigm),
                              nreads = colData(fluidigm)$NREADS,
                              ralign = colData(fluidigm)$RALIGN,
                              gene_filter = is_common,
                              pos_controls = rownames(fluidigm) %in% hk,
 
                              zcut = 3, mixture = FALSE,
                              plot = TRUE)
```

<img src="/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M4B/metric_sample_filter-1.png" title="plot of chunk metric_sample_filter" alt="plot of chunk metric_sample_filter" height="1000px" />

```r
 # Simplify to a single logical
 mfilt = !apply(simplify2array(mfilt[!is.na(mfilt)]),1,any)
```
 
 In the call above, we have set the following parameters:
 
 * zcut = 3. Filter leniency (see below).
 * mixture = FALSE. Mixture modeling will not be used (see below).
 * plot = TRUE. Plot distributions of metrics before and after filtering.
 
### On Threshold Selection
 
 Let's take a closer look at the computation behind selecting the ralign filter.
 In choosing a threshold value 67.7, `metric_sample_filter` is taking 4 values
 into account:
 
 1) `hard_ralign`, the default "hard" threshold at 15 - rather forgiving... 2) 3
 (`zcut`) times the standard deviation below the mean `ralign` value. 3) 3
 (`zcut`) times the MAD below the median `ralign` value. 4) `suff_ralign`, the
 sufficient threshold set to NULL by default.
 

```r
 hist(qc$RALIGN, breaks = 0:100)
 # Hard threshold
 abline(v = 15, col = "yellow", lwd = 2)
 # 3 (zcut) standard deviations below the mean ralign value
 abline(v = mean(qc$RALIGN) - 3*sd(qc$RALIGN), col = "green", lwd = 2)
 # 3 (zcut) MADs below the median ralign value
 abline(v = median(qc$RALIGN) - 3*mad(qc$RALIGN), col = "red", lwd = 2)
 # Sufficient threshold
 abline(v = NULL, col = "grey", lwd = 2)
 
 # Final threshold is the minimum of 
 # 1) the sufficient threshold and 
 # 2) the max of all others
 thresh = min(NULL,
              max(c(15,mean(qc$RALIGN) - 3*sd(qc$RALIGN),
                    median(qc$RALIGN) - 3*mad(qc$RALIGN))))
 abline(v = thresh, col = "blue", lwd = 2, lty = 2)
 
 legend("topleft",legend = c("Hard","SD","MAD","Sufficient","Final"),
        lwd = 2, col = c("yellow","green","red","grey","blue"),
        lty = c(1,1,1,1,2), cex = .5)
```

<img src="/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M4B/thresh-1.png" title="plot of chunk thresh" alt="plot of chunk thresh" width="600px" height="400px" />
 
 We see here that the 3rd "MAD" threshold exceeds the first two thresholds
 ("Hard" and "SD"), and as the "Sufficient" threshold is NULL
 `metric_sample_filter` settles for the the third threshold. If the "Sufficient"
 threshold was not NULL and was exceeded by any of the other three thresholds
 ("Hard","SD","MAD"), `metric_sample_filter` would settle for the "Sufficient"
 threshold. Note also that if `mixture=TRUE` an additional criterion is
 considered: distributions may be fit to a two-component mixture model, and a
 threshold is defined with respect to the mean and standard deviation of the
 "best" component.
 
 As `metric_sample_filter` relies on a maximum of candidate thresholds, we
 recommend users treat this function as a stringent sample filter.
 
### Applying the sample filter
 
 With the `metric_sample_filter` output in hand, it is fairly straightforward to
 remove the one "poor" sample from our study:
 

```r
 goodDat = fluidigm[,mfilt]
 
 # Final Gene Filtering: Highly expressed in at least 5 cells
 num_reads = quantile(assay(fluidigm)[assay(fluidigm) > 0])[4]
 num_cells = 5
 is_quality = rowSums(assay(fluidigm) >= num_reads ) >= num_cells
```
 
### Running and Scoring Normalization Workflows with `scone`
 
 Not only does `scone` normalize expression data, but it also provides a
 framework for evaluating the performance of normalization workflows.
 
#### Creating a SconeExperiment Object
 
 Prior to running main `scone` function we will want to define a
 `SconeExperiment` object that contains the primary expression data,
 experimental metadata, and control gene sets.
 

```r
 # Expression Data (Required)
 expr = assay(goodDat)[is_quality,]
 
 # Biological Origin - Variation to be preserved (Optional)
 bio = factor(colData(goodDat)$Biological_Condition)
 
 # Processed Alignment Metrics - Variation to be removed (Optional)
 qc = colData(goodDat)[,metadata(goodDat)$which_qc]
 ppq = scale(qc[,apply(qc,2,sd) > 0],center = TRUE,scale = TRUE)
 
 # Positive Control Genes - Prior knowledge of DE (Optional)
 poscon = intersect(rownames(expr),strsplit(paste0("ALS2, CDK5R1, CYFIP1,",
                                                   " DPYSL5, FEZ1, FEZ2, ",
                                                   "MAPT, MDGA1, NRCAM, ",
                                                   "NRP1, NRXN1, OPHN1, ",
                                                   "OTX2, PARD6B, PPT1, ",
                                                   "ROBO1, ROBO2, RTN1, ",
                                                   "RTN4, SEMA4F, SIAH1, ",
                                                   "SLIT2, SMARCA1, THY1, ",
                                                   "TRAPPC4, UBB, YWHAG, ",
                                                   "YWHAH"),split = ", ")[[1]])
 
 # Negative Control Genes - Uniformly expressed transcripts (Optional)
 negcon = intersect(rownames(expr),hk)
 
 # Creating a SconeExperiment Object
 my_scone <- SconeExperiment(expr,
                 qc=ppq, bio = bio,
                 negcon_ruv = rownames(expr) %in% negcon,
                 poscon = rownames(expr) %in% poscon
 )
```
 
### Defining Normalization Modules
 
 Before we can decide which workflows (normalizations) we will want to compare,
 we will also need to define the types of scaling functions we will consider in
 the comparison of normalizations:
 

```r
 ## ----- User-defined function: Dividing by number of detected genes -----
 
 EFF_FN = function (ei)
 {
   sums = colSums(ei > 0)
   eo = t(t(ei)*sums/mean(sums))
   return(eo)
 }
 
 ## ----- Scaling Argument -----
 
 scaling=list(none=identity, # Identity - do nothing
 
              eff = EFF_FN, # User-defined function
 
              sum = SUM_FN, # SCONE library wrappers...
              tmm = TMM_FN, 
              uq = UQ_FN,
              fq = FQT_FN,
              deseq = DESEQ_FN)
```
 
 If imputation is to be included in the comparison, imputation arguments must
 also be provided by the user:
 

```r
 # Simple FNR model estimation with SCONE::estimate_ziber
 fnr_out = estimate_ziber(x = expr, bulk_model = TRUE,
                          pos_controls = rownames(expr) %in% hk,
                          maxiter = 10000)
 
 ## ----- Imputation List Argument -----
 imputation=list(none=impute_null, # No imputation
                 expect=impute_expectation) # Replace zeroes
 
 ## ----- Imputation Function Arguments -----
 # accessible by functions in imputation list argument
 impute_args = list(p_nodrop = fnr_out$p_nodrop, mu = exp(fnr_out$Alpha[1,]))
 
 my_scone <- scone(my_scone,
                 imputation = imputation, impute_args = impute_args,
                 scaling=scaling,
                 k_qc=3, k_ruv = 3,
                 adjust_bio="no",
                 run=FALSE)
```
 
 Note, that because the imputation step is quite slow, we do not run it here,
 but will run scone without imputation.
 
### Selecting SCONE Workflows
 
 The main `scone` method arguments allow for a lot of flexibility, but a user
 may choose to run very specific combinations of modules. For this purpose,
 `scone` can be run in `run=FALSE` mode, generating a list of workflows to be
 performed and storing this list within a `SconeExperiment` object. After
 running this command the list can be extracted using the `get_params` method.
 

```r
 my_scone <- scone(my_scone,
                 scaling=scaling,
                 k_qc=3, k_ruv = 3,
                 adjust_bio="no",
                 run=FALSE)
 
 head(get_params(my_scone))
```

```
##                                 imputation_method scaling_method uv_factors adjust_biology
## none,none,no_uv,no_bio,no_batch              none           none      no_uv         no_bio
## none,eff,no_uv,no_bio,no_batch               none            eff      no_uv         no_bio
## none,sum,no_uv,no_bio,no_batch               none            sum      no_uv         no_bio
## none,tmm,no_uv,no_bio,no_batch               none            tmm      no_uv         no_bio
## none,uq,no_uv,no_bio,no_batch                none             uq      no_uv         no_bio
## none,fq,no_uv,no_bio,no_batch                none             fq      no_uv         no_bio
##                                 adjust_batch
## none,none,no_uv,no_bio,no_batch     no_batch
## none,eff,no_uv,no_bio,no_batch      no_batch
## none,sum,no_uv,no_bio,no_batch      no_batch
## none,tmm,no_uv,no_bio,no_batch      no_batch
## none,uq,no_uv,no_bio,no_batch       no_batch
## none,fq,no_uv,no_bio,no_batch       no_batch
```
 
 In the call above, we have set the following parameter arguments:
 
 * k_ruv = 3. 
   The maximum number of RUVg factors to consider.
 * k_qc = 3. 
   The maximum number of quality PCs (QPCs) to be included in a linear model,
   analogous to RUVg normalization. The qc argument must be provided.
 * adjust_bio = "no." Biological origin will NOT be included in RUVg or QPC
   regression models. The bio argument will be provided for evaluation purposes.
 
 These arguments translate to the following set of options:
 

```r
 apply(get_params(my_scone),2,unique)
```

```
## $imputation_method
## [1] "none"
## 
## $scaling_method
## [1] "none"  "eff"   "sum"   "tmm"   "uq"    "fq"    "deseq"
## 
## $uv_factors
## [1] "no_uv"   "ruv_k=1" "ruv_k=2" "ruv_k=3" "qc_k=1"  "qc_k=2"  "qc_k=3" 
## 
## $adjust_biology
## [1] "no_bio"
## 
## $adjust_batch
## [1] "no_batch"
```
 
 Some scaling methods, such as scaling by gene detection rate (`EFF_FN()`), will
 not make sense within the context of imputed data, as imputation replaces
 zeroes with non-zero values. We can use the `select_methods` method to produce
 a `SconeExperiment` object initialized to run only meaningful normalization
 workflows.
 

```r
 is_screened = ((get_params(my_scone)$imputation_method == "expect") &
                  (get_params(my_scone)$scaling_method %in% c("none",
                                                              "eff")))
 
 my_scone = select_methods(my_scone,
                           rownames(get_params(my_scone))[!is_screened ])
```
 
### Calling `scone` with `run=TRUE`
 
 Now that we have selected our workflows, we can run `scone` in `run=TRUE` mode.
 As well as arguments used in `run=FALSE` mode, this mode relies on a few
 additional arguments. In order to understand these arguments, we must first
 understand the 8 metrics used to evaluate each normalization. The first 6
 metrics rely on a reduction of the normalized data down to 3 dimensions via PCA
 (default). Each metric is taken to have a positive (higher is better) or
 negative (lower is better) signature.
 
 * BIO_SIL: Preservation of Biological Difference. 
   The average silhouette width of clusters defined by `bio`, defined with
   respect to a Euclidean distance metric over the first 3 expression PCs.
   Positive signature.
 * BATCH_SIL: Removal of Batch Structure.
   The average silhouette width of clusters defined by `batch`, defined with
   respect to a Euclidean distance metric over the first 3 expression PCs.
   Negative signature.
 * PAM_SIL: Preservation of Single-Cell Heterogeneity.
   The maximum average silhouette width of clusters defined by PAM clustering,
   defined with respect to a Euclidean distance metric over the first 3 
   expression PCs. 
   Positive signature.
 * EXP_QC_COR: Removal of Alignment Artifacts. 
   Maximum squared Spearman correlation between first 3 expression PCs and 
   first `k_qc` QPCs. 
   Negative signature.
 * EXP_UV_COR: Removal of Expression Artifacts. 
   Maximum squared Spearman correlation between first 3 expression PCs and 
   first 3 PCs of the negative control (specified by `eval_negcon` or 
   `ruv_negcon` by default) sub-matrix of the original (raw) data. 
   Negative signature.
 * EXP_WV_COR: Preservation of Biological Variance.
   Maximum squared Spearman correlation between first 3 expression PCs and
   first 3 PCs of the positive control (specified by `eval_poscon`) 
   sub-matrix of the original (raw) data. 
   Positive signature.
 * RLE_MED: Reduction of Global Differential Expression.
   The mean squared-median Relative Log Expression (RLE). 
   Negative signature.
 * RLE_IQR: Reduction of Global Differential Variability.
   The variance of the inter-quartile range (IQR) of the RLE. 
   Negative signature.
 

```r
 BiocParallel::register(
   BiocParallel::SerialParam()
 ) # Register BiocParallel Serial Execution
 
 my_scone <- scone(my_scone,
                   scaling=scaling,
                   run=TRUE,
                   eval_kclust = 2:6,stratified_pam = TRUE,
                   return_norm = "in_memory",
                   zero = "preadjust")
```
 
 In the call above, we have set the following parameter arguments:
 
 * eval_kclust = 2:6. 
   For PAM_SIL, range of k (# of clusters) to use when computing 
   maximum average silhouette width of PAM clusterings.
 * stratified_pam = TRUE.
   For PAM_SIL, apply separate PAM clusterings to each biological 
   batch rather than across all batches. Average is weighted by 
   batch group size.
 * return_norm = "in_memory".
   Store all normalized matrices in addition to evaluation data. 
   Otherwise normalized data is not returned in the resulting 
   object.
 * zero = "preadjust". 
   Restore data entries that are originally zeroes back to zero
   after the scaling step.
 
 The output will contain various updated elements:
 

```r
 # View Metric Scores
 head(get_scores(my_scone))
```

```
##                                     BIO_SIL   PAM_SIL EXP_QC_COR EXP_UV_COR EXP_WV_COR     RLE_MED
## none,fq,no_uv,no_bio,no_batch    0.30191470 0.5558012 -0.7282716 -0.7792265  0.5828036 -0.01333402
## none,fq,qc_k=2,no_bio,no_batch   0.17493117 0.4127906 -0.4402420 -0.2854710  0.4211321 -0.01711373
## none,fq,qc_k=1,no_bio,no_batch   0.16855305 0.4686958 -0.5893234 -0.3024496  0.4254814 -0.01835348
## none,uq,ruv_k=1,no_bio,no_batch  0.26899478 0.5595772 -0.5381136 -0.3694012  0.3341094 -0.06235379
## none,fq,qc_k=3,no_bio,no_batch   0.04594019 0.4336716 -0.2257120 -0.2631728  0.3205924 -0.02772356
## none,tmm,ruv_k=1,no_bio,no_batch 0.19590108 0.4902965 -0.3196080 -0.3979276  0.6215272 -0.07162959
##                                      RLE_IQR
## none,fq,no_uv,no_bio,no_batch    -0.04482276
## none,fq,qc_k=2,no_bio,no_batch   -0.06635901
## none,fq,qc_k=1,no_bio,no_batch   -0.08912365
## none,uq,ruv_k=1,no_bio,no_batch  -0.11459679
## none,fq,qc_k=3,no_bio,no_batch   -0.07460980
## none,tmm,ruv_k=1,no_bio,no_batch -0.45831856
```

```r
 # View Mean Score Rank
 head(get_score_ranks(my_scone))
```

```
##    none,fq,no_uv,no_bio,no_batch   none,fq,qc_k=2,no_bio,no_batch   none,fq,qc_k=1,no_bio,no_batch 
##                         34.42857                         32.42857                         32.00000 
##  none,uq,ruv_k=1,no_bio,no_batch   none,fq,qc_k=3,no_bio,no_batch none,tmm,ruv_k=1,no_bio,no_batch 
##                         31.71429                         31.42857                         30.85714
```

```r
 # Extract normalized data from top method
 out_norm = get_normalized(my_scone,
                           method = rownames(get_params(my_scone))[1])
```
 
 `get_scores` returns the 8 raw metrics for each normalization multiplied by
 their signature - or "scores." `get_score_ranks` returns the mean score rank
 for each normalization. Both of these are sorted in decreasing order by mean
 score rank. Finally `get_normalized` returns the normalized expression data for
 the requested method. If the normalized data isn't stored in the object it will
 be recomputed.
 
### Step 3: Selecting a normalization for downstream analysis
 
 Based on our sorting criteria, it would appear that
 `none,fq,no_uv,no_bio,no_batch` performs well compared to other normalization
 workflows. A useful way to visualize this method with respect to others is the
 `biplot_color` function
 

```r
 pc_obj = prcomp(apply(t(get_scores(my_scone)),1,rank),
                 center = TRUE,scale = FALSE)
 bp_obj = biplot_color(pc_obj,y = -get_score_ranks(my_scone),expand = .6)
```

![plot of chunk biplot_color](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M4B/biplot_color-1.png)
 
 We have colored each point above according the corresponding method's mean
 score rank (yellow vs blue ~ good vs bad), and we can see that workflows span a
 continuum of metric performance. Most importantly - and perhaps to no surprise
 - there is evidence of strong trade-offs between i) Preserving clustering and
 wanted variation and ii) removing unwanted variation. At roughly 90 degrees to
 this axis is a direction in which distributional properties of relative
 log-expression (RLE_MED and RLE_IQR) improve. Let's visualize the
 top-performing method and it's relation to un-normalized data ("no-op"
 normalization):
 

```r
 bp_obj = biplot_color(pc_obj,y = -get_score_ranks(my_scone),expand = .6)
 
 points(t(bp_obj[1,]), pch = 1, col = "red", cex = 1)
 points(t(bp_obj[1,]), pch = 1, col = "red", cex = 1.5)
 
 points(t(bp_obj[rownames(bp_obj) == "none,none,no_uv,no_bio,no_batch",]),
        pch = 1, col = "blue", cex = 1)
 points(t(bp_obj[rownames(bp_obj) == "none,none,no_uv,no_bio,no_batch",]),
        pch = 1, col = "blue", cex = 1.5)
 
 arrows(bp_obj[rownames(bp_obj) == "none,none,no_uv,no_bio,no_batch",][1],
        bp_obj[rownames(bp_obj) == "none,none,no_uv,no_bio,no_batch",][2],
        bp_obj[1,][1],
        bp_obj[1,][2],
        lty = 2, lwd = 2)
```

![plot of chunk biplot_color4](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M4B/biplot_color4-1.png)
 
 The arrow we've added to the plot traces a line from the "no-op" normalization
 to the top-ranked normalization in SCONE. We see that SCONE has selected a
 method in-between the two extremes, reducing the signal of unwanted variation
 while preserving biological signal.
 
 Finally, another useful function for browsing results is `sconeReport`. This
 function launches a shiny app for evaluating performance of specific
 normalization workflows.
 

```r
 # Methods to consider
 scone_methods = c(rownames(get_params(my_scone))[1:12],
                   "none,none,no_uv,no_bio,no_batch")
 
 # Shiny app
 sconeReport(my_scone,methods = scone_methods,
             qc = ppq,
             bio = bio,
             negcon = negcon, poscon = poscon)
```
 
 
 
 ### Parameter settings:
   * WRKDIR = /mnt100/home/Dropbox/SingleCell/Jun2017/R
   * FN = M4B-sconeTutorial
   * Scripts = Scripts
   * RUN DATE = Sun Jun 25 20:31:55 2017
 
 

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
##  [1] scRNAseq_1.2.0             RColorBrewer_1.1-2         scone_1.0.0               
##  [4] SummarizedExperiment_1.6.3 DelayedArray_0.2.7         matrixStats_0.52.2        
##  [7] Biobase_2.36.2             GenomicRanges_1.28.3       GenomeInfoDb_1.12.2       
## [10] IRanges_2.10.2             S4Vectors_0.14.3           BiocGenerics_0.22.0       
## [13] BiocParallel_1.10.1        knitr_1.16                 rmarkdown_1.6             
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
##  [79] stringr_1.2.0            robustbase_0.92-7        caTools_1.17.1          
##  [82] purrr_0.2.2.2            EDASeq_2.10.0            mime_0.5                
##  [85] scran_1.4.4              R.oo_1.21.0              biomaRt_2.32.1          
##  [88] compiler_3.4.0           beeswarm_0.2.3           plotly_4.7.0            
##  [91] tibble_1.3.3             statmod_1.4.30           geneplotter_1.54.0      
##  [94] stringi_1.1.5            highr_0.6                GenomicFeatures_1.28.3  
##  [97] RSpectra_0.12-0          lattice_0.20-35          trimcluster_0.1-2       
## [100] Matrix_1.2-10            tensorA_0.36             data.table_1.10.4       
## [103] bitops_1.0-6             httpuv_1.3.3             rtracklayer_1.36.3      
## [106] R6_2.2.2                 latticeExtra_0.6-28      hwriter_1.3.2           
## [109] ShortRead_1.34.0         gridExtra_2.2.1          KernSmooth_2.23-15      
## [112] vipor_0.4.5              codetools_0.2-15         boot_1.3-19             
## [115] energy_1.7-0             MASS_7.3-47              gtools_3.5.0            
## [118] assertthat_0.2.0         rhdf5_2.20.0             rjson_0.2.15            
## [121] pkgmaker_0.22            rprojroot_1.2            RUVSeq_1.10.0           
## [124] GenomicAlignments_1.12.1 Rsamtools_1.28.0         GenomeInfoDbData_0.99.0 
## [127] diptest_0.75-7           grid_3.4.0               tidyr_0.6.3             
## [130] class_7.3-14             segmented_0.5-2.1        shiny_1.0.3             
## [133] ggbeeswarm_0.5.3
```

