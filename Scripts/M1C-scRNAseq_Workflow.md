### A step-by-step workflow for low-level analysis of single-cell RNA-seq data 

<!-- Part III - Normalizing based on spike-in coverage -->

<style type="text/css"> body, td { font-size: 14px; } code.r{ font-size: 12px; } pre { font-size: 12px } </style>

[bioconductor workflows simpleSingleCell](https://www.bioconductor.org/help/workflows/simpleSingleCell/)



<!-- ***************************************************** -->


<!-- ***************************************************** -->


<!-- ***************************************************** -->



<!-- ***************************************************** -->


## Normalizing based on spike-in coverage


Scaling normalization strategies for scRNA-seq data can be broadly divided into two classes. The first class assumes that there exists a subset of genes that are not DE between samples, as previously described. The second class uses the fact that the same amount of spike-in RNA was added to each cell. Differences in the coverage of the spike-in transcripts can only be due to cell-specific biases, e.g., in capture efficiency or sequencing depth. Scaling normalization is then applied to equalize spike-in coverage across cells.
 
The choice between these two normalization strategies depends on the biology of the cells and the features of interest. If the majority of genes are expected to be DE and there is no reliable house-keeping set, spike-in normalization may be the only option for removing cell-specific biases. Spike-in normalization should also be used if differences in the total RNA content of individual cells are of interest. In any particular cell, an increase in the amount of endogenous RNA will not increase spike-in coverage (with or without library quantification). Thus, the former will not be represented as part of the bias in the latter, which means that the effects of total RNA content on expression will not be removed upon scaling. With non-DE normalization, an increase in RNA content will systematically increase the expression of all genes in the non-DE subset, such that it will be treated as bias and removed.

We demonstrate the use of spike-in normalization on a dataset involving different cell types – namely, mouse embryonic stem cells (mESCs) and mouse embryonic fibroblasts (MEFs) (Islam et al. 2011). The count table was obtained from NCBI GEO as a supplementary file under the accession [GSE29087](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE29087). We load the counts into R and specify the rows corresponding to spike-in transcripts. The negative control wells do not contain any cells and are useful for quality control but need to be removed prior to downstream analysis.



```r
 counts <- read.table(file.path(EXT_DATA,"GSE29087_L139_expression_tab.txt.gz"), 
      colClasses=c(list("character", NULL, NULL, NULL, NULL, NULL, NULL), 
      rep("integer", 96)), skip=6, sep='\t', row.names=1)
 print(dim(counts))

 GSE29087.sce <- newSCESet(countData=counts)
 GSE29087.sce$grouping <- rep(c("mESC", "MEF", "Neg"), c(48, 44, 4))
 GSE29087.sce <- GSE29087.sce[,GSE29087.sce$grouping!="Neg"] # Removing negative control wells.
 GSE29087.sce <- calculateQCMetrics(GSE29087.sce, feature_controls=list(spike=grep("SPIKE", rownames(counts))))
 setSpike(GSE29087.sce) <- "spike"


 # save
 save(GSE29087.sce, file=file.path(WRKDIR, 'Data', 'GSE29087.sce'))
```
<!-- *************************************************************** -->
 
We then apply the **computeSpikeFactors** method to estimate size factors for all cells. This method computes the total count over all spike-in transcripts in each cell, and calculates size factors to equalize the total spike-in count across cells. Here, we set **general.use=TRUE** as we intend to apply the spike-in factors to all counts.


```r
 load(file.path(WRKDIR, 'Data', 'GSE29087.sce'))

 print(ncol(GSE29087.sce))
```

```
## Samples 
##      92
```

```r
 print(nrow(GSE29087.sce))
```

```
## Features 
##    22936
```

```r
 GSE29087.sce <- computeSpikeFactors(GSE29087.sce, general.use=TRUE)
```
<!-- *************************************************************** -->

Applying **normalize** will use the spike-in-based size factors to compute normalized log-expression values. Unlike in the previous analyses, we do not have to set separate size factors for the spike-in transcripts. This is because the relevant factors are already being used for all genes and spike-in transcripts when **general.use=TRUE**. (The exception is if the experiment uses multiple spike-in sets that behave differently and need to be normalized separately.)


```r
 GSE29087.sce <- normalize(GSE29087.sce)
```
<!-- *************************************************************** -->
 
For comparison, we also compute the deconvolution size factors and plot them against the spike-in factors. We observe a negative correlation between the two sets of values (Figure 27). This is because MEFs contain more endogenous RNA, which reduces the relative spike-in coverage in each library (thereby decreasing the spike-in size factors) but increases the coverage of endogenous genes (thus increasing the deconvolution size factors). If the spike-in size factors were applied to the counts, the expression values in MEFs would be scaled up while expression in mESCs would be scaled down. However, the opposite would occur if deconvolution size factors were used.

Figure 27:
* Size factors from spike-in normalization, plotted against the size factors from deconvolution for all cells in the mESC/MEF dataset. Axes are shown on a log-scale, and cells are coloured according to their identity. Deconvolution size factors were computed with small pool sizes owing to the low number of cells of each type.


```r
 colours <- c(mESC="red", MEF="grey")
 deconv.sf <- computeSumFactors(GSE29087.sce, sf.out=TRUE, 
             cluster=GSE29087.sce$grouping, sizes=1:4*10)
 plot(sizeFactors(GSE29087.sce), deconv.sf, 
      col=colours[GSE29087.sce$grouping], pch=16, log="xy", 
      xlab="Size factor (spike-in)", ylab="Size factor (deconvolution)")
 legend("bottomleft", col=colours, legend=names(colours), pch=16)
```

![plot of chunk normplotspike](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M1C/normplotspike-1.png)
<!-- *************************************************************** -->

## Blocking on the cell cycle phase


Cell cycle phase is usually uninteresting in studies focusing on other aspects of biology. However, the effects of cell cycle on the expression profile can mask other effects and interfere with the interpretation of the results. This cannot be avoided by simply removing cell cycle marker genes, as the cell cycle can affect a substantial number of other transcripts (Buettner et al. 2015). Rather, more sophisticated strategies are required, one of which is demonstrated below using data from a study of T Helper 2 ($T_H2$) cells (Mahata et al. 2014). Buettner et al. (2015) have already applied quality control and normalized the data, so we can use them directly as log-expression values (accessible as Supplementary Data 1 of https://dx.doi.org/10.1038/nbt.3102).



```r
 incoming <- as.data.frame(read_excel(file.path(EXT_DATA,"nbt.3102-S7.xlsx"), sheet=1))
 print(dim(incoming))
```

```
## [1]   81 7074
```

```r
 rownames(incoming) <- incoming[,1]
 incoming <- incoming[,-1]
 incoming <- incoming[,!duplicated(colnames(incoming))] # Remove duplicated genes.
 TH2.sce <- newSCESet(exprsData=t(incoming))
```
<!-- *************************************************************** -->

We empirically identify the cell cycle phase using the pair-based classifier in **cyclone**. The majority of cells in Figure 28 seem to lie in G1 phase, with small numbers of cells in the other phases.


```r
 mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
 library(org.Mm.eg.db)
 fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))
 set.seed(100)
```
<!-- *************************************************************** -->

#### Figure 28:
* Cell cycle phase scores from applying the pair-based classifier on the T~H~2 dataset, where each point represents a cell.


```r
 anno <- select(org.Mm.eg.db, keys=rownames(TH2.sce), keytype="SYMBOL", column="ENSEMBL")
 ensembl <- anno$ENSEMBL[match(rownames(TH2.sce), anno$SYMBOL)]
 assignments <- cyclone(TH2.sce, mm.pairs, gene.names=ensembl, assay="exprs")
 plot(assignments$score$G1, assignments$score$G2M, xlab="G1 score", ylab="G2/M score", pch=16)
 abline(0,1, col='red')
 abline(h=0.5, v=0.5, col='blue', lty=2)
```

![plot of chunk phaseplot](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M1C/phaseplot-1.png)
<!-- *************************************************************** -->
 
We can block directly on the phase scores in downstream analyses. This is more graduated than using a strict assignment of each cell to a specific phase, as the magnitude of the score considers the uncertainty of the assignment. The phase covariates in the design matrix will absorb any phase-related effects on expression such that they will not affect estimation of the effects of other experimental factors. Users should also ensure that the phase score is not confounded with other factors of interest. For example, model fitting is not possible if all cells in one experimental condition are in one phase, and all cells in another condition are in a different phase.



```r
 design <- model.matrix(~ G1 + G2M, assignments$score)
 fit.block <- trendVar(TH2.sce, use.spikes=NA, trend="loess", design=design)
 dec.block <- decomposeVar(TH2.sce, fit.block)
```
<!-- *************************************************************** -->
 
For analyses that do not use design matrices, we remove the cell cycle effect directly from the expression values using **removeBatchEffect**. The result of this procedure is visualized with some PCA plots in Figure 29. Before removal, the distribution of cells along the first two principal components is strongly associated with their G1 and G2/M scores. This is no longer the case after removal, which suggests that the cell cycle effect has been mitigated.

#### Figure 29:
* PCA plots before (left) and after (right) removal of the cell cycle effect in the T~H~2 dataset. Each cell is represented by a point with colour and size determined by the G1 and G2/M scores, respectively. Only HVGs were used to construct each plot."----


```r
 # Finding HVGs without blocking on phase score.
 fit <- trendVar(TH2.sce, use.spikes=NA, trend="loess") 
 dec <- decomposeVar(TH2.sce, fit)
 top.hvgs <- which(dec$FDR <= 0.05 & dec$bio >= 0.5)
 TH2.sce$G1score <- assignments$score$G1
 TH2.sce$G2Mscore <- assignments$score$G2M
 out <- plotPCA(TH2.sce, feature_set=top.hvgs, colour_by="G1score", size_by="G2Mscore") + 
 fontsize + ggtitle("Before removal")
 
 # Using HVGs after blocking on the phase score.
 top.hvgs2 <- which(dec.block$FDR <= 0.05 & dec.block$bio >= 0.5)
 norm_exprs(TH2.sce) <- removeBatchEffect(exprs(TH2.sce), covariates=assignments$score[,c("G1", "G2M")])
 out2 <- plotPCA(TH2.sce, exprs_values="norm_exprs", feature_set=top.hvgs2, colour_by="G1score", 
 size_by="G2Mscore") + fontsize + ggtitle("After removal")
 multiplot(out, out2, cols=2)
```

![plot of chunk pcaplot](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M1C/pcaplot-1.png)
<!-- **************************************************************** -->

As an aside, this dataset contains cells at various stages of differentiation (Mahata et al. 2014). This is an ideal use case for diffusion maps which perform dimensionality reduction along a continuous process. In Figure 30, cells are arranged along a trajectory in the low-dimensional space. The first diffusion component is likely to correspond to $T_H2$ differentiation, given that a key regulator Gata3 (Zhu et al. 2006) changes in expression from left to right.


#### Figure 30:
* A diffusion map for the T~H~2 dataset, where each cell is coloured by its expression of _Gata3_.


```r
 plotDiffusionMap(TH2.sce, exprs_values="norm_exprs", colour_by="Gata3") + fontsize
```

![plot of chunk diffusionPlot](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M1C/diffusionPlot-1.png)


## Extracting annotation from Ensembl identifiers

Feature-counting tools typically report genes in terms of standard identifiers from Ensembl or Entrez. These identifiers are used as they are unambiguous and highly stable. However, they are difficult to interpret compared to the gene symbols which are more commonly used in the literature. We can easily convert from one to the other using annotation packages like **org.Mm.eg.db**. This is demonstrated below for Ensembl identifiers in a mESC dataset (Kołodziejczyk et al. 2015) obtained from http://www.ebi.ac.uk/teichmann-srv/espresso. The **select** call extracts the specified data from the annotation object, and the **match** call ensures that the first gene symbol is used if multiple symbols correspond to a single Ensembl identifier.



```r
 incoming <- read.table(file.path(EXT_DATA,"counttable_es.csv"), header=TRUE, row.names=1)
 my.ids <- rownames(incoming)
 anno <- select(org.Mm.eg.db, keys=my.ids, keytype="ENSEMBL", column="SYMBOL")
 anno <- anno[match(my.ids, anno$ENSEMBL),]
 head(anno)
```

```
##              ENSEMBL SYMBOL
## 1 ENSMUSG00000000001  Gnai3
## 2 ENSMUSG00000000003   Pbsn
## 3 ENSMUSG00000000028  Cdc45
## 4 ENSMUSG00000000031    H19
## 5 ENSMUSG00000000037  Scml2
## 6 ENSMUSG00000000049   Apoh
```


To identify which rows correspond to mitochondrial genes, we need to use extra annotation describing the genomic location of each gene. For Ensembl, this involves using the TxDb.Mmusculus.UCSC.mm10.ensGene package.


```r
 suppressMessages(require(TxDb.Mmusculus.UCSC.mm10.ensGene))

 location <- select(TxDb.Mmusculus.UCSC.mm10.ensGene, keys=my.ids, 
     column="CDSCHROM", keytype="GENEID")
 location <- location[match(my.ids, location$GENEID),]
 is.mito <- location$CDSCHROM == "chrM" & !is.na(location$CDSCHROM)
 sum(is.mito)
```

```
## [1] 13
```

Identification of rows that correspond to spike-in transcripts is much easier, given that the ERCC spike-ins were used.


```r
 is.spike <- grepl("^ERCC", my.ids)
 sum(is.spike)
```

```
## [1] 92
```
All of this information can be consolidated into a **SCESet** object for further manipulation. Alternatively, annotation from BioMart resources can be directly added to the object using the **getBMFeatureAnnos** function from **scater**.



```r
 anno <- anno[,-1,drop=FALSE]
 rownames(anno) <- my.ids
 mESC.sce <- newSCESet(countData=incoming, featureData=AnnotatedDataFrame(anno))
 mESC.sce <- calculateQCMetrics(mESC.sce, feature_controls=list(ERCC=is.spike))
 setSpike(mESC.sce) <- "ERCC"
```

We filter out rows that do not correspond to endogenous genes or spike-in transcripts. This will remove rows containing mapping statistics such as the number of unaligned or unassigned reads, which would be misleading if treated as gene expression values. The object is then ready for downstream analyses as previously described.

 

```r
 mESC.sce <- mESC.sce[grepl("ENSMUS", rownames(mESC.sce)) | isSpike(mESC.sce),]
 dim(mESC.sce)
```

```
## Features  Samples 
##    38653      704
```
 

 
## Conclusions

This workflow provides a step-by-step guide for performing basic analyses of single-cell RNA-seq data in R. It provides instructions for a number of low-level steps such as quality control, normalization, cell cycle phase assignment, data exploration, HVG and marker gene detection, and clustering. This is done with a number of different datasets to provide a range of usage examples. The workflow is modular so individual steps can be substituted with alternative methods according to user preferences. In addition, the processed data can be easily used for higher-level analyses with other Bioconductor packages. We anticipate that this workflow will assist readers in assembling analyses of their own scRNA-seq data.

## Software availability

All software packages used in this workflow are publicly available from the Comprehensive R Archive Network (https://cran.r-project.org) or the Bioconductor project (http://bioconductor.org). The specific version numbers of the packages used are shown below, along with the version of the R installation. Version numbers of all Bioconductor packages correspond to release version 3.5 of the Bioconductor project. The workflow takes less than an hour to run on a desktop computer with 8 GB of memory.

 
### Parameter settings:
   * WRKDIR = /mnt100/home/Dropbox/SingleCell/Jun2017/R
   * FN = M1C-scRNAseq_Workflow
   * Scripts = Scripts
   * EXT_DATA = /mnt100/home/Dropbox/SingleCell//Data
   * TABLE_DIR = /mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/tables/M1C 
   * RUN DATE = Wed Jun 21 23:22:12 2017
 
 

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
##  [1] TxDb.Mmusculus.UCSC.mm10.ensGene_3.4.0 GenomicFeatures_1.28.3                
##  [3] org.Mm.eg.db_3.4.1                     AnnotationDbi_1.38.1                  
##  [5] scran_1.4.4                            BiocParallel_1.10.1                   
##  [7] destiny_2.4.0                          mvoutlier_2.0.8                       
##  [9] sgeostat_1.0-27                        Rtsne_0.13                            
## [11] readxl_1.0.0                           R.utils_2.5.0                         
## [13] R.oo_1.21.0                            R.methodsS3_1.7.1                     
## [15] BiocStyle_2.4.0                        scater_1.4.0                          
## [17] ggplot2_2.2.1                          readr_1.1.1                           
## [19] irlba_2.2.1                            Matrix_1.2-10                         
## [21] plyr_1.8.4                             tidyr_0.6.3                           
## [23] reshape2_1.4.2                         stringr_1.2.0                         
## [25] vcd_1.4-3                              GenomicRanges_1.28.3                  
## [27] GenomeInfoDb_1.12.2                    IRanges_2.10.2                        
## [29] S4Vectors_0.14.3                       Biobase_2.36.2                        
## [31] BiocGenerics_0.22.0                    edgeR_3.18.1                          
## [33] limma_3.32.2                           RColorBrewer_1.1-2                    
## [35] data.table_1.10.4                      magrittr_1.5                          
## [37] dplyr_0.7.0                            knitr_1.16                            
## [39] rmarkdown_1.6                         
## 
## loaded via a namespace (and not attached):
##   [1] backports_1.1.0            Hmisc_4.0-3                RcppEigen_0.3.3.3.0       
##   [4] igraph_1.0.1               lazyeval_0.2.0             sp_1.2-4                  
##   [7] shinydashboard_0.6.1       splines_3.4.0              digest_0.6.12             
##  [10] htmltools_0.3.6            viridis_0.4.0              checkmate_1.8.2           
##  [13] memoise_1.1.0              cluster_2.0.6              Biostrings_2.44.1         
##  [16] matrixStats_0.52.2         xts_0.9-7                  colorspace_1.3-2          
##  [19] rrcov_1.4-3                RCurl_1.95-4.8             tximport_1.4.0            
##  [22] lme4_1.1-13                survival_2.41-3            zoo_1.8-0                 
##  [25] glue_1.1.0                 gtable_0.2.0               zlibbioc_1.22.0           
##  [28] XVector_0.16.0             MatrixModels_0.4-1         DelayedArray_0.2.7        
##  [31] car_2.1-4                  kernlab_0.9-25             prabclus_2.2-6            
##  [34] DEoptimR_1.0-8             SparseM_1.77               VIM_4.7.0                 
##  [37] scales_0.4.1               mvtnorm_1.0-6              DBI_0.7                   
##  [40] GGally_1.3.1               Rcpp_0.12.11               sROC_0.1-2                
##  [43] htmlTable_1.9              viridisLite_0.2.0          xtable_1.8-2              
##  [46] laeken_0.4.6               proxy_0.4-17               foreign_0.8-68            
##  [49] mclust_5.3                 Formula_1.2-1              DT_0.2                    
##  [52] htmlwidgets_0.8            FNN_1.1                    fpc_2.1-10                
##  [55] acepack_1.4.1              modeltools_0.2-21          reshape_0.8.6             
##  [58] XML_3.98-1.8               flexmix_2.3-14             nnet_7.3-12               
##  [61] dynamicTreeCut_1.63-1      locfit_1.5-9.1             labeling_0.3              
##  [64] rlang_0.1.1                munsell_0.4.3              cellranger_1.1.0          
##  [67] tools_3.4.0                RSQLite_1.1-2              pls_2.6-0                 
##  [70] evaluate_0.10              cvTools_0.3.2              yaml_2.1.14               
##  [73] robustbase_0.92-7          nlme_3.1-131               mime_0.5                  
##  [76] quantreg_5.33              biomaRt_2.32.1             compiler_3.4.0            
##  [79] pbkrtest_0.4-7             beeswarm_0.2.3             e1071_1.6-8               
##  [82] statmod_1.4.30             smoother_1.1               tibble_1.3.3              
##  [85] robCompositions_2.0.3      pcaPP_1.9-61               stringi_1.1.5             
##  [88] highr_0.6                  lattice_0.20-35            trimcluster_0.1-2         
##  [91] nloptr_1.0.4               lmtest_0.9-35              cowplot_0.7.0             
##  [94] bitops_1.0-6               rtracklayer_1.36.3         httpuv_1.3.3              
##  [97] latticeExtra_0.6-28        R6_2.2.2                   gridExtra_2.2.1           
## [100] vipor_0.4.5                boot_1.3-19                MASS_7.3-47               
## [103] assertthat_0.2.0           SummarizedExperiment_1.6.3 rhdf5_2.20.0              
## [106] rprojroot_1.2              rjson_0.2.15               GenomicAlignments_1.12.1  
## [109] Rsamtools_1.28.0           GenomeInfoDbData_0.99.0    diptest_0.75-7            
## [112] mgcv_1.8-17                hms_0.3                    rpart_4.1-11              
## [115] class_7.3-14               minqa_1.2.4                TTR_0.23-1                
## [118] shiny_1.0.3                base64enc_0.1-3            ggbeeswarm_0.5.3
```

