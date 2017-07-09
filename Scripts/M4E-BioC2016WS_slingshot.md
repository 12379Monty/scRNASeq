### BioC 2016 workshop - Analysis of single-cell RNA-seq data with R and Bioconductor

### III - Lineage Reconstruction

<style type="text/css"> body, td { font-size: 14px; } code.r{ font-size: 12px; } pre { font-size: 12px } </style>

[BioC 2016 workshop](https://github.com/drisso/bioc2016singlecell)
[slingshot](https://github.com/kstreet13/slingshot)


<!-- ***************************************************** -->


<!-- ***************************************************** -->



<!-- ***************************************************** -->



<!-- ***************************************************** -->


<!--
%\VignetteEngine{knitr::rmarkdown}
%\VignetteIndexEntry{slingshot Vignette}
-->



# Introduction

This is the third part of the Bioc2016 workshop "Analysis of single-cell RNA-seq data with R and Bioconductor."

In this part we will cover lineage reconstruction with the `Githubpkg("kstreet13/slingshot")` package.

The goal of `slingshot` is to use clusters to uncover global structure and convert this structure into smooth lineages represented by one-dimensional variables, often called ``pseudotime.'' We give tools for learning cluster relationships in an unsupervised or semi-supervised manner and constructing smooth curves representing each lineage, with visualization methods for each step.

### Basic `slingshot` analysis

The minimal input to `slingshot` is a matrix representing the cells in a reduced-dimensional space and a vector of clustering results. The analysis then procedes:

* Find connections between clusters with the `get_lineages` function, optionally specifying known start and end points.
* Construct smooth curves and pseudotime variables with the `get_curves` function.
* Assess the cluster connectivity with `plot_tree` and the curve stability with `plot_curves`.

Using connections between clusters to define global structure improves the stability of inferred lineages. And our use of smooth curves in place of piecewise-linear ones reduces variability in the inferred pseudotime vectors.

### An example dataset

We will take our inputs from the previous sections: the normalized counts matrix obtained with `scone` and the cluster assignments obtained with `clusterExperiment`, both of which can be loaded directly from the workshop package.


```r
 data('full_pca')
 
 ### Examine dimensionality reduction
 ###plot3d(pcaX, aspect = 'iso')
 pairs(pcaX[,1:3], asp = 1)
```

![plot of chunk datain](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M4E/datain-1.png)
 
### Step 1: Assign clusters to lineages with `get_lineages`
 
 The `get_lineages` function takes as input an `n x p` matrix and a vector of clustering results of length `n`. It then maps connections between adjacent clusters using a minimum spanning tree (MST) and identifies paths through these connections that represent potential lineages.
 
 This analysis can be performed in an entirely unsupervised manner or in a semi-supervised manner by specifying known beginning and end point clusters. We recommend that you specify a root cluster; this will have no effect on how the clusters are connected, but it will allow for nicer curves in datasets with a branching structure. Pre-specified end point clusters will be constrained to only one connection.
 

```r
 l1 <- get_lineages(pcaX, clus)
###plot_tree(pcaX, clus, l1, threeD = T)
 plot_tree(pcaX, clus, l1, dim = 3)
```

![plot of chunk lines_unsup](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M4E/lines_unsup-1.png)
 
 Running `get_lineages` with no supervision produces the connections shown above. Since no root cluster was specified, `slingshot` picked one of the leaf-node clusters to be the beginning, based on a simple parsimony rule. The root cluster is the leaf-node cluster connected by a green line.
 

```r
 l2 <- get_lineages(pcaX, clus, start.clus = 'm10')
###plot_tree(pcaX, clus, l2, threeD = T)
 plot_tree(pcaX, clus, l2, dim = 3)
```

![plot of chunk lines_sup_start](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M4E/lines_sup_start-1.png)
 
 When we specify a root cluster we get the same connections and the only difference is which line is drawn in green.
 

```r
 l3 <- get_lineages(pcaX, clus, start.clus = 'm10', end.clus = 'm17')
###plot_tree(pcaX, clus, l3, threeD = T)
 plot_tree(pcaX, clus, l3, dim = 3)
```

![plot of chunk lines_sup_end](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M4E/lines_sup_end-1.png)
 
 Here we demonstrate the ability to specify end point clusters, which puts a constraint on the connections. We now draw the MST subject to the constraint that given end point clusters must be leaves. Pre-specified end point clusters are connected by red lines.
 
 There are a few additional arguments we could have passed to `get_lineages` for more greater control:
 
 * `dist.fun` is a function for computing distances between clusters. The default is squared distance between cluster centers normalized by their joint covariance matrix.
 * `omega` is a granularity parameter, allowing the user to set an upper limit on connection distances. It takes values between 0 and 1 (or `Inf`), representing a percentage of the largest observed distance.
 * `distout` is a logical value, indicating whether the user wants the pairwise cluster distance matrix to be returned with the output.
 
 After constructing the MST, `get_lineages` identifies paths through the tree to designate as lineages. At this stage, a lineage will consist of an ordered set of cluster names, starting with the root cluster and ending with a leaf. The output of `get_lineages` is a list of these vectors, along with some additional information on how they were constructed.
 
### Step 2: Construct smooth lineages and order cells with `get_curves`
 
 In order to model development along these various lineages, we will construct smooth curves with the function `get_curves`. Using smooth curves based on all the cells eliminates the problem of cells projecting onto vertices of piece-wise linear trajectories and makes `slingshot` more robust to noise in the clustering results.
 
 In order to construct smooth lineages, `get_curves` follows an iterative process similar to that of principal curves presented in [@princurve]. When there is only a single lineage, the resulting curve is simply the principal curve through the center of the data, with one adjustment: the initial curve is constructed with the linear connections between cluster centers rather than the first prinicpal component of the data. This adjustment adds stability and typically hastens the algorithm's convergence.
 
 When there are two or more lineages, we add an additional step to the algorithm: averaging curves near shared cells. Both lineages should agree fairly well on cells that have yet to differentiate, so at each iteration we average the curves in the neighborhood of these cells. This increases the stability of the algorithm and produces smooth branching lineages.
 

```r
 crv <- get_curves(pcaX, clus, l2)
###plot_curves(pcaX, clus, crv, threeD = T)
 plot_curves(pcaX, clus, crv, dim = 3)
```

![plot of chunk curves](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M4E/curves-1.png)
 
 The output of `get_curves` is a list with one element per curve. Each element is an object of the `principal.curve` class, with the following slots:
 
 * `s`: the matrix of points that make up the curve. These correspond to the orthogonal projections of the data points, but oredered such the `lines(s)` will produce a smooth curve.
 * `tag`: the indices of the original data points in `s`.
 * `lambda`: arclengths of the points in `s` along the curve.
 * `dist`: the total squared distance between data points and their projections onto the curve.
 * `pseudotime`: the vector of pseudotime values along this lineage.
 
### Step 3: Find temporally expressed genes
 
 Typically, the next step will be to find genes that change their expression as a function of developmental time. This can be done using the full genes-by-samples data matrix, but we will use the subset consisting of the 1,000 most variable genes.
 

```r
 data('var_genes')
```
 
 For a quick analysis, we will regress each gene on the two pseudotime vectors we have generated, using a general additive model (GAM). This allows us to detect non-linear patterns in gene expression over developmental time.
 

```r
 gam.pval <- vector("list",length(crv))
 for(l in 1:length(crv)){
   t <- crv[[l]]$pseudotime
   y <- vargenes[, match(rownames(crv[[l]]$s[!is.na(t),]),rownames(pcaX))]
   t <- t[! is.na(t)]
   gam.pval[[l]] <- apply(y,1,function(z){
     d <- data.frame(z=z, t=t)
     tmp <- gam(z ~ lo(t), data=d)
     p <- summary(tmp)[4][[1]][1,5]
     p
   })
 }
```
 
 We can then pick out the top genes for each lineage and visualize their expression over developmental time with a heatmap.
 

```r
 topgenes1 <- names(sort(gam.pval[[1]], decreasing = F))[1:100]
 heatdata1 <- vargenes[rownames(vargenes) %in% topgenes1, 
                       order(crv[[1]]$pseudotime, na.last = NA)]
 heatclus1 <- clus[clus!='-1'][order(crv[[1]]$pseudotime, na.last = NA)]
 ce1 <- clusterExperiment(heatdata1, heatclus1, transformation=identity)
 plotHeatmap(ce1, clusterSamplesData="orderSamplesValue")
```

![plot of chunk heatmaps](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M4E/heatmaps-1.png)

```r
 topgenes2 <- names(sort(gam.pval[[2]], decreasing = F))[1:100]
 heatdata2 <- vargenes[rownames(vargenes) %in% topgenes2, 
                       order(crv[[2]]$pseudotime, na.last = NA)]
 heatclus2 <- clus[clus!='-1'][order(crv[[2]]$pseudotime, na.last = NA)]
 ce2 <- clusterExperiment(heatdata2, heatclus2, transformation=identity)
 plotHeatmap(ce2, clusterSamplesData="orderSamplesValue")
```

![plot of chunk heatmaps](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M4E/heatmaps-2.png)
 
 
### References
 
  
#### Parameter settings:
   * WRKDIR = /mnt100/home/Dropbox/SingleCell/Jun2017/R
   * FN = M4E-BioC2016WS_slingshot
   * Scripts = Scripts
   * RUN DATE = Wed Jun 28 01:34:54 2017
 
 

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
##  [1] splines   stats4    parallel  methods   stats     graphics  grDevices utils     datasets 
## [10] base     
## 
## other attached packages:
##  [1] gam_1.14-4                 foreach_1.4.3              clusterExperiment_1.2.0   
##  [4] SummarizedExperiment_1.6.3 DelayedArray_0.2.7         matrixStats_0.52.2        
##  [7] GenomicRanges_1.28.3       GenomeInfoDb_1.12.2        IRanges_2.10.2            
## [10] S4Vectors_0.14.3           slingshot_0.0.3-5          princurve_1.1-12          
## [13] Biobase_2.36.2             BiocGenerics_0.22.0        bioc2016singlecell_0.0.2  
## [16] devtools_1.13.2            BiocParallel_1.10.1        knitr_1.16                
## [19] rmarkdown_1.6             
## 
## loaded via a namespace (and not attached):
##  [1] nlme_3.1-131            bitops_1.0-6            bold_0.4.0              progress_1.1.2         
##  [5] doParallel_1.0.10       RColorBrewer_1.1-2      httr_1.2.1              rprojroot_1.2          
##  [9] prabclus_2.2-6          tools_3.4.0             backports_1.1.0         R6_2.2.2               
## [13] lazyeval_0.2.0          colorspace_1.3-2        ade4_1.7-6              trimcluster_0.1-2      
## [17] nnet_7.3-12             withr_1.0.2             prettyunits_1.0.2       gridExtra_2.2.1        
## [21] compiler_3.4.0          xml2_1.1.1              pkgmaker_0.22           diptest_0.75-7         
## [25] scales_0.4.1            DEoptimR_1.0-8          mvtnorm_1.0-6           robustbase_0.92-7      
## [29] NMF_0.20.6              stringr_1.2.0           digest_0.6.12           XVector_0.16.0         
## [33] htmltools_0.3.6         highr_0.6               limma_3.32.2            rlang_0.1.1            
## [37] howmany_0.3-1           jsonlite_1.5            mclust_5.3              dplyr_0.7.0            
## [41] dendextend_1.5.2        RCurl_1.95-4.8          magrittr_1.5            modeltools_0.2-21      
## [45] GenomeInfoDbData_0.99.0 Matrix_1.2-10           Rcpp_0.12.11            munsell_0.4.3          
## [49] ape_4.1                 abind_1.4-5             viridis_0.4.0           stringi_1.1.5          
## [53] whisker_0.3-2           MASS_7.3-47             zlibbioc_1.22.0         flexmix_2.3-14         
## [57] MAST_1.2.1              plyr_1.8.4              grid_3.4.0              rncl_0.8.2             
## [61] lattice_0.20-35         uuid_0.1-2              igraph_1.0.1            taxize_0.8.4           
## [65] fpc_2.1-10              rngtools_1.2.4          reshape2_1.4.2          codetools_0.2-15       
## [69] glue_1.1.0              XML_3.98-1.8            evaluate_0.10           RNeXML_2.0.7           
## [73] data.table_1.10.4       locfdr_1.1-8            tidyr_0.6.3             gtable_0.2.0           
## [77] reshape_0.8.6           assertthat_0.2.0        kernlab_0.9-25          ggplot2_2.2.1          
## [81] gridBase_0.4-7          phylobase_0.8.4         xtable_1.8-2            class_7.3-14           
## [85] viridisLite_0.2.0       tibble_1.3.3            iterators_1.0.8         registry_0.3           
## [89] memoise_1.1.0           cluster_2.0.6           rgl_0.95.1441
```

