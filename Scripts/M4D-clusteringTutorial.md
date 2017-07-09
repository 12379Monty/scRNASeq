### Clustering Tutorial


<style type="text/css"> body, td { font-size: 14px; } code.r{ font-size: 12px; } pre { font-size: 12px } </style>

[Clustering tutorial](https://github.com/epurdom/clusterExperiment/tree/master/vignettes)



<!-- ***************************************************** -->


<!-- ***************************************************** -->



<!-- ***************************************************** -->



<!-- ***************************************************** -->

 
<!--
%\VignetteEngine{knitr::rmarkdown}
%\VignetteIndexEntry{clusterExperiment Vignette}
-->




### Introduction {#Intro}

The goal of this package is to encourage the user to try many different clustering algorithms in one package structure. We give tools for running many different clusterings and choices of parameters. We also provide visualization to compare many different clusterings and algorithm tools to find common shared clustering patterns. We implement common post-processing steps unrelated to the specific clustering algorithm (e.g. subsampling the data for stability, finding cluster-specific markers via differential expression, etc). 

The other main goal of this package is to implement strategies that we have developed for finding a single robust clustering based on the many clusterings that the user might create by perturbing various parameters of a clustering algorithm. There are several steps to these strategies that we call our standard clustering workflow. Our RSEC algorithm (Resampling-based Sequential Ensemble Clustering) is our preferred realization of this workflow that depends on subsampling on and other ensembl methods to provide robust clusterings, particularly for single-cell sequencing experiments and other large mRNA-Seq experiments.

We also provide a class `clusterExperiment` that inherits from `SummarizedExperiment` to store the many clusterings and related information.  

**All of our methods also have a barebones version that allows input of matrices and greater control. This comes at the expense of the user having to manage and keep track of the clusters, input data, transformation of the data, etc**. We do not discuss these barebone versions in this tutorial. Instead, we focus on using `SummarizedExperiment` object as the input and working with the resulting `ClusterExperiment` object. See the help pages of each method for more on how to allow for matrix input.

Although this package was developed with (single-cell) RNA-seq data in mind, its use is not limited to RNA-seq or even to gene expression data. Any dataset characterized by high dimensionality could benefit from the methods implemented here.

## The clustering workflow 

The package encodes many common practices that are shared across clustering algorithms, like subsampling the data, computing silhouette width, sequential clustering procedures, and so forth. It also provides novel strategies that we developed as part of the RSEC algorithm. 

As mentioned above, RSEC is a specific algorithm for creating a clustering that follows these basic steps:

* Implement many different clusterings using different choices of parameters using the function `clusterMany`. This results in a large collection of clusterings, where each clustering is based on different parameters. 
* Find a unifying clustering across these many clusterings using the `combineMany` function. 
* Determine whether some clusters should be merged together into larger clusters. This involves two steps:
    - Find a hierarchical clustering of the clusters found by `combineMany` using `makeDendrogram`
    - Merge together clusters of this hierarchy based on the percentage of differential expression, using `mergeClusters`.

The basic premise of RSEC is to find small, robust clusters of samples, and then merge them into larger clusters as relevant. We find that many algorithmic methods for choosing the appropriate number of clusters for methods  err on the side of too few clusters. However, we find in practice that we tend to prefer to err on finding many clusters and then merging them based on examining the data.

RSEC makes many specific choices in this basic workflow, and many steps of this workflow are useful separately from RSEC. For this reason, the `clusterExperiment` package generalizes this workflow so that the user can use these tools with their own choices. We call this generalization the clustering workflow, as oppose to the specific choices set in RSEC. We will introduce this workflow in its general setting, and only ing specific sections devoted to the `RSEC` function will discussed the RSEC algorithm. For this reason, the examples shown here use simple clustering choices in the workflow to save on computational overhead; RSEC uses extensive subsampling techniques and takes much longer to run than is practical for a vignette.

Users can also run or add their own clusters to this workflow at different stages. Additional functionality for clustering is available in the `clusterSingle` function, and the user should see the documentation in the help page of that function. <!-- The more technical documentation in XXX also describes in more detail additional functionality that is available that allows the user to specify the underlying clustering algorithm and greater control of the underlying parameters. -->  However, there is special functionality to ease in quickly visualizing and managing the results of this workflow.

### The RSEC Routine

The RSEC algorithm (Resampling-based Sequential Ensemble Clustering) was developed specifically for finding robust clusters in single cell sequencing data. It follows our above suggested workflow, but makes specific choices along the way that we find to be advantageous in clustering large mRNA expression datasets, like those of single cell sequencing. It is implemented in the function `RSEC`. This vignette serves to show how the pieces of the workflow operate, and use choices that are shown are not those recommended by RSEC. This is both to show the flexibility of the functions and because RSEC is computationally more time-intensive, since it is based on both resampling, and sequential discovery of clusters. 

### Finding related features/genes

A common practice after determining clusters is to perform differential gene expression analysis in order to find genes that show the greatest differences amongst the clusters. We would stress that this is purely an exploratory technique, and any p-values that result from this analysis are not valid, in the sense that they are likely to be inflated. This is because the same data was used to define the clusters and to perform differential expression analysis.

Since this is a common task, we provide the function `getBestFeatures` to perform various kinds of differential expression analysis between the clusters. A common F-statistic between groups can be chosen. **However, we find that it is far more informative to do pairwise comparisons between clusters, or one cluster against all, in order to find genes that are specific to a particular cluster**. An option for all of these choices is provided in the `getBestFeatures` function. The `getBestFeatures` function uses the DE analysis provided by the `limma` package. In addition, the `getBestFeatures` function provides an option to do use the "voom" correction in the `limma` package to account for the mean-variance relationship that is common in count data. 

### Visualization Tools

We provide a visualization to compare many clusterings of the same data in the function `plotClusters`.  We also provide some functions that are wrappers for common visualization tasks in clustering gene expression data. We provide a heatmap function, `plotHeatmap` that is an interface to the `aheatmap` function in the `NMF` package. 

# Quickstart {#quickstart}

We will go quickly through the standard steps of clustering using the  `clusterExperiment` package, before turning to more details. The standard workflow we envision is the following:

* `clusterMany` -- run desired clusterings
* `combineMany` -- get a unified clustering
* `makeDendrogram` -- get a hierarchical relationship between the clusters
* `mergeClusters` -- merge together clusters with little DE between their genes.
* `getBestFeatures` -- Find Features that are differential between the final clustering.


## Data example {#data}

We will make use of a single cell RNA sequencing experiment made available in the `scRNAseq` package.


```r
 set.seed(14456) ## for reproducibility, just in case
 suppressMessages(require(scRNAseq))
 data("fluidigm")
```

We will use the `fluidigm` dataset (see `help("fluidigm")`). This dataset is stored as a SummarizedExperiment object. We can access the data with `assay` and metadata on the samples with `colData`. 


```r
 assay(fluidigm)[1:5,1:10]
```

```
##          SRR1275356 SRR1274090 SRR1275251 SRR1275287 SRR1275364 SRR1275269 SRR1275263 SRR1275242
## A1BG              0          0          0          0          0          0          0          0
## A1BG-AS1          0          0          0          0          0          0          0          0
## A1CF              0          0          0          0          0          0          0          0
## A2M               0          0          0         31          0         46          0          0
## A2M-AS1           0          0          0          0          0          0          0          0
##          SRR1275338 SRR1274117
## A1BG              0          0
## A1BG-AS1          0          0
## A1CF              0          0
## A2M               0         29
## A2M-AS1           0        133
```

```r
 colData(fluidigm)[,1:5]
```

```
## DataFrame with 130 rows and 5 columns
##               NREADS  NALIGNED    RALIGN TOTAL_DUP    PRIMER
##            <numeric> <numeric> <numeric> <numeric> <numeric>
## SRR1275356  10554900   7555880   71.5862   58.4931 0.0217638
## SRR1274090    196162    182494   93.0323   14.5122 0.0366826
## SRR1275251   8524470   5858130   68.7213   65.0428 0.0351827
## SRR1275287   7229920   5891540   81.4884   49.7609 0.0208685
## SRR1275364   5403640   4482910   82.9609   66.5788 0.0298284
## ...              ...       ...       ...       ...       ...
## SRR1275259   5949930   4181040   70.2705   52.5975 0.0205253
## SRR1275253  10319900   7458710   72.2747   54.9637 0.0205342
## SRR1275285   5300270   4276650   80.6873   41.6394 0.0227383
## SRR1275366   7701320   6373600   82.7600   68.9431 0.0266275
## SRR1275261  13425000   9554960   71.1727   62.0001 0.0200522
```

```r
 NCOL(fluidigm) #number of samples
```

```
## [1] 130
```


### Step 0: Filtering and normalization {#step0}

While there are 130 samples, there are only 65 cells, because each cell is sequenced twice at different sequencing depth. We will limit the analysis to the samples corresponding to high sequencing depth. 


```r
 se <- fluidigm[,colData(fluidigm)[,"Coverage_Type"]=="High"]
```

We also filter out lowly expressed genes: we retain only those genes with at least 10 reads in at least 10 cells.


```r
 wh_zero <- which(rowSums(assay(se))==0)
 pass_filter <- apply(assay(se), 1, function(x) length(x[x >= 10]) >= 10)
 se <- se[pass_filter,]
 dim(se)
```

```
## [1] 7069   65
```

This removed 19186 genes out of 26255. We now have 7069 genes (or features) remaining. Notice that it is important to remove genes with zero counts in all samples (we had 9673 genes which were zero in all samples here). Otherwise, PCA dimensionality reductions and other implementations may have a problem. 

Normalization is an important step in any RNA-seq data analysis and many different normalization methods have been proposed in the literature. Comparing normalization methods or finding the best performing normalization in this dataset is outside of the scope of this vignette. Instead, we will use a simple quantile normalization that will at least make our clustering reflect the biology rather than the difference in sequencing depth among the different samples.


```r
 fq <- round(limma::normalizeQuantiles(assay(se)))
 assays(se) <- list(normalized_counts=fq)
```

### Step 1: Clustering with `clusterMany` {#step1}

`clusterMany` lets the user quickly pick between many clustering options and run all of the clusterings in one single command. In the quick start we pick a simple set of clusterings based on varying the dimensionality reduction options. The way to designate which options to vary is to give multiple values to an argument. Due to a bug in R, we need to set `getClass.msg=FALSE` or otherwise a slew of annoying warnings will spit out at every call; this should be fixed in the next patch to R. 

Here is our call to `clusterMany`:


```r
 suppressMessages(require(clusterExperiment))
 ce<-clusterMany(se, clusterFunction="pam",ks=5:10,
       isCount=TRUE,dimReduce=c("PCA","var"),nVarDims=c(100,500,1000),
       nPCADims=c(5,15,50),run=TRUE)
```

In this call to `clusterMany` we made the follow choices about what to vary: 

* set `dimReduce=c("PCA", "var")` meaning run the clustering algorithm after both a dimensionality reduction to the top principal componetns, and also after filtering to the top most variable genes
* For PCA reduction, choose 5,15, and 50 principal components for the reduced data set (set `nPCADims=c(5,15,50)`)
* For most variable reduction, choose 100, 500, and 1000 most variable genes (set `nVarDims=c(100,500,1000)`)
* For the number of clusters, vary from $k=5,\ldots,10$ (set `ks=5:10`)

By giving only a single value to the relative argument, we keep the other possible options fixed, for example: 

* we used 'pam' for all clustering (`clusterFunction`) 
* we left the default of `minSize` which requires clusters to be of size at least 5; smaller clusters will be ignored and samples assigned to them given the unassigned value of -1.

We also set `isCount=TRUE` to indicate that our input data are counts. This means that operations for clustering and visualizations will internally transform the data as $log_2(x+1)$ (We could have alternatively explicitly set a transformation by giving a function to the  `transFun` argument, for example if we wanted $\sqrt(x)$ or $log(x+\epsilon)$ or just `identity`). 

We can visualize the resulting clusterings using the `plotClusters` command. It is also useful to change the amount of space to allow for the labels of the clusterings, so we will reset the `mar` option in `par` (we also change `axisLine` argument for this reason).


```r
 defaultMar<-par("mar")
 plotCMar<-c(1.1,8.1,4.1,1.1)
 par(mar=plotCMar)
 plotClusters(ce,main="Clusters from clusterMany", 
              whichClusters="workflow", 
              sampleData=c("Biological_Condition","Cluster2"), 
              axisLine=-1)
```

![plot of chunk plotClusterEx1](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M4D/plotClusterEx1-1.png)

This plot shows the samples in the columns, and different clusterings on the rows. Each sample is color coded based on its clustering for that row, where the colors have been chosen to try to match up clusters across different clusterings that show large overlap. Moreover, the samples have been ordered so that each subsequent clustering (starting at the top and going down) will try to order the samples to keep the clusters together, without rearranging the clustering blocks of the previous clustering/row.

Notice that we also added the `sampleData` argument in our call, indicating that we also want to visualize some information about the samples saved in the `colData` slot (inherited from our original `fluidigm` object). We chose the columns "Biological_Condition" and "Cluster2" from `colData`, which correspond to the original biological condition of the experiment, and the clusters reported in the original paper, respectively. These are shown at the bottom of the plot.

Notice that some samples are white. This indicates that they have the value -1, meaning they were not clustered. This is from our choices to require at least 5 samples to make a cluster. 

We can see that some clusters are fairly stable across different choices of dimensions while others can vary dramatically. 

We have set `whichClusters="workflow"` to only plot clusters created from the workflow. Right now that's all there are anyway, but as commands get rerun with different options, other clusterings can build up in the object (see discussion in [this section](#rerun) about how multiple calls to workflow are stored). So setting `whichClusters="workflow"` means that we are only going to see our most recent calls to the workflow.

The labels shown are those given by `clusterMany` but can be a bit much for plotting. We can assign new labels if we prefer, for example to be more succinct, by changing the `clusterLabels` of the object (note we are permanently changing the labels here within the `ClusterExperiment` object). We choose to remove "Features" as being too wordy:


```r
 cl<-clusterLabels(ce)
 cl<-gsub("Features","",cl)
 clusterLabels(ce)<-cl
```

We will also show the clusters in a different order, which corresponds to varying the number of dimensions, rather than k. 


```r
 cl<-clusterLabels(ce)
 ndims<-sapply(sapply(strsplit(cl,","),
        function(x){strsplit(x[1],"=")}),function(x){x[2]})
 ord<-order(as.numeric(ndims))
 par(mar=plotCMar)
 plotClusters(ce,main="Clusters from clusterMany", 
              whichClusters=ord, 
              sampleData=c("Biological_Condition","Cluster2"), 
              axisLine=-1)
```

![plot of chunk plotCluster_newOrder](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M4D/plotCluster_newOrder-1.png)

We see that the order in which the clusters are given to `plotClusters` changes the plot greatly. There are many different options for how to run `plotClusters` discussed in in the detailed section on  [plotClusters](#plotClusters), but for now, this plot is good enough for a quick visualization. 

### The output 

The output of `clusterMany` is a `ClusterExperiment` object. This is a class built for this package and explained in the section on [ClusterExperiment Objects](#ceobjects). In the object `ce` the clusters are stored, names and colors for each cluster within a clustering are assigned, and other information about the clusterings is recorded. Furthermore, all of the information in the original `SummarizedExperiment` is retained. 

There are many accessor functions that help you get at the information in a `ClusterExperiment` object and some of the most relevant are described in the section on [ClusterExperiment Objects](#ceobjects). (`ClusterExperiment` objects are S4 objects, and are not lists). 

For right now we will only mention the most basic such function that retrieves the actual cluster assignments. This is the `clusterMatrix` function, that returns a matrix where the columns are the different clusterings and the rows are samples. Within a clustering, the clusters are encoded by consecutive integers. 


```r
 head(clusterMatrix(ce)[,1:3])
```

```
##      nVAR=100,k=5 nVAR=500,k=5 nVAR=1000,k=5
## [1,]            1            1             1
## [2,]            2            1             1
## [3,]            3            2             2
## [4,]            1            1             1
## [5,]            4            3             2
## [6,]            1            4             2
```

Because we have changed `clusterLabels` above, these new cluster labels are shown here. Notice that some of the samples are assigned the value of -1. -1 has the significance of encoding samples that are not assigned to any cluster. Why certain samples are not clustered depends on the underlying choices of the clustering routine. In this case, the default in `clusterMany` set the minimum size of a cluster to be 5, which resulted in -1 assignments.

Another special value is -2 discussed in the section on [ClusterExperiment objects](#ceobjects)

### Step 2: Find a consensus with `combineMany` {#step2}

To find a consensus clustering across the many different clusterings created by `clusterMany` the function `combineMany` can be used next. 


```r
 ce<-combineMany(ce)
```

```
## Note: no clusters specified to combine, using results from clusterMany
```

Notice we get a warning that we did not specify any clusters to combine, so it is using the default -- those from the `clusterMany` call. 

If we look at the `clusterMatrix` of the returned `ce` object, we see that the new cluster from `combineMany` has been added to the existing clusterings. This is the basic strategy of the functions in this package. Any clustering that is created is added to existing clusterings, so the user does not need to keep track of past clusterings and can easily compare what has changed. 

We can again run `plotClusters`, which will now also show the result of `combineMany`:


```r
 head(clusterMatrix(ce)[,1:3])
```

```
##      combineMany nVAR=100,k=5 nVAR=500,k=5
## [1,]          -1            1            1
## [2,]          -1            2            1
## [3,]          -1            3            2
## [4,]          -1            1            1
## [5,]          -1            4            3
## [6,]          -1            1            4
```

```r
 par(mar=plotCMar)
 plotClusters(ce,whichClusters="workflow")
```

![plot of chunk lookAtCombineMany](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M4D/lookAtCombineMany-1.png)

The default result of `combineMany` is not usually a great choice, and certainly isn't helpful here. The clustering from the default `combineMany` leaves most samples unassigned (white in the above plot). This is because the default way of combining is very conservative -- it requires samples to be in the same cluster in *every clustering* to be assigned a cluster. This is quite stringent. We can vary this by setting the `proportion` argument to indicate the minimum proportion of times they should be together with other samples in the cluster they are assigned to.  Explicit details on how `combineMany` makes these clusters are discussed  in the section on [combineMany](#combineMany).

So let's label the one we found "combineMany, default" and then create a new one. (Making an informative label will make it easier to keep track of this particular clustering later, particularly if we make multiple calls to the workflow). 


```r
 wh<-which(clusterLabels(ce)=="combineMany")
 if(length(wh)!=1) stop() else clusterLabels(ce)[wh]<-"combineMany,default"
```

Now we'll rerun `combineMany` with `proportion=0.7`. This time, we will give it an informative label upfront in our call to `combineMany`.


```r
 ce<-combineMany(ce,proportion=0.7,clusterLabel="combineMany,0.7")
```

```
## Note: no clusters specified to combine, using results from clusterMany
```

```r
 par(mar=plotCMar)
 plotClusters(ce,whichClusters="workflow")
```

![plot of chunk combineMany_newCombineMany](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M4D/combineMany_newCombineMany-1.png)

We see that more clusters are detected. Those that are still not assigned a cluster from `combineMany` clearly vary across the clusterings as to whether the samples are clustered together or not. Varying the `proportion` argument will adjust whether some of the unclustered samples get added to a cluster. There is also a `minSize` parameter for `combineMany`, with the default of `minSize=5`. We could reduce that requirement as well and more of the unclustered samples would be grouped into a cluster. Here, we reduce it to `minSize=3` (we'll call this "combineMany,final"):


```r
 ce<-combineMany(ce,proportion=0.7,minSize=3,clusterLabel="combineMany,final")
```

```
## Note: no clusters specified to combine, using results from clusterMany
```

```r
 par(mar=plotCMar)
 plotClusters(ce,whichClusters="workflow",main="Min. Size=3")
```

![plot of chunk combineMany_changeMinSize](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M4D/combineMany_changeMinSize-1.png)


We can also visualize the proportion of times these clusters were together across these clusterings (this information was made and stored in the ClusterExperiment object when we called `combineMany` as long as `proportion` value is <1):


```r
 plotCoClustering(ce)
```

![plot of chunk plotCoClustering_quickstart](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M4D/plotCoClustering_quickstart-1.png)

This visualization can help in determining whether to change the value of `proportion` (though see [combineMany](#combineMany) for how -1 assignments affect `combineMany`).

### Step 3: Merge clusters together with `makeDendrogram` and `mergeClusters` {#step3}

It is not uncommon in practice to create forty or more clusterings with `clusterMany`, in which case the results of `combineMany` can often still result in too many small clusters. We might wonder if they are necessary or could be logically combined together. We could change the value of `proportion` in our call to `combineMany`. But we have found that it is often after looking at the clusters and how different they look on individual genes that we best make this determination, rather than the proportion of times they are together in different clustering routines. 

For this reason, we often find the need for an additional clustering step that merges clusters together that are not different, based on running tests of differential expression between the clusters found in `combineMany`. We often display and use both sets of clusters side-by-side (that from `combineMany` and that from `mergeClusters`).

`mergeClusters` needs a hierarchical clustering of the clusters; it then goes progressively up that hierarchy, deciding whether two adjacent clusters can be merged. The function `makeDendrogram` makes such a hierarchy between clusters (by applying `hclust` to the medoids of the clusters). Because the results of `mergeClusters` are so dependent on that hierarchy, we require the user to call `makeDendrogram` rather than calling it internally. This is because different options in `makeDendrogram` can affect how the clusters are hierarchically ordered, and we want to encourage the user make these choices.

As an example, here we use the 500 most variable genes to make the cluster hierarchy.


```r
 ce<-makeDendrogram(ce,dimReduce="var",ndims=500)
 plotDendrogram(ce)
```

![plot of chunk makeDendrogram](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M4D/makeDendrogram-1.png)

We can see that clusters 1 and 3 are most closely related, at least in the top 500 most variable genes. 

If we look at the summary of `ce`, it now has 'makeDendrogram' marked as 'Yes'. 


```r
 ce
```

```
## class: ClusterExperiment 
## dim: 7069 65 
## Primary cluster type: combineMany 
## Primary cluster label: combineMany,final 
## Table of clusters (of primary clustering):
##  c1 c-1  c2  c3  c4  c5  c6 
##   9  12   7  15  14   4   4 
## Total number of clusterings: 39 
## Dendrogram run on 'combineMany,final' (cluster index: 1)
## -----------
## Workflow progress:
## clusterMany run? Yes 
## combineMany run? Yes 
## makeDendrogram run? Yes 
## mergeClusters run? No
```


Now we are ready to actually merge clusters together. We now run `mergeClusters` that will go up this hierarchy and compare the level of differential expression (DE) in each pair. In other words, if we focus on the left side of the tree, DE tests are run, between 1 and 3, and between 6 and 8. If there is not enough DE between each of these (based on a cutoff that can be set by the user), then clusters 1 and 3 and/or 6 and 8 will be merged. And so on up the tree.

It is useful to first run `mergeClusters` without actually creating any merged clusters so as to preview what the final clustering will be (and perhaps to help in setting the cutoff).


```r
 mergeClusters(ce,mergeMethod="adjP",plotInfo="mergeMethod")
```

```
## Note: Merging will be done on ' combineMany,final ', with clustering index 1
```

![plot of chunk mergeClustersPlot](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M4D/mergeClustersPlot-1.png)

Then we can decide on a cutoff and visualize the resulting clustering. 

<!---We'll add some more information to the `plotCoClustering` making the clustering of the samples be via the dendrogram on which `mergeClusters` was made:
-->


```r
 ce<-mergeClusters(ce,mergeMethod="adjP",cutoff=0.02)
```

```
## Note: Merging will be done on ' combineMany,final ', with clustering index 1
```

![plot of chunk mergeClusters](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M4D/mergeClusters-1.png)

```r
 par(mar=plotCMar)
 plotClusters(ce,whichClusters="workflow", 
              sampleData=c("Biological_Condition","Cluster2"))
 plotCoClustering(ce,whichClusters=c("mergeClusters","combineMany"),
                 sampleData=c("Biological_Condition","Cluster2"),annLegend=FALSE)
```

![plot of chunk mergeClusters](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M4D/mergeClusters-2.png)![plot of chunk mergeClusters](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M4D/mergeClusters-3.png)

Notice that `mergeClusters` combines clusters based on the actual values of the features, while the `coClustering` plot shows how often the samples clustered together. It is not uncommon that `mergeClusters` will merge clusters that don't look "close" on the `coClustering` plot. This can be due to just the choices of the hierarchical clustering, but can also be because the two merged clusters are not often confused for each other across the clustering algorithms, yet don't have strong differences on individual genes. This can be the case especially when the clustering is done on reduced PCA space, where an accumulation of small differences might consistently separate the samples (so the two clusters are not "confused" as to the samples), but because the differences are not strong on individual genes, `mergeClusters` combines them. These are ultimately different criteria.

Finally, we can do a heatmap visualizing this final step of clustering.


```r
 plotHeatmap(ce,clusterSamplesData="dendrogramValue",breaks=.99,
            sampleData=c("Biological_Condition", "Cluster1", "Cluster2"))
```

![plot of chunk plotHeatmap](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M4D/plotHeatmap-1.png)

By choosing "dendrogramValue" for the clustering of the samples, we will be showing the clusters according to the hierarchical ordering of the clusters found by `makeDendrogram`. The argument `breaks=0.99` means that the last color of the heatmap spectrum will be forced to be the top 1% of the data (rather than evenly spaced through the entire range of values). This can be helpful in making sure that rare extreme values in the upper range do not absorb too much space in the color spectrum. There are many more options for `plotHeatmap`, some of which are discussed in the section on [plotHeatmap](#plotHeatmap).

### Step 4: Finding Features related to the clusters {#step4}

The last step is to then find features that are different between these clusters, as a start to understanding biological differences in the samples. The function `getBestFeatures` performs tests for differential expression (i.e. different mean levels) between the clusters for each feature. It relies on `limma` to run the differential expression analysis, with `voom` correction if the data are indicated by the user to be counts.

There are several types of tests that can be performed to identify features that are different between the clusters. Here we perform all pairwise tests between the clusters. 


```r
 pairsAll<-getBestFeatures(ce,contrastType="Pairs",p.value=0.05,
                           number=nrow(ce),isCount=TRUE)
```

```
## Note: If `isCount=TRUE` the data will be transformed with voom() rather than
## with the transformation function in the slot `transformation`.
## This makes sense only for counts.
```

```r
 head(pairsAll)
```

```
##   IndexInOriginal Contrast Feature     logFC    AveExpr         t      P.Value    adj.P.Val
## 1            2210    X1-X2    GLI3 -8.794644  3.3201396 -8.479338 9.915358e-13 3.285547e-10
## 2             552    X1-X2  BCL11B  9.201065  2.2054809  8.155363 4.262018e-12 1.205128e-09
## 3            1483    X1-X2  DNAJB1 -9.369208  1.2717802 -7.954335 1.051694e-11 2.671051e-09
## 4            4501    X1-X2 PRTFDC1 -8.573379 -0.1792585 -7.206006 2.969915e-10 4.998650e-08
## 5            4590    X1-X2   PTPRD  8.377537  4.3540717  6.972703 8.326707e-10 1.197183e-07
## 6            1866    X1-X2 FAM107A -8.753312  0.6885497 -6.973648 8.292109e-10 1.196264e-07
##          B
## 1 17.15910
## 2 16.16650
## 3 15.27657
## 4 12.03648
## 5 11.48524
## 6 11.34457
```

We can visualize only these significantly different pair-wise features with `plotHeatmap` by using the column "IndexInOriginal" in the result of `getBestFeatures` to quickly identify the genes to be used in the heatmap. Notice that the same genes can be replicated across different contrasts, so we will not always have unique genes:


```r
 length(pairsAll$Feature)==length(unique(pairsAll$Feature))
```

```
## [1] FALSE
```

In this case they are not unique. Hence, we will make sure we take only unique gene values so that they are not plotted multiple times in our heatmap. (This is a good practice even if in a particular case the genes are unique).


```r
 plotHeatmap(ce, clusterSamplesData="dendrogramValue",
             clusterFeaturesData=unique(pairsAll[,"IndexInOriginal"]),
             main="Heatmap of features w/ significant pairwise differences",
             breaks=.99)
```

![plot of chunk getBestFeatures_heatmap](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M4D/getBestFeatures_heatmap-1.png)

Notice that the samples clustered into the -1 cluster (i.e. not assigned) are clustered as an outgroup. They can also be mixed into the dendrogram (see [makeDendrogram](#makeDendrogram))

# ClusterExperiment Objects {#ceobjects}

The `ce` object that we created by calling `clusterMany` is a `ClusterExperiment` object. The `ClusterExperiment` class is used by this package to keep the data and the clusterings together. It inherits from `SummarizedExperiment`, which means the data and `colData` and other information orginally in the `fluidigm` object are retained and can be accessed with the same functions as before. The `ClusterExperiment` object additionally stores clusterings and information about the clusterings along side the data. This helps keep everything together, and like the `SummarizedExperiment` object, allows for simple things like subsetting to a reduced set of subjects and being confident that the corresponding clusterings, colData, and so forth are similarly subset.

Typing the name at the control prompt results in a quick summary of the object. 


```r
 ce
```

```
## class: ClusterExperiment 
## dim: 7069 65 
## Primary cluster type: mergeClusters 
## Primary cluster label: mergeClusters 
## Table of clusters (of primary clustering):
##  m1 m-1  m2  m3  m4 
##  17  12   7  15  14 
## Total number of clusterings: 40 
## Dendrogram run on 'combineMany,final' (cluster index: 2)
## -----------
## Workflow progress:
## clusterMany run? Yes 
## combineMany run? Yes 
## makeDendrogram run? Yes 
## mergeClusters run? Yes
```

This summary tells us the total number of clusterings (40), and gives some indication as to what parts of the standard workflow have been completed and stored in this object. It also gives information regarding the `primaryCluster` of the object. The `primaryCluster` is just one of the clusterings that has been chosen to be the "primary" clustering, meaning that by default various functions will turn to this clustering as the desired clustering to use. The "primaryCluster" can be reset by the user (see `primaryClusterIndex`). `clusterMany` arbitrarily sets the 'primaryCluster' to the first one, and each later step of the workflow sets the primary index to the most recent, but the user can set a specific clustering to be the primaryCluster with `primaryClusterIndex`.

There are also additional commands to access the clusterings and their related information (type `help("ClusterExperiment-methods")` for more).

The cluster assignments are stored in the `clusterMatrix` slot of `ce`, with samples on the rows and different clusterings on the columns. We can look at the cluster matrix and the primary cluster with the commands `clusterMatrix` and `primaryCluster`


```r
 head(clusterMatrix(ce))[,1:5]
```

```
##      mergeClusters combineMany,final combineMany,0.7 combineMany,default nVAR=100,k=5
## [1,]            -1                -1              -1                  -1            1
## [2,]            -1                -1              -1                  -1            2
## [3,]            -1                -1              -1                  -1            3
## [4,]             1                 1               1                  -1            1
## [5,]            -1                -1              -1                  -1            4
## [6,]            -1                -1              -1                  -1            1
```

```r
 primaryCluster(ce)
```

```
##  [1] -1 -1 -1  1 -1 -1  2  3  3  3  4  4  1  4 -1  1  2  3  2 -1  4  4  4  3 -1  3  3  3  1  4  4  1
## [33]  1  2  2  1 -1  4 -1  1  2  3  3  1  3  3  1 -1  2  4  4  4  3  4 -1  1  1  3  3  1  1  1  4  1
## [65]  1
```

Remember that we made multiple calls to `combineMany`: only the last such call will be shown when we use `whichClusters="workflow"` in our plotting (see this [section](#rerun) for a discussion of how these repeated calls are handled.)

**Negative Valued Cluster Assignments** The different clusters are stored as consecutive integers, with '-1' and '-2' having special meaning. '-1' refers to samples that were not clustered by the clustering algorithm. In our example, we removed clusters that didn't meet specific size criterion, so they were assigned '-1'. '-2' is for samples that were not included in the original input to the clustering. This is useful if, for example, you cluster on a subset of the samples, and then want to store this clustering with the clusterings done on all the data. You can create a vector of clusterings that give '-2' to the samples not originally used and then add these clusterings to the `ce` object manually with `addClusters`. 

`clusterLabels` gives the column names of the `clusterMatrix`; `clusterMany` has given column names based on the parameter choices, and later steps in the workflow also give a name (or allow the user to set them). Clusterings might also have no specific label if the user created them. As we've seen, the user can also change these labels.


```r
 head(clusterLabels(ce),10)
```

```
##  [1] "mergeClusters"       "combineMany,final"   "combineMany,0.7"     "combineMany,default"
##  [5] "nVAR=100,k=5"        "nVAR=500,k=5"        "nVAR=1000,k=5"       "nPCA=5,k=5"         
##  [9] "nPCA=15,k=5"         "nPCA=50,k=5"
```

`clusterTypes` on the other hand indicates what call made the clustering. Unlike the labels, it is wise to not change the values of `clusterTypes` lightly. 


```r
 head(clusterTypes(ce),10)
```

```
##  [1] "mergeClusters" "combineMany"   "combineMany.2" "combineMany.1" "clusterMany"   "clusterMany"  
##  [7] "clusterMany"   "clusterMany"   "clusterMany"   "clusterMany"
```

The information that was in the original `fluidigm` object has also been preserved, like `colData` that contains information on each sample.


```r
 colData(ce)[,1:5]
```

```
## DataFrame with 65 rows and 5 columns
##               NREADS  NALIGNED    RALIGN TOTAL_DUP    PRIMER
##            <numeric> <numeric> <numeric> <numeric> <numeric>
## SRR1275356  10554900   7555880   71.5862   58.4931 0.0217638
## SRR1275251   8524470   5858130   68.7213   65.0428 0.0351827
## SRR1275287   7229920   5891540   81.4884   49.7609 0.0208685
## SRR1275364   5403640   4482910   82.9609   66.5788 0.0298284
## SRR1275269  10729700   7806230   72.7536   50.4285 0.0204349
## ...              ...       ...       ...       ...       ...
## SRR1275259   5949930   4181040   70.2705   52.5975 0.0205253
## SRR1275253  10319900   7458710   72.2747   54.9637 0.0205342
## SRR1275285   5300270   4276650   80.6873   41.6394 0.0227383
## SRR1275366   7701320   6373600   82.7600   68.9431 0.0266275
## SRR1275261  13425000   9554960   71.1727   62.0001 0.0200522
```

Another important slot in the `ClusterExperiment` object is the `clusterLegend` slot. This consists of a list, one element per column or clustering of `clusterMatrix`. 


```r
 length(clusterLegend(ce))
```

```
## [1] 40
```

```r
 clusterLegend(ce)[1:2]
```

```
## $mergeClusters
##    clusterIds color     name 
## -1 "-1"       "white"   "m-1"
## 1  "1"        "#1F78B4" "m1" 
## 2  "2"        "#33A02C" "m2" 
## 3  "3"        "#E31A1C" "m3" 
## 4  "4"        "#FF7F00" "m4" 
## 
## $`combineMany,final`
##    clusterIds color     name 
## -1 "-1"       "white"   "c-1"
## 1  "1"        "#1F78B4" "c1" 
## 2  "2"        "#33A02C" "c2" 
## 3  "3"        "#E31A1C" "c3" 
## 4  "4"        "#FF7F00" "c4" 
## 5  "5"        "#6A3D9A" "c5" 
## 6  "6"        "#B15928" "c6"
```

We can see that each element of `clusterLegend` consists of a matrix, with number of rows equal to the number of clusters in the clustering. The columns store information about that cluster. `clusterIds` is the internal id used in `clusterMatrix` to identify the cluster, `name` is a name for the cluster, and `color` is a color for that cluster. `color` is used in plotting and visualizing the clusters, and `name` is an arbitrary character string for a cluster. They are automatically given default values when the `ClusterExperiment` object is created, but we will see under the description of visualization methods how the user might want to manipulate these for better plotting results. 


# Visualizing the data {#visual}

## Plotting the clusters {#plotClusters}
We demonstrated during the quick start that we can visualize multiple clusterings using the `plotClusters` command.


```r
 par(mar=plotCMar)
 plotClusters(ce,main="Clusters from clusterMany", whichClusters="workflow", 
              axisLine=-1)
```

![plot of chunk plotClusterEx1_redo](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M4D/plotClusterEx1_redo-1.png)



We have seen that we can get very different plots depending on how we order the clusterings, and what clusterings are included. The argument `whichClusters` allows the user to choose different clusterings or provide an explicit ordering of the clusterings. `whichClusters` can take either a single character value, or a vector of either characters  or indices. If `whichClusters` matches either "all" or "workflow", then the clusterings chosen are either all, or only those from the most recent calls to the workflow functions. Choosing "workflow" removes from the visualization both user-defined clusterings and also previous calls to the workflow that have since been rerun. Setting `whichClusters="workflow"` can be a useful if you have called a method like `combineMany` several times, as we did, only with different parameters. All of those runs are saved (unless `eraseOld=TRUE`), but you may not want to plot them.

If `whichClusters` is a character that is not one of these designated values, the entries should match a clusterType value (like `clusterMany`) or a clusterLabel value (with exact matching). Alternatively, the user can specify numeric indices corresponding to the columns of `clusterMatrix` that provide the order of the clusters.


```r
 par(mar=plotCMar)
 plotClusters(ce,whichClusters="clusterMany",
                main="Only Clusters from clusterMany",axisLine=-1)
```

![plot of chunk plotClusterEx1_newOrder](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M4D/plotClusterEx1_newOrder-1.png)

We can also add to our plot (categorical) information on each of our subjects from the `colData` of our SummarizedExperiment object (which is also retained in our ClusterExperiment object). This can be helpful to see if the clusters correspond to other features of the samples, such as sample batches. Here we add the values from the columns "Biological_Condition" and "Cluster2" that were present in the `fluidigm` object and given with the published data.


```r
 par(mar=plotCMar)
 plotClusters(ce,whichClusters="workflow", 
              sampleData=c("Biological_Condition","Cluster2"), 
              main="Workflow clusters plus other data",axisLine=-1)
```

![plot of chunk plotClusterEx1_addData](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M4D/plotClusterEx1_addData-1.png)

### Saving the alignment of plotClusters {#plotClustersAlign}

`plotClusters` invisibly returns a `ClusterExperiment` object. In our earlier calls to `plotCluster`, this would be the same as the input object and so there is no reason to save it. However, the alignment and color assignments created by `plotClusters` can be requested to be saved via the `resetNames`, `resetColors` and `resetOrderSamples` arguments. If any of these are set to TRUE, then the object returned will be different than those of the input. Specifically, if `resetColors=TRUE` the `colorLegend` of the returned object will be changed so that the colors assigned to each cluster will be as were shown in the plot. Similarly, if `resetNames=TRUE` the names of the clusters will be changed to be integer values, but now the integers will be aligned to try to be the same across clusters (and therefore not consecutive integers, which is why these are saved as names for the clusters and not `clusterIds`). If `resetOrderSamples=TRUE`, then the order of the samples shown in the plot will be similarly saved in the slot `orderSamples`. 

Here we make a call to plotClusters, but now ask to reset everything to match this ordering. We'll save this as a different object so that this is not a permanent change for the rest of the vignette.


```r
 par(mar=plotCMar)
 ce_temp<-plotClusters(ce,whichClusters="workflow", 
                sampleData=c("Biological_Condition","Cluster2"), 
                main="Clusters from clusterMany, 
                different order",axisLine=-1,resetNames=TRUE,
                resetColors=TRUE,resetOrderSamples=TRUE)
```

![plot of chunk plotClusterEx1_setColors](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M4D/plotClusterEx1_setColors-1.png)

```r
 clusterLegend(ce_temp)[c("mergeClusters","combineMany,final")]
```

```
## $mergeClusters
##    clusterIds color     name
## -1 "-1"       "white"   "-1"
## 1  "1"        "#1F78B4" "1" 
## 2  "2"        "#E31A1C" "3" 
## 3  "3"        "#FF7F00" "4" 
## 4  "4"        "#6A3D9A" "5" 
## 
## $`combineMany,final`
##    clusterIds color     name
## -1 "-1"       "white"   "-1"
## 1  "1"        "#1F78B4" "1" 
## 2  "2"        "#E31A1C" "3" 
## 3  "3"        "#FF7F00" "4" 
## 4  "4"        "#6A3D9A" "5" 
## 5  "5"        "#B15928" "6" 
## 6  "6"        "#2ef4ca" "7"
```

Now, the `clusterLegend` slot of the object no longer has the default color/name assignments, but it has names and colors that match across the clusters. Notice, that this means the prefix "m" or "c" that was previously given to distinguish the `combineMany` result from the `mergeClusters` result is now gone (the user could manually add them by changing the clusterLegend values).

We can also force `plotClusters` to use the existing color definitions, rather than create its own. This makes particular sense if you want to have continuity between plots -- i.e. be sure that a particular cluster always has a certain color -- but would like to do different variations of plotClusters to get a sense of how similar the clusters are.

For example, we set the colors above based on the cluster order from `plotClusters` where the clusters were ordered according to the workflow. But now we want to plot only the clusters from `clusterMany`, yet keep the same colors as before so we can compare them. We do this by setting the argument `existingColors="all"`, meaning use all of the existing colors (currently this is the only option available for how to use the existing colors).


```r
 par(mar=plotCMar)
 par(mfrow=c(1,2))
 plotClusters(ce_temp, sampleData=c("Biological_Condition","Cluster2"),
              whichClusters="workflow", existingColors="all",
              main="Workflow Clusters, fix the color of clusters",axisLine=-1)
 
 plotClusters(ce_temp, sampleData=c("Biological_Condition","Cluster2"), 
              existingColors="all", whichClusters="clusterMany",
              main="clusterMany Clusters, fix the color of clusters",
              axisLine=-1)
```

![plot of chunk plotClusterEx1_forceColors](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M4D/plotClusterEx1_forceColors-1.png)

## Heatmap with the clusters {#plotHeatmap}
There is also a default heatmap command for a ClusterExperiment object that we used in the Quick Start. By default it clusters on the most variable features (after transforming the data) and shows the primaryCluster alongside the data. The primaryCluster now that we've run the workflow will be set as that from the mergeClusters step. 


```r
 par(mfrow=c(1,1))
 par(mar=defaultMar)
 plotHeatmap(ce,main="Heatmap with clusterMany")
```

![plot of chunk plotHeatmap_Ex1](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M4D/plotHeatmap_Ex1-1.png)

The `plotHeatmap` command has numerous options, in addition to those of `aheatmap`. `plotHeatmap` mainly provides additional functionality in the following areas: 

* Easy inclusion of clustering information or sample information, based on the ClusterExperiment object.
* Additional methods for ordering/clustering the samples that makes use of the clustering information.
* Use of separate input data for clustering and for visualization.

### Displaying clustering or sample information

Like `plotClusters`, `plotHeatmap` has a `whichClusters` option that behaves similarly to that of `plotClusters`.  In addition to the options "all" and "workflow" that we saw with `plotClusters`, `plotHeatmap` also takes the option "none"" (no clusters shown) and "primary" (only the primaryCluster). The user can also request a subset of the clusters by giving specific indices to `whichClusters` like in `plotClusters`.

Here we create a heatmap that shows the clusters from the workflow. Notice that we choose only the last 2 -- from `combineMany` and `mergeClusters`. If we chose all "workflow" clusters it would be too many.


```r
 whClusterPlot<-1:2
 plotHeatmap(ce,whichClusters=whClusterPlot, annLegend=FALSE)
```

![plot of chunk plotHeatmap_Ex1.1](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M4D/plotHeatmap_Ex1.1-1.png)

Notice we also passed the option 'annLegend=FALSE' to the underlying `aheatmap` command (with many clusterings shown, it is often not useful to have a legend for all the clusters because the legend doesn't fit on the page!). The many detailed commands of `aheatmap` that are not set internally by `plotHeatmap` can be passed along as well. 

Like `plotClusters`, `plotHeatmap` takes an argument `sampleData`, which refers to columns of the `colData` of that object and can be included.  

### Additional options for clustering/ordering samples

We can choose to not cluster the samples, but order the samples by cluster. This time we'll just show the primary cluster (the `mergeCluster` result) by setting `whichClusters="primaryCluster"`:


```r
 plotHeatmap(ce,clusterSamplesData="primaryCluster",
             whichClusters="primaryCluster",
             main="Heatmap with clusterMany",annLegend=FALSE)
```

![plot of chunk plotHeatmap_primaryCluster](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M4D/plotHeatmap_primaryCluster-1.png)

As an improvement upon this, we can cluster the clusters into a dendrogram so that the most similar clusters will be near each other. We already did this before with our call to `makeDendrogram`. We haven't done anything to change that, so the dendrogram from that call is still stored in the object. We can check this in the information shown in our object:


```r
 show(ce)
```

```
## class: ClusterExperiment 
## dim: 7069 65 
## Primary cluster type: mergeClusters 
## Primary cluster label: mergeClusters 
## Table of clusters (of primary clustering):
##  m1 m-1  m2  m3  m4 
##  17  12   7  15  14 
## Total number of clusterings: 40 
## Dendrogram run on 'combineMany,final' (cluster index: 2)
## -----------
## Workflow progress:
## clusterMany run? Yes 
## combineMany run? Yes 
## makeDendrogram run? Yes 
## mergeClusters run? Yes
```


We can see, as we expected, that the dendrogram we made from "combineMany,final" is still the active dendrogram saved in the `ce` object. Now we will call `plotHeatmap` choosing to display the "mergeClusters" and "combineMany,final" clustering, and ordering the samples by the dendrogram in `ce`. We will also display some data about the samples from the `colData` slot using the `sampleData` argument. Notice that instead of giving the index, we can also give the clusterLabels of the clusters we want to show.


```r
 plotHeatmap(ce,clusterSamplesData="dendrogramValue",
             whichClusters=c("mergeClusters","combineMany"),
             main="Heatmap with clusterMany",
             sampleData=c("Biological_Condition","Cluster2"),annLegend=FALSE)
```

![plot of chunk plotHeatmap_dendro](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M4D/plotHeatmap_dendro-1.png)

If there is not a dendrogram stored, `plotHeatmap` will call `makeDendrogram` based on the primary cluster (with the default settings of `makeDendrogram`); calling `makeDendrogram` on `ce` is preferred so that the user can control the choices in how it is done (which we will discuss below). For visualization purposes, the dendrogram for the `combineMany` cluster is preferred to that of the `mergeCluster` cluster, since "combineMany,final" is just a finer partition of the "mergeClusters" clustering. 

### Using separate input data for clustering and for visualization

While count data is a common type of data, it is also common that the input data in the SummarizedExperiment object might be normalized data from a normalization package such as `RUVSeq`. In this case, the clustering and all numerical calculations should be done on the normalized data (which may or may not need a log transform). However, these normalized data might not be on a logical count scale (for example, in `RUVSeq`, the normalize data are residuals after subtracting out gene-specific batch effects). 

In this case, it can be convenient to have the *visualization* of the data (i.e. the color scale), be based on a count scale that is interpretable, even while the clustering is done based on the normalized data. This is possible by giving a new matrix of values to the argument `visualizeData`. In this case, the color scale (and clustering of the features) is based on the input `visualizeData` matrix, but all clustering of the samples is done on the internal data in the `ClusterExperiment` object. 


# The clustering workflow {#workflow}

We will now go into more detail about important options for the main parts of the clustering workflow.

## clusterMany {#clusterMany}

### Overview of the implemented clustering procedures
In the quick start section we picked some simple and familiar clustering options that would run quickly and needed little explanation. However, our workflow generally assumes more complex options and more parameter variations are tried. Before getting into the specific options of `clusterMany`, let us first describe some of these more complicated setups, since many of the arguments of `clusterMany` depend on understanding them. 

**Clustering Algorithms (`clusterD`):** Clustering algorithms generally start with a particular predefined distance or dissimilarity between the samples, usually encoded in a $n x n$ matrix, $D$. In our package, we consider only such clustering algorithms. The input could also be similarities, though we will continue to call such a matrix $D$. 

The simplest scenario is a simple calculation of $D$ and then a clustering of $D$, usually dependent on particular parameters.  In this package we try to group together algorithms that cluster $D$ based on common parameters and operations that can be done. Currently there are two "types" of algorithms we consider, which we call type "K" and "01". The help page of `clusterD` documents these choices more fully, but we give an overview here. Many of the parameters that are allowed to vary in `clusterMany` refer to parameters for the clustering of the $D$ matrix, either for "K" or "01" type algorithms.

 The "K" algorithms are so called because their main parameter requirement is that the user specifies the number of clusters ($K$) to be created. They assume the input $D$ are dissimilarities, and depending on the algorithm, may have additional expectations. "pam" and "kmeans" are examples of such types of algorithms. 

The "01" algorithms are so named because the algorithm assumes that the input $D$ consists of *similarities* between samples and that the similarities encoded in $D$ are on a scale of 0-1. They use this fact to make the primary user-specified parameter be not the number of final clusters, but a measure $\alpha$ of how dissimilar samples in the same cluster can be (on a scale of 0-1). Given $\alpha$, the algorithm then implements a method to then determine the clusters (so $\alpha$ implicitly determines $K$). These methods rely on the assumption that because the 0-1 scale has special significance, the user will be able to make an determination more easily as to the level of dissimilarity allowed in a true cluster, rather than predetermine the number of clusters $K$. The current 01 methods are "tight", "hierarchical01" and "pam". 



**Subsampling** In addition to the basic clustering algorithms on a matrix $D$, we also implement many other common cluster processing steps that are relevant to the result of the clustering. We have already seen such an example with dimensionality reduction, where the input $D$ is determined based on different input data. A more significant processing that we perform is calculation of a dissimilarity matrix $D$ not by a distance on the vector of data, but by subsampling and clustering of the subsampled data. The resulting $D$ matrix is a matrix of co-clustering percentages. Each entry is a percentage of subsamples where the two samples shared a clustering over the many subsamplings of the data (there are slight variations as how this  can be calculated, see help pages of `subsampleClustering` ). Note that this implies there actually two different clustering algorithms (and sets of corresponding parameters) -- one for the clustering on the subsampled data, and one for the clustering of the resulting $D$ of the percentage of coClustering of samples. `clusterMany` focuses on varying parameters related to the clustering of $D$, and generally assumes that the underlying clustering on the subsampled data is simple (e.g. "pam"). (The underlying clustering machinery in our package is performed by a function `clusterSingle` which is not discussed in this tutorial, but can be called explicitly for fine-grained control over all of the features ). The subsampling option is computationally expensive, and when coupled with comparing many parameters, does result in a lengthy evaluation of `clusterMany`. However, we recommend it as one of the most useful methods for getting stable clustering results.

**Sequential Detection of Clusters** Another complicated addition to the clustering that requires additional explanation is the implementation of sequential clustering. This refers to clustering of the data, then removing the "best" cluster, and then re-clustering the remaining samples, and then continuing this iteration until all samples are clustered (or the algorithm in some other way calls a stop). Such sequential clustering can often be convenient when there is very dominant cluster, for example, that is far away from the other mass of data. Removing samples in these clusters and resampling can sometimes be more productive and result in results more robust to the choice of samples. A particular implementation of such a sequential method, based upon [@tseng2005], is implemented in the `clusterExperiment` package. Because of the iterative nature of this process, there are many possible parameters (see `help(seqCluster)`). `clusterMany` does not allow variation of very many of these parameters, but instead just has the choice of running the sequential clustering or not. Sequential clustering can also be quite computationally expensive, particularly when paired with subsampling to determine $D$ at each step of the iteration.

### Arguments of `clusterMany`
Now that we've explained the underlying architecture of the clustering provided in the package, we discuss the parameters that can be varied in `clusterMany`. There are additional arguments available for `clusterMany` but right now we focus on only the ones that can be given multiple options. Recall that arguments in `clusterMany` that take on multiple values mean that the combinations of all the multiple valued arguments will be given as input for a clustering routine. 

* `sequential` This parameter consists of logical values, TRUE and/or FALSE, indicating whether the sequential strategy should be implemented or not. 
* `subsample` This parameter consists of logical values, TRUE and/or FALSE, indicating whether the subsampling strategy for determining $D$ should be implemented or not. 
* `clusterFunction` The clustering functions to be tried. If `subsample=TRUE` is part of the combination, then `clusterFunction` the method that will be used on the co-clustering matrix $D$ created from subsampling the data (where pam clustering is used on the subsampled data). Otherwise, `clusterFunction` is the clustering method that will be used on  `dist(t(x))`
* `ks` The argument 'ks' is interpreted differently for different choices of the other parameters. If `sequential=TRUE` is part of the combination, `ks` defines the argument `k0` of sequential clustering (see `help(seqCluster)`), which is approximately like the initial starting point for the number of clusters in the sequential process.  Otherwise, `ks` is passed to set $K$ of both the clustering of subsampled data and the actual clustering of the data (if a $K$ needs to be set, i.e. a type "K" algorithm). When/if `findBestK=TRUE` is part of the combination, `ks` also defines the range of values to search for the best k (see the details in `help(clusterMany)` for more).
* `dimReduce` These are character strings indicating what choices of dimensionality reduction should be tried. The choices are "PCA", indicating clustering on the top principal components, "var", indicating clustering on the top most variable features, and "none", indicating the whole data set should be used. If either "PCA" or "var" are chosen, the following parameters indicate the number of such features to be used (and can be a vector of values as we have seen):
    * `nVarDims`
    * `nPCADims`
* `distFunction` These are character values giving functions that provide a distance matrix between the samples, when applied to the data. These functions should be accessible in the global environment (`clusterMany` applies `get` to the global environment to access these functions). To make them compatible with the `dist` function, these functions should assume the samples are in the rows, i.e. they should work when applied to t(assay(ce)). We give an example in the next subsection below.
* `minSizes` these are integer values determining the minimum size required for a cluster (passed to the `clusterD` part of clustering). 
* `alphas` These are the $\alpha$ parameters for "01" clustering techniques; these values are only relevant if one of the `clusterFunction` values is a "01" clustering algorithm. The values given to `alphas` should be between 0 and 1, with smaller values indicating greater similarity required between the clusters. 
* `betas` These are the $\beta$ parameters for sequential clustering; these values are only relevant if `sequential=TRUE` and determine the level of stability required between changes in the parameters to determine that a cluster is stable.
* `findBestK` This option is for "K" clustering techniques, and indicates that $K$ should be chosen automatically as the $K$ that gives the largest silhouette distance between clusters.
* `removeSil` A logical value as to whether samples with small silhouette distance to their assigned cluster are "removed", in the sense that they are not given their original cluster assignment but instead assigned -1. This option is for "K" clustering techniques as a method of removing poorly clustered samples (the "01" techniques used by `clusterMany` generally do this intrinsically as part of the algorithm). 
  * `silCutoff` If `removeSil` is TRUE, then `silCutoff` determines the cutoff on silhouette distance for unassigning the sample. 

`clusterMany` tries to have generally simple interface, and for this reason makes choices about what is meant by certain combinations. For example, in combinations where `findBestK=TRUE`, `ks=2:10` is taken to mean that the clustering should find the best $k$ out of the range of 2-10. However, in other combinations `ks` might indicate the specific number of clusters, $k$, that should be found. For parameter combinations that are not what is desired, the user should consider making direct calls to `clusterSingle` where all of these options (and many more) can be explicitly called. 


Other parameters for the clustering are kept fixed. As described above, there are many more possible parameters in play than are considered in `clusterMany`. These parameters can be set via the arguments `clusterDArgs`, `subsampleArgs` and `seqArgs`. These arguments correspond to the different processes described above (the clustering of $D$, the creation of $D$ via subsampling, and the sequential clustering process, respectively). These arguments take a list of arguments that are sent directly to `clusterSingle`. However, these arguments may be overridden by the interpretation of `clusterMany` of how different combinations interact; again for complete control direct calls to `clusterSingle` are necessary.


  
Argument| Dependencies |  Passed to | Argument passed to
---------------|-----------------|:-------------:|------:|
ks             | sequential=TRUE |  seqCluster   |  k0 
-    | sequential=FALSE, findBestK=FALSE, clusterFunction of type 'K' | clusterD | k
-              | sequential=FALSE, findBestK=FALSE, subsample=TRUE | subsampleClustering | k
-               | sequential=FALSE, findBestK=TRUE, clusterFunction of type 'K' | clusterD | kRange
dimReduce      | none            | transform     | dimReduce
nVarDims       | dimReduce in 'mad','cv','var' | transform | nVarDims
nPCADims       | dimReduce='PCA' | transform     | nPCADims
clusterFunction| none            | clusterD      | clusterFunction    
minSizes       | none            | clusterD      | minSize
distFunction   | subsample=FALSE | clusterD      | distFunction
alphas         | clusterFunction of type '01'| clusterD | alpha
findBestK      | clusterFunction of type 'K' | clusterD | findBestK
removeSil      | clusterFunction of type 'K' | clusterD | removeSil
silCutoff      | clusterFunction of type 'K' | clusterD | silCutoff
betas          | sequential=TRUE | seqCluster    | beta

### Example changing the distance function and clustering algorithm

Providing different distance functions is slightly more involved than the other parameters, so we give an example here. 

First we define distances that we would like to compare. We are going to define two distances that take values between 0-1 based on different choices of correlation. 


```r
 corDist<-function(x){(1-cor(t(x),method="pearson"))/2}
 spearDist<-function(x){(1-cor(t(x),method="spearman"))/2}
```

These distances are defined so as to give distance of 0 between samples with correlation 1, and distance of 1 for correlation -1.

We will also compare using different algorithms for clustering. Since we chose distances between 0-1, we can use any algorithm. Currently, `clusterMany` requires that the distances work with all of the `clusterFunction` choices given. Since some of the `clusterFunction` algorithms require a distance matrix between 0-1, this means we can only compare all of the algorithms when the distance is a 0-1 distance. (Future versions will probably try to create a work around so that the algorithm just skips algorithms that don't match the distance).

**Note on 0-1 clustering when `subsample=FALSE`** We would note that the default values $\alpha$ for the 0-1 clustering were set with the distance $D$ the result of subsampling or other concensus summary in mind. In generally, subsampling creates a $D$ matrix with  high similarity for many samples who share a cluster (the proportion of times samples are seen together for well clustered samples can easily be in the .8-.95 range, or even exactly 1). For this reason the default $\alpha$ is 0.1 which requires distances between samples in the 0.1 range or less (i.e. a similarity in the range of 0.9 or more). We show an example of the $D$ matrix from subsampling; we make use of the `clusterSingle` which is the workhorse mentioned above that runs a single clustering command directly, which gives the output $D$ from the sampling in the "coClustering" slot of `ce`. Note that the result is $1-p_{ij}$ where $p_{ij}$ is the proportion of times sample $i$ and $j$ clustered together.


```r
 ceSub<-clusterSingle(ce,dimReduce="mad",ndims=1000,
         subsample=TRUE,clusterFunction="hierarchical01",
         subsampleArgs=list(k=8),clusterLabel="subsamplingCluster",
         clusterDArgs=list(minSize=5))
 plotCoClustering(ceSub,colorScale=rev(seqPal5))
```

![plot of chunk visualizeSubsamplingD](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M4D/visualizeSubsamplingD-1.png)

We see even here, the default of $\alpha=0.1$ was perhaps too conservative since only two clusters came out (with size greater than 5).

The distances based on correlation calculated directly on the data, such as we created above, are often used for clustering expression data. But they are unlikely to have distances as low as seen in subsampling, even for well clustered samples. Here's a visualization of the correlation distance matrix we defined above (using Spearman's correlation) on the top 1000 most variable features:


```r
 dSp<-spearDist(t(transform(ce,dimReduce="mad",nVarDims=1000)))
 plotHeatmap(dSp,isSymmetric=TRUE,colorScale=rev(seqPal5))
```

![plot of chunk visualizeSpearmanDist](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M4D/visualizeSpearmanDist-1.png)

We can see that the choice of $\alpha$ must be much higher (and we are likely to be more sensitive to it).

Notice to calculate the distance in the above plot, we made use of the `transform` function applied to our `ce` object to get the results of dimensionality reduction. The `transform` function gave us a data matrix back that has been transformed, and also reduced in dimensions, like would be done in our clustering routines. `transform` has similar parameters as seen in `clusterMany`,`makeDendrogram` or `clusterSingle` and is useful when you want to manually apply something to transformed and/or dimensionality reduced data; and you can be sure you are getting the same matrix of data back that the clustering algorithms are using.


**Comparing distance functions with `clusterMany`** Now that we have defined the distances we want to compare in our global environment, we can give these to the argument "distFunction" in `clusterMany`. They should be given as character strings giving the names of the functions. For computational ease for this vignette, we will just choose the dimensionality reduction to be the top 1000 features based on MAD and set K=8 or $\alpha=0.45$. We will save these results as a separate object so as to not disrupt the earlier workflow.


```r
 ceDist<-clusterMany(ce, k=7:9, alpha=c(0.35,0.4,0.45), 
              clusterFunction=c("tight","hierarchical01","pam","hierarchicalK"),
              findBestK=FALSE,
              removeSil=c(FALSE),dist=c("corDist","spearDist"),
              isCount=TRUE,dimReduce=c("mad"),nPCADims=1000,run=TRUE)
 clusterLabels(ceDist)<-gsub("clusterFunction","alg",clusterLabels(ceDist))
 clusterLabels(ceDist)<-gsub("Dist","",clusterLabels(ceDist))
 clusterLabels(ceDist)<-gsub("distFunction","dist",clusterLabels(ceDist))
 clusterLabels(ceDist)<-gsub("hierarchical","hier",clusterLabels(ceDist))
 par(mar=c(1.1,15.1,1.1,1.1))
 plotClusters(ceDist,axisLine=-2,sampleData=c("Biological_Condition"))
```

![plot of chunk clusterManyDiffDist](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M4D/clusterManyDiffDist-1.png)

Notice that using the 01 methods did not give relevant results

### Dealing with large numbers of clusterings

A good first check before running `clusterMany` is to determine how many  clusterings you are asking for. `clusterMany` has some limited internal checks to not do unnecessary duplicates (e.g. `removeSil` only works with some clusterFunctions so `clusterMany` would detect that), but generally takes all combinations. This can take a while for more complicated clustering techniques, so it is a good idea to check what you are getting into. You can do this by running `clusterMany` with `run=FALSE`.

In the following we consider expanding our original clustering choices to consider individual choices of $K$ (rather than just `findBestK=TRUE`). 


```r
 checkParam<-clusterMany(se, clusterFunction="pam", ks=2:10,
                         removeSil=c(TRUE,FALSE), isCount=TRUE,
                         dimReduce=c("PCA","var"),
                         nVarDims=c(100,500,1000),nPCADims=c(5,15,50),run=FALSE)
 dim(checkParam$paramMatrix) #number of rows is the number of clusterings
```

```
## [1] 108  12
```

Each row of the matrix `checkParam$paramMatrix` is a requested clustering (the columns indicate the value of a possible parameter). Our selections indicate 108 different clusterings (!). 

We can set `ncores` argument to have these clusterings done in parallel. If `ncores>1`, the parallelization is done via `mclapply` and should not be done in the Rgui interface (see help pages for `mclapply`). 

<!---The output can also be fed back into clusterMany [**Check this**], which allows the user to modify the input (though do so carefully, there are not robust checks as to whether the input is of the right form). The best reason for modifying this would be to remove some of the clusterings. For example, if you only want to do some of the combinations, you could remove the unwanted ones by removing those rows of the `paramMatrix`.


We can now run this, either by giving the information in `checkParam$paramMatrix` to clusterMany argument `paramMatrix`, or by recalling the function from scratch. If the user has not changed `paramMatrix`, there's no advantage in giving `paramMatrix` to `clusterMany` rather than just recalling `clusterMany`, but we'll do it here just to show how it is done. 


```r
 # ce<-clusterMany(se,  paramMatrix=checkParam$paramMatrix, clusterDArgs=checkParam$clusterDArgs, seqArgs=checkParam$seqArgs,subsampleArgs=checkParam$subsampleArgs)
 ce<-clusterMany(ce, clusterFunction="pam",ks=2:10,findBestK=TRUE,
     removeSil=c(TRUE), isCount=TRUE,dimReduce=c("PCA","var"),
     nVarDims=c(100,500,1000),nPCADims=c(5,15,50),run=TRUE)
```

Note that we also provided in the above call the additional arguments `clusterDArgs`, `seqArgs` and `subsampleArgs` which normally we might neglect with a direct call to `clusterMany`. This is because in creating the `paramMatrix`, `clusterMany` may internally change these default values, and we want to make sure we exactly replicate what we would get from a direct call.
-->

## Create a unified cluster from many clusters with `combineMany` {#combineMany}

After creating many clusterings, `combineMany` finds a single cluster based on what samples were in the same clusters throughout the many clusters found by `clusterMany`. While subsampling the data helps be robust to outlying samples, combining across many clustering parameters can help be robust to choice in parameters, particularly when the parameters give roughly similar numbers of clusters. 

As mentioned in the Quick Start section, the default option for `combineMany` is to only define a cluster when *all* of the samples are in the same clusters across *all* clusterings. However, this is generally too conservative and just results in most samples not being assigned to a cluster. 

Instead `combineMany` has a parameter `proportion` that governs in what proportion of clusterings the samples should be together. Internally, `combineMany` makes a coClustering matrix $D$. Like the $D$ created by subsampling in `clusterMany`, the coClustering matrix takes on values 0-1 for the proportion of times the samples are together in the clustering. This $D$ matrix is saved in the `ce` object and can be visualized with `plotCoClustering` (which is just a call to `plotHeatmap`). Recall the one we last made in the QuickStart, with our last call to `combineMany` (`proportion=0.7` and `minSize=3`). 


```r
 plotCoClustering(ce)
```

![plot of chunk combineMany_detailed](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M4D/combineMany_detailed-1.png)

`combineMany` performs the clustering by running a "01" clustering algorithm on the $D$ matrix of percentage co-clustering (the default being "hierarchical01"). The `alpha` argument to the 01 clustering is `1-proportion`. Also passed to the clustering algorithm is the parameter `minSize` which sets the minimum size of a cluster.


We can also manually choose the set of clusters to use in `combineMany` with the argument `whichClusters`. Here we choose only the clusters that correspond to using dimensionality reduction using the most variable features. We also set `minSize` to be lower than the default of 5 to allow for smaller clusters


```r
 wh<-grep("nVAR",clusterLabels(ce))
 ce<-combineMany(ce,whichCluster=wh,proportion=0.7,minSize=3,
                 clusterLabel="combineMany,nVAR")
 plotCoClustering(ce)
```

![plot of chunk combineMany_chooseClusters](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M4D/combineMany_chooseClusters-1.png)

We can compare to all of our other versions of `combineMany`. While they do not all have `clusterTypes` equal to "combineMany" (only the most recent call has clusterType exactly equal to "combineMany"), they all have "combineMany" as part of their clusterType, even though they have different clusterLabels (and now we'll see that it was useful to give them different labels!)


```r
 wh<-grep("combineMany",clusterTypes(ce))
 par(mar=plotCMar)
 plotClusters(ce,whichClusters=rev(wh),axisLine=-1)
```

![plot of chunk combineMany_showDifferent](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M4D/combineMany_showDifferent-1.png)


**Treatment of Unclustered assignments** -1 values are treated separately in the calculation. In particular, they are not considered in the calculation of percentage co-clustering -- the percent co-clustering is taken only with respect to  those clusterings where both samples were assigned. However, a post-processing is done to the clusters found from running the clustering on the $D$ matrix. For each sample, the percentage of times that they were marked -1 in the clusterings is calculated. If this percentage is greater than the argument `propUnassigned` then the sample is marked as -1 in the clustering returned by `combineMany`.


**Good scenarios for running `combineMany`** Varying certain parameters result in clusterings better for `combineMany` than other sets of parameters. In particular, if there are huge discrepancies in the set of clusterings given to `combineMany`, the results will be a shattering of the samples into many small clusters. Similarly, if the number of clusters $K$ is very different, the end result will likely be like that of the large $K$, and how much value that is (rather than just picking the clustering with the largest $K$), is debatable. However, for "01" clustering algorithms or clusterings using the sequential algorithm, varying the underlying parameters $\alpha$ or $k_0$ often results in roughly similar clusterings across the parameters so that creating a consensus across them is highly informative.

## Creating a Hierarchy of Clusters and Merging clusters {#hierarchy}

As mentioned above, we find that merging clusters together based on the extent of differential expression between the features to be a useful method for combining many small clusters.

We provide a method for doing this that consists of two steps. Making a hierarchy between the clusterings and then estimating the amount of differential expression at each branch of the hierarchy.

### makeDendrogram {#makeDendrogram}

`makeDendrogram` creates a hierarchical clustering of the clusters as determined by the primaryCluster of the `ClusterExperiment` object. In addition to being used for merging clusters, the dendrograms created by `makeDendrogram` are also useful for ordering the clusters in `plotHeatmap` as has been shown above.

`makeDendrogam` performs hierarchical clustering of the cluster medoids (after transformation of the data) and provides a dendrogram that will order the samples according to this clustering of the clusters. The hierarchical ordering of the dendrograms is saved internally in the `ClusterExperiment` object.

Like clustering, the dendrogram can depend on what features are included from the data. The same options for clustering are available for the hierarchical clustering of the clusters, namely choices of dimensionality reduction via `dimReduce` and the number of dimensions via `ndims`.


```r
 ce<-makeDendrogram(ce,dimReduce="var",ndims=500)
 plotDendrogram(ce)
```

![plot of chunk makeDendrogram_reducedFeatures](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M4D/makeDendrogram_reducedFeatures-1.png)
 Notice that the plot of the dendrogram shows the hierarchy of the clusters (and color codes them according to the colors stored in colorLegend slot).
 

Recall that the most recent clustering made is from our call to `combineMany`, where we experimented with using on some of the clusterings from `clusterMany`, so that is our current primaryCluster:


```r
 show(ce)
```

```
## class: ClusterExperiment 
## dim: 7069 65 
## Primary cluster type: combineMany 
## Primary cluster label: combineMany,nVAR 
## Table of clusters (of primary clustering):
##  c1 c-1  c2  c3  c4  c5  c6  c7 
##   5  10   5   9  15  12   4   5 
## Total number of clusterings: 41 
## Dendrogram run on 'combineMany,nVAR' (cluster index: 1)
## -----------
## Workflow progress:
## clusterMany run? Yes 
## combineMany run? Yes 
## makeDendrogram run? Yes 
## mergeClusters run? No
```

This is the clustering from combining only the clusterings from `clusterMany` that use the top  most variable genes. Because it is the primaryCluster, it was the clustering that was used by default to make the dendrogram. 

We might prefer to get back to the dendrogram based on our `combineMany` in quick start (the "combineMany, final" clustering). We've lost that dendrogram when we called `makeDendrogram` again. However, we can rerun `makeDendrogram` and choose a different clustering from which to make the dendrogram. 



```r
 ce<-makeDendrogram(ce,dimReduce="var",ndims=500,
                    whichCluster="combineMany,final")
 plotDendrogram(ce)
```

![plot of chunk remakeMakeDendrogram](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M4D/remakeMakeDendrogram-1.png)

Note that the clusterType of this clustering is not "combineMany", but "combineMany.x", where "x" indicates what iteration it was:


```r
 clusterTypes(ce)[which(clusterLabels(ce)=="combineMany,final")]
```

```
## [1] "combineMany.3"
```

We can choose reset this past call to `combineMany` to the current 'combineMany' output (which will also set this clustering to be the primaryCluster). We don't need to recall `makeDendrogram`, since we already used this clustering to make the dendrogram by our explicit call of the argument `whichCluster`. 



```r
 ce<-setToCurrent(ce,whichCluster="combineMany,final")
 show(ce)
```

```
## class: ClusterExperiment 
## dim: 7069 65 
## Primary cluster type: combineMany 
## Primary cluster label: combineMany,final 
## Table of clusters (of primary clustering):
##  c1 c-1  c2  c3  c4  c5  c6 
##   9  12   7  15  14   4   4 
## Total number of clusterings: 41 
## Dendrogram run on 'combineMany,final' (cluster index: 3)
## -----------
## Workflow progress:
## clusterMany run? Yes 
## combineMany run? Yes 
## makeDendrogram run? Yes 
## mergeClusters run? No
```

### Merging clusters with little differential expression {#mergeClusters}
We then can use this hierarchy of clusters to merge clusters that show little difference in expression. We do this by testing, for each node of the dendrogram, for which features is the mean of the set of clusters to the right split of the node is equal to the mean on the left split. This is done via the `getBestFeatures` (see section on [getBestFeatures](#getBestFeatures)), where the `type` argument is set to "Dendro".

Starting at the bottom of the tree, those clusters that have the percentage of features with differential expression below a certain value (determined by the argument `cutoff`) are merged into a larger cluster. This testing of differences and merging continues until the estimated percentage of non-null DE features is above `cutoff`. This means lower values of `cutoff` result in less merging of clusters. There are multiple methods of estimation of the percentage of non-null features implemented. The option `mergeMethod="adjP"` which we showed earlier is the simplest: the proportion found significant by calculating the proportion of DE genes at a given False Discovery Rate threshold (using the Benjamini-Hochberg procedure). However, other methods are also implemented (see the help of `mergeClusters`). 

Notice that `mergeClusters` will always run based on the clustering that made the currently existing dendrogram. So it is always good to check that it is what we expect. 


```r
 ce
```

```
## class: ClusterExperiment 
## dim: 7069 65 
## Primary cluster type: combineMany 
## Primary cluster label: combineMany,final 
## Table of clusters (of primary clustering):
##  c1 c-1  c2  c3  c4  c5  c6 
##   9  12   7  15  14   4   4 
## Total number of clusterings: 41 
## Dendrogram run on 'combineMany,final' (cluster index: 3)
## -----------
## Workflow progress:
## clusterMany run? Yes 
## combineMany run? Yes 
## makeDendrogram run? Yes 
## mergeClusters run? No
```

`mergeClusters` can also be run without merging the cluster, and simply drawing a plot showing the dendrogram along with the estimates of the percentage of non-null features to aid in deciding a cutoff and method. By setting `plotType="all"`, all of the estimates of the different methods are displayed simultaneously, while before in the QuickStart, we only showed the default values.


```r
 mergeClusters(ce,mergeMethod="none",plotType="all")
```

```
## Note: Merging will be done on ' combineMany,final ', with clustering index 3
```

![plot of chunk mergeClusters_plot](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M4D/mergeClusters_plot-1.png)

Now we can pick a cutoff. We'll give it a label to keep it separate from the previous run we had made.


```r
 ce<-mergeClusters(ce,cutoff=0.05,mergeMethod="adjP",clusterLabel="mergeClusters,v2")
```

```
## Note: Merging will be done on ' combineMany,final ', with clustering index 3
```

![plot of chunk mergeClusters_ex](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M4D/mergeClusters_ex-1.png)

```r
 ce
```

```
## class: ClusterExperiment 
## dim: 7069 65 
## Primary cluster type: mergeClusters 
## Primary cluster label: mergeClusters,v2 
## Table of clusters (of primary clustering):
##  m1 m-1  m2  m3 
##  31  12   7  15 
## Total number of clusterings: 42 
## Dendrogram run on 'combineMany,final' (cluster index: 4)
## -----------
## Workflow progress:
## clusterMany run? Yes 
## combineMany run? Yes 
## makeDendrogram run? Yes 
## mergeClusters run? Yes
```


If we want to rerun `mergeClusters` with a different cutoff, we can go ahead and do that (notice in the summary above, that the dendrogram is still there, and associated with the same cluster, so we can rerun `mergeClusters` and know that it will still pull from the same input clustering from `combineMany` as we were before)


```r
 ce<-mergeClusters(ce,cutoff=0.15,mergeMethod="MB",
                   clusterLabel="mergeClusters,v3")
```

```
## Note: Merging will be done on ' combineMany,final ', with clustering index 4
```

![plot of chunk mergeClusters_redo](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M4D/mergeClusters_redo-1.png)

```r
 ce
```

```
## class: ClusterExperiment 
## dim: 7069 65 
## Primary cluster type: mergeClusters 
## Primary cluster label: mergeClusters,v3 
## Table of clusters (of primary clustering):
##  m1 m-1  m2  m3  m4  m5  m6 
##   9  12   7  15  14   4   4 
## Total number of clusterings: 43 
## Dendrogram run on 'combineMany,final' (cluster index: 5)
## -----------
## Workflow progress:
## clusterMany run? Yes 
## combineMany run? Yes 
## makeDendrogram run? Yes 
## mergeClusters run? Yes
```

## Keeping track of and rerunning elements of the workflow {#rerun}

The commands we have shown above show a workflow which continually saves the results over the previous object, so that additional information just gets added to the existing object. 

What happens if some parts of the clustering workflow are re-run? For example, in the above we reran parts of the workflow when we talked about them in more detail, or to experiment with parameter settings. 

The workflow commands check for existing clusters of the workflow (based on the `clusterTypes` of the clusterings). If there exist clusterings from previous runs *and* such clusterings came from calls that are "downstream" of the requested clustering, then the method will  change their clusterTypes value by adding a ".i", where $i$ is a numerical index keeping track of replicate calls.

For example, suppose we rerun 'combineMany', say with a different parameter choice of the proportion similarity to require. Then `combineMany` searches the existing clusterings in the input object. Any existing `combineMany` results will have their `clusterTypes` changed from `combineMany` to `combineMany.x`, where $x$ is the largest such number needed to be greater than any existing `combineMany.x` (after all, you might do this many times!). Their labels will also be updated if they just have the default label, but if the user has given different labels to the clusters those will be preserved. 

Moreover, since `mergeClusters` is downstream of `combineMany` in the workflow, currently existing `mergeClusters` will also get bumped to `mergeClusters.x`. However, `clusterMany` is upstream of `combineMany` (i.e. you expect there to be existing `clusterMany` before you run `combineMany`) so nothing will happen to `clusterMany`. 

This is handled internally, and may never be apparent to the user unless they choose `whichClusters="all"` in a plotting command. Indeed this is one reason to always pick `whichClusters="workflow"`, so that these saved previous versions are not displayed. 

However, if the user wants to "go back" to previous versions and make them the current iteration, we have seen that the `setToCurrent` command will do this. `setToCurrent` follows the same process as described above, only with an existing cluster set to the current part of the pipeline. 

Note that there is nothing that governs or protects the `clusterTypes` values to be of a certain kind. This means that if the user decides to name a clusterTypes of a clustering one of these protected names, that is allowed. However, it could clearly create some havoc if done poorly.

**Erasing old clusters** You can also choose to have **all** old versions erased by choosing the options `eraseOld=TRUE` in the call to `clusterMany`, `combineMany`,`mergeClusters` and/or `setToCurrent`. `eraseOld=TRUE` in any of these functions will delete ALL past workflow results except for those that are both in the current workflow *and* "upstream" of the requested command. You can also manually remove clusters with `removeClusters`.


**Finding workflow iterations** Sometimes which numbered iteration a particular call is in will not be obvious if there are many calls to the workflow. You may have a `mergeClusters.2` cluster but no `mergeClusters.1` because of an upstream workflow call in the middle that bumped the iteration value up to 2 without ever making a `mergeClusters.1`. If you really want to, you can see more about the existing iterations and where they are in the `clusterMatrix`. "0" refers to the current iteration; otherwise the smaller the iteration number, the earlier it was run.


```r
 workflowClusterTable(ce)
```

```
##                Iteration
## Type             0  1  2  3  4
##   final          0  0  0  0  0
##   mergeClusters  1  0  0  1  1
##   combineMany    1  1  1  0  1
##   clusterMany   36  0  0  0  0
```

Explicit details about every workflow cluster and their index in `clusterMatrix` is given by `workflowClusterDetails`:


```r
 head(workflowClusterDetails(ce),8)
```

```
##   index          type iteration               label
## 1     1 mergeClusters         0    mergeClusters,v3
## 2     2 mergeClusters         4    mergeClusters,v2
## 3     3   combineMany         4    combineMany,nVAR
## 4     4 mergeClusters         3     mergeClusters.3
## 5     5   combineMany         0   combineMany,final
## 6     6   combineMany         2     combineMany,0.7
## 7     7   combineMany         1 combineMany,default
## 8     8   clusterMany         0        nVAR=100,k=5
```

**A note on the `whichCluster` argument** Many functions take the `whichCluster` argument for identifying a clustering or clusterings on which to perform an action. These arguments all act similarly across functions, and allow the user to give character arguments. As described above, these can be shortcuts like "workflow", or they can match either clusterTypes or clusterLabels of the object. It is important to note that matching is first done to clusterTypes, and then if not successful to clusterLabels. Since neither clusterTypes nor clusterLabels is guaranteed to be unique, the user should be careful in how they make the call. And, of course, `whichCluster` arguments can also take explicit numeric integers that identify the column(s) of the clusterMatrix that should be used.

### Designate a Final Clustering

A final protected clusterTypes is "final". This is not created by any method, but can be set to be the clusterType of a clustering by the user (via the `clusterTypes` command). Any clustering marked `final` will be considered one of the workflow values for commands like `whichClusters="workflow"`. However, they will NOT be renamed with ".x" or removed if `eraseOld=TRUE`. This is a way for a user to 'save' a clustering as important/final so it will not be changed internally by any method, yet still have it show up with the "workflow" clustering results. There is no limit to the number of such clusters that are so marked, but the utility of doing so will drop if too many such clusters are chosen. 

For best functionality, particularly if a user has determined a single final clustering after completing clustering, a user will probably want to set the primaryClusterIndex to be that of the final cluster and rerun makeDendrogram. This will help in plotting and visualizing. The `setToFinal` command does this.

Here we will demonstrate marking a cluster as final. We go back to our previous mergeClusters that we found with `cutoff=0.05` and mark it as our final clustering. First we need to find which cluster it is. We see from our above call to the workflow functions above, that it is clusterType equal to "mergeClusters.4" and  label equal to "mergeClusters,v2". In our call to `setToFinal` we will decide to change it's label as well.


```r
 ce<-setToFinal(ce,whichCluster="mergeClusters,v2",
                clusterLabel="Final Clustering")
 par(mar=plotCMar)
 plotClusters(ce,whichClusters="workflow")
```

![plot of chunk markFinal](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M4D/markFinal-1.png)

Note that because it is labeled as "final" it shows up automatically in "workflow" clusters in our `plotClusters` plot. It has also been set as our primaryCluster and has the new clusterLabel we gave it in the call to `setToFinal`.

This didn't get rid of our undesired `mergeClusters` result that is most recent. It still shows up as "the" mergeClusters result. This might be undesired. We could remove that "mergeClusters" result with `removeClusters`. Alternatively, we could manually change the clusterTypes to `mergeClusters.x` so that it doesn't show up as current. 
A cleaner way to do this would have been to first set the desired cluster ("mergeClusters.4") to the most current iteration with `setToCurrent`, which would have bumped up the existing `mergeClusters` result to be no longer current. 

## RSEC

The clustering workflow described so far is a generalization of our RSEC algorithm for single-cell sequencing data. The RSEC algorithm is particularly oriented around using subsampling and sequential discovery of clusters to find very robust signals. 

In particular, `RSEC` is a single function that follows the entire workflow described above, but makes the choices to set `subsample=TRUE` and `sequential=TRUE`. Furthermore, the only clustering functions that are allowed are the "01" types ("hierarchical01" and "tight"). This removes a number of options from clusterMany, making for a slightly reduced set of commands. `RSEC` also implements the `combineMany`, `makeDendrogram` and `mergeClusters` steps, again with not all arguments available to be set. Furthermore, the defaults set in `RSEC` are those we choose for our algorithm, and occassionally vary from stand-alone method. The final output is a `clusterExperiment` object as you would get from following the workflow. 

We give the following correspondence to help see what arguments of each component are fixed by RSEC, and which are allowed to be set by the user (as well as their correspondence to arguments in the workflow functions).

  
| | Arguments Internally Fixed | Arguments for the user |  Correspondence | Notes |
|:-----------------|:-----------------|:-------------|:------|:------|
*clusterMany*     | |  | |
-| sequential=TRUE | k0s              | ks | only sets 'k0', no other k 
- | distFunction=NA | clusterFunction | | only tight or hierarchical01 
- | removeSil=FALSE | dimReduce | | 
- | subsample=TRUE | nVarDims | | 
- |  silCutoff=0 | nPCADims | | 
- |  | alphas | |
- |  | betas | |
- |  | minSizes | |
-  | | clusterDArgs | |
-  | | subsampleArgs | |
- |  | seqArgs | |
- |  | run | |
- |  | ncores | |
- |  | random.seed | |
- |  | isCount | |
- |  | transFun | |
*combineMany* | |  | |
-  | propUnassigned = *(default)* | combineProportion | proportion 
-  | combineMinSize    | minSize | |
*makeDendrogram* | |  | |
-  | ignoreUnassignedVar=TRUE | dendroReduce      | dimReduce | 
-  | unassignedSamples= *(default)* | dendroNDims       | ndims |  
*mergeClusters*  | |  | |
-  | plotType='none' | mergeMethod     | | 
-  | | mergeCutoff | cutoff | 
-  |  | isCount | | used for both mergeMethod and clusterMany 


# Finding Features related to a Clustering {#getBestFeatures}

The function `getBestFeatures` finds features in the data that are strongly differentiated between the clusters of a given clustering. Finding the best features is generally the last step in the workflow, once a final clustering has been decided upon, though as we have seen it is also called internally in `mergeClusters` to decide between which clusters to merge together. 

The function `getBestFeatures` calls `limma` on input data to determine the gene features most associated with a particular clustering.  `getBestFeatures` picks the `primaryCluster` of a `ClusterExperiment` object as the clustering to use to find features. If the standard workflow is followed, this will be the last completed step (usually the result of `mergeClusters` or manually choosing a final cluster via `setToFinal`). The primaryCluster can of course be changed by setting `primaryClusterIndex` to point to a different clustering. 

Since our most recent clustering (the one determined by our `setToFinal`) only has 2 clusterings, we are going to reset the primary clustering to be our result from `combineMany`, with the label "mergeClusters,v3". This will be better for explaining the functionality of the `getBestFeatures` method.


```r
 wh<-which(clusterLabels(ce)=="mergeClusters,v3")
 if(length(wh)==1) primaryClusterIndex(ce)<-wh else stop()
```

The basic implementation is that of `limma` which fits a linear model per feature and tests for the significance of parameters of that linear model. The main contribution of `getBestFeatures` is to interface with `limma` so as to pick appropriate parameters or tests for comparing clusters. Naturally, `getBestFeatures` also seamlessly works with `ClusterExperiment` objects to minimize the burden on the user. The output is in the form of `topTable` in `limma`, i.e. a data.frame giving the relevant features, the p-value, etc. 


## Types of Significance Tests

There are several choices of what is the most appropriate test to determine whether a feature is differentially expressed across the clusterings. All of these methods first fit a linear model where the clusters categories of the clustering is the explanatory factor in the model (samples with -1 or -2 are ignored). The methods differ only in what significance tests they then perform, which is controlled by the argument `type`. By default,  `getBestFeatures` finds significant genes based on a F-test between the clusters (`type="F"`). This is a very standard test to compare clusters, which is why it is the default, however it may not be the one that gives the best or most specific results. Indeed, in our "Quick Start", we did not use the $F$ test, but rather all pair-wise comparisons between the clusters. 

The $F$ test is a test for whether there are *any* differences in expression between the clusters for a feature. Three other options are available that try to detect instead specific kinds of differences between clusters that might be of greater interest. Specifically, these differences are encoded as "contrasts", meaning specific types of differences between the means of clusters. 

Note that for all of these contrasts, we are making use of all of the data, not just the samples in the particular cluster pairs being compared. This means the variance is estimated with all the samples. Indeed, the same linear model is being used for all of these comparisons. 

### All Pairwise
The option `type="Pairs"`, which we saw earlier, performs all pair-wise tests between the clusters for each feature, testing for each pair of clusters whether the mean of the feature is different between the two clusters. Here is the example from above using all pairwise comparisons:


```r
 pairsAllTop<-getBestFeatures(ce,contrastType="Pairs",p.value=0.05)
 dim(pairsAllTop)
```

```
## [1] 150   9
```

```r
 head(pairsAllTop)
```

```
##   IndexInOriginal Contrast Feature     logFC   AveExpr         t      P.Value    adj.P.Val
## 1            6565    X1-X2     VIM -8.608470  7.366498 -7.374382 5.322547e-10 3.865591e-07
## 2            3309    X1-X2    MIAT  7.875230  7.987637  6.479908 1.815383e-08 7.103104e-06
## 3             552    X1-X2  BCL11B  9.186117  4.120686  6.468115 1.901249e-08 7.357625e-06
## 4            4501    X1-X2 PRTFDC1 -7.635334  1.913974 -6.235150 4.724549e-08 1.522698e-05
## 5            2210    X1-X2    GLI3 -7.460028  5.121621 -6.189493 5.643513e-08 1.765221e-05
## 6            1284    X1-X2   CXADR  6.796402 10.516764  5.723578 3.408406e-07 7.608638e-05
##           B
## 1 10.834597
## 2  8.032926
## 3  7.995912
## 4  7.265297
## 5  7.122301
## 6  5.669800
```

Notice that compared to the quick start guide, We didn't set the parameter `number` which is passed to topTable, so we can get out *at most* 10 significant features for each contrast/comparison (because the default value of `number` in `topTable` is 10). Similarly, if we didn't set a value for `p.value`, `topTable` would return the top `number` genes per contrast, regardless of whether they were all significant or not. These are the defaults of `topTable`, which we purposefully do not modify, but we urge the user to read the documentation of `topTable` carefully to understand what is being asked for. In the QuickStart, we set `number=NROW(ce)` to make sure we got *all* significant genes.

In addition to the columns provided by `topTable`, the column "Contrast" tells us what pairwise contrast the result is from. "X1-X2" means a comparison of cluster 1 and cluster 2. The column "IndexInOriginal" gives the index of the gene to the original input data matrix, namely `assay(ce)`. The other columns are given by `topTable` (with the column "Feature" renamed -- it is usually "ProbeID" in `limma`).  

### One Against All
The choice `type="OneAgainsAll"` performs a comparison of a cluster against the mean of all of the other clusters. 


```r
 best1vsAll<-getBestFeatures(ce,contrastType="OneAgainstAll",p.value=0.05)
 head(best1vsAll)
```

```
##   IndexInOriginal ContrastName              Contrast Feature     logFC  AveExpr         t
## 1            6565            1 X1-(X2+X3+X4+X5+X6)/5     VIM -5.569268 7.366498 -6.372501
## 2            5753            1 X1-(X2+X3+X4+X5+X6)/5   STRAP -5.808820 7.975761 -6.358567
## 3            2960            1 X1-(X2+X3+X4+X5+X6)/5    LDHA -6.447274 7.200328 -5.436302
## 4            6090            1 X1-(X2+X3+X4+X5+X6)/5 TMEM258 -4.881551 6.563927 -5.423463
## 5             552            1 X1-(X2+X3+X4+X5+X6)/5  BCL11B  5.678299 4.120686  5.340430
## 6            2350            1 X1-(X2+X3+X4+X5+X6)/5   GSTP1 -5.005114 8.755384 -5.125553
##        P.Value    adj.P.Val        B
## 1 2.764171e-08 6.776853e-06 8.270257
## 2 2.918907e-08 6.877918e-06 8.223424
## 3 1.015365e-06 1.469819e-04 5.159886
## 4 1.065741e-06 1.516858e-04 5.117965
## 5 1.456449e-06 1.936990e-04 4.847522
## 6 3.245503e-06 3.781724e-04 4.153399
```

Notice that now there is both a "Contrast" and a "ContrastName" column. Like before, "Contrast" gives an explicit definition of what is the comparisons, in the form of "(X2+X3+X4+X5+X6)/5-X1", meaning the mean of the means of clusters 2-6 is compared to the mean of cluster1. "ContrastName" interprets this into a more usable name, namely that this contrast can be easily identified as a test of "cluster 1". 


### Dendrogram
The option `type="Dendro"` is more complex; it assumes that there is a hierarchy of the clusters (created by `makeDendrogram` and stored in the `ClusterExperiment` object).  Then for each *node* of the dendrogram, `getBestFeatures` defines a contrast or comparison of the mean expression between the daughter nodes. 

Notice the problem if we try to run the command, however.


```r
 bestDendroTest1<-getBestFeatures(ce,contrastType="Dendro",p.value=0.05)
```

```
## Error in .local(x, ...): Primary cluster does not match the cluster on which the dendrogram was made. Either replace existing dendrogram with on using the primary cluster (via 'makeDendrogram'), or reset primaryCluster with 'primaryClusterIndex' to be equal to index of 'dendo_index' slot
```

This is because the dendrogram we have does not correspond to the primaryCluster. 


```r
 show(ce)
```

```
## class: ClusterExperiment 
## dim: 7069 65 
## Primary cluster type: mergeClusters 
## Primary cluster label: mergeClusters,v3 
## Table of clusters (of primary clustering):
##  m1 m-1  m2  m3  m4  m5  m6 
##   9  12   7  15  14   4   4 
## Total number of clusterings: 43 
## Dendrogram run on 'combineMany,final' (cluster index: 5)
## -----------
## Workflow progress:
## clusterMany run? Yes 
## combineMany run? Yes 
## makeDendrogram run? Yes 
## mergeClusters run? Yes
```

`getBestFeatures` only works on the primaryCluster, and has noticed that the dendrogram in our `ce` object was not made from the primaryCluster. So we need to rerun `makeDendrogram` with our primaryCluster. It's regardless generally a good idea to make sure the dendrogram you have is the one you want. 


```r
 ce<-makeDendrogram(ce,dimReduce="var",ndims=500)
 bestDendro<-getBestFeatures(ce,contrastType="Dendro",p.value=0.05)
 head(bestDendro)
```

```
##   IndexInOriginal ContrastName Contrast Feature     logFC  AveExpr         t      P.Value
## 1            1817        Node5    X1-X6    ETV5 -9.971069 4.512155 -4.905724 7.282393e-06
## 2            6990        Node5    X1-X6  ZNF644  8.883001 5.897892  4.681077 1.640846e-05
## 3             886        Node5    X1-X6    CDK6 -8.794753 3.470607 -4.549650 2.620739e-05
## 4            1046        Node5    X1-X6    CLK3 -7.932513 2.126972 -4.498168 3.143613e-05
## 5            4581        Node5    X1-X6   PTPN1 -6.495897 1.794408 -4.485792 3.283714e-05
## 6             449        Node5    X1-X6  ATP5G1 -7.700498 4.695179 -4.435686 3.915636e-05
##      adj.P.Val        B
## 1 0.0008416443 2.771380
## 2 0.0015590244 2.156401
## 3 0.0021846701 1.801288
## 4 0.0024746326 1.663223
## 5 0.0025564506 1.630123
## 6 0.0029259648 1.496488
```

Again, there is both a "ContrastName" and "Contrast" column. The "Contrast" column identifies which clusters where on each side of the node (and hence commpared) and "ContrastName" is a name for the node.


```r
 levels((bestDendro)$Contrast)
```

```
## NULL
```

We can plot the dendrogram showing the node names to help make sense of which contrasts go with which nodes (plotDendrogram calls `plot.phylo` from the `ape` package and can take those arguments).



```r
 plotDendrogram(ce,show.node.label=TRUE)
```

![plot of chunk dendroWithNodeNames](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M4D/dendroWithNodeNames-1.png)


## Analysis for count and other RNASeq data

The `getBestFeatures` method for `ClusterExperiment` objects has an argument `isCount`. If this is marked `TRUE` then the data in `assay(x)` is assumed to be counts, and the call to `limma` uses the `voom` correction. This correction deals with the mean-variance relationship that is found with count data. This means that the differential expression analysis is done on $log_2(x+0.5)$. This is *regardless of what transformation is stored in the `ClusterExperiment` object*! The `voom` call within `getBestFeatures` however, is by default set to `normalize.method = "none"` in the call to `voom` (though the user can set `normalize.method` in the call to `getBestFeatures`).  

If instead `isCount=FALSE`, then `limma` is performed on `transform(x)`, i.e. the data after transformation of the data with the transformation stored in the `ClusterExperiment` object. In this case, there is no `voom` correction. 

Unlike edgeR or DESeq, the voom correction does not explicitly require a count matrix, and therefore it has been proposed that it can be used on FPKM or TPM entries, or data normalized via RUV. Setting `isCount=TRUE` even if the data in the assay slot is not count will have this effect. However, the authors of the package do not recommend using voom on anything other than counts, see e.g. [this discussion](https://support.bioconductor.org/p/45749/). 

## Piping into other DE routines

Ultimately, for many settings, the user may prefer to use other techniques for differential expression analysis or have more control over certain aspects of it. The function `clusterContrasts` may be called by the user to get the contrasts that are defined within `getBestFeatures` (e.g. dendrogram contrasts or pairwise contrasts). These contrasts, which are in the format needed for `limma` can be piped into programs that allow for contrasts in their linear models like edgeR [@Robinson:2010cw] for mRNA-Seq or MAST [@Finak:2015id] for single-cell sequencing. 

Similarly, more complicated normalizations, like RUV [@GagnonBartsch:2011jv], adjust each gene individually for unwanted batch or other variation within the linear model. In this case, a matrix $W$ that describes this variation should be included in the linear model. Again, this can be done in other programs, using the contrasts provided by `clusterContrasts`


## Additional considerations

The user should be careful about questions of multiple comparisons when all of these multiple contrasts are being performed on each feature; the default is to correct across all of these tests (see the help of `getBestFeatures` and the argument `contrastAdj` for more). As noted in the introduction, p-values created in this way are reusing the data (since the data was also used for creating the clusters) and hence should not be considered valid p-values regardless.


As mentioned, `getBestFeatures` accepts arguments to `limma`'s function `topTable` to decide which genes should be returned (and in what order). In particular, we can set an adjusted p-value cutoff for each contrast, and set `number` to control the number of genes returned *for each contrast*. By setting `number` to be the length of all genes, and `p.value=0.05`, we can return all genes for each contrast that have adjusted p-values less than 0.05. All of the arguments to `topTable` regarding what results are returned and in what order can be given by the user at the call to `getBestFeatures`.


# References
 
 ### Parameter settings:
   * WRKDIR = /mnt100/home/Dropbox/SingleCell/Jun2017/R
   * FN = M4D-clusteringTutorial
   * Scripts = Scripts
   * RUN DATE = Tue Jun 27 02:26:46 2017
 
 

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
##  [1] clusterExperiment_1.2.0    scRNAseq_1.2.0             SummarizedExperiment_1.6.3
##  [4] DelayedArray_0.2.7         matrixStats_0.52.2         Biobase_2.36.2            
##  [7] GenomicRanges_1.28.3       GenomeInfoDb_1.12.2        IRanges_2.10.2            
## [10] S4Vectors_0.14.3           BiocGenerics_0.22.0        BiocParallel_1.10.1       
## [13] knitr_1.16                 rmarkdown_1.6             
## 
## loaded via a namespace (and not attached):
##  [1] nlme_3.1-131            bitops_1.0-6            bold_0.4.0              progress_1.1.2         
##  [5] doParallel_1.0.10       RColorBrewer_1.1-2      httr_1.2.1              rprojroot_1.2          
##  [9] prabclus_2.2-6          tools_3.4.0             backports_1.1.0         R6_2.2.2               
## [13] lazyeval_0.2.0          colorspace_1.3-2        ade4_1.7-6              trimcluster_0.1-2      
## [17] nnet_7.3-12             prettyunits_1.0.2       gridExtra_2.2.1         compiler_3.4.0         
## [21] xml2_1.1.1              pkgmaker_0.22           diptest_0.75-7          scales_0.4.1           
## [25] DEoptimR_1.0-8          mvtnorm_1.0-6           robustbase_0.92-7       NMF_0.20.6             
## [29] stringr_1.2.0           digest_0.6.12           XVector_0.16.0          htmltools_0.3.6        
## [33] highr_0.6               limma_3.32.2            rlang_0.1.1             howmany_0.3-1          
## [37] jsonlite_1.5            mclust_5.3              dendextend_1.5.2        dplyr_0.7.0            
## [41] RCurl_1.95-4.8          magrittr_1.5            modeltools_0.2-21       GenomeInfoDbData_0.99.0
## [45] Matrix_1.2-10           Rcpp_0.12.11            munsell_0.4.3           ape_4.1                
## [49] abind_1.4-5             viridis_0.4.0           stringi_1.1.5           whisker_0.3-2          
## [53] MASS_7.3-47             zlibbioc_1.22.0         flexmix_2.3-14          MAST_1.2.1             
## [57] plyr_1.8.4              grid_3.4.0              rncl_0.8.2              lattice_0.20-35        
## [61] splines_3.4.0           uuid_0.1-2              taxize_0.8.4            fpc_2.1-10             
## [65] rngtools_1.2.4          reshape2_1.4.2          codetools_0.2-15        XML_3.98-1.8           
## [69] glue_1.1.0              evaluate_0.10           RNeXML_2.0.7            data.table_1.10.4      
## [73] foreach_1.4.3           locfdr_1.1-8            gtable_0.2.0            tidyr_0.6.3            
## [77] reshape_0.8.6           kernlab_0.9-25          assertthat_0.2.0        ggplot2_2.2.1          
## [81] gridBase_0.4-7          phylobase_0.8.4         xtable_1.8-2            class_7.3-14           
## [85] viridisLite_0.2.0       tibble_1.3.3            iterators_1.0.8         registry_0.3           
## [89] cluster_2.0.6
```

