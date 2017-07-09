### Seurat - Guided Clustering Tutorial


<style type="text/css"> body, td { font-size: 14px; } code.r{ font-size: 12px; } pre { font-size: 12px } </style>

[Seurat - Guided Clustering Tutorial](https://http://satijalab.org/seurat/pbmc-tutorial.html)



<!-- ***************************************************** -->


<!-- ***************************************************** -->



<!-- ***************************************************** -->



<!-- ***************************************************** -->


### Setup the Seurat Object

In this tutorial, we will be analyzing the a dataset of Peripheral Blood Mononuclear Cells (PBMC) freely available from 10X Genomics. There are 2,700 single cells that were sequenced on the Illumina NextSeq 500. The raw data can be found [here](https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz).

We start by reading in the data. All features in Seurat have been configured to work with both regular and sparse matrices, but sparse matrices result in significant memory and speed savings for Drop-seq/inDrop/10x data.


```r
 suppressMessages(require(Seurat))
 suppressMessages(require(dplyr))
 suppressMessages(require(Matrix))

 # Load the PBMC dataset
 pbmc.data <- Read10X(EXT_DATA)

 #Examine the memory savings between regular and sparse matrices
 dense.size <- object.size(as.matrix(pbmc.data))
 dense.size
```

```
## 709264728 bytes
```

```r
 sparse.size <- object.size(pbmc.data)
 sparse.size
```

```
## 38715120 bytes
```

```r
 dense.size/sparse.size
```

```
## 18.3 bytes
```

```r
 # Initialize the Seurat object with the raw (non-normalized data)
 # Note that this is slightly different than the older Seurat workflow, 
 # where log-normalized values were passed in directly.
 # You can continue to pass in log-normalized values, just set do.logNormalize=F 
 # in the next step.
 pbmc <- new("seurat", raw.data = pbmc.data)

 # Keep all genes expressed in >= 3 cells, keep all cells with >= 200 genes
 # Perform log-normalization, first scaling each cell to a total of 1e4 molecules 
 # (as in Macosko et al. Cell 2015)
 pbmc <- Setup(pbmc, min.cells = 3, min.genes = 200, 
              do.logNormalize = T, total.expr = 1e4, project = "10X_PBMC")
```
<!-- ***************************************************** -->


### Basic QC and selecting cells for further analysis

While the setup function imposes a basic minimum gene-cutoff, you may want to filter out cells at this stage based on technical or biological parameters. Seurat allows you to easily explore QC metrics and filter cells based on any user-defined criteria. In the example below, we visualize gene and molecule counts, plot their relationship, and exclude cells with a clear outlier number of genes detected as potential **multiplets**. Of course this is not a guarenteed method to exclude **cell doublets**, but we include this as an example of filtering user-defined outlier cells. We also filter cells based on the percentage of mitochondrial genes present.


```r
 # The number of genes and UMIs (nGene and nUMI) are automatically calculated for 
 # every object by Seurat.  For non-UMI data, nUMI represents the sum of the non-normalized 
 # values within a cell.  We calculate the percentage of mitochondrial genes here 
 # and store it in percent.mito using the AddMetaData. 
 # The % of UMI mapping to MT-genes is a common scRNA-seq QC metric.
 # NOTE: You must have the Matrix package loaded to calculate the percent.mito values. 
 mito.genes <- grep("^MT-", rownames(pbmc@data), value = T)
 percent.mito <- colSums(expm1(pbmc@data[mito.genes, ]))/colSums(expm1(pbmc@data))

 #n AddMetaData adds columns to object@data.info, and is a great place to stash QC stats
 pbmc <- AddMetaData(pbmc, percent.mito, "percent.mito")
 VlnPlot(pbmc, c("nGene", "nUMI", "percent.mito"), nCol = 3)
```

![plot of chunk basicQC](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M2A/basicQC-1.png)

```r
 # GenePlot is typically used to visualize gene-gene relationships, 
 # but can be used for anything calculated by the object, 
 # i.e. columns in object@data.info, PC scores etc.
 # Since there is a rare subset of cells with an outlier level of high 
 # mitochondrial percentage, and also low UMI content, we filter these as well
 par(mfrow = c(1, 2))
 GenePlot(pbmc, "nUMI", "percent.mito")
 GenePlot(pbmc, "nUMI", "nGene") 
```

![plot of chunk basicQC](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M2A/basicQC-2.png)

```r
 # We filter out cells that have unique gene counts over 2,500
 # Note that accept.high and accept.low can be used to define a 'gate', 
 # and can filter cells not only based on nGene but on anything in the object
 # (as in GenePlot above)
 pbmc <- SubsetData(pbmc, subset.name = "nGene", accept.high = 2500)
 pbmc <- SubsetData(pbmc, subset.name = "percent.mito", accept.high = 0.05)

 # Redo gene plots

 par(mfrow = c(1, 2))
 GenePlot(pbmc, "nUMI", "percent.mito")
 GenePlot(pbmc, "nUMI", "nGene") 
```

![plot of chunk basicQC](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M2A/basicQC-3.png)

### Regress out unwanted sources of variation

Your single cell dataset likely contains ‘uninteresting’ sources of variation. This could include not only technical noise, but batch effects, or even biological sources of variation (cell cycle stage). As suggested in Buettner et al, NBT, 2015, regressing these signals out of the analysis can improve downstream dimensionality reduction and clustering. Seurat implements a basic version of this by constructing linear models to predict gene expression based on user-defined variables. Seurat stores the z-scored residuals of these models in the scale.data slot, and they are used for dimensionality reduction and clustering.

We typically regress out cell-cell variation in gene expression driven by batch (if applicable), cell alignment rate (as provided by Drop-seq tools for Drop-seq data), the number of detected molecules, and mitochondrial gene expression. For cycling cells, we can also learn a ‘cell-cycle’ score (as in Macosko et al) and Regress this out as well. In this simple example here for post-mitotic blood cells, we simply regress on the number of detected molecules per cell as well as the percentage mitochondrial gene content an example.


```r
 # Note that this overwrites pbmc@scale.data. Therefore, if you intend to use RegressOut, 
 # you can set do.scale=F and do.center=F in the original object to save some time.
 pbmc <- RegressOut(pbmc, latent.vars = c("nUMI", "percent.mito"))
```

```
## [1] "Regressing out nUMI"         "Regressing out percent.mito"
```

### Detection of variable genes across the single cells

Seurat calculates highly variable genes and focuses on these for downstream analysis. **MeanVarPlot**(), which works by calculating the average expression and dispersion for each gene, placing these genes into bins, and then calculating a z-score for dispersion within each bin. This helps control for the relationship between variability and average expression. This function is unchanged from (Macosko et al.), but new methods for variable gene expression identification are coming soon. We suggest that users set these parameters to mark visual outliers on the dispersion plot, but the exact parameter settings may vary based on the data type, heterogeneity in the sample, and normalization strategy. The parameters here identify ~2,000 variable genes, and represent typical parameter settings for UMI data that is normalized to a total of 1e4 molecules.


```r
 pbmc <- MeanVarPlot(pbmc ,fxn.x = expMean, fxn.y = logVarDivMean, 
         x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, do.contour = F)
```

![plot of chunk detectVG](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M2A/detectVG-1.png)

```r
 length(pbmc@var.genes)
```

```
## [1] 1838
```

### Perform linear dimensional reduction

Perform PCA on the scaled data. By default, the genes in object@var.genes are used as input, but can be defined using pc.genes. We have typically found that running dimensionality reduction on genes with high-dispersion can improve performance. However, with UMI data - particularly after using RegressOut, we often see that PCA returns similar (albeit slower) results when run on much larger subsets of genes, including the whole transcriptome.


```r
 pbmc <- PCA(pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, 
             pcs.print = 5, genes.print = 5)
```

```
## [1] "PC1"
## [1] "CST3"   "TYROBP" "FCN1"   "LST1"   "AIF1"  
## [1] ""
## [1] "PTPRCAP" "IL32"    "LTB"     "CD2"     "CTSW"   
## [1] ""
## [1] ""
## [1] "PC2"
## [1] "NKG7" "GZMB" "PRF1" "CST7" "GZMA"
## [1] ""
## [1] "CD79A"    "MS4A1"    "HLA-DQA1" "TCL1A"    "HLA-DQB1"
## [1] ""
## [1] ""
## [1] "PC3"
## [1] "CYBA"     "HLA-DPA1" "HLA-DPB1" "HLA-DRB1" "CD37"    
## [1] ""
## [1] "PF4"   "PPBP"  "SDPR"  "SPARC" "GNG11"
## [1] ""
## [1] ""
## [1] "PC4"
## [1] "IL32"   "GIMAP7" "AQP3"   "FYB"    "MAL"   
## [1] ""
## [1] "CD79A"    "HLA-DQA1" "CD79B"    "MS4A1"    "HLA-DQB1"
## [1] ""
## [1] ""
## [1] "PC5"
## [1] "FCER1A"  "LGALS2"  "MS4A6A"  "S100A8"  "CLEC10A"
## [1] ""
## [1] "FCGR3A"        "CTD-2006K23.1" "IFITM3"        "ABI3"          "CEBPB"        
## [1] ""
## [1] ""
```

ProjectPCA scores each gene in the dataset (including genes not included in the PCA) based on their correlation with the calculated components. Though we don’t use this further here, it can be used to identify markers that are strongly correlated with cellular heterogeneity, but may not have passed through variable gene selection. The results of the projected PCA can be explored by setting use.full=T in the functions below.


```r
 pbmc <- ProjectPCA(pbmc)
```

```
## [1] "PC1"
##  [1] "CST3"     "S100A9"   "FTL"      "TYROBP"   "FCN1"     "LYZ"      "LST1"     "AIF1"    
##  [9] "S100A8"   "FTH1"     "TYMP"     "LGALS2"   "CFD"      "FCER1G"   "LGALS1"   "CD68"    
## [17] "CTSS"     "SERPINA1" "SAT1"     "NPC2"     "GRN"      "CFP"      "IFITM3"   "COTL1"   
## [25] "IFI30"    "PSAP"     "SPI1"     "CD14"     "GPX1"     "MS4A6A"  
## [1] ""
##  [1] "MALAT1"  "RPS27A"  "PTPRCAP" "RPSA"    "RPS3A"   "RPL23A"  "IL32"    "RPL3"    "RPL9"   
## [10] "RPS27"   "CD3D"    "RPL21"   "LTB"     "RPS3"    "RPS6"    "RPS15A"  "RPL31"   "CD3E"   
## [19] "RPL13A"  "LDHB"    "RPS25"   "RPL30"   "RPS12"   "RPS18"   "RPS23"   "RPS29"   "RPL27A" 
## [28] "RPL5"    "EEF1A1"  "AES"    
## [1] ""
## [1] ""
## [1] "PC2"
##  [1] "NKG7"   "GZMB"   "PRF1"   "CST7"   "GZMA"   "FGFBP2" "GNLY"   "CTSW"   "SPON2"  "CCL4"  
## [11] "GZMH"   "CCL5"   "FCGR3A" "KLRD1"  "XCL2"   "SRGN"   "CLIC3"  "GZMM"   "B2M"    "CD247" 
## [21] "AKR1C3" "S100A4" "PRSS23" "TTC38"  "S1PR5"  "HCST"   "IGFBP7" "ITGB2"  "HOPX"   "GPR56" 
## [1] ""
##  [1] "CD79A"     "MS4A1"     "RPL18A"    "HLA-DQA1"  "TCL1A"     "HLA-DQB1"  "CD79B"     "LINC00926"
##  [9] "LTB"       "RPL32"     "VPREB3"    "RPL13A"    "HLA-DRA"   "RPL13"     "RPS23"     "RPS27"    
## [17] "RPL11"     "RPS18"     "RPS5"      "RPL8"      "HLA-DQA2"  "RPLP2"     "RPL12"     "RPS2"     
## [25] "CD37"      "FCER2"     "CD74"      "BANK1"     "RPS12"     "HLA-DRB1" 
## [1] ""
## [1] ""
## [1] "PC3"
##  [1] "RPL10"  "RPS2"   "RPL11"  "RPL28"  "RPL32"  "RPL18A" "RPL12"  "RPL19"  "RPS6"   "TMSB10"
## [11] "RPS14"  "RPL13"  "RPS19"  "RPS15"  "RPL6"   "RPL29"  "RPLP1"  "RPS3"   "RPL15"  "RPL26" 
## [21] "RPS4X"  "EEF1A1" "RPS16"  "RPS7"   "RPL14"  "GNB2L1" "RPL13A" "RPL3"   "RPS18"  "RPL8"  
## [1] ""
##  [1] "PF4"          "PPBP"         "SDPR"         "SPARC"        "GNG11"        "HIST1H2AC"   
##  [7] "GP9"          "NRGN"         "TUBB1"        "RGS18"        "CLU"          "AP001189.4"  
## [13] "ITGA2B"       "CD9"          "PTCRA"        "TMEM40"       "CA2"          "ACRBP"       
## [19] "MMD"          "TREML1"       "F13A1"        "SEPT5"        "PGRMC1"       "MYL9"        
## [25] "TSC22D1"      "MPP1"         "CMTM5"        "PTGS1"        "SNCA"         "RP11-367G6.3"
## [1] ""
## [1] ""
## [1] "PC4"
##  [1] "CD3D"    "LDHB"    "IL7R"    "RPS14"   "CD3E"    "VIM"     "IL32"    "RPL32"   "RPS12"  
## [10] "NOSIP"   "RPL28"   "RPS25"   "GIMAP7"  "RPL11"   "JUNB"    "RPL13"   "RPS3"    "AQP3"   
## [19] "ZFP36L2" "FYB"     "RPL10"   "RGCC"    "MAL"     "FOS"     "LEF1"    "RPLP1"   "CD2"    
## [28] "RPL35A"  "RPL36"   "RPS28"  
## [1] ""
##  [1] "CD79A"     "HLA-DQA1"  "CD79B"     "MS4A1"     "CD74"      "HLA-DQB1"  "HLA-DPB1"  "HLA-DPA1" 
##  [9] "HLA-DRB1"  "TCL1A"     "HLA-DRA"   "LINC00926" "HLA-DRB5"  "HLA-DQA2"  "VPREB3"    "HLA-DMA"  
## [17] "GZMB"      "HLA-DMB"   "HVCN1"     "FCER2"     "BANK1"     "FGFBP2"    "HLA-DOB"   "PDLIM1"   
## [25] "FCRLA"     "TSPAN13"   "PRF1"      "GNLY"      "CD72"      "EAF2"     
## [1] ""
## [1] ""
## [1] "PC5"
##  [1] "FCER1A"   "LGALS2"   "MS4A6A"   "S100A8"   "CLEC10A"  "FOLR3"    "GPX1"     "CD14"    
##  [9] "GSTP1"    "ALDH2"    "SERPINF1" "S100A12"  "ID1"      "CD1C"     "GRN"      "RNASE6"  
## [17] "GSN"      "IER3"     "CSF3R"    "BLVRB"    "RPL13"    "ASGR1"    "S100A9"   "SAT2"    
## [25] "LYZ"      "RNASE2"   "VCAN"     "QPCT"     "CEBPD"    "FCGR1A"  
## [1] ""
##  [1] "FCGR3A"        "CDKN1C"        "MS4A7"         "HES4"          "CKB"           "RP11-290F20.3"
##  [7] "RHOC"          "CTSL"          "MS4A4A"        "LILRA3"        "SIGLEC10"      "CTD-2006K23.1"
## [13] "IFITM3"        "HMOX1"         "ABI3"          "IFITM2"        "LRRC25"        "LILRB1"       
## [19] "BATF3"         "PTP4A3"        "CEBPB"         "PILRA"         "CSF1R"         "HCK"          
## [25] "CXCL16"        "VMO1"          "C1QA"          "TPPP3"         "TNFSF10"       "TCF7L2"       
## [1] ""
## [1] ""
```

Seurat provides several useful ways of visualizing both cells and genes that define the PCA, including PrintPCA(), VizPCA(), PCAPlot(), and PCHeatmap()


```r
 PrintPCA(pbmc, pcs.print = 1:5, genes.print = 5, use.full = TRUE)
```

```
## [1] "PC1"
## [1] "CST3"   "S100A9" "FTL"    "TYROBP" "FCN1"  
## [1] ""
## [1] "MALAT1"  "RPS27A"  "PTPRCAP" "RPSA"    "RPS3A"  
## [1] ""
## [1] ""
## [1] "PC2"
## [1] "NKG7" "GZMB" "PRF1" "CST7" "GZMA"
## [1] ""
## [1] "CD79A"    "MS4A1"    "RPL18A"   "HLA-DQA1" "TCL1A"   
## [1] ""
## [1] ""
## [1] "PC3"
## [1] "RPL10" "RPS2"  "RPL11" "RPL28" "RPL32"
## [1] ""
## [1] "PF4"   "PPBP"  "SDPR"  "SPARC" "GNG11"
## [1] ""
## [1] ""
## [1] "PC4"
## [1] "CD3D"  "LDHB"  "IL7R"  "RPS14" "CD3E" 
## [1] ""
## [1] "CD79A"    "HLA-DQA1" "CD79B"    "MS4A1"    "CD74"    
## [1] ""
## [1] ""
## [1] "PC5"
## [1] "FCER1A"  "LGALS2"  "MS4A6A"  "S100A8"  "CLEC10A"
## [1] ""
## [1] "FCGR3A" "CDKN1C" "MS4A7"  "HES4"   "CKB"   
## [1] ""
## [1] ""
```

```r
 VizPCA(pbmc, 1:2)
```

![plot of chunk explorePCs](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M2A/explorePCs-1.png)

```r
 PCAPlot(pbmc, 1, 2)
```

![plot of chunk explorePCs](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M2A/explorePCs-2.png)

In particular PCHeatmap() allows for easy exploration of the primary sources of heterogeneity in a dataset, and can be useful when trying to decide which PCs to include for further downstream analyses. Both cells and genes are ordered according to their PCA scores. Setting cells.use to a number plots the ‘extreme’ cells on both ends of the spectrum, which dramatically speeds plotting for large datasets. Though clearly a supervised analysis, we find this to be a valuable tool for exploring correlated gene sets.



```r
 PCHeatmap(pbmc, pc.use = 1, cells.use = 100, do.balanced = TRUE)
```

![plot of chunk PCHeatmap](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M2A/PCHeatmap-1.png)

```r
 PCHeatmap(pbmc, pc.use = 1:12, cells.use = 500, do.balanced = TRUE,
           label.columns = FALSE, use.full = FALSE)
```

![plot of chunk PCHeatmap](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M2A/PCHeatmap-2.png)

### Determine statistically significant principal components

To overcome the extensive technical noise in any single gene for scRNA-seq data, Seurat clusters cells based on their PCA scores, with each PC essentially representing a ‘metagene’ that combines information across a correlated gene set. Determining how many PCs to include downstream is therefore an important step.

In Macosko et al, we implemented a resampling test inspired by the jackStraw procedure. We randomly permute a subset of the data (1% by default) and rerun PCA, constructing a ‘null distribution’ of gene scores, and repeat this procedure. We identify ‘significant’ PCs as those who have a strong enrichment of low p-value genes.


```r
 # NOTE: This process can take a long time for big datasets, comment out for expediency.
 # More approximate techniques such as those implemented in PCElbowPlot() 
 # can be used to reduce computation time
 pbmc <- JackStraw(pbmc, num.replicate = 100, do.print = FALSE)
```

The JackStrawPlot() function provides a visualization tool for comparing the distribution of p-values for each PC with a uniform distribution (dashed line). ‘Significant’ PCs will show a strong enrichment of genes with low p-values (solid curve above the dashed line). In this case it appears that PCs 1-10 are significant.


```r
JackStrawPlot(pbmc, PCs = 1:12)
```

![plot of chunk JackStrawPlot](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M2A/JackStrawPlot-1.png)
A more ad hoc method for determining which PCs to use is to look at a plot of the standard deviations of the principle components and draw your cutoff where there is a clear elbow in the graph. This can be done with PCElbowPlot(). In this example, it looks like the elbow would fall around PC 9.


```r
 PCElbowPlot(pbmc)
```

![plot of chunk PCElbowPlot](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M2A/PCElbowPlot-1.png)

PC selection – identifying the true dimensionality of a dataset – is an important step for Seurat, but can be challenging/uncertain for the user. We therefore suggest these three approaches to consider. The first is more supervised, exploring PCs to determine relevant sources of heterogeneity, and could be used in conjunction with GSEA for example. The second implements a statistical test based on a random null model, but is time-consuming for large datasets, and may not return a clear PC cutoff. The third is a heuristic that is commonly used, and can be calculated instantly. In this example all three approaches yielded similar results, but we might have been justified in choosing anything between PC 7-10 as a cutoff. We followed the jackStraw here, admittedly buoyed by seeing the PCHeatmap returning interpretable signals (including canonical dendritic cell markers) throughout these PCs. Though the results are only subtly affected by small shifts in this cutoff (you can test below), we strongly suggest always explore the PCs they choose to include downstream.


### Cluster the cells

Seurat now includes an graph-based clustering approach compared to (Macosko et al.). Importantly, the *distance metric* which drives the clustering analysis (based on previously identified PCs) remains the same. However, our approach to partioning the cellular distance matrix into clusters has dramatically improved. Our approach was heavily inspired by recent manuscripts which applied graph-based clustering approaches to scRNA-seq data **[SNN-Cliq, Xu and Su, Bioinformatics, 2015]** and CyTOF data **[PhenoGraph, Levine et al., Cell, 2015]**. Briefly, these methods embed cells in a graph structure - for example a K-nearest neighbor (KNN) graph, with edges drawn between cells with similar gene expression patterns, and then attempt to partition this graph into highly interconnected ‘quasi-cliques’ or ‘communities’. As in PhenoGraph, we first construct a KNN graph based on the euclidean distance in PCA space, and refine the edge weights between any two cells based on the shared overlap in their local neighborhoods (Jaccard distance). To cluster the cells, we apply the smart local moving algorithm **[SLM, Blondel et al., Journal of Statistical Mechanics]**, to iteratively group cell groupings together with the goal of optimizing the standard modularity function.

The **FindClusters**() function implements the procedure, and contains a resolution parameter that sets the ‘granularity’ of the downstream clustering, with increased values leading to a greater number of clusters. We find that setting this parameter between 0.6-1.2 typically returns good results for single cell datasets of around 3K cells. Optimal resolution often increases for larger datasets. The clusters are saved in the object@ident slot.

#### HACK:
* To fix error:   cannot open file '/usr/local/lib/R/site-library/output39550.txt': No such file or directory
* sudo cp /usr/local/lib/R/site-library/edge_39550.txt /usr/local/lib/R/site-library/output39550.txt























