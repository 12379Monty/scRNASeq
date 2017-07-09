### Monocle:  Cell counting, differential expression, and trajectoryanalysis for single-cell RNA-Seq experiments


<style type="text/css"> body, td { font-size: 14px; } code.r{ font-size: 12px; } pre { font-size: 12px } </style>

[Monocle Vignette](https://bioconductor.org/packages/release/bioc/html/monocle.html)



<!-- ***************************************************** -->


<!-- ***************************************************** -->



<!-- ***************************************************** -->



<!-- ***************************************************** -->


```r
 ###################################################
 ### code chunk number 1: package_loads
 ###################################################
 suppressMessages(require(Biobase))
 suppressMessages(require(knitr))
 suppressMessages(require(reshape2))
 suppressMessages(require(ggplot2))
 
 knitr::opts_chunk$set(autodep=TRUE, cache=FALSE, warning=FALSE)
 set.seed(0)
 
 
 ###################################################
 ### code chunk number 2: init_monocle
 ###################################################
 suppressMessages(require(HSMMSingleCell))
 suppressMessages(require(monocle))
```

### Introduction

The monocle package provides a toolkit for analyzing single cell gene expression experiments.  This vignette provides an  overview  of  a  single  cell  RNA-Seq  analysis  workflow  with  Monocle.  Monocle  was  originally  developed  to  analyze dynamic biological processes such as cell differentiation, although it also supports simpler experimental settings.

**Monocle 2** includes new, improved algorithms for classifying and counting cells, performing differential expression analysis between subpopulations of cells, and reconstructing cellular trajectories.  Monocle 2 has also been re-engineered to work well with very large single-cell RNA-Seq experiments containing tens of thousands of cells or even more.  Monocle can help you perform three main types of analysis:

1. **Clustering, classifying, and counting cells**.  Single-cell RNA-Seq experiments allow you to discover new (andpossibly rare) subtypes of cells.  Monocle helps you identify them
2. **Constructing  single-cell  trajectories**.  In  development,  disease,  and  throughout  life,  cells  transition  from  onestate to another.  Monocle helps you discover these transitions
3. **Differential expression analysis**.  Characterizing new cell types and states begins with comparing them to other,better understood cells.  Monocle includes a sophisticated but easy to use system for differential expression

Before we look at Monocle’s functions for each of these common analysis tasks, let’s see how to load up single-celldatasets in Monocle.

### Getting started with Monocle

The *monocle* package takes a matrix of gene expression values as calculated by Cufflinks [2] or another gene expression estimation program.  Monocle can work with relative expression values (e.g.  FPKM or TPM units) or absolute transcript counts (e.g.  from UMI experiments).  Monocle also works “out-of-the-box” with the transcript count matrices producedby CellRanger, the software pipeline for analyzing experiments from the 10X Genomics Chromium instrument.  Monocle also works well with data from other RNA-Seq workflows such as sci-RNA-Seq and instruments like the Biorad ddSEQ.  Although  Monocle  can  be  used  with  raw  read  counts,  these  are  not  directly  proportional  to  expression  values  unless you  normalize  them  by  length,  so  some  Monocle  functions  could  produce  nonsense  results.   If  you  don’t  have  UMI counts, We recommend you load up FPKM or TPM values instead of raw read counts.


#### The CellDataSet class

monocle holds single cell expression data in objects of the **CellDataSetclass**.  The class is derived from the Bioconductor **ExpressionSet** class, which provides a common interface familiar to those who have analyzed microarray experiments with Bioconductor.  The class requires three input files:

1. **exprs**, a numeric matrix of expression values, where rows are genes, and columns are cells

2. **phenoData**, an AnnotatedDataFrame object, where rows are cells, and columns are cell attributes (such as celltype, culture condition, day captured, etc.)

3. **featureData**,  an AnnotatedDataFrame object,  where  rows  are  features  (e.g.   genes),  and  columns  are  gene attributes, such as biotype, gc content, etc.

The expression value matrix must have the same number of columns as the phenoData has rows, and it must have the  same  number  of  rows  as  the featureData data  frame  has  rows.   Row  names  of  the phenoDataobject  should match the column names of the expression matrix.  Row names of the feature Dataobject should match row names of the expression matrix.  Also, one of the columns of the featureData must be named ”geneshortname”.

You can create a new CellDataSet object as follows:


```r
 # Load data 
 data(HSMM_expr_matrix)
 data(HSMM_gene_annotation)
 data(HSMM_sample_sheet)


 # Once these tables are loaded, create the CellDataSet object like this:
 pd <- new("AnnotatedDataFrame", data = HSMM_sample_sheet)
 fd <- new("AnnotatedDataFrame", data = HSMM_gene_annotation)

 HSMM <- newCellDataSet(as(as.matrix(HSMM_expr_matrix), "sparseMatrix"),
                        phenoData = pd, 
                        featureData = fd,
                        lowerDetectionLimit=0.5,
                        expressionFamily=negbinomial.size())
```


This will create a CellDataSet object with expression values measured in FPKM, a measure of relative expression reported by Cufflinks.  By default, Monocle assumes that your expression data is in units of transcript counts and usesa  negative  binomial  model  to  test  for  differential  expression  in  downstream  steps.   However,  if  you’re  using  relative expression  values  such  as  TPM  or  FPKM  data,  see  below  for  how  to  tell  Monocle  how  to  model  it  in  downstream steps.

NOTE: if you do have UMI data, you should not normalize it yourself prior to creating your CellDataSet.  You should  also not try  to  convert  the  UMI  counts  to  relative  abundances  (by  converting  it  to  FPKM/TPM  data).  You should not use relative2abs() as  discussed  below  in  section  2.4.   Monocle  will  do  all  needed  normalization  steps internally.  Normalizing it yourself risks breaking some of Monocle’s key steps. 
 
#### Choosing a distribution for expression data

Monocle works well with both relative expression data and count-based measures (e.g.  UMIs).  In general, it works best with  transcript  count  data,  especially  UMI  data.  Whatever  your  data  type,  it  is critical to  specify  the  appropriate distribution  for  it.   FPKM/TPM  values  are  generally  log-normally  distributed,  while  UMIs  or  read  counts  are  better modeled  with  the  negative  binomial.   To  work  with  count  data,  specify  the  negative  binomial  distribution  as  the **expressionFamily** argument to **newCellDataSet**:


```r
 # NOT RUN
 HSMM <- newCellDataSet(count_matrix,
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily=negbinomial.size())
```

There are several allowed values for expressionFamily, which expects a“family function” from the **VGAM** package.

**Using  the  wrong  expressionFamily  for  your  data  will  lead  to  bad  results,  errors  from  Monocle,  or  both**.  However, if you have FPKM/TPM data, you can still use negative binomial if you first convert your relative expression values to transcript counts using **relative2abs.**  This often leads to much more accurate results than using **tobit**().  See 2.4 for details.
 
#### Working with large data sets

Some single-cell RNA-Seq experiments report measurements from tens of thousands of cells or more.  As instrumentation improves and costs drop, experiments will become ever larger and more complex, with many conditions, controls,and replicates.  A matrix of expression data with 50,000 cells and a measurement for each of the 25,000+ genes in the human  genome  can  take  up  a  lot  of  memory.  However,  because  current  protocols  typically  don’t  capture  all  or  even most of the mRNA molecules in each cell, many of the entries of expression matrices are zero.  Using sparse matrices can help you work with huge datasets on a typical computer.  We generally recommend the use of **sparseMatrices** for most users, as it speeds up many computations even for more modestly sized datasets.  To  work  with  your  data  in  a  sparse  format,  simply  provide  it  to  Monocle  as  a  sparse  matrix  from  the Matrix package:



```r
 # NOT RUN
 HSMM <- newCellDataSet(as(as.matrix(rpc_matrix), "sparseMatrix"),
                        phenoData = pd,
                        featureData = fd,
                        lowerDetectionLimit=0.5,
                        expressionFamily=negbinomial.size())
```

**NOTE:** The output from a number of RNA-Seq pipelines, including CellRanger, is already in a sparseMatrix format(e.g.   MTX).  If  so,  you  should  just  pass  it  directly  to  newCellDataSet  without  first  converting  it  to  a  dense  matrix(via(as.matrix)),  because  that  may  exceed  your  available  memory.   If  you  have  10X  Genomics  data  and  are  using cellrangerRkit, you can use it to load your data and then pass that to Monocle as follows:


```r
 # NOT RUN
  cellranger_pipestance_path <- "/path/to/your/pipeline/output/directory"
  gbm <- load_cellranger_matrix(cellranger_pipestance_path)
  gbm_cds <- newCellDataSet(exprs(gbm),
                            phenoData = pData(gbm),
                            featureData = fData(gbm),
                            lowerDetectionLimit=0.5,
                            expressionFamily=negbinomial.size())
```

Monocle’s sparse matrix support is provided by thei *Matrix* package.  Other sparse matrix packages,  such as *slam* or *SparseM* are not supported.

#### Converting relative expression values into mRNA counts

If you performed your single-cell RNA-Seq experiment using spike-in standards, you can convert these measurements into mRNAs per cell (RPC). RPC values are often easier to analyze than FPKM or TPM values, because have better statistical  tools  to  model  them.   In  fact,  it’s  possible  to  convert  FPKM  or  TPM  values  to  RPC  values  even  if  there were no spike-in standards included in the experiment.  Monocle 2 includes an algorithm called *Census* which performs this conversion (Qiu et al, submitted).  You can convert to RPC values before creating your CellDataSet object using the *relative2abs* function, as follows:


```r
 # NO RUN
 pd <- new("AnnotatedDataFrame", data = HSMM_sample_sheet)
 fd <- new("AnnotatedDataFrame", data = HSMM_gene_annotation)

 # First create a CellDataSet from the relative expression levels
 HSMM <- newCellDataSet(as.matrix(HSMM_expr_matrix),
                        phenoData = pd,
                        featureData = fd,
                        lowerDetectionLimit=0.1,
                        expressionFamily=tobit(Lower=0.1))

 # Next, use it to estimate RNA counts
 rpc_matrix <- relative2abs(HSMM)

 # Now, make a new CellDataSet using the RNA counts
 HSMM <- newCellDataSet(as(as.matrix(rpc_matrix), "sparseMatrix"),
                        phenoData = pd,
                        featureData = fd,
                        lowerDetectionLimit=0.5,
                        expressionFamily=negbinomial.size())
```

Note that since we are using RPC values, we have changed the value of lowerDetectionLimit to reflect the new scale of expression.  Importantly, we have also set the expressionFamily to negbinomial(), which tells Monocle touse the negative binomial distribution in certain downstream statistical tests.  Failing to change these two options can create problems later on, so make sure not to forget them when using RPC values. Finally,  we’ll  also  call  two  functions  that  pre-calculate  some  information  about  the  data.   Size  factors  help  us normalize  for  differences  in  mRNA  recovered  across  cells,  and “dispersion” values  will  help  us  perform  differential expression analysis later.


```r
 HSMM <- estimateSizeFactors(HSMM)
 HSMM <- estimateDispersions(HSMM)
```

We’re now ready to start using the HSMM object in our analysis.


#### Filtering low-quality cells

The  first  step  in  any  single-cell  RNA-Seq  analysis  is  identifying  poor-quality  libraries  from  further  analysis.   Most single-cell workflows will include at least some libraries made from dead cells or empty wells in a plate.  It’s also crucial to remove doublets:  libraries that were made from two or more cells accidentally.  These cells can disrupt downstream steps  like  pseudo time  ordering  or  clustering.  This  section  walks  through  typical  quality  control  steps  that  should  be performed as part of all analyses with Monocle. It  is  often  convenient  to  know  how  many  express  a  particular  gene,  or  how  many  genes  are  expressed  by  a  given cell.  Monocle provides a simple function to compute those statistics:


```r
 HSMM <- detectGenes(HSMM, min_expr = 0.1)
 print(head(fData(HSMM)))
```

```
##                    gene_short_name        biotype num_cells_expressed use_for_ordering
## ENSG00000000003.10          TSPAN6 protein_coding                 184            FALSE
## ENSG00000000005.5             TNMD protein_coding                   0            FALSE
## ENSG00000000419.8             DPM1 protein_coding                 211            FALSE
## ENSG00000000457.8            SCYL3 protein_coding                  18            FALSE
## ENSG00000000460.12        C1orf112 protein_coding                  47             TRUE
## ENSG00000000938.8              FGR protein_coding                   0            FALSE
```

```r
 expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 10))
```

The  vector **expressed_genes** now  holds  the  identifiers  for  genes  expressed  in  at  least  50  cells  of  the  data  set.We  will  use  this  list  later  when  we  put  the  cells  in  order  of  biological  progress.   It  is  also  sometimes  convenient  to exclude genes expressed in few if any cells from thenCellDataSet object so as not to waste CPU time analyzing them for differential expression.

Let’s  start  trying  to  remove  unwanted,  poor  quality  libraries  from  the  CellDataSet.   Your  single  cell  RNA-Seq protocol may have given you the opportunity to image individual cells after capture but prior to lysis.  This image data allows you to score your cells, confirming that none of your libraries were made from empty wells or wells with excess cell debris.  With some protocols and instruments, you may get more than one cell captured instead just a single cell.  You  should  exclude  libraries  that  you  believe  did  not  come  from  a  single  cell,  if  possible.   Empty  well  or  debris  well libraries can be especially problematic for Monocle.  It’s also a good idea to check that each cell’s RNA-seq library was sequenced to an acceptible degree.  While there is no widely accepted minimum level for what constitutes seequencing “deeply enough”, use your judgement:  a cell sequenced with only a few thousand reads is unlikely to yield meaningful measurements.

CellDataSet objects provide a convenient place to store per-cell scoring data:  the phenoData slot.  Simply include scoring  attributes  as  columns  in  the  data  frame  you  used  to  create  your CellDataSet container.  You  can  then  easily filter  out  cells  that  don’t  pass  quality  control.   You  might  also  filter  cells  based  on  metrics  from  high  throughput sequencing quality assessment packages such as FastQC. Such tools can often identify RNA-Seq libraries made from heavily degraded RNA, or where the library contains an abnormally large amount of ribosomal, mitochondrial, or other RNA type that you might not be interested in.

The HSMM dataset included with this package has scoring columns built in:


```r
 print(head(pData(HSMM)))
```

```
##                Library Well Hours Media Mapped.Fragments Pseudotime State Size_Factor
## T0_CT_A01 SCC10013_A01  A01     0    GM          1958074  23.916673     1    1.392811
## T0_CT_A03 SCC10013_A03  A03     0    GM          1930722   9.022265     1    1.311607
## T0_CT_A05 SCC10013_A05  A05     0    GM          1452623   7.546608     1    1.218922
## T0_CT_A06 SCC10013_A06  A06     0    GM          2566325  21.463948     1    1.013981
## T0_CT_A07 SCC10013_A07  A07     0    GM          2383438  11.299806     1    1.085580
## T0_CT_A08 SCC10013_A08  A08     0    GM          1472238  67.436042     2    1.099878
##           num_genes_expressed
## T0_CT_A01                6850
## T0_CT_A03                6947
## T0_CT_A05                7019
## T0_CT_A06                5560
## T0_CT_A07                5998
## T0_CT_A08                6055
```

**This dataset has already been filtered using the following commands:**


```r
 # NOT RUN
   valid_cells <- row.names(subset(pData(HSMM), 
          Cells.in.Well == 1 & 
          Control == FALSE & 
          Clump == FALSE & 
          Debris == FALSE & 
          Mapped.Fragments > 1000000))
   HSMM <- HSMM[,valid_cells]
```

If you are using RPC values to measure expression, as we are in this vignette, it’s also good to look at the distribution of mRNA totals across the cells:


```r
 pData(HSMM)$Total_mRNAs <- Matrix::colSums(exprs(HSMM))
 HSMM <- HSMM[,pData(HSMM)$Total_mRNAs < 1e6]
 upper_bound <- 10^(mean(log10(pData(HSMM)$Total_mRNAs)) + 2*sd(log10(pData(HSMM)$Total_mRNAs)))
 lower_bound <- 10^(mean(log10(pData(HSMM)$Total_mRNAs)) - 2*sd(log10(pData(HSMM)$Total_mRNAs)))
 qplot(Total_mRNAs, data=pData(HSMM), color=Hours, geom="density") +
 geom_vline(xintercept=lower_bound) +
 geom_vline(xintercept=upper_bound)
```

![plot of chunk lkmRNATotal](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M3A/lkmRNATotal-1.png)

```r
 HSMM <- HSMM[,pData(HSMM)$Total_mRNAs > lower_bound &
               pData(HSMM)$Total_mRNAs < upper_bound]
 HSMM <- detectGenes(HSMM, min_expr = 0.1)
```

We’ve  gone  ahead  and  removed  the  few  cells  with  either  very  low  mRNA  recovery  or  far  more  mRNA  that  the typical cell.  Often, doublets or triplets have roughly twice the mRNA recovered as true single cells, so the latter filter is another means of excluding all but single cells from the analysis.  Such filtering is handy if your protocol doesn’t allow direct visualization of cell after they’ve been captured.  Note that these thresholds of 10,000 and 40,000 mRNAs are specific to this dataset.  You may need to adjust filter thresholds based on your experimental protocol.

Once you’ve excluded cells that do not pass your quality control filters, you should verify that the expression values stored in your CellDataSet follow a distribution that is roughly lognormal:


```r
 # Log-transform each value in the expression matrix.
 L <- log(exprs(HSMM[expressed_genes,]))

 # Standardize each gene, so that they are all on the same scale,
 # Then melt the data with plyr so we can plot it easily"
 melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))

 # Plot the distribution of the standardized gene expression values.
 qplot(value, geom="density", data=melted_dens_df) +  
 stat_function(fun = dnorm, size=0.5, color='red') +
 xlab("Standardized log(FPKM)") +
 ylab("Density")
```

![plot of chunk lkLogNormal](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M3A/lkLogNormal-1.png)

### Classifying and counting cells of different types

Single  cell  experiments  are  often  performed  on  complex  mixtures  of  multiple  cell  types.   Dissociated  tissue  samples might  contain  two,  three,  or  even  many  different  cells  types.  In  such  cases,  it’s  often  nice  to  classify  cells  based  on type  using  known  markers.  In  the  myoblast  experiment,  the  culture  contains  fibroblasts  that  came  from  the  original muscle  biopsy  used  to  establish  the  primary  cell  culture.   Myoblasts  express  some  key  genes  that  fibroblasts  don’t. Selecting only the genes that express, for example, sufficiently high levels ofMYF5 excludes the fibroblasts.  Likewise, fibroblasts express high levels of ANPEP(CD13), while myoblasts tend to express few if any transcripts of this gene.

#### Classifying cells with CellTypeHierarchy

Monocle  provides  a  simple  system  for  tagging  cells  based  on  the  expression  of  marker  genes  of  your  choosing.   You simply provide a set of functions that Monocle can use to annotate each cell.  For example, you could provide a function for each of several cell types.  These functions accept as input the expression data for each cell, and return TRUE to tell Monocle that a cell meets the criteria defined by the function.  So you could have one function that returns TRUE for cells that express myoblast-specific genes, another function for fibroblast-specific genes, etc.  Here’s an example of such a set of “gating” functions:


```r
 MYF5_id <- row.names(subset(fData(HSMM), gene_short_name == "MYF5"))
 ANPEP_id <- row.names(subset(fData(HSMM), gene_short_name == "ANPEP"))

 cth <- newCellTypeHierarchy()
 cth <- addCellType(cth, "Myoblast", classify_func=function(x) {x[MYF5_id,] >= 1})
 cth <- addCellType(cth, "Fibroblast", classify_func=function(x)
        {x[MYF5_id,] < 1 & x[ANPEP_id,] > 1})
```

The functions are organized into a small data structure called a CellTypeHierarchy, that Monocle uses to classify the  cells.  You  first  initialize  a  new  CellTypeHierarchy  object,  then  register  your  gating  functions  within  it.  Once  the data structure is set up, you can use it to classify all the cells in the experiment:


```r
 HSMM <- classifyCells(HSMM, cth, 0.1)
```

The function **classifyCells** applies each gating function to each cell, classifies the cells according to the gating functions,  and  returns  the CellDataSetwith  a  new  column, CellType in  its pData table.  We  can  now  count  how many cells of each type there are in the experiment.


```r
 table(pData(HSMM)$CellType)
```

```
## 
## Fibroblast   Myoblast    Unknown 
##         56         85        121
```

```r
 pie <- ggplot(pData(HSMM), aes(x = factor(1), fill = factor(CellType))) +
     geom_bar(width = 1)
 
 pie + coord_polar(theta = "y") +
   theme(axis.title.x=element_blank(), axis.title.y=element_blank())
```

![plot of chunk lkCellType](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M3A/lkCellType-1.png)

Note  that  many  cells  are  marked “Unknown”.  This  is  common,  largely  because  of  the  low  rate  of  mRNA  capture in most single-cell RNA-Seq experiments.  A cell might express a few MYF5 mRNAs, but we weren’t lucky enough to capture one of them.  When a cell doesn’t meet any of the criteria specified in your classification functions, it’s marked “Unknown”.  If  it  meets  multiple  functions’  criteria,  it’s  marked “Ambiguous”.  You  could  just  exclude  such  cells,  but you’d be throwing out a lot of your data.  In this case, we’d lose more than half of the cells!

#### Unsupervised cell clustering
 
Monocle provides an algorithm you can use to impute the types of the “Unknown” cells.  This algorithm, implemented in  the  function **clusterCells**,  groups  cells  together  according  to  global  expression  profile.   That  way,  if  your  cell expresses lots of genes specific to myoblasts, but just happens to lack MYF5, we can still recognize it as a myoblast.  **clusterCells** can be used in an unsupervised manner, as well as in a “semi-supervised” mode, which allows to assist the algorithm with some expert knowledge.  Let’s look at the unsupervised mode first.

The first step is to decide which genes to use in clustering the cells.  We could use all genes, but we’d be including a lot of genes that are not expressed at a high enough level to provide a meaningful signal.  Including them would just add noise to the system.  We can filter genes based on average expression level, and we can additionally select genes that are unusually variable across cells.  These genes tend to be highly informative about cell state.


```r
 disp_table <- dispersionTable(HSMM)
 unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
 HSMM <- setOrderingFilter(HSMM, unsup_clustering_genes$gene_id)
 plot_ordering_genes(HSMM)
```

![plot of chunk unsup_clust](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M3A/unsup_clust-1.png)

The **setOrderingFilter** function marks genes that will be used for clustering in subsequent calls to **clusterCells**,although we will be able to provide other lists of genes if we want.  **plot_ordering_genes** shows how variability (dispersion) in a gene’s expression depends on the average expression across cells.  The red line shows Monocle’s expectation of the dispersion based on this relationship.  The genes we marked for use in clustering are shown as black dots, while the others are shown as grey dots.


Now we’re ready to try clustering the cells:


```r
 # HSMM@auxClusteringData[["tSNE"]]£variance_explained <- NULL
 plot_pc_variance_explained(HSMM, return_all = F)#norm_method ='log',
```

![plot of chunk plot_pc_var_exp](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M3A/plot_pc_var_exp-1.png)

```r
 HSMM <- reduceDimension(HSMM, max_components=2, num_dim = 6,
                         reduction_method ='tSNE', verbose = T)

 ## Remove noise by PCA ...
 ## Reduce dimension by tSNE ...

 HSMM <- clusterCells(HSMM,
                      num_clusters=2)
```

```
## Distance cutoff calculated to 0.9229564
```

```r
 ## Distance cutoff calculated to 1.243471
 plot_cell_clusters(HSMM, 1, 2, color="CellType", markers=c("MYF5", "ANPEP"))
```

![plot of chunk plot_pc_var_exp](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M3A/plot_pc_var_exp-2.png)

Monocle  uses  t-SNE  [3]  to  cluster  cells,  using  an  approach  that’s  very  similar  to  and  inspired  by  Rahul  Satija’s excellent Seurat package, which itself was inspired by viSNE from Dana Pe’er’s lab.

The cells tagged as myoblasts by our gating functions are marked in green, while the fibroblasts are tagged in red.  The cells that don’t express either marker are blue.  In many experiments,  cells of different types are clearly separate from one another.  Unfortunately, in this experiment, the cells don’t simply cluster by type - there’s not a clear space between the green cells and the red cells.  This isn’t all that surprising, because myoblasts and contaminating interstitial fibroblasts  express  many  of  the  same  genes  in  these  culture  conditions,  and  there  are  multiple  culture  conditions  in the  experiment.  That  is,  there  are  other  sources  of  variation  in  the  experiment  that  might  be  driving  the  clustering.  One source of variation in the experiment stems from the experimental design.  To initiate myoblast differentiation, we switch media from a high-mitogen growth medium (GM) to a low-mitogen differentiation medium (DM). Perhaps the cells are clustering based on the media they’re cultured in?


```r
 plot_cell_clusters(HSMM, 1, 2, color="Media")
```

![plot of chunk plot_cell_clust](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M3A/plot_cell_clust-1.png)

Monocle  allows  us  to  subtract  the  effects  of “uninteresting” sources  of  variation  to  reduce  their  impact  on  the clustering.    You  can  do  this  with  the **residualModelFormulaStr** argument  to **clusterCells** and  several  other Monocle functions.  This argument accepts an R model formula string specifying the effects you want to subtract prior to clustering.


```r
 HSMM <- reduceDimension(HSMM, max_components=2, num_dim = 2, reduction_method ='tSNE',
         residualModelFormulaStr="~Media + num_genes_expressed", verbose = T)

 ## Removing batch effects
 ## Remove noise by PCA ...
 ## Reduce dimension by tSNE ...

 HSMM <- clusterCells(HSMM, num_clusters=2)
```

```
## Distance cutoff calculated to 0.9967294
```

```r
 ## Distance cutoff calculated to 1.277005`

 plot_cell_clusters(HSMM, 1, 2, color="CellType")
```

![plot of chunk RUV](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M3A/RUV-1.png)

Now that we’ve accounted for some unwanted sources of variation, we’re ready to take another crack at classifying the cells by unsupervised clustering:


```r
 HSMM <- clusterCells(HSMM, num_clusters=2)
```

```
## Distance cutoff calculated to 0.9967294
```

```r
 ## Distance cutoff calculated to 1.277005

 plot_cell_clusters(HSMM, 1, 2, color="Cluster") + facet_wrap(~CellType)
```

![plot of chunk reClust](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M3A/reClust-1.png)

**In the vignette example** Now, most of the myoblasts are in one cluster, most of the fibroblasts are in the other, and the unknowns are spread across both.  However, we still see some cells of both types in each cluster.  This could be due to lack of specificity in our marker genes and our **CellTypeHierarchy** functions, but it could also be due to suboptimal clustering.  To help rule out the latter, let’s try running **clusterCells** in its semi-supervised mode.


#### Semi-supervised cell clustering with known marker genes

First, we’ll select a different set of genes to use for clustering the cells.  Before we just picked genes that were highly expressed  and  highly  variable.   Now,  we’ll  pick  genes  that co-vary with  our  markers.   In  a  sense,  we’ll  be  building  a large list of genes to use as markers, so that even if a cell doesn’t have MYF5, it might be recognizable as a myoblast based on other genes.


```r
 Dummy <- 1
 marker_diff <- markerDiffTable(HSMM[expressed_genes,],
                                cth,
                                residualModelFormulaStr="~Media + num_genes_expressed",
                                cores=1)
```

The function **markerDiffTable** takes a *CellDataSet*  and a *CellTypeHierarchy*  and classifies all the cells into types according to your provided functions.  It then removes all the “Unknown” and “Ambiguous” functions before identifying genes that are differentially expressed between the types.  Again, you can provide a residual model of effects to exclude from  this  test.   The  function  then  returns  a  data  frame  of  test  results,  and  you  can  use  this  to  pick  the  genes  you want  to  use  for  clustering.   Often  it’s  best  to  pick  the  top  10  or  20  genes  that  are  most  specific  for  each  cell  type.  This ensures that the clustering genes aren’t dominated by markers for one cell type.  You generally want a balanced panel of markers for each type if possible.  Monocle provides a handy function for ranking genes by how restricted their expression is for each type.


```r
 candidate_clustering_genes <- row.names(subset(marker_diff, qval < 0.01))
 marker_spec <- calculateMarkerSpecificity(HSMM[candidate_clustering_genes,], cth)
 head(selectTopMarkers(marker_spec, 3))
```

```
##              gene_id   CellType specificity
## 1 ENSG00000019991.11 Fibroblast   0.9894231
## 2 ENSG00000128340.10 Fibroblast   0.9999611
## 3  ENSG00000163710.3 Fibroblast   0.9718761
## 4  ENSG00000156298.8   Myoblast   0.9884418
## 5  ENSG00000233494.1   Myoblast   1.0000000
## 6  ENSG00000270123.1   Myoblast   1.0000000
```

The  last  line  above  shows  the  top  three  marker  genes  for  myoblasts  and  fibroblasts.   The  ”specificity” score  is calculated  using  the  metric  described  in  Cabili  et  al  [?]  and  can  range  from  zero  to  one.  The  closer  it  is  to  one,  the more restricted it is to the cell type in question.  You can use this feature to define new markers for known cell types,or  pick  out  genes  you  can  use  in  purifying  newly  discovered  cell  types.   This  can  be  highly  valuable  for  downstream follow up experiments.

To cluster the cells, we’ll choose the top 500 markers for each of these cell types:

```{ semiSupClust}
 semisup_clustering_genes <- unique(selectTopMarkers(marker_spec, 500)$gene_id)
 HSMM <- setOrderingFilter(HSMM, semisup_clustering_genes)
 plot_ordering_genes(HSMM)
```

Note that we’ve got smaller set of genes, and some of them are not especially highly expressed or variable across the  experiment.   However,  they  are  great  for  distinguishing  cells  that  express MYF5 from  those  that  have ANPEP.  We’ve  already  marked  them  for  use  in  clustering,  but  even  if  we  hadn’t,  we  could  still  use  them  by  providing  them directly to clusterCells.


```r
 plot_pc_variance_explained(HSMM, return_all = F)
```

![plot of chunk plot_pc_variance](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M3A/plot_pc_variance-1.png)


```r
 HSMM <- reduceDimension(HSMM, max_components=2, num_dim = 2, reduction_method = 'tSNE',
                         residualModelFormulaStr="~Media + num_genes_expressed", verbose = T)
 HSMM <- clusterCells(HSMM, num_clusters=2)
```

```
## Distance cutoff calculated to 0.9967294
```

```r
 plot_cell_clusters(HSMM, 1, 2, color="CellType")
```

![plot of chunk RUV2](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M3A/RUV2-1.png)

#### Imputing cell type

Note  that  we’ve  reduce  the  number  of “contaminating” fibroblasts  in  the  myoblast  cluster,  and  vice  versa.  But  what about  the “Unknown” cells?   If  you  provide clusterCells with  a  the  CellTypeHierarcy,  Monocle  will  use  it  classify whole clusters, rather than just individual cells.  Essentially, cluserCells works exactly as before, except after the clusters are built, it counts the frequency of each cell type in each cluster.  When a cluster is composed of more than a certain percentage (in this case, 10%) of a certain type, all the cells in the cluster are set to that type.  If a cluster is composed of more than one cell type, the whole thing is marked “Ambiguous”.  If there’s no cell type thats above the threshold, the cluster is marked “Unknown”.  Thus, Monocle helps you impute the type of each cell even in the presence of missing marker data.


```r
 HSMM <- clusterCells(HSMM,
                      num_clusters=2,
                      frequency_thresh=0.1,
                      cell_type_hierarchy=cth)
```

```
## Distance cutoff calculated to 0.9967294
```

```r
 plot_cell_clusters(HSMM, 1, 2, color="CellType", markers = c("MYF5", "ANPEP"))
```

![plot of chunk imputeCellType](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M3A/imputeCellType-1.png)

As you can see, the clusters are fairly pure in terms of MYF5 expression.  There are some cells expressing ANPEP in both clusters, but those in the myoblast cluster also express MYF5.  This is not surprising, as ANPEP isn’t a very specific marker of fibroblasts.  Overall, we’ve successfully classified all the cells:


```r
 pie <- ggplot(pData(HSMM), aes(x = factor(1), fill = factor(CellType))) +
    geom_bar(width = 1)
 pie + coord_polar(theta = "y") +
    theme(axis.title.x=element_blank(), axis.title.y=element_blank())
```

![plot of chunk pie2](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M3A/pie2-1.png)

Finally, we subset the CellDataSet object to create HSMM_myo, which includes only myoblasts.  We’ll use this in the rest of the analysis.



```r
 ### cell classification did not work as in vignette!!!
 HSMM_myo <- HSMM[,pData(HSMM)$CellType == "Myoblast"]
 HSMM_myo <- estimateDispersions(HSMM_myo)
```


### Constructing single cell trajectories

During development, in response to stimuli, and througout life, cells transition from one functional “state” to another. Cells in different states express different sets of genes, producing a dynamic repetoire of proteins and metabolites that carry out their work. As cells move between states, undergo a process of transcriptional re-configuration, with some genes being silenced and others newly activated. These transient states are often hard to characterize because purifying cells in between more stable endpoint states can be difficult or impossible. Single-cell RNA-Seq can enable you to see these states without the need for purification. However, to do so, we must determine where each cell is the range of possible states.

Monocle introduced the strategy of using RNA-Seq for *single cell trajectory analysis*. Rather than purifying cells into discrete states experimentally, Monocle uses an algorithm to learn the sequence of gene expression changes each cell must go through as part of a dynamic biological process. Once it has learned the overall “trajectory” of gene expression changes, Monocle can place each cell at its proper position in the trajectory. You can then use Monocle’s differential analysis toolkit to find genes regulated over the course of the trajectory, as described in section 5.3. If there are multiple outcome for the process, Monocle will reconstruct a “branched” trajectory. These branches correspond to cellular“decisions”, and Monocle provides powerful tools for identifying the genes affected by them and involved in making them. You can see how to analyze branches in section 8.6. Monocle relies on a machine learning technique called reversed graph embedding to construct single-cell trajectories. You can read more about the theoretical foundations of Monocle’s approach in section 8, or consult the references shown below in section ??.


#### “Pseudotime”: a measure of progress through a biological process

In many biological processes, cells do not progress in perfect synchrony. In single-cell expression studies of processes such as cell differentiation, captured cells might be widely distributed in terms of progress. That is, in a population of cells captured at exactly the same time, some cells might be far along, while others might not yet even have begun the process. This asynchrony creates major problems when you want to understand the sequence of regulatory changes that occur as cells transition from one state to the next. Tracking the expression across cells captured at the same time produces a very compressed sense of a gene’s kinetics, and the apparent variability of that gene’s expression will be very high.

By ordering each cell according to its progress along a learned trajectory, Monocle alleviates the problems that arise due to asynchrony. Instead of tracking changes in expression as a function of time, Monocle tracks changes as a function of progress along the trajectory, which we term“pseudotime”. Pseudotime is an abstract unit of progress: it’s simply the distance between a cell and the start of the trajectory, measured along the shortest path. The trajectory’s total length is defined in terms of the total amount of transcriptional change that a cell undergoes as it moves from the starting state to the end state. For further details, see section 8.

#### The ordering algorithm
##### Choosing genes for ordering
Inferring a single-cell trajectory is a hard machine learning problem. The first step is to select the genes Monocle will use as input for its machine learning approach. This is called feature selection, and it has a major impact in the shape of the trajectory. In single-cell RNA-Seq, genes expressed at low levels are often very noisy, but some may contain important information regarding the state of the cell. If we simply provide all the input data, Monocle might get confused by the noise, or fix on a feature of the data that isn’t biologically meaningful, such as batch effects arising from collecting data on different days. Monocle provides you with a variety of tools to select genes that will yield a robust, accurate, and biologically meaningful trajectory. You can use these tools to either perform a completely“unsupervised”analysis, in which Monocle has no forehand knowledge of which gene you consider important. Alternatively, you can make use of expert knowledge in the form of genes that are already known to define biolgical progress to shape Monocle’s trajectory. We consider this mode “semi-supervised”, because Monocle will augment the markers you provide with other, related genes. Monocle then uses these genes to produce trajectories consistent with known biology but that often reveal new regulatory structure. We return to the muscle data to illustrate both of these modes.

##### Reducing the dimensionality of the data

Once we have selected the genes we will use to order the cells, Monocle applies a dimensionality reduction to the data, which will drastically improves the quality of the trajectory. Monocle reduces the dimensionality of the data with the reduceDimension function. This function has a number of options, so you should familiarize yourself with them by consulting the its manual page. You can choose from two algorithms in reduceDimension. The first, termed Independent Component Analysis, is a classic linear technique for decomposing data that powered the original version of Monocle. The second, called DDRTree, is a much more powerful nonlinear technique that is the default for Monocle 2. For more on how these both work, see section 8.

##### Ordering the cells in pseudotime

With the expression data projected into a lower dimensional space, Monocle is ready to learn the trajectory that describes how cells transition from one state into another. Monocle assumes that the trajectory has a tree structure, with one end of it the “root”, and the others the “leaves”. A cell at the beginning of the biological process starts at the root and progresses along the trunk until it reaches the first branch, if there is one. It then chooses a branch, and moves further and further along the tree until it reaches a leaf. These mathematical assumptions translate into some important biological ones. First, that the data includes all the major stages of the biological process. If your experiment failed to capture any cells at a key developmental transition, Monocle won’t know its there. Second, that gene expression changes are smooth as a cell moves from one stage to the next. This assumption is realistic: major discontinuities in the trajectory would amount to a cell almost instantaneously turning over its transcriptome, which probably doesn’t happen in most biolgical processes.

To order your cells, Monocle uses the orderCells function. This routine has one important argument, which allows you to set the root of the tree, and thus the beginning of the process. See the manual page for orderCells for more details.

### Unsupervised ordering
In this section, we discuss ordering cells in a completely unsupervised fashion. First, we must decide which genes we will use to define a cell’s progress through myogenesis. Monocle orders cells by examining the pattern of expression of these genes across the cell population. Monocle looks for genes that vary in“interesting”(i.e. not just noisy) ways, and uses these to structure the data. We ultimately want a set of genes that increase (or decrease) in expression as a function of progress through the process we’re studying.
Ideally, we’d like to use as little prior knowledge of the biology of the system under study as possible. We’d like to discover the important ordering genes from the data, rather than relying on literature and textbooks, because that might introduce bias in the ordering. One effective way to isolate a set of ordering genes is to simply compare the cells collected at the beginning of the process to those at the end and find the differentially expressed genes, as described above. The command below will find all genes that are differentially expressed in response to the switch from growth medium to differentiation medium:


```r
 Dummy <- 1

 diff_test_res <- differentialGeneTest(HSMM_myo[expressed_genes,],
                                      fullModelFormulaStr="~Media")
 ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
```

Choosing genes based on differential analysis of time points is often highly effective, but what if we don’t have time series data? If the cells are asynchronously moving through our biological process (as is usually the case), Monocle can often reconstruct their trajectory from a single population captured all at the same time. Below are two methods to select genes that require no knowledge of the design of the experiment at all.

#### Selecting genes with high dispersion across cells

Genes that vary a lot are often highly informative for identifying cell subpopulations or ordering cells along a trajectory. In RNA-Seq, a gene’s variance typically depends on its mean, so we have to be a bit careful about how we select genes based on their variance.


```r
disp_table <- dispersionTable(HSMM_myo)
ordering_genes <- subset(disp_table,
                         mean_expression >= 0.5 &
                         dispersion_empirical >= 1 * dispersion_fit)$gene_id
```

Once we have a list of gene ids to be used for ordering, we need to set them in the HSMM object, because the next several functions will depend on them.


```r
HSMM_myo <- setOrderingFilter(HSMM_myo, ordering_genes)
plot_ordering_genes(HSMM_myo)
```

![plot of chunk setOrderingGenes](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M3A/setOrderingGenes-1.png)

The genes we’ve chosen to use for ordering define the state space of the cells in our data set. Each cell is a point in this space, which has dimensionality equal to the number of genes we’ve chosen. So if there are 500 genes used forordering, each cell is a point in a 500-dimensional space. For a number of reasons, Monocle works better if we can reduce the dimensionality of that space before we try to put the cells in order. In this case, we will reduce the space down to one with two dimensions, which we will be able to easily visualize and interpret while Monocle is ordering the cells.


```r
 HSMM_myo <- reduceDimension(HSMM_myo, max_components=2)
 HSMM_myo <- orderCells(HSMM_myo)
```

Once the cells are ordered, we can visualize the trajectory in the reduced dimensional space.


```r
 plot_cell_trajectory(HSMM_myo, color_by="Hours")
```

![plot of chunk unnamed-chunk-1](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M3A/unnamed-chunk-1-1.png)

The trajectory has a tree-like structure. We can see that the cells collected at time zero are located near one of the tips of the tree, while the others are distributed amongst the two“branches”. Monocle doesn’t know a priori which of the trajectory of the tree to call the“beginning”, so we often have to call orderCells again using the root_state argument to specify the beginning. First, we plot the trajectory, this time coloring the cells by“State”:


```r
plot_cell_trajectory(HSMM_myo, color_by="State")
```

![plot of chunk unnamed-chunk-2](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M3A/unnamed-chunk-2-1.png)

“State” is just Monocle’s term for the segment of the tree. The function below is handy for identifying the State which contains most of the cells from time zero. We can then pass that to orderCells:


```r
GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$Hours)[,"0"]
return(as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))])) }else { return (1) }
}
HSMM_myo <- orderCells(HSMM_myo, root_state=GM_state(HSMM_myo))
plot_cell_trajectory(HSMM_myo, color_by="Pseudotime")
```

![plot of chunk unnamed-chunk-3](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M3A/unnamed-chunk-3-1.png)
If there are a ton of states in your tree, it can be a little hard to make out where each one falls on the tree. Sometimes it can be handy to“facet”the trajectory plot so it’s easier to see where each of the states are located:


```r
plot_cell_trajectory(HSMM_myo, color_by="State") + facet_wrap(~State, nrow=1)
```

![plot of chunk unnamed-chunk-4](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M3A/unnamed-chunk-4-1.png)

And if you don’t have a timeseries, you might need to set the root based on where certain marker genes are expressed, using your biological knowledge of the system. For example, in this experiment, a highly proliferative population of progenitor cells are generating two types of post-mitotic cells. So the root should have cells that express high levels of proliferation markers. We can use the jitter plot to pick figure out which state corresponds to rapid proliferation:


```r
 blast_genes <- row.names(subset(fData(HSMM_myo), gene_short_name %in% c("CCNB2", "MYOD1", "MYOG")))
 plot_genes_jitter(HSMM_myo[blast_genes,], grouping="State", min_expr=0.1)
```

![plot of chunk unnamed-chunk-5](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M3A/unnamed-chunk-5-1.png)

To confirm that the ordering is correct we can select a couple of markers of myogenic progress. Plotting these genes demonstrates that ordering looks good:


```r
 HSMM_expressed_genes <- row.names(subset(fData(HSMM_myo), num_cells_expressed >= 10))
 HSMM_filtered <- HSMM_myo[HSMM_expressed_genes,]

 my_genes <- row.names(subset(fData(HSMM_filtered),
                              gene_short_name %in% c("CDK1", "MEF2C", "MYH3")))
cds_subset <- HSMM_filtered[my_genes,]
plot_genes_in_pseudotime(cds_subset, color_by="Hours")
```

![plot of chunk unnamed-chunk-6](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M3A/unnamed-chunk-6-1.png)

**Note**: in the event that Monocle produces a linear trajectory, there will only be one state. In this case, you can tell it which end is the beginning of the trajectory by passing orderCells an optional argument argument: the reverse flag. The reverse flag tells Monocle to reverse the orientation of the entire process as it’s being discovered from the data, so that the cells that would have been assigned to the end are instead assigned to the beginning, and so on.


#### Selecting genes based on PCA loading

A number of single-cell clustering studies have found that principal component analysis (PCA) is an effective way of finding genes that vary widely across cells. You can use PCA to pick out a good set of ordering genes in a completely unsupervised manner. The procedure below first normalizes the expression data by adding a pseudocount, log-transforming it, and then standardizing it so that all genes have roughly the same dynamic range. Then, it uses the irlba package to perform the PCA. The code below is a little more complicated than just calling prcomp(), but it has the advantage of performing the entire computation with sparse matrix operations. This is important for large, sparse CellDataSet objects like the ones common in droplet-based single-cell RNA-Seq.


```r
HSMM_myo <- HSMM_myo[HSMM_expressed_genes,]
exprs_filtered <- t(t(exprs(HSMM_myo)/pData(HSMM_myo)$Size_Factor)) #nz_genes <- which(exprs_filtered != 0, arr.ind = T) exprs_filtered@x <- log(exprs_filtered@x + 1)
# Calculate the variance across genes without converting to a dense
# matrix:
expression_means <- Matrix::rowMeans(exprs_filtered)
expression_vars <- Matrix::rowMeans((exprs_filtered - expression_means)^2)
# Filter out genes that are constant across all cells:
genes_to_keep <- expression_vars > 0
exprs_filtered <- exprs_filtered[genes_to_keep,]
expression_means <- expression_means[genes_to_keep]

expression_vars <- expression_vars[genes_to_keep]
# Here's how to take the top PCA loading genes, but using
# sparseMatrix operations the whole time, using irlba. Note # that the v matrix from irlba is the loading matrix set.seed(0)
irlba_pca_res <- irlba(t(exprs_filtered),
                      nu=0,
                      center=expression_means,
                      scale=sqrt(expression_vars),
                      right_only=TRUE)$v
row.names(irlba_pca_res) <- row.names(exprs_filtered)
# Here, we will just
# take the top 200 genes from components 2 and 3.
# Component 1 usually is driven by technical noise.
# We could also use a more principled approach,
# similar to what dpFeature does below
PC2_genes <- names(sort(abs(irlba_pca_res[, 2]), decreasing = T))[1:200] 
PC3_genes <- names(sort(abs(irlba_pca_res[, 3]), decreasing = T))[1:200]
ordering_genes <- union(PC2_genes, PC3_genes)
```

Using these to order the cells as above yields the following trajectory:


```r
HSMM_myo <- setOrderingFilter(HSMM_myo, ordering_genes)
HSMM_myo <- reduceDimension(HSMM_myo, max_components=2)
HSMM_myo <- orderCells(HSMM_myo)
HSMM_myo <- orderCells(HSMM_myo, root_state=GM_state(HSMM_myo))
plot_cell_trajectory(HSMM_myo, color_by="Hours")
```

![plot of chunk unnamed-chunk-8](/mnt100/home/Dropbox/SingleCell/Jun2017/R/Scripts/figure/M3A/unnamed-chunk-8-1.png)


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
##  [1] bindrcpp_0.2           HSMMSingleCell_0.110.0 reshape2_1.4.2         BiocParallel_1.10.1   
##  [5] monocle_2.4.0          DDRTree_0.1.5          irlba_2.2.1            VGAM_1.0-3            
##  [9] ggplot2_2.2.1          Biobase_2.36.2         BiocGenerics_0.22.0    Matrix_1.2-10         
## [13] knitr_1.16             rmarkdown_1.6         
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_0.12.11       bindr_0.1          highr_0.6          compiler_3.4.0     RColorBrewer_1.1-2
##  [6] plyr_1.8.4         tools_3.4.0        densityClust_0.2.1 digest_0.6.12      evaluate_0.10     
## [11] Rtsne_0.13         tibble_1.3.3       gtable_0.2.0       lattice_0.20-35    pkgconfig_2.0.1   
## [16] rlang_0.1.1        igraph_1.0.1       fastICA_1.2-1      dplyr_0.7.0        stringr_1.2.0     
## [21] cluster_2.0.6      combinat_0.0-8     rprojroot_1.2      grid_3.4.0         glue_1.1.0        
## [26] R6_2.2.2           qlcMatrix_0.9.5    pheatmap_1.0.8     limma_3.32.2       magrittr_1.5      
## [31] codetools_0.2-15   matrixStats_0.52.2 backports_1.1.0    scales_0.4.1       htmltools_0.3.6   
## [36] assertthat_0.2.0   colorspace_1.3-2   labeling_0.3       proxy_0.4-17       stringi_1.1.5     
## [41] lazyeval_0.2.0     munsell_0.4.3      slam_0.1-40
```



