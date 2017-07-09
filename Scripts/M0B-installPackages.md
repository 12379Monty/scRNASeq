###   Install R Packages
 
<style type="text/css"> body, td { font-size: 14px; } code.r{ font-size: 12px; } pre { font-size: 12px } </style>


<!-- ***************************************************** -->

### Install additional packages

Install all packages used in
[Hicks](https://github.com/stephaniehicks/scBatchPaper)

#### NOTE: could not install 'mvoutlier', 'destiny' 

Cause 'http://ab-initio.mit.edu/nlopt/nlopt-2.4.2.tar.gz' is non-responsive.
Can download nlopt-2.4.2.tar.gz from 
https://osdn.net/projects/sfnet_daetools/downloads/gnu-linux-libs/nlopt-2.4.2.tar.gz/

To install:
* tar xvfz nlopt-2.4.2.tar.gz
* cd nlopt-2.4.2
* ./configure
* makes
* make install

Then try installing  'mvoutlier' and 'destiny'  again.  **Still Fails!!!**

#### FIXED after ab-initio.mit.edu cma eback online:
* had to "make uninstall" nlopt-2.4.2 first.
* then biocLite(c(mvoutlier', 'destiny'))
* Done!


```r
 RequiredLibs.vec <- unique(
  c('R.utils', 'readxl', 'Rtsne', 'mvoutlier', 'destiny', 'scran', 'BiocParallel',
    'Biobase', 'GenomicRanges', 'vcd', 'stringr', 'reshape2',
    'tidyr', 'plyr', 'dplyr', 'irlba', 'readr', 'scater',
    'knitr', 'ggplot2', 'grid', 'cowplot', 'Biobase', 'GenomicRanges', 
    'stringr', 'reshape2', 'tidyr', 'RColorBrewer', 'plyr', 'dplyr', 
    'irlba', 'rafalib' , 'rafalib', 'SummarizedExperiment', 'MultiAssayExperiment', 
    'scRNASeqMouseShalekBMDCs', 'scRNASeqHCT116CellLinesWu', 'scRNASeqMouseShalekBMDCs', 
    'scRNASeqHCT116CellLinesWu', 'scRNASeqMouseShalekBMDCs', 'scRNASeqHCT116CellLinesWu', 
    'scRNASeqMouseShalekBMDCs', 'scRNASeqHCT116CellLinesWu', 'scRNASeqMouseShalekBMDCs', 
    'scRNASeqHCT116CellLinesWu', 'scRNASeqMouseDengMonoAllelic', 'scRNASeqMouseJaitinSpleen',
    'scRNASeqMouseKumarPSC', 'SummarizedExperiment', 'scRNASeqHumanPatelGlioblastoma', 
    'scRNASeqMouseShalekDendritic', 'scRNASeqMouseTreutleinLineage', 'scRNASeqHumanBosePrinting', 
    'scRNASeqMouseBurnsInnerEar', 'scRNASeqHumanGuoGermCells', 'scRNASeqMouseKowalczykAging', 
    'scRNASeqHumanLengOscillatoryGenes', 'scRNASeqDaniSajitaSeurat', 'scRNASeqMouseZeiselCortex', 
    'scRNASeqHumanPatelGlioblastoma'))


 InstalledLibs.vec <- installed.packages()[,1]
 ToIntall.vec <- setdiff(RequiredLibs.vec,InstalledLibs.vec)


 #ToIntall.vec <- c('nloptr', 'lme4', 'pbkrtest', 'car', 'VIM', 'robCompositions')
 

 if(length(ToIntall.vec)){
  source("http://www.bioconductor.org/biocLite.R")
  biocLite(ToIntall.vec)
 }
```

```
## Bioconductor version 3.5 (BiocInstaller 1.26.0), ?biocLite for help
```

```
## BioC_mirror: https://bioconductor.org
```

```
## Using Bioconductor 3.5 (BiocInstaller 1.26.0), R 3.4.0 (2017-04-21).
```

```
## Installing package(s) 'scRNASeqMouseShalekBMDCs',
##   'scRNASeqHCT116CellLinesWu', 'scRNASeqMouseDengMonoAllelic',
##   'scRNASeqMouseJaitinSpleen', 'scRNASeqMouseKumarPSC',
##   'scRNASeqHumanPatelGlioblastoma', 'scRNASeqMouseShalekDendritic',
##   'scRNASeqMouseTreutleinLineage', 'scRNASeqHumanBosePrinting',
##   'scRNASeqMouseBurnsInnerEar', 'scRNASeqHumanGuoGermCells',
##   'scRNASeqMouseKowalczykAging', 'scRNASeqHumanLengOscillatoryGenes',
##   'scRNASeqDaniSajitaSeurat', 'scRNASeqMouseZeiselCortex'
```

```
## Warning: packages 'scRNASeqMouseShalekBMDCs', 'scRNASeqHCT116CellLinesWu',
## 'scRNASeqMouseDengMonoAllelic', 'scRNASeqMouseJaitinSpleen',
## 'scRNASeqMouseKumarPSC', 'scRNASeqHumanPatelGlioblastoma',
## 'scRNASeqMouseShalekDendritic', 'scRNASeqMouseTreutleinLineage',
## 'scRNASeqHumanBosePrinting', 'scRNASeqMouseBurnsInnerEar',
## 'scRNASeqHumanGuoGermCells', 'scRNASeqMouseKowalczykAging',
## 'scRNASeqHumanLengOscillatoryGenes', 'scRNASeqDaniSajitaSeurat',
## 'scRNASeqMouseZeiselCortex' are not available (for R version 3.4.0)
```

```
## Old packages: 'bindrcpp', 'DBI', 'RSQLite', 'statmod', 'XML'
```


```r
 require(devtools)
 install_github("willtownes/patel2014gliohuman")
 install_github("stephaniehicks/trapnell2014myoblasthuman")
```

### Parameter settings:
  * WRKDIR = /Users/francois/Dropbox/SingleCell/Jun2017/R/
  * FN = M0B-installPackages
  * Scripts = Scripts
  * RUN DATE = Mon Jun 19 19:27:44 2017


```
## R version 3.4.0 (2017-04-21)
## Platform: x86_64-apple-darwin15.6.0 (64-bit)
## Running under: macOS Sierra 10.12.3
## 
## Matrix products: default
## BLAS: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRblas.0.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  base     
## 
## other attached packages:
## [1] BiocInstaller_1.26.0 knitr_1.16           rmarkdown_1.6       
## 
## loaded via a namespace (and not attached):
##  [1] compiler_3.4.0  backports_1.1.0 magrittr_1.5    rprojroot_1.2  
##  [5] htmltools_0.3.6 tools_3.4.0     Rcpp_0.12.11    stringi_1.1.5  
##  [9] methods_3.4.0   digest_0.6.12   stringr_1.2.0   evaluate_0.10
```

