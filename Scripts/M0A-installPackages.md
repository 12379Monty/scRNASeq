###   Install R Packages
 
<style type="text/css"> body, td { font-size: 14px; } code.r{ font-size: 12px; } pre { font-size: 12px } </style>


<!-- ***************************************************** -->

### Update R

[updateR](https://andreacirilloblog.wordpress.com/2015/10/22/updater-package-update-r-version-with-a-function-on-mac-osx/)


```r
  # DO ONCE
  require(devtools)
  install_github('andreacirilloac/updateR')
 
  suppressMessages(require(updateR))

  updateR(admin_password = "***os_admin_user_password***")
```
#### Everything went smoothly, R was updated to R version 3.4.0 (2017-04-21) -- "You Stupid Darkness"


### Upgrade bioconductor


```r
  source("http://www.bioconductor.org/biocLite.R")
  biocLite("BiocUpgrade")

  # check versions
```
#### Error: Bioconductor version 3.5 can be upgraded, but only to 'devel'; see ?useDevel.

#### A worfklow for low-level analyses of single-cell RNA-seq data

[simpleSingleCell](https://www.bioconductor.org/help/workflows/simpleSingleCell/)

#### Be sure to sudo R!

#### Doesn't work!


```r
  # simpleSingleCell workflow
  source("http://bioconductor.org/workflows.R")
  workflowInstall("simpleSingleCell")
```

### Install the rest of the packages

Install all packages used in
[Hicks](https://github.com/stephaniehicks/scBatchPaper)


```r
 RequiredLibs.vec <- unique(
  c('Biobase', 'GenomicRanges', 'vcd', 'stringr', 'reshape2',
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
 
 # Force install of ALL using sudo
 ToIntall.vec <- RequiredLibs.vec

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
## Installing package(s) 'Biobase', 'GenomicRanges', 'vcd', 'stringr',
##   'reshape2', 'tidyr', 'plyr', 'dplyr', 'irlba', 'readr', 'scater',
##   'knitr', 'ggplot2', 'grid', 'cowplot', 'RColorBrewer', 'rafalib',
##   'SummarizedExperiment', 'MultiAssayExperiment',
##   'scRNASeqMouseShalekBMDCs', 'scRNASeqHCT116CellLinesWu',
##   'scRNASeqMouseDengMonoAllelic', 'scRNASeqMouseJaitinSpleen',
##   'scRNASeqMouseKumarPSC', 'scRNASeqHumanPatelGlioblastoma',
##   'scRNASeqMouseShalekDendritic', 'scRNASeqMouseTreutleinLineage',
##   'scRNASeqHumanBosePrinting', 'scRNASeqMouseBurnsInnerEar',
##   'scRNASeqHumanGuoGermCells', 'scRNASeqMouseKowalczykAging',
##   'scRNASeqHumanLengOscillatoryGenes', 'scRNASeqDaniSajitaSeurat',
##   'scRNASeqMouseZeiselCortex'
```

```
## Warning: packages 'grid', 'scRNASeqMouseShalekBMDCs',
## 'scRNASeqHCT116CellLinesWu', 'scRNASeqMouseDengMonoAllelic',
## 'scRNASeqMouseJaitinSpleen', 'scRNASeqMouseKumarPSC',
## 'scRNASeqHumanPatelGlioblastoma', 'scRNASeqMouseShalekDendritic',
## 'scRNASeqMouseTreutleinLineage', 'scRNASeqHumanBosePrinting',
## 'scRNASeqMouseBurnsInnerEar', 'scRNASeqHumanGuoGermCells',
## 'scRNASeqMouseKowalczykAging', 'scRNASeqHumanLengOscillatoryGenes',
## 'scRNASeqDaniSajitaSeurat', 'scRNASeqMouseZeiselCortex' are not available
## (for R version 3.4.0)
```

```
## Warning: package 'grid' is a base package, and should not be updated
```

```
## 
##   There is a binary version available but the source version is
##   later:
##       binary source needs_compilation
## dplyr  0.5.0  0.7.0              TRUE
## 
## 
## The downloaded binary packages are in
## 	/var/folders/27/s0bmn0055zg162hrbt5nhxpm0000gq/T//RtmpMOXvhe/downloaded_packages
```

```
## installing the source package 'dplyr'
```

```
## Old packages: 'bindrcpp', 'DBI', 'RSQLite', 'statmod', 'XML'
```


```r
 require(devtools)
```

```
## Loading required package: devtools
```

```r
 install_github("willtownes/patel2014gliohuman")
```

```
## Downloading GitHub repo willtownes/patel2014gliohuman@master
## from URL https://api.github.com/repos/willtownes/patel2014gliohuman/zipball/master
```

```
## Installing patel2014gliohuman
```

```
## '/Library/Frameworks/R.framework/Resources/bin/R' --no-site-file  \
##   --no-environ --no-save --no-restore --quiet CMD INSTALL  \
##   '/private/var/folders/27/s0bmn0055zg162hrbt5nhxpm0000gq/T/RtmpMOXvhe/devtools2a9f45c5e185/willtownes-patel2014gliohuman-143f352'  \
##   --library='/Library/Frameworks/R.framework/Versions/3.4/Resources/library'  \
##   --install-tests
```

```
## 
```

```
## Installation failed: Command failed (3)
```

```r
 install_github("stephaniehicks/trapnell2014myoblasthuman")
```

```
## Skipping install of 'trapnell2014myoblasthuman' from a github remote, the SHA1 (a4f4dde2) has not changed since last install.
##   Use `force = TRUE` to force installation
```

### Parameter settings:
  * WRKDIR = /Users/francois/Dropbox/SingleCell/Jun2017/R/
  * FN = M0A-installPackages
  * Scripts = Scripts
  * RUN DATE = Mon Jun 19 19:23:26 2017


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
## [1] devtools_1.13.2      BiocInstaller_1.26.0 knitr_1.16          
## [4] rmarkdown_1.6       
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_0.12.11    digest_0.6.12   withr_1.0.2     rprojroot_1.2  
##  [5] R6_2.2.2        backports_1.1.0 git2r_0.18.0    magrittr_1.5   
##  [9] evaluate_0.10   httr_1.2.1      stringi_1.1.5   curl_2.6       
## [13] tools_3.4.0     stringr_1.2.0   compiler_3.4.0  memoise_1.1.0  
## [17] htmltools_0.3.6 methods_3.4.0
```

