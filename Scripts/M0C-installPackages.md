###   Install R Packages - conitued
 
Here add installs required for [seurat](//http://satijalab.org/seurat)

<style type="text/css"> body, td { font-size: 14px; } code.r{ font-size: 12px; } pre { font-size: 12px } </style>


<!-- ***************************************************** -->

#### Install seurat From Source


```r
 require(devtools)
```

```
## Loading required package: devtools
```

```r
 install_github("satijalab/seurat")
```

```
## Downloading GitHub repo satijalab/seurat@master
## from URL https://api.github.com/repos/satijalab/seurat/zipball/master
```

```
## Installing Seurat
```

```
## Installing ape
```

```
## '/usr/lib/R/bin/R' --no-site-file --no-environ --no-save --no-restore  \
##   --quiet CMD INSTALL '/tmp/RtmpY1facf/devtoolsea412e9308/ape'  \
##   --library='/usr/local/lib/R/site-library' --install-tests
```

```
## 
```

```
## Installing caret
```

```
## Installing ModelMetrics
```

```
## '/usr/lib/R/bin/R' --no-site-file --no-environ --no-save --no-restore  \
##   --quiet CMD INSTALL '/tmp/RtmpY1facf/devtoolsea47d5523e2/ModelMetrics'  \
##   --library='/usr/local/lib/R/site-library' --install-tests
```

```
## 
```

```
## '/usr/lib/R/bin/R' --no-site-file --no-environ --no-save --no-restore  \
##   --quiet CMD INSTALL '/tmp/RtmpY1facf/devtoolsea45a784c2d/caret'  \
##   --library='/usr/local/lib/R/site-library' --install-tests
```

```
## 
```

```
## Installing fastICA
```

```
## '/usr/lib/R/bin/R' --no-site-file --no-environ --no-save --no-restore  \
##   --quiet CMD INSTALL '/tmp/RtmpY1facf/devtoolsea45d1dc487/fastICA'  \
##   --library='/usr/local/lib/R/site-library' --install-tests
```

```
## 
```

```
## Installing lars
```

```
## '/usr/lib/R/bin/R' --no-site-file --no-environ --no-save --no-restore  \
##   --quiet CMD INSTALL '/tmp/RtmpY1facf/devtoolsea41052a02f/lars'  \
##   --library='/usr/local/lib/R/site-library' --install-tests
```

```
## 
```

```
## Installing mixtools
```

```
## Installing segmented
```

```
## '/usr/lib/R/bin/R' --no-site-file --no-environ --no-save --no-restore  \
##   --quiet CMD INSTALL '/tmp/RtmpY1facf/devtoolsea44983311e/segmented'  \
##   --library='/usr/local/lib/R/site-library' --install-tests
```

```
## 
```

```
## '/usr/lib/R/bin/R' --no-site-file --no-environ --no-save --no-restore  \
##   --quiet CMD INSTALL '/tmp/RtmpY1facf/devtoolsea424d48b85/mixtools'  \
##   --library='/usr/local/lib/R/site-library' --install-tests
```

```
## 
```

```
## Installing pbapply
```

```
## '/usr/lib/R/bin/R' --no-site-file --no-environ --no-save --no-restore  \
##   --quiet CMD INSTALL '/tmp/RtmpY1facf/devtoolsea43e5fd636/pbapply'  \
##   --library='/usr/local/lib/R/site-library' --install-tests
```

```
## 
```

```
## Installing ranger
```

```
## '/usr/lib/R/bin/R' --no-site-file --no-environ --no-save --no-restore  \
##   --quiet CMD INSTALL '/tmp/RtmpY1facf/devtoolsea44d680a7d/ranger'  \
##   --library='/usr/local/lib/R/site-library' --install-tests
```

```
## 
```

```
## Installing RcppProgress
```

```
## '/usr/lib/R/bin/R' --no-site-file --no-environ --no-save --no-restore  \
##   --quiet CMD INSTALL '/tmp/RtmpY1facf/devtoolsea46597c79/RcppProgress'  \
##   --library='/usr/local/lib/R/site-library' --install-tests
```

```
## 
```

```
## Installing ROCR
```

```
## '/usr/lib/R/bin/R' --no-site-file --no-environ --no-save --no-restore  \
##   --quiet CMD INSTALL '/tmp/RtmpY1facf/devtoolsea41b9e3d60/ROCR'  \
##   --library='/usr/local/lib/R/site-library' --install-tests
```

```
## 
```

```
## Installing tclust
```

```
## Installing sn
```

```
## Installing mnormt
```

```
## '/usr/lib/R/bin/R' --no-site-file --no-environ --no-save --no-restore  \
##   --quiet CMD INSTALL '/tmp/RtmpY1facf/devtoolsea42adc696c/mnormt'  \
##   --library='/usr/local/lib/R/site-library' --install-tests
```

```
## 
```

```
## Installing numDeriv
```

```
## '/usr/lib/R/bin/R' --no-site-file --no-environ --no-save --no-restore  \
##   --quiet CMD INSTALL '/tmp/RtmpY1facf/devtoolsea42527fd06/numDeriv'  \
##   --library='/usr/local/lib/R/site-library' --install-tests
```

```
## 
```

```
## '/usr/lib/R/bin/R' --no-site-file --no-environ --no-save --no-restore  \
##   --quiet CMD INSTALL '/tmp/RtmpY1facf/devtoolsea4479708e1/sn'  \
##   --library='/usr/local/lib/R/site-library' --install-tests
```

```
## 
```

```
## '/usr/lib/R/bin/R' --no-site-file --no-environ --no-save --no-restore  \
##   --quiet CMD INSTALL '/tmp/RtmpY1facf/devtoolsea43f63ce69/tclust'  \
##   --library='/usr/local/lib/R/site-library' --install-tests
```

```
## 
```

```
## Installing tsne
```

```
## '/usr/lib/R/bin/R' --no-site-file --no-environ --no-save --no-restore  \
##   --quiet CMD INSTALL '/tmp/RtmpY1facf/devtoolsea424251451/tsne'  \
##   --library='/usr/local/lib/R/site-library' --install-tests
```

```
## 
```

```
## Installing VGAM
```

```
## '/usr/lib/R/bin/R' --no-site-file --no-environ --no-save --no-restore  \
##   --quiet CMD INSTALL '/tmp/RtmpY1facf/devtoolsea4404826ac/VGAM'  \
##   --library='/usr/local/lib/R/site-library' --install-tests
```

```
## 
```

```
## '/usr/lib/R/bin/R' --no-site-file --no-environ --no-save --no-restore  \
##   --quiet CMD INSTALL  \
##   '/tmp/RtmpY1facf/devtoolsea41544c0e7/satijalab-seurat-3bd092a'  \
##   --library='/usr/local/lib/R/site-library' --install-tests
```

```
## 
```

### Install additional packages required


```r
 RequiredLibs.vec <- unique(
  c(
    ))

 InstalledLibs.vec <- installed.packages()[,1]
 ToIntall.vec <- setdiff(RequiredLibs.vec,InstalledLibs.vec)

 if(length(ToIntall.vec)){
  source("http://www.bioconductor.org/biocLite.R")
  biocLite(ToIntall.vec)
 }
```

### Parameter settings:
  * WRKDIR = /mnt100/home/Dropbox/SingleCell/Jun2017/R
  * FN = M0C-installPackages
  * Scripts = Scripts
  * RUN DATE = Wed Jun 21 23:47:52 2017


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
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  base     
## 
## other attached packages:
## [1] devtools_1.13.2 knitr_1.16      rmarkdown_1.6  
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_0.12.11    digest_0.6.12   withr_1.0.2     rprojroot_1.2  
##  [5] R6_2.2.2        backports_1.1.0 git2r_0.18.0    magrittr_1.5   
##  [9] evaluate_0.10   httr_1.2.1      stringi_1.1.5   curl_2.6       
## [13] tools_3.4.0     stringr_1.2.0   compiler_3.4.0  memoise_1.1.0  
## [17] htmltools_0.3.6 methods_3.4.0
```

