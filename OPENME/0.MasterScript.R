## Author: Víctor Taracena
## email: vtaracena@gmail.com

## This script, and the consequent ones, are based on Maggie Wagner's and collaborators' scripts for her research on microbial 
## communities in Boecheria stricta  (Dryad: http://dx.doi.org/10.5061/dryad.g60r3; 
## Article: http://dx.doi.org/10.1038/ncomms12151) and the bioconductor workflow by Callaghan et al. 
## (https://f1000research.com/articles/5-1492/v2). 

                        

                          ## Master Script ##

## This script will create all the directories and move all the files to their place that you need
## to run this workflow. Also it will install and load all the libraries you need, set a function for rarefaction curves
## and the session information of my PC when I ran these scripts. 

#### Note: This commands works in Linux / MacOS, if you are using Windows you'll have to change
#### to the commands that Windows use or do it manually. 

#### IMPORTANT: set working directory to Source file location. 
#### Every file and directory will be created inside this directory.

#### Create the directories ####
system('mkdir RawData')
system('mkdir bin')
system('mkdir Data')
system('mkdir Plots')

#### Move files to their new homes ####
system('mv Initial/1.PreProcessing.R bin/')
system('mv Initial/2.AssignTaxonomy.R bin/')
system('mv Initial/3.DiversityAnalyses.R bin/')
system('mv Initial/4.MultivariateAnalyses.R bin/')
system('mv Initial/OTUtable.txt.bz2 RawData/')
system('mv Initial/contaminants.fasta RawData/')
system('mv Initial/phylogeny.tre.bz2 RawData/')
system('mv Initial/SampleMD.txt RawData/')
system('mv Initial/TaxAssignment.txt.bz2 RawData/')
system('mv Initial/SiteCoords.txt RawData/')

#### Unzip the .bz2 files ####
system('bzip2 -dk RawData/OTUtable.txt.bz2')
system('bzip2 -dk RawData/phylogeny.tre.bz2')
system('bzip2 -dk RawData/TaxAssignment.txt.bz2')

#### Remove unnecessary files ####
system('rm -r Initial')
system('rm RawData/OTUtable.txt.bz2')
system('rm RawData/phylogeny.tre.bz2')
system('rm RawData/TaxAssignment.txt.bz2')

#### Clear workspace ####
rm(list=ls())

#### Install required packages ####
.cran_packages <- c("grid","gridExtra","plyr","dplyr", "tidyr", "rmarkdown",
                    "ggplot2", "gridExtra", "magrittr", "ape", "ggplot2", 
                    "vegan", "reshape2", "VennDiagram", "maps")
.bioc_packages <- c("BiocStyle","dada2", "phyloseq", "DECIPHER", "phangorn",
                    "ShortRead")

.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {
  install.packages(.cran_packages[!.inst])
}

.inst <- .bioc_packages %in% installed.packages()
if(any(!.inst)) {
  source("http://bioconductor.org/biocLite.R")
  biocLite(.bioc_packages[!.inst], ask = F)
}

#### Load packages into session ####
sapply(c(.cran_packages, .bioc_packages), require, character.only = TRUE)
set.seed(100)

#### Set some functions we will need later ####
calculate_rarefaction_curves <- function(psdata, measures, depths) {
  require('plyr') # ldply
  require('reshape2') # melt
  
  estimate_rarified_richness <- function(psdata, measures, depth) {
    if(max(sample_sums(psdata)) < depth) return()
    psdata <- prune_samples(sample_sums(psdata) >= depth, psdata)
    
    rarified_psdata <- rarefy_even_depth(psdata, depth, verbose = FALSE)
    
    alpha_diversity <- estimate_richness(rarified_psdata, measures = measures)
    
    # as.matrix forces the use of melt.array, which includes the Sample names (rownames)
    molten_alpha_diversity <- melt(as.matrix(alpha_diversity), varnames = c('Sample', 'Measure'), value.name = 'Alpha_diversity')
    
    molten_alpha_diversity
  }
  
  names(depths) <- depths # this enables automatic addition of the Depth to the output by ldply
  rarefaction_curve_data <- ldply(depths, estimate_rarified_richness, psdata = psdata, 
                                  measures = measures, .id = 'Depth', .progress = ifelse(interactive(), 'text', 'none'))
  
  # convert Depth from factor to numeric
  rarefaction_curve_data$Depth <- as.numeric(levels(rarefaction_curve_data$Depth))[rarefaction_curve_data$Depth]
  
  rarefaction_curve_data
}

#### Session Info ####
sessionInfo()

"R version 3.4.0 (2017-04-21)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 16.04.2 LTS

Matrix products: default
BLAS: /usr/lib/libblas/libblas.so.3.6.0
LAPACK: /usr/lib/lapack/liblapack.so.3.6.0

locale:
[1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=es_MX.UTF-8       
[4] LC_COLLATE=en_US.UTF-8     LC_MONETARY=es_MX.UTF-8    LC_MESSAGES=en_US.UTF-8   
[7] LC_PAPER=es_MX.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
[10] LC_TELEPHONE=C             LC_MEASUREMENT=es_MX.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods  
[9] base     

other attached packages:
[1] BiocStyle_2.4.0            rmarkdown_1.5              biomformat_1.4.0          
[4] BiocInstaller_1.26.0       shiny_1.0.3                phangorn_2.2.0            
[7] ape_4.1                    DECIPHER_2.4.0             RSQLite_1.1-2             
[10] dada2_1.4.0                Rcpp_0.12.10               gridExtra_2.2.1           
[13] ggplot2_2.2.1              magrittr_1.5               phyloseq_1.20.0           
[16] ShortRead_1.34.0           GenomicAlignments_1.12.1   SummarizedExperiment_1.6.1
[19] DelayedArray_0.2.2         matrixStats_0.52.2         Biobase_2.36.2            
[22] Rsamtools_1.28.0           GenomicRanges_1.28.2       GenomeInfoDb_1.12.0       
[25] BiocParallel_1.10.1        Biostrings_2.44.0          XVector_0.16.0            
[28] IRanges_2.10.1             S4Vectors_0.14.1           BiocGenerics_0.22.0       

loaded via a namespace (and not attached):
[1] nlme_3.1-131            bitops_1.0-6            RColorBrewer_1.1-2     
[4] rprojroot_1.2           tools_3.4.0             backports_1.0.5        
[7] R6_2.2.1                vegan_2.4-3             rpart_4.1-10           
[10] Hmisc_4.0-3             DBI_0.6-1               lazyeval_0.2.0         
[13] mgcv_1.8-16             colorspace_1.3-2        permute_0.9-4          
[16] ade4_1.7-6              nnet_7.3-12             DESeq2_1.16.1          
[19] compiler_3.4.0          htmlTable_1.9           scales_0.4.1           
[22] checkmate_1.8.2         genefilter_1.58.1       quadprog_1.5-5         
[25] stringr_1.2.0           digest_0.6.12           foreign_0.8-67         
[28] base64enc_0.1-3         htmltools_0.3.6         htmlwidgets_0.8        
[31] rlang_0.1               hwriter_1.3.2           jsonlite_1.4           
[34] acepack_1.4.1           RCurl_1.95-4.8          GenomeInfoDbData_0.99.0
[37] Formula_1.2-1           Matrix_1.2-8            munsell_0.4.3          
[40] yaml_2.1.14             stringi_1.1.5           MASS_7.3-45            
[43] zlibbioc_1.22.0         rhdf5_2.20.0            plyr_1.8.4             
[46] grid_3.4.0              lattice_0.20-35         splines_3.4.0          
[49] multtest_2.32.0         annotate_1.54.0         locfit_1.5-9.1         
[52] knitr_1.15.1            igraph_1.0.1            geneplotter_1.54.0     
[55] reshape2_1.4.2          codetools_0.2-15        fastmatch_1.1-0        
[58] XML_3.98-1.7            evaluate_0.10           latticeExtra_0.6-28    
[61] data.table_1.10.4       RcppParallel_4.3.20     httpuv_1.3.3           
[64] foreach_1.4.3           gtable_0.2.0            mime_0.5               
[67] xtable_1.8-2            survival_2.41-3         tibble_1.3.1           
[70] iterators_1.0.8         AnnotationDbi_1.38.0    memoise_1.1.0          
[73] cluster_2.0.6"


          ## Please continue to script 1.PreProcessing.R ##

