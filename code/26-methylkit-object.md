---
title: "26- methylKit object construction & filtering"
author: Steven Roberts
date: "07 July, 2026" 
output: 
  github_document:
    toc: true
    toc_depth: 3
    number_sections: true
    html_preview: true
  html_document:
    theme: readable
    highlight: zenburn
    toc: true
    toc_float: true
    number_sections: true
    code_folding: show
    code_download: true
---



------------------------------------------------------------------------

# Overview

Builds the shared methylKit objects consumed by every downstream analysis from
the merged-CpG `cov` files produced in `25-meth-extract.Rmd`. Two sample sets
are constructed and saved as `.rds` so later steps are resumable
(per the plan's reproducibility notes):

| Object | Set | Samples | Use |
|--------|-----|---------|-----|
| `myobj_24` | 24-sample | all incl. **93M** | Goal 1 landscape / QC (`27`–`29`) |
| `myobj_23` | 23-sample | 24 minus **93M** | Goal 2 differential / biomarkers (`30`–`33`) |

Treatment coding (per `plan.md`, `code/10`/`code/11`): low-PAH **HC, AP, BS = 0**;
high-PAH **EB, SA, SC = 1**.

Coverage filtering follows the plan: **10x primary** (`lo.count = 10`) with a
**5x sensitivity** pass, and `hi.perc = 98` to cap PCR/duplication artefacts.

- Input: `../output/25-meth-extract/<sample>.CpG_report.merged_CpG_evidence.cov`
- Output: `../output/26-methylkit-object/` (`.rds` objects, filtering tables)


``` bash
mkdir -p ../output/26-methylkit-object ../output/25-meth-extract
```

# Fetch the per-sample cov files

The 23 published merged-CpG cov files live on gannet; **93M** is generated
locally in `25.1-93M-cov-recovery` and already sits in `../output/25-meth-extract`.
Pull anything not already present (~770 MB each, idempotent).


``` bash
cov_dir="../output/25-meth-extract"
gannet="https://gannet.fish.washington.edu/seashell/bu-github/project-mytilus-methylation/data"

# every sample except 93M (which is built locally)
for s in 105M 106M 107M 109M 239M 240M 241M 242M 269M 270M 271M 272M \
         69M 70M 71M 72M 78M 79M 80M 81M 92M 94M 95M; do
  f="$cov_dir/${s}.CpG_report.merged_CpG_evidence.cov"
  [ -s "$f" ] || { echo "downloading ${s}"; wget -q -O "$f" "$gannet/${s}.CpG_report.merged_CpG_evidence.cov"; }
done

echo "cov files present: $(ls -1 "$cov_dir"/*.CpG_report.merged_CpG_evidence.cov | wc -l) (expect 24)"
```

```
## cov files present: 24 (expect 24)
```


``` r
# methylKit is expected to be preinstalled; install only if missing.
if (!requireNamespace("methylKit", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("methylKit")
}
```

# Sample sheet


``` r
cov_dir <- "../output/25-meth-extract"

# 24-sample set (landscape). Order groups low (0) then high (1); 93M retained.
samples_24 <- c("105M","106M","107M","109M",   # HC  (low, 0)
                "239M","240M","241M","242M",   # AP  (low, 0)
                "269M","270M","271M","272M",   # BS  (low, 0)
                "69M","70M","71M","72M",       # EB  (high, 1)
                "78M","79M","80M","81M",       # SA  (high, 1)
                "92M","93M","94M","95M")       # SC  (high, 1)

site_24 <- c(rep("HC",4), rep("AP",4), rep("BS",4),
             rep("EB",4), rep("SA",4), rep("SC",4))
pah_24  <- ifelse(site_24 %in% c("HC","AP","BS"), 0L, 1L)

# 23-sample set (differential): drop 93M (coverage/quality outlier)
keep_23   <- samples_24 != "93M"
samples_23 <- samples_24[keep_23]
site_23    <- site_24[keep_23]
pah_23     <- pah_24[keep_23]

meta <- data.frame(sample = samples_24, site = site_24, pah = pah_24)
write.csv(meta, "../output/26-methylkit-object/sample_metadata.csv", row.names = FALSE)
meta
```

```
##    sample site pah
## 1    105M   HC   0
## 2    106M   HC   0
## 3    107M   HC   0
## 4    109M   HC   0
## 5    239M   AP   0
## 6    240M   AP   0
## 7    241M   AP   0
## 8    242M   AP   0
## 9    269M   BS   0
## 10   270M   BS   0
## 11   271M   BS   0
## 12   272M   BS   0
## 13    69M   EB   1
## 14    70M   EB   1
## 15    71M   EB   1
## 16    72M   EB   1
## 17    78M   SA   1
## 18    79M   SA   1
## 19    80M   SA   1
## 20    81M   SA   1
## 21    92M   SC   1
## 22    93M   SC   1
## 23    94M   SC   1
## 24    95M   SC   1
```


``` r
cov_path <- function(s) file.path(cov_dir, paste0(s, ".CpG_report.merged_CpG_evidence.cov"))

files_24 <- as.list(cov_path(samples_24))
files_23 <- as.list(cov_path(samples_23))

stopifnot(all(file.exists(unlist(files_24))))
```

# Read 24-sample set (landscape, descriptive mincov = 1)

For the descriptive landscape pass we read with `mincov = 1` so low-coverage
sites remain visible before any filtering (per plan 1.1).


``` r
myobj_24 <- methRead(
  location   = files_24,
  sample.id  = as.list(samples_24),
  assembly   = "GCF_036588685.1",
  pipeline   = "bismarkCoverage",
  context    = "CpG",
  mincov     = 1,
  treatment  = pah_24
)

saveRDS(myobj_24, "../output/26-methylkit-object/myobj_24.rds")
```

# Read 23-sample set (differential, mincov = 2)


``` r
myobj_23 <- methRead(
  location   = files_23,
  sample.id  = as.list(samples_23),
  assembly   = "GCF_036588685.1",
  pipeline   = "bismarkCoverage",
  context    = "CpG",
  mincov     = 2,
  treatment  = pah_23
)

saveRDS(myobj_23, "../output/26-methylkit-object/myobj_23.rds")
```

# Coverage filtering — 10x primary + 5x sensitivity


``` r
myobj_23 <- readRDS("../output/26-methylkit-object/myobj_23.rds")

# Primary: 10x minimum, cap top 2% of coverage
filt_23_10x <- filterByCoverage(myobj_23,
                                lo.count = 10, lo.perc = NULL,
                                hi.count = NULL, hi.perc = 98)

# Sensitivity: 5x minimum
filt_23_5x  <- filterByCoverage(myobj_23,
                                lo.count = 5,  lo.perc = NULL,
                                hi.count = NULL, hi.perc = 98)

saveRDS(filt_23_10x, "../output/26-methylkit-object/filt_23_10x.rds")
saveRDS(filt_23_5x,  "../output/26-methylkit-object/filt_23_5x.rds")
```

# unite (require coverage in >=5 individuals/group)

cov files are already CpG-merged, so `destrand = FALSE` is correct (destranding
only applies to per-strand cytosine reports).


``` r
meth_23_10x <- unite(filt_23_10x, min.per.group = 5L, destrand = FALSE)
meth_23_5x  <- unite(filt_23_5x,  min.per.group = 5L, destrand = FALSE)

saveRDS(meth_23_10x, "../output/26-methylkit-object/meth_23_10x.rds")
saveRDS(meth_23_5x,  "../output/26-methylkit-object/meth_23_5x.rds")

c(sites_10x = nrow(meth_23_10x), sites_5x = nrow(meth_23_5x))
```

```
## sites_10x  sites_5x 
##   6281686  12765151
```

# CpG-retention trade-off table (1x / 5x / 10x)

Transparent reporting of how many CpGs survive each threshold (plan 1.5 / 2.x).


``` r
retention <- function(obj, thresholds = c(1, 5, 10)) {
  do.call(rbind, lapply(seq_along(obj), function(i) {
    cov <- getData(obj[[i]])$coverage
    data.frame(sample = obj[[i]]@sample.id,
               t(sapply(thresholds, function(t) sum(cov >= t))))
  })) |>
    setNames(c("sample", paste0("CpG_ge_", thresholds, "x")))
}

ret <- retention(myobj_23)
write.csv(ret, "../output/26-methylkit-object/cpg_retention_by_threshold.csv", row.names = FALSE)
ret
```

```
##    sample CpG_ge_1x CpG_ge_5x CpG_ge_10x
## 1    105M  16460075  12551409    6878087
## 2    106M  16224662  12018400    6021521
## 3    107M  16493720  12597887    6963756
## 4    109M  16546232  12678732    7092892
## 5    239M  16856479  13279109    8061758
## 6    240M  16497009  12241169    6486599
## 7    241M  16299956  12206985    6371095
## 8    242M  15909468  11487293    6648088
## 9    269M  16599996  12837018    7335962
## 10   270M  16697049  12963599    7464001
## 11   271M  16474926  12465459    6651761
## 12   272M  16423173  12364749    6534830
## 13    69M  16675373  12682545    7017737
## 14    70M  16700013  12630774    7158347
## 15    71M  16618575  12646397    6946220
## 16    72M  16633582  12709293    6916750
## 17    78M  16763254  12996234    7441908
## 18    79M  16407468  12443337    6703649
## 19    80M  16568264  12722427    7104547
## 20    81M  16513522  12580826    6856760
## 21    92M  16381469  12440635    6700329
## 22    94M  16674446  12810561    7174576
## 23    95M  16259956  12082393    6089099
```

# Session info


``` r
sessionInfo()
```

```
## R version 4.2.3 (2023-03-15)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 24.04.4 LTS
## 
## Matrix products: default
## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so
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
## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
##  [1] methylKit_1.24.0     GenomicRanges_1.50.2 GenomeInfoDb_1.34.9 
##  [4] IRanges_2.32.0       S4Vectors_0.36.2     BiocGenerics_0.44.0 
##  [7] lubridate_1.9.4      forcats_1.0.0        stringr_1.5.1       
## [10] dplyr_1.1.4          purrr_1.1.0          readr_2.1.5         
## [13] tidyr_1.3.1          tibble_3.3.0         ggplot2_4.0.0       
## [16] tidyverse_2.0.0      knitr_1.50          
## 
## loaded via a namespace (and not attached):
##  [1] MatrixGenerics_1.10.0       Biobase_2.58.0             
##  [3] splines_4.2.3               R.utils_2.13.0             
##  [5] gtools_3.9.5                S7_0.2.0                   
##  [7] GenomeInfoDbData_1.2.9      Rsamtools_2.14.0           
##  [9] yaml_2.3.10                 numDeriv_2016.8-1.1        
## [11] pillar_1.11.0               lattice_0.22-7             
## [13] glue_1.8.0                  limma_3.54.2               
## [15] bbmle_1.0.25.1              RColorBrewer_1.1-3         
## [17] XVector_0.38.0              qvalue_2.30.0              
## [19] colorspace_2.1-1            Matrix_1.5-3               
## [21] R.oo_1.27.1                 plyr_1.8.9                 
## [23] XML_3.99-0.18               pkgconfig_2.0.3            
## [25] emdbook_1.3.14              zlibbioc_1.44.0            
## [27] mvtnorm_1.3-3               scales_1.4.0               
## [29] tzdb_0.5.0                  BiocParallel_1.32.6        
## [31] timechange_0.3.0            generics_0.1.4             
## [33] farver_2.1.2                withr_3.0.2                
## [35] SummarizedExperiment_1.28.0 cli_3.6.5                  
## [37] magrittr_2.0.4              crayon_1.5.3               
## [39] mclust_6.1.1                evaluate_1.0.4             
## [41] R.methodsS3_1.8.2           MASS_7.3-60                
## [43] tools_4.2.3                 data.table_1.14.10         
## [45] hms_1.1.3                   matrixStats_1.5.0          
## [47] BiocIO_1.8.0                lifecycle_1.0.4            
## [49] DelayedArray_0.24.0         Biostrings_2.66.0          
## [51] compiler_4.2.3              fastseg_1.44.0             
## [53] rlang_1.1.6                 grid_4.2.3                 
## [55] RCurl_1.98-1.17             rjson_0.2.23               
## [57] bitops_1.0-9                restfulr_0.0.16            
## [59] gtable_0.3.6                codetools_0.2-20           
## [61] reshape2_1.4.4              R6_2.6.1                   
## [63] GenomicAlignments_1.34.1    rtracklayer_1.58.0         
## [65] bdsmatrix_1.3-7             stringi_1.8.7              
## [67] parallel_4.2.3              Rcpp_1.1.0                 
## [69] vctrs_0.6.5                 tidyselect_1.2.1           
## [71] xfun_0.52                   coda_0.19-4.1
```
