---
title: "29- CpG O/E mosaic structure (Goal 1)"
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

# Overview — Goal 1.4

Computes **CpG observed/expected (CpG O/E)** per gene/CDS from the reference
FASTA to characterize the mosaic genome architecture, then relates CpG O/E to
observed methylation. Bivalve gene sets typically show a **bimodal** CpG O/E
distribution; the low-O/E mode corresponds to historically germline-methylated,
constitutively expressed, more heavily methylated gene bodies.

\[ \text{CpG O/E} = \frac{N_{CpG}}{N_C \times N_G} \times \frac{L^2}{L-1} \]

- Genome FASTA: `../output/22-genome-prep/GCF_036588685.1_PNRI_Mtr1.1.1.hap1_genomic.fa`
- Gene GFF: `../output/22-genome-prep/genomic.gff`
- Observed methylation: `../output/26-methylkit-object/myobj_24.rds`
- Output dir: `../output/29-cpg-oe-mosaic/`


``` bash
mkdir -p ../output/29-cpg-oe-mosaic
```

# Fetch annotation GFF if missing


``` bash
gff="../output/22-genome-prep/genomic.gff"
if [ ! -s "$gff" ]; then
  url="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/036/588/685/GCF_036588685.1_PNRI_Mtr1.1.1.hap1/GCF_036588685.1_PNRI_Mtr1.1.1.hap1_genomic.gff.gz"
  wget -q -O "$gff.gz" "$url" && gunzip -f "$gff.gz"
fi
ls -lh "$gff"
```

```
## -rw-rw-r-- 1 sr320 sr320 316M Mar  6  2024 ../output/22-genome-prep/genomic.gff
```

# Load genome + gene coordinates


``` r
genome <- readDNAStringSet("../output/22-genome-prep/GCF_036588685.1_PNRI_Mtr1.1.1.hap1_genomic.fa")
# Normalize names to first token (matches GFF seqid)
names(genome) <- sub("\\s.*$", "", names(genome))

gff   <- import("../output/22-genome-prep/genomic.gff")
genes <- gff[gff$type == "gene"]
genes <- genes[as.character(seqnames(genes)) %in% names(genome)]
length(genes)
```

```
## [1] 47041
```

# CpG O/E per gene


``` r
cpg_oe <- function(seq) {
  s  <- as.character(seq)
  L  <- nchar(s)
  if (L < 2) return(NA_real_)
  nC   <- lengths(gregexpr("C", s, fixed = TRUE)); nC <- ifelse(grepl("C", s), nC, 0)
  nG   <- lengths(gregexpr("G", s, fixed = TRUE)); nG <- ifelse(grepl("G", s), nG, 0)
  nCpG <- lengths(gregexpr("CG", s, fixed = TRUE)); nCpG <- ifelse(grepl("CG", s), nCpG, 0)
  if (nC == 0 || nG == 0) return(NA_real_)
  (nCpG / (nC * nG)) * (L^2 / (L - 1))
}

# Extract gene sequences from the DNAStringSet by coordinate. (getSeq has no
# DNAStringSet method; subseq is the vectorized equivalent.) CpG O/E is
# strand-symmetric — nC*nG, nCpG and L are invariant under reverse-complement —
# so plus-strand extraction is correct regardless of gene strand.
gi        <- match(as.character(seqnames(genes)), names(genome))
gene_seqs <- subseq(genome[gi], start = start(genes), end = end(genes))
oe <- vapply(seq_along(gene_seqs), function(i) cpg_oe(gene_seqs[[i]]), numeric(1))

gene_oe <- tibble(
  gene_id = genes$ID,
  chr     = as.character(seqnames(genes)),
  start   = start(genes),
  end     = end(genes),
  width   = width(genes),
  cpg_oe  = oe
) |> filter(!is.na(cpg_oe))

write.csv(gene_oe, "../output/29-cpg-oe-mosaic/gene_cpg_oe.csv", row.names = FALSE)
summary(gene_oe$cpg_oe)
```

```
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  0.0000  0.4808  0.7753  0.7299  0.9001  2.7257
```

# Bimodal CpG O/E distribution


``` r
p <- ggplot(gene_oe, aes(cpg_oe)) +
  geom_histogram(bins = 60, fill = "grey30") +
  geom_vline(xintercept = median(gene_oe$cpg_oe), linetype = 2, colour = "red") +
  labs(x = "CpG O/E (per gene)", y = "gene count",
       title = "Per-gene CpG O/E distribution (expect bimodality)") +
  theme_minimal()
ggsave("../output/29-cpg-oe-mosaic/cpg_oe_distribution.png", p, width = 6, height = 4)
```


``` r
# Optional 2-component mixture to formalize the two modes
if (requireNamespace("mixtools", quietly = TRUE)) {
  set.seed(1)
  mix <- mixtools::normalmixEM(gene_oe$cpg_oe, k = 2)
  saveRDS(mix, "../output/29-cpg-oe-mosaic/cpg_oe_mixture.rds")
  data.frame(mode = 1:2, mean = mix$mu, sd = mix$sigma, lambda = mix$lambda)
}
```

# Relate CpG O/E to observed methylation

Aggregate observed per-gene methylation (mean over covered CpGs, 24-sample
pooled >=10x) and join to CpG O/E. Expect low-O/E genes to be more methylated.


``` r
library(methylKit)
myobj_24 <- readRDS("../output/26-methylkit-object/myobj_24.rds")

cpg_meth <- lapply(seq_along(myobj_24), function(i) {
  d <- getData(myobj_24[[i]]); d <- d[d$coverage >= 10, ]
  data.frame(chr = d$chr, start = d$start, pct = 100 * d$numCs / d$coverage)
}) |> bind_rows()

cpg_gr  <- GRanges(cpg_meth$chr, IRanges(cpg_meth$start, width = 1))
gene_gr <- GRanges(gene_oe$chr, IRanges(gene_oe$start, gene_oe$end), gene_id = gene_oe$gene_id)

hits <- findOverlaps(cpg_gr, gene_gr)
gene_meth <- tibble(gene_id = gene_gr$gene_id[subjectHits(hits)],
                    pct     = cpg_meth$pct[queryHits(hits)]) |>
  group_by(gene_id) |>
  summarise(gene_mean_meth = mean(pct), n_cpg = n(), .groups = "drop")

oe_meth <- inner_join(gene_oe, gene_meth, by = "gene_id")
write.csv(oe_meth, "../output/29-cpg-oe-mosaic/cpg_oe_vs_methylation.csv", row.names = FALSE)

p <- ggplot(oe_meth, aes(cpg_oe, gene_mean_meth)) +
  geom_point(alpha = 0.2, size = 0.6) +
  geom_smooth(method = "loess") +
  labs(x = "CpG O/E", y = "mean gene-body % methylation",
       title = "Low-O/E genes are the methylation-prone gene bodies") +
  theme_minimal()
ggsave("../output/29-cpg-oe-mosaic/oe_vs_methylation.png", p, width = 6, height = 4)

cor.test(oe_meth$cpg_oe, oe_meth$gene_mean_meth, method = "spearman")
```

```
## 
## 	Spearman's rank correlation rho
## 
## data:  oe_meth$cpg_oe and oe_meth$gene_mean_meth
## S = 1.5107e+13, p-value < 2.2e-16
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##        rho 
## -0.6823926
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
##  [1] methylKit_1.24.0     rtracklayer_1.58.0   GenomicRanges_1.50.2
##  [4] Biostrings_2.66.0    GenomeInfoDb_1.34.9  XVector_0.38.0      
##  [7] IRanges_2.32.0       S4Vectors_0.36.2     BiocGenerics_0.44.0 
## [10] lubridate_1.9.4      forcats_1.0.0        stringr_1.5.1       
## [13] dplyr_1.1.4          purrr_1.1.0          readr_2.1.5         
## [16] tidyr_1.3.1          tibble_3.3.0         ggplot2_4.0.0       
## [19] tidyverse_2.0.0      knitr_1.50          
## 
## loaded via a namespace (and not attached):
##  [1] MatrixGenerics_1.10.0       Biobase_2.58.0             
##  [3] splines_4.2.3               R.utils_2.13.0             
##  [5] gtools_3.9.5                S7_0.2.0                   
##  [7] GenomeInfoDbData_1.2.9      Rsamtools_2.14.0           
##  [9] yaml_2.3.10                 numDeriv_2016.8-1.1        
## [11] pillar_1.11.0               lattice_0.22-7             
## [13] limma_3.54.2                glue_1.8.0                 
## [15] bbmle_1.0.25.1              RColorBrewer_1.1-3         
## [17] qvalue_2.30.0               colorspace_2.1-1           
## [19] Matrix_1.5-3                R.oo_1.27.1                
## [21] plyr_1.8.9                  XML_3.99-0.18              
## [23] pkgconfig_2.0.3             emdbook_1.3.14             
## [25] zlibbioc_1.44.0             mvtnorm_1.3-3              
## [27] scales_1.4.0                tzdb_0.5.0                 
## [29] BiocParallel_1.32.6         timechange_0.3.0           
## [31] mgcv_1.9-3                  generics_0.1.4             
## [33] farver_2.1.2                withr_3.0.2                
## [35] SummarizedExperiment_1.28.0 cli_3.6.5                  
## [37] mclust_6.1.1                magrittr_2.0.4             
## [39] crayon_1.5.3                evaluate_1.0.4             
## [41] R.methodsS3_1.8.2           nlme_3.1-168               
## [43] MASS_7.3-60                 textshaping_1.0.3          
## [45] data.table_1.14.10          tools_4.2.3                
## [47] hms_1.1.3                   BiocIO_1.8.0               
## [49] lifecycle_1.0.4             matrixStats_1.5.0          
## [51] DelayedArray_0.24.0         compiler_4.2.3             
## [53] fastseg_1.44.0              systemfonts_1.2.3          
## [55] rlang_1.1.6                 grid_4.2.3                 
## [57] RCurl_1.98-1.17             rjson_0.2.23               
## [59] bitops_1.0-9                labeling_0.4.3             
## [61] restfulr_0.0.16             gtable_0.3.6               
## [63] codetools_0.2-20            reshape2_1.4.4             
## [65] R6_2.6.1                    GenomicAlignments_1.34.1   
## [67] bdsmatrix_1.3-7             ragg_1.4.0                 
## [69] stringi_1.8.7               parallel_4.2.3             
## [71] Rcpp_1.1.0                  vctrs_0.6.5                
## [73] tidyselect_1.2.1            xfun_0.52                  
## [75] coda_0.19-4.1
```
