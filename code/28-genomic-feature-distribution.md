---
title: "28- Genomic-feature distribution (Goal 1)"
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

# Overview — Goal 1.3

Distributes covered CpGs (and their methylation levels) across genomic features
to demonstrate the canonical molluscan **gene-body–enriched mosaic** pattern.
Builds feature BEDs from the `GCF_036588685.1` GFF (`22-genome-prep.Rmd`) with
BEDTools, intersects covered CpG positions against each feature class, and
summarizes CpG fraction + mean methylation per feature.

- GFF: `../output/22-genome-prep/genomic.gff`
- Genome index: `../output/22-genome-prep/GCF_036588685.1_PNRI_Mtr1.1.1.hap1_genomic.fa.fai`
- 24-sample object: `../output/26-methylkit-object/myobj_24.rds`
- Output dir: `../output/28-genomic-feature-distribution/`


``` bash
mkdir -p ../output/28-genomic-feature-distribution/features
```

# Annotation GFF

Fetch the `GCF_036588685.1` GFF from NCBI if it is not already staged in
`22-genome-prep` (the gannet genome mirror serves the FASTA 403; NCBI is used
for both). Idempotent.


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

# Per-feature GFFs from the annotation


``` bash
input_gff="../output/22-genome-prep/genomic.gff"
out_dir="../output/28-genomic-feature-distribution/features"

features=("gene" "mRNA" "exon" "CDS" "lnc_RNA" "ncRNA" "tRNA" "rRNA" "pseudogene")

for feature in "${features[@]}"; do
  awk -v ft="$feature" '$3 == ft && $1 !~ /^#/ {print}' "$input_gff" \
    > "$out_dir/${feature}.gff"
  echo "Created $out_dir/${feature}.gff ($(wc -l < "$out_dir/${feature}.gff") records)"
done
```

```
## Created ../output/28-genomic-feature-distribution/features/gene.gff (47041 records)
## Created ../output/28-genomic-feature-distribution/features/mRNA.gff (53269 records)
## Created ../output/28-genomic-feature-distribution/features/exon.gff (558468 records)
## Created ../output/28-genomic-feature-distribution/features/CDS.gff (495644 records)
## Created ../output/28-genomic-feature-distribution/features/lnc_RNA.gff (4877 records)
## Created ../output/28-genomic-feature-distribution/features/ncRNA.gff (1 records)
## Created ../output/28-genomic-feature-distribution/features/tRNA.gff (5782 records)
## Created ../output/28-genomic-feature-distribution/features/rRNA.gff (1615 records)
## Created ../output/28-genomic-feature-distribution/features/pseudogene.gff (739 records)
```

# Derived contexts (intron / promoter / intergenic)


``` bash
# BEDTools genome file from the FASTA index (chrom <tab> length).
# Index the FASTA first if it was not indexed in 22-genome-prep.
fa=../output/22-genome-prep/GCF_036588685.1_PNRI_Mtr1.1.1.hap1_genomic.fa
[ -s "$fa.fai" ] || /home/shared/samtools-1.12/samtools faidx "$fa"

cut -f1,2 "$fa.fai" \
  | sort -k1,1 > ../output/28-genomic-feature-distribution/genome.txt
wc -l ../output/28-genomic-feature-distribution/genome.txt
head ../output/28-genomic-feature-distribution/genome.txt
```

```
## 541 ../output/28-genomic-feature-distribution/genome.txt
## NC_007687.1	18652
## NC_086373.1	113486924
## NC_086374.1	103707636
## NC_086375.1	94312406
## NC_086376.1	96324366
## NC_086377.1	88129055
## NC_086378.1	91445414
## NC_086379.1	85406108
## NC_086380.1	80675551
## NC_086381.1	82295107
```


``` bash
feat="../output/28-genomic-feature-distribution/features"
genome="../output/28-genomic-feature-distribution/genome.txt"

# GFF -> sorted BED helper (cols 1,4-1,5 -> 0-based start)
gff2bed () { awk 'BEGIN{OFS="\t"} $1!~/^#/ {print $1, $4-1, $5}' "$1" | sort -k1,1 -k2,2n; }

gff2bed "$feat/gene.gff" > "$feat/gene.bed"
gff2bed "$feat/exon.gff" > "$feat/exon.bed"
gff2bed "$feat/mRNA.gff" > "$feat/mRNA.bed"

# Introns = gene minus exon
bedtools subtract -a "$feat/gene.bed" -b "$feat/exon.bed" > "$feat/intron.bed"

# Promoters/TSS = 1 kb upstream flank of genes (strand-aware via gene GFF)
awk 'BEGIN{OFS="\t"} $1!~/^#/ {print $1,$4-1,$5,".",".",$7}' "$feat/gene.gff" \
  | sort -k1,1 -k2,2n > "$feat/gene_stranded.bed"
bedtools flank -i "$feat/gene_stranded.bed" -g "$genome" -l 1000 -r 0 -s \
  | cut -f1-3 | sort -k1,1 -k2,2n > "$feat/promoter_1kb.bed"

# Intergenic = genome complement of genes
bedtools complement -i "$feat/gene.bed" -g "$genome" > "$feat/intergenic.bed"

wc -l "$feat"/{gene,exon,intron,promoter_1kb,intergenic}.bed
```

```
##    47041 ../output/28-genomic-feature-distribution/features/gene.bed
##   558468 ../output/28-genomic-feature-distribution/features/exon.bed
##   255783 ../output/28-genomic-feature-distribution/features/intron.bed
##    47039 ../output/28-genomic-feature-distribution/features/promoter_1kb.bed
##    43995 ../output/28-genomic-feature-distribution/features/intergenic.bed
##   952326 total
```

> Repeats: if a RepeatMasker `.out`/`.bed` track becomes available, add a
> `repeat.bed` here and include it below; otherwise note as an annotation gap.

# Covered CpG positions -> BED (>=10x, pooled across 24 samples)


``` r
myobj_24 <- readRDS("../output/26-methylkit-object/myobj_24.rds")

# Pool: a CpG is "covered" if >=10x in any sample; record mean % meth across
# samples meeting that threshold.
cpg_tab <- lapply(seq_along(myobj_24), function(i) {
  d <- getData(myobj_24[[i]]); d <- d[d$coverage >= 10, ]
  data.frame(chr = d$chr, start = d$start, end = d$end,
             pct = 100 * d$numCs / d$coverage)
}) |> bind_rows() |>
  group_by(chr, start, end) |>
  summarise(pct_meth = mean(pct), .groups = "drop")

options(scipen = 999)   # never write BED coordinates in scientific notation
bed <- cpg_tab |>
  transmute(chr,
            start = as.integer(start - 1),   # integers can't render as 1.06e+08
            end   = as.integer(end),
            pct_meth) |>
  arrange(chr, start)

write.table(bed, "../output/28-genomic-feature-distribution/covered_cpg_10x.bed",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
nrow(bed)
```

```
## [1] 17159645
```

# Intersect covered CpGs against each feature


``` bash
feat="../output/28-genomic-feature-distribution/features"
cpg="../output/28-genomic-feature-distribution/covered_cpg_10x.bed"
out="../output/28-genomic-feature-distribution"

sort -k1,1 -k2,2n "$cpg" > "$out/covered_cpg_10x.sorted.bed"

for f in gene exon intron mRNA promoter_1kb intergenic; do
  bedtools intersect -a "$out/covered_cpg_10x.sorted.bed" -b "$feat/${f}.bed" -u \
    > "$out/cpg_in_${f}.bed"
  echo "${f} $(wc -l < "$out/cpg_in_${f}.bed")"
done
```

```
## gene 7452241
## exon 1242316
## intron 6235238
## mRNA 7125665
## promoter_1kb 614337
## intergenic 9709347
```

# Feature x methylation-level summary + barplot


``` r
out <- "../output/28-genomic-feature-distribution"
features <- c("gene","exon","intron","mRNA","promoter_1kb","intergenic")

total_cpg <- nrow(read.table(file.path(out, "covered_cpg_10x.sorted.bed")))

feat_summary <- map_dfr(features, function(f) {
  fp <- file.path(out, paste0("cpg_in_", f, ".bed"))
  if (file.size(fp) == 0) return(tibble(feature = f, n_cpg = 0,
                                         frac_cpg = 0, mean_pct_meth = NA))
  d <- read.table(fp)
  tibble(feature = f, n_cpg = nrow(d),
         frac_cpg = nrow(d) / total_cpg,
         mean_pct_meth = mean(d$V4))
})

write.csv(feat_summary, file.path(out, "feature_methylation_summary.csv"),
          row.names = FALSE)

p <- ggplot(feat_summary, aes(reorder(feature, -mean_pct_meth), mean_pct_meth)) +
  geom_col(fill = "steelblue") +
  labs(x = "feature", y = "mean % CpG methylation",
       title = "Methylation by genomic feature (expect gene-body enrichment)") +
  theme_minimal()
ggsave(file.path(out, "feature_methylation_barplot.png"), p, width = 6, height = 4)
feat_summary
```

```
## # A tibble: 6 × 4
##   feature        n_cpg frac_cpg mean_pct_meth
##   <chr>          <int>    <dbl>         <dbl>
## 1 gene         7452241   0.434          19.7 
## 2 exon         1242316   0.0724         19.6 
## 3 intron       6235238   0.363          19.7 
## 4 mRNA         7125665   0.415          20.0 
## 5 promoter_1kb  614337   0.0358          7.30
## 6 intergenic   9709347   0.566           5.44
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
## [43] textshaping_1.0.3           tools_4.2.3                
## [45] data.table_1.14.10          hms_1.1.3                  
## [47] matrixStats_1.5.0           BiocIO_1.8.0               
## [49] lifecycle_1.0.4             DelayedArray_0.24.0        
## [51] Biostrings_2.66.0           compiler_4.2.3             
## [53] systemfonts_1.2.3           fastseg_1.44.0             
## [55] rlang_1.1.6                 grid_4.2.3                 
## [57] RCurl_1.98-1.17             rjson_0.2.23               
## [59] labeling_0.4.3              bitops_1.0-9               
## [61] restfulr_0.0.16             gtable_0.3.6               
## [63] codetools_0.2-20            reshape2_1.4.4             
## [65] R6_2.6.1                    GenomicAlignments_1.34.1   
## [67] rtracklayer_1.58.0          bdsmatrix_1.3-7            
## [69] utf8_1.2.6                  ragg_1.4.0                 
## [71] stringi_1.8.7               parallel_4.2.3             
## [73] Rcpp_1.1.0                  vctrs_0.6.5                
## [75] tidyselect_1.2.1            xfun_0.52                  
## [77] coda_0.19-4.1
```
