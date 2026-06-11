# project-mytilus-methylation

Whole-genome bisulfite sequencing (WGBS) of mantle tissue from **24 individual *Mytilus trossulus*** collected across six Puget Sound sites that contrast **high** versus **low Σ16 PAH** (sum of 16 priority polycyclic aromatic hydrocarbons) exposure. Samples come from the **2021–2022 Washington Department of Fish and Wildlife (WDFW) Nearshore Monitoring deployment**. This repository pairs (1) a descriptive characterization of the *M. trossulus* CpG methylation landscape with (2) discovery of candidate DNA-methylation biomarkers that distinguish high- from low-PAH mussels.

See [`plan.md`](plan.md) for the authoritative study design and full analysis plan, and [`paper.md`](paper.md) for manuscript context.

## Background & motivation

Legacy PAH contamination persists in nearshore sediments and bioaccumulates in sessile filter feeders, yet traditional enzyme-activity biomarkers respond only within narrow windows that are hard to align with regional monitoring logistics. DNA methylation is a more stable molecular record of exposure. In bivalves, methylation occurs almost exclusively at CpG sites in a **mosaic pattern** enriched in the gene bodies of constitutively expressed genes, so differentially methylated loci can reveal cellular responses to environmental stressors. This project asks whether methylation in *M. trossulus* mantle tissue carries a PAH-associated signal that could seed cost-effective, targeted biomarker assays for the broader Nearshore Monitoring Program.

## Study design & samples

- **Species / tissue:** *Mytilus trossulus*, mantle
- **Sequencing:** paired-end 150 bp WGBS (Illumina NovaSeq)
- **Design:** 6 sites (3 low PAH, 3 high PAH), 4 individuals per site = 24 mussels
- **Exposure metric:** Σ16 PAH (ng/g); low range 97.4–127.26, high range 3,304.96–12,078.42
- In methylKit, low-PAH sites are coded as treatment `0` and high-PAH sites as treatment `1`.

| Code | Site | PAH group | Treatment | n |
|------|------|-----------|-----------|---|
| HC | Hood Canal | Low | 0 | 4 |
| AP | Holly Aiston Preserve | Low | 0 | 4 |
| BS | Broad Spit | Low | 0 | 4 |
| EB | Elliot Bay | High | 1 | 4 |
| SA | Seattle Aquarium Pier 59 | High | 1 | 4 |
| SC | Smith Cove (Terminal 91) | High | 1 | 4 |

All 24 individuals are sequenced and processed through Bismark. Sample **93M** (Smith Cove) is a coverage/quality outlier: it is retained in the **24-sample** landscape characterization (Goal 1) but excluded from the **23-sample** differential set used for biomarker discovery (Goal 2). Per-sample IDs and descriptions are listed in full in [`plan.md`](plan.md).

## Analysis goals

- **Goal 1 — Methylation landscape.** A genome-wide, treatment-independent description of CpG methylation: per-sample global %CpG methylation and coverage, the bivalve bimodal/mosaic distribution, CpG class proportions, genomic-feature enrichment (gene-body bias), CpG observed/expected structure, and sample correlation/clustering/PCA (including the formal 93M outlier decision).
- **Goal 2 — High-vs-low-PAH biomarkers.** Per-CpG differential methylation (DML) and tiled-window DMRs between low- (HC/AP/BS) and high-PAH (EB/SA/SC) mussels, with overdispersion correction, site covariate, and per-site concordance filtering; annotation against the reference genome; functional/xenobiotic enrichment; and a ranked shortlist of biomarker candidates for targeted assays.

## Repository structure

| Path | Contents |
|------|----------|
| [`code/`](code/) | R Markdown analysis scripts; see [`code/README.md`](code/README.md) for naming conventions |
| [`output/`](output/) | One directory per script (matching name), per [`output/README.md`](output/README.md) |
| [`data/`](data/) | Reference genome / NCBI dataset and raw-data download instructions |
| [`plan.md`](plan.md) | Authoritative study design, sample table, and full analysis plan |
| [`paper.md`](paper.md) | Manuscript draft (title, abstract, methods, results, discussion) |

## Pipeline / workflow

The `20`–`34` series is a complete, self-contained end-to-end pipeline that runs from raw WGBS reads through differential methylation and biomarker ranking. It does **not** depend on the legacy `00.*`–`14` scripts; each step regenerates its own inputs and writes to a matching `output/<name>/` directory.

### Upstream processing (all 24 samples)

- [`code/20-raw-data-qc.Rmd`](code/20-raw-data-qc.Rmd): Download raw paired-end FASTQs from `nightingales/M_trossulus`, verify MD5 checksums, raw-read FastQC + MultiQC.
- [`code/21-fastp-trimming.Rmd`](code/21-fastp-trimming.Rmd): fastp trimming (Illumina adapter auto-detect + 20 bp hard trim on both ends) and post-trim FastQC/MultiQC.
- [`code/22-genome-prep.Rmd`](code/22-genome-prep.Rmd): Download NCBI `GCF_036588685.1` genome + GFF; `samtools faidx`; `bismark_genome_preparation` (bisulfite Bowtie2 index).
- [`code/23-bismark-align.Rmd`](code/23-bismark-align.Rmd): Bismark/Bowtie2 alignment of trimmed reads to the converted genome.
- [`code/24-dedup-sort.Rmd`](code/24-dedup-sort.Rmd): `deduplicate_bismark` then SAMtools coordinate sort + index.
- [`code/25-meth-extract.Rmd`](code/25-meth-extract.Rmd): `bismark_methylation_extractor` + `coverage2cytosine --merge_CpG --zero_based` to produce merged-CpG `.cov` reports.

### methylKit object construction

- [`code/26-methylkit-object.Rmd`](code/26-methylkit-object.Rmd): `methRead` the `cov` files into the 24-sample (landscape, incl. 93M) and 23-sample (differential) sets; `filterByCoverage` 10× primary + 5× sensitivity (`hi.perc = 98`); `unite`; save `.rds`.

### Goal 1 — methylation landscape

- [`code/27-methylation-landscape.Rmd`](code/27-methylation-landscape.Rmd): Per-sample global %CpG methylation + coverage stats, bimodality plot, CpG class proportions, correlation/clustering/PCA, 93M outlier decision.
- [`code/28-genomic-feature-distribution.Rmd`](code/28-genomic-feature-distribution.Rmd): Build gene/exon/intron/promoter/intergenic BEDs from the GFF; `bedtools intersect` covered CpGs × features; feature × methylation barplot (gene-body enrichment).
- [`code/29-cpg-oe-mosaic.Rmd`](code/29-cpg-oe-mosaic.Rmd): Per-gene CpG O/E from the reference FASTA; bimodal O/E distribution; O/E vs observed methylation.

### Goal 2 — high vs low PAH biomarkers

- [`code/30-dml-dmr-methylkit.Rmd`](code/30-dml-dmr-methylkit.Rmd): Per-CpG DML (`calculateDiffMeth`, overdispersion `"MN"` + site covariate); 1 kb tiled DMRs (optional `dmrseq`); thresholds q < 0.01 / ≥55% (DML), ≥25% (DMR); BH-FDR + SLIM; per-site concordance; 5× sensitivity.
- [`code/31-dml-dmr-annotation.Rmd`](code/31-dml-dmr-annotation.Rmd): DML/DMR → bedGraph; `bedtools intersect`/`closest` vs `genomic.gff`; recover LOC gene IDs, features, NCBI descriptions, genomic context.
- [`code/32-functional-enrichment.Rmd`](code/32-functional-enrichment.Rmd): Feature & GO enrichment vs covered-CpG/gene background; explicit CYP450/GST/sulfotransferase/ABC xenobiotic check.
- [`code/33-biomarker-candidates.Rmd`](code/33-biomarker-candidates.Rmd): Apply biomarker criteria; site-as-replicate (n = 3 vs 3) sensitivity; rank + shortlist candidates for targeted assay.

### Visualization

- [`code/34-igv-visualization.Rmd`](code/34-igv-visualization.Rmd): Per-sample + DML/DMR bedGraphs; IGV session XML / `igvR` for top loci; summary figures.

## Data availability

- **Raw FASTQs:** [owl.fish.washington.edu/nightingales/M_trossulus/](https://owl.fish.washington.edu/nightingales/M_trossulus/) — download instructions (recursive `wget` and per-file commands) are in [`data/raw-wgbs/README.md`](data/raw-wgbs/README.md). Large intermediates (BAMs, `.cov` files) are kept off-repo and pulled in-script.
- **Reference genome:** *M. trossulus* NCBI assembly [`GCF_036588685.1`](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_036588685.1/) (PNRI_Mtr1.1.1.hap1) with the corresponding `genomic.gff` annotation; see [`data/README.md`](data/README.md).

## Reproducibility & tools

Key tools (pin versions in each script): **FastQC** v0.12.1, **MultiQC** v1.12, **fastp** v1.3.2, **Bismark** v0.24–0.25.1 with **Bowtie2** v2.5.5, **SAMtools** v1.x, **methylKit** (R ≥ 4.4.2), and **BEDTools** v2.30.0. Each analysis is an R Markdown template following the `NN.NN-topic.Rmd` convention, with outputs written to a 1:1 matching `output/<name>/` directory. Scripts record `sessionInfo()`, set seeds where stochastic, and `saveRDS` key objects so downstream steps are resumable. See [`code/README.md`](code/README.md) and [`output/README.md`](output/README.md) for the full conventions.
