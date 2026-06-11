`project-mytilus-methylation/code`

---

Code directory.

All file names should adhere to this format:

00.00-<code_topic>.Rmd

- `00.00`: Numbers to left of decimal should be incremented for each new type of analysis, which is unrelated to other existing analyses. Numbers to the right of the decimal should be incremented for analyses related to existing analyses.

	- E.g. `01.00-exon-counts.Rmd`
	- E.g. `01.01-exon-counts-per-gene.Rmd`
    - E.g. `02.00-hisat2-genome-indexing.Rmd`
    - E.g. `02.01-hisat2-alignments.Rmd`

- `<code_topic>`: Provide a brief description of code functionality. E.g. `exon-counts`.

- `.Rmd`: Most code is anticipated to be run using R Markdown. However, other formats are fine (e.g. Jupyter Notebooks). 

All outputs from a script should be directed to a corresponding directory with the same name (excluding the suffix) in the [`../ouput directory`](../output/) .

---

- [`00.00-WGBSseq-reads-FastQC-MultiQC.Rmd`](./00.00-WGBSseq-reads-FastQC-MultiQC.Rmd): Initial quality assessment of sequencing reads using FastQC and MultiQC.

- [`00.10-trimming-fastp-FastQC-MultiQC.Rmd`](./00.10-trimming-fastp-FastQC-MultiQC.Rmd): Trimming of reads using fastp, followed by assessment with FastQC and MultiQC.

---

## `20+` self-contained end-to-end pipeline

The `20`-series is a complete, standalone workflow that begins from the raw
WGBS reads and runs through differential methylation and biomarker ranking. It
does **not** depend on the outputs of the `00.*`–`14` scripts; each step regenerates
its own inputs. Outputs go to a matching `../output/<name>/` directory.

### Upstream processing (all 24 samples)

- [`20-raw-data-qc.Rmd`](./20-raw-data-qc.Rmd): Download raw paired-end WGBS FastQs from `nightingales/M_trossulus`, verify MD5 checksums, raw-read FastQC + MultiQC.
- [`21-fastp-trimming.Rmd`](./21-fastp-trimming.Rmd): fastp trimming (Illumina adapter auto-detect + 20 bp hard trim on both ends) and post-trim FastQC/MultiQC.
- [`22-genome-prep.Rmd`](./22-genome-prep.Rmd): Download NCBI `GCF_036588685.1` genome + GFF; `samtools faidx`; `bismark_genome_preparation` (bisulfite Bowtie2 index).
- [`23-bismark-align.Rmd`](./23-bismark-align.Rmd): Bismark/Bowtie2 alignment of trimmed reads to the converted genome.
- [`24-dedup-sort.Rmd`](./24-dedup-sort.Rmd): `deduplicate_bismark` then SAMtools coordinate sort + index.
- [`25-meth-extract.Rmd`](./25-meth-extract.Rmd): `bismark_methylation_extractor` + `coverage2cytosine --merge_CpG --zero_based` to produce merged-CpG `.cov` reports.

### methylKit object construction

- [`26-methylkit-object.Rmd`](./26-methylkit-object.Rmd): `methRead` the `cov` files into the 24-sample (landscape, incl. 93M) and 23-sample (differential) sets; `filterByCoverage` 10x primary + 5x sensitivity (`hi.perc = 98`); `unite`; save `.rds`.

### Goal 1 — methylation landscape

- [`27-methylation-landscape.Rmd`](./27-methylation-landscape.Rmd): Per-sample global %CpG methylation + coverage stats, bimodality plot, CpG class proportions, correlation/clustering/PCA, 93M outlier decision.
- [`28-genomic-feature-distribution.Rmd`](./28-genomic-feature-distribution.Rmd): Build gene/exon/intron/promoter/intergenic BEDs from the GFF; `bedtools intersect` covered CpGs × features; feature × methylation barplot (gene-body enrichment).
- [`29-cpg-oe-mosaic.Rmd`](./29-cpg-oe-mosaic.Rmd): Per-gene CpG O/E from the reference FASTA; bimodal O/E distribution; O/E vs observed methylation.

### Goal 2 — high vs low PAH biomarkers

- [`30-dml-dmr-methylkit.Rmd`](./30-dml-dmr-methylkit.Rmd): Per-CpG DML (`calculateDiffMeth`, overdispersion `"MN"` + site covariate); 1 kb tiled DMRs (optional `dmrseq`); thresholds q<0.01 / ≥55% (DML), ≥25% (DMR); BH-FDR + SLIM; per-site concordance; 5x sensitivity.
- [`31-dml-dmr-annotation.Rmd`](./31-dml-dmr-annotation.Rmd): DML/DMR → bedGraph; `bedtools intersect`/`closest` vs `genomic.gff`; recover LOC gene IDs, features, NCBI descriptions, genomic context.
- [`32-functional-enrichment.Rmd`](./32-functional-enrichment.Rmd): Feature & GO enrichment vs covered-CpG/gene background; explicit CYP450/GST/sulfotransferase/ABC xenobiotic check.
- [`33-biomarker-candidates.Rmd`](./33-biomarker-candidates.Rmd): Apply biomarker criteria; site-as-replicate (n=3 vs 3) sensitivity; rank + shortlist candidates for targeted assay.
- [`34-igv-visualization.Rmd`](./34-igv-visualization.Rmd): Per-sample + DML/DMR bedGraphs; IGV session XML / `igvR` for top loci; summary figures.