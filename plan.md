# Project plan: *Mytilus trossulus* DNA methylation

## Overview

This repository analyzes whole-genome bisulfite sequencing (WGBS) of mantle tissue from **24 individual *M. trossulus*** collected during the **2021–2022 Washington Department of Fish and Wildlife (WDFW) Nearshore Monitoring deployment** in Puget Sound. Samples were chosen to contrast **high-** and **low-Σ16 PAH** (sum of 16 priority polycyclic aromatic hydrocarbons) exposure sites: three sites per group, four mussels per site.

The primary analysis compares DNA methylation between exposure groups using Bismark alignment, methylKit differential methylation, and gene annotation (see `paper.md`, `code/10-methylkit.Rmd`, `code/11-methylkit-klone.Rmd`).

## Study design

| Factor | Detail |
|--------|--------|
| Species | *Mytilus trossulus* |
| Tissue | Mantle |
| Sequencing | Paired-end 150 bp WGBS (Illumina NovaSeq) |
| Sites | 6 (3 low PAH, 3 high PAH) |
| Individuals | 24 (4 per site) |
| Exposure metric | Σ16 PAH (ng/g) from WDFW deployment |
| Low PAH range | 97.4–127.26 ng/g |
| High PAH range | 3,304.96–12,078.42 ng/g |

Two additional high-PAH Commencement Bay sites were excluded before extraction due to tissue degradation (`paper.md`). Myrtle Edwards is listed among candidate high-PAH sites in the manuscript but is not represented in the sequenced set.

## Sampling sites

| Code | Site | PAH group | Samples |
|------|------|-----------|---------|
| HC | Hood Canal | Low | 105M, 106M, 107M, 109M |
| AP | Holly Aiston Preserve | Low | 239M, 240M, 241M, 242M |
| BS | Broad Spit | Low | 269M, 270M, 271M, 272M |
| EB | Elliot Bay | High | 69M, 70M, 71M, 72M |
| SA | Seattle Aquarium Pier 59 | High | 78M, 79M, 80M, 81M |
| SC | Smith Cove (Terminal 91) | High | 92M, 93M, 94M, 95M |

Site codes are defined in `code/10-methylkit.Rmd` and `code/11-methylkit-klone.Rmd`. In methylKit, low-PAH sites (HC, AP, BS) are coded as treatment `0` and high-PAH sites (EB, SA, SC) as treatment `1`.

## The 24 samples

Each sample is one mussel (`*M` = mussel ID from the WDFW deployment). All underwent QC (FastQC/MultiQC), trimming (fastp), Bismark alignment, and methylation calling. Raw reads are hosted at [owl.fish.washington.edu/nightingales/M_trossulus/](https://owl.fish.washington.edu/nightingales/M_trossulus/).

### Low-PAH sites (n = 12)

| Sample | Site | Description |
|--------|------|-------------|
| **105M** | Hood Canal (HC) | Individual *M. trossulus* mantle WGBS library from a low-Σ16 PAH Hood Canal deployment site. |
| **106M** | Hood Canal (HC) | Second replicate mussel from Hood Canal; low-PAH exposure group. |
| **107M** | Hood Canal (HC) | Third replicate mussel from Hood Canal; low-PAH exposure group. |
| **109M** | Hood Canal (HC) | Fourth replicate mussel from Hood Canal; low-PAH exposure group. |
| **239M** | Holly Aiston Preserve (AP) | Individual mussel from Holly Aiston Preserve on Hood Canal; low-Σ16 PAH site. |
| **240M** | Holly Aiston Preserve (AP) | Second replicate from Holly Aiston Preserve; low-PAH exposure group. |
| **241M** | Holly Aiston Preserve (AP) | Third replicate from Holly Aiston Preserve; low-PAH exposure group. |
| **242M** | Holly Aiston Preserve (AP) | Fourth replicate from Holly Aiston Preserve; low-PAH exposure group. |
| **269M** | Broad Spit (BS) | Individual mussel from Broad Spit; low-Σ16 PAH reference site in Hood Canal. |
| **270M** | Broad Spit (BS) | Second replicate from Broad Spit; low-PAH exposure group. |
| **271M** | Broad Spit (BS) | Third replicate from Broad Spit; low-PAH exposure group. |
| **272M** | Broad Spit (BS) | Fourth replicate from Broad Spit; low-PAH exposure group. |

### High-PAH sites (n = 12)

| Sample | Site | Description |
|--------|------|-------------|
| **69M** | Elliot Bay (EB) | Individual mussel from Elliot Bay, Seattle; high-Σ16 PAH urban site. |
| **70M** | Elliot Bay (EB) | Second replicate from Elliot Bay; high-PAH exposure group. |
| **71M** | Elliot Bay (EB) | Third replicate from Elliot Bay; high-PAH exposure group. |
| **72M** | Elliot Bay (EB) | Fourth replicate from Elliot Bay; high-PAH exposure group. |
| **78M** | Seattle Aquarium Pier 59 (SA) | Individual mussel from Pier 59 at Seattle Aquarium; high-Σ16 PAH site. |
| **79M** | Seattle Aquarium Pier 59 (SA) | Second replicate from Seattle Aquarium Pier 59; high-PAH exposure group. |
| **80M** | Seattle Aquarium Pier 59 (SA) | Third replicate from Seattle Aquarium Pier 59; high-PAH exposure group. |
| **81M** | Seattle Aquarium Pier 59 (SA) | Fourth replicate from Seattle Aquarium Pier 59; high-PAH exposure group. |
| **92M** | Smith Cove / Terminal 91 (SC) | Individual mussel from Smith Cove (Terminal 91); high-Σ16 PAH industrial harbor site. |
| **93M** | Smith Cove / Terminal 91 (SC) | Second replicate from Smith Cove; high-PAH exposure group. Sequenced and processed through Bismark; **excluded from methylKit differential analysis** (23-sample set in `code/10-methylkit.Rmd`). |
| **94M** | Smith Cove / Terminal 91 (SC) | Third replicate from Smith Cove; high-PAH exposure group. |
| **95M** | Smith Cove / Terminal 91 (SC) | Fourth replicate from Smith Cove; high-PAH exposure group. |

## Pipeline notes

- **QC & trimming:** `code/00.00-WGBSseq-reads-FastQC-MultiQC.Rmd`, `code/00.10-trimming-fastp-FastQC-MultiQC.Rmd`
- **Alignment:** `code/04-bismark-pipeline` — all 24 samples in `output/04-bismark-pipeline/bismark_summary_report.txt`
- **Differential methylation:** `code/11-methylkit-klone.Rmd` — 23 samples (93M omitted); compares high vs low PAH at q < 0.01 and ≥55% methylation difference
- **Reference genome:** *M. trossulus* NCBI GCF_036588685.1

## Data sources

- Raw FASTQs: [nightingales/M_trossulus](https://owl.fish.washington.edu/nightingales/M_trossulus/)
- Sample list in repo: `data/raw-wgbs/README.md`
- Manuscript context: `paper.md`

---

# Analysis plan

This section lays out the analysis from the current intermediate files forward. It assumes the upstream pipeline is complete: FastQC/MultiQC QC (`00.00`, `00.10`), fastp trimming (20 bp hard-trim both ends), Bismark/Bowtie2 alignment + deduplication, methylation extraction, and `coverage2cytosine --merge_CpG --zero_based` (`code/05-bismark-cov.Rmd`, `code/09-meth-quant.Rmd`). The key per-sample inputs going forward are:

- **`*.CpG_report.merged_CpG_evidence.cov`** — per-CpG merged (both strands) coverage files in Bismark `cov` format (`chr  start  end  %meth  count_methylated  count_unmethylated`). These are the primary input for methylKit (`methRead(..., pipeline = "bismarkCoverage")`) and for genome-wide landscape summaries.
- **`*_pe.deduplicated.sorted.bam`** — used by `processBismarkAln()` as an alternative entry point, but the `cov` files are preferred (faster, deterministic, and already CpG-merged).
- **Reference + annotation:** *M. trossulus* `GCF_036588685.1` (`data/GCF_036588685.1_PNRI_Mtr1.1.1.hap1_genomic.fa` and `data/ncbi_dataset/data/GCF_036588685.1/genomic.gff`).

Two analysis goals are addressed below, each as a self-contained workflow that builds on the existing methylKit/BEDTools framework, with justified improvements.

## Sample-set policy (applies to both goals)

| Set | Samples | Use |
|-----|---------|-----|
| **24-sample set** | all individuals incl. **93M** | **Landscape characterization (Goal 1)** — descriptive, per-sample summaries do not require balanced groups, so 93M is retained to maximize the descriptive picture and to *diagnose* why it was dropped. |
| **23-sample set** | 24 minus **93M** | **Differential methylation / biomarker discovery (Goal 2)** — matches `code/11-methylkit-klone.Rmd`. 93M is excluded because it behaved as a coverage/quality outlier; keeping it out preserves a balanced-ish design (HC/AP/BS = 4/4/4 low; EB/SA = 4/4, SC = 3 high). |

**Recommendation:** run the global landscape on all 24 first, formally evaluate 93M against the QC criteria below (coverage, global %CpG methylation, correlation/PCA position), and document the exclusion decision with a figure rather than carrying it as an unstated convention. If 93M passes QC, reconsider adding it back to the 23-sample comparison as a sensitivity analysis.

## Goal 1 — Describing the *M. trossulus* methylation landscape

Objective: a genome-wide, descriptive characterization of CpG methylation independent of treatment, establishing the expected bivalve mosaic / gene-body–enriched pattern and the sample-level QC baseline.

### 1.1 Per-sample global CpG methylation summary
- Read all 24 `*.CpG_report.merged_CpG_evidence.cov` files into methylKit (`methRead`, `pipeline = "bismarkCoverage"`, `context = "CpG"`, `mincov = 1` for the descriptive pass so low-coverage sites are visible before filtering).
- Report per-sample: number of CpGs covered, mean/median coverage, and **global weighted % CpG methylation** (Σ methylated counts / Σ total counts). The manuscript reports a global mean of ~12% (range 10.3–13%); confirm and tabulate per sample.
- Produce per-sample methylation histograms (`getMethylationStats(..., plot = TRUE)`) and coverage histograms (`getCoverageStats(..., plot = TRUE)`).

### 1.2 Methylation-level distribution & CpG classification
- Bivalve methylomes are **bimodal/mosaic**: most CpGs are near-0% or near-100% methylated. Plot the genome-wide distribution of per-CpG % methylation (pooled and per sample) to show bimodality.
- Classify CpGs (on sufficiently covered sites, ≥10x; see 1.5) into three bins, a standard invertebrate-methylome scheme:

| Class | % methylation | Expected (bivalve) |
|-------|---------------|--------------------|
| Methylated | ≥ 50% (or ≥ 75% strict) | Minority; concentrated in gene bodies |
| Intermediate | 10–50% | Small fraction |
| Unmethylated/sparse | ≤ 10% | Majority of genome |

- Report the proportion of CpGs in each class genome-wide and per sample; expect a small methylated fraction dominated by gene bodies.

### 1.3 Genomic-context distribution (feature enrichment)
- Build feature BED/GFF subsets from `GCF_036588685.1` `genomic.gff` (the `code/12-IGV-prep.Rmd` feature-splitting loop already produces per-feature GFFs: `gene`, `mRNA`, `exon`, `CDS`, `lnc_RNA`, `ncRNA`, etc.). Derive additional contexts not directly in the GFF:
  - **Introns** = gene minus exon (BEDTools `subtract`).
  - **Promoters/TSS** = ±1–2 kb flanks of gene/mRNA starts (BEDTools `flank`/`slop` against a genome `.fai`; the FASTA is already indexed in `12-IGV-prep`).
  - **Intergenic** = genome complement of genes (BEDTools `complement`).
  - **Repeats** = if a RepeatMasker/`.out` track is available; otherwise note as a gap.
- Convert covered CpG positions (≥10x) to BED and **`bedtools intersect`** against each feature class. Report the fraction of CpGs, and the mean methylation level, per feature.
- **Expected result to highlight:** gene-body (exon/CDS/intron) enrichment of methylation with low intergenic and low promoter methylation — the canonical molluscan mosaic pattern (Männer et al. 2021; Roberts & Gavery 2012). Present as a feature × methylation-level barplot.

### 1.4 CpG density / CpG O/E (mosaic structure)
- Compute **CpG observed/expected (CpG O/E = (CpG / (C × G)) × (L² / (L−1)))** per gene/CDS from the reference FASTA to characterize the mosaic genome architecture. Bivalve gene sets typically show a **bimodal CpG O/E** distribution: a low-O/E mode (historically germline-methylated, depleted CpG via deamination) and a high-O/E mode (sparsely methylated).
- Relate CpG O/E to observed methylation (genes in the low-O/E mode should be the more heavily methylated, constitutively expressed gene bodies). This links the static genome signature to the empirical methylome and frames which genes are "methylation-prone."

### 1.5 Sample-level QC of the landscape
- **Coverage threshold — recommendation: 10× minimum per sample as primary, with a 5× sensitivity pass.** The manuscript and `code/10` used 5×; `code/11-methylkit-klone.Rmd` actually filters at `lo.count = 10`. Rationale for 10×: per-CpG methylation is a binomial proportion whose variance scales with 1/coverage, so 5× estimates are noisy (e.g., 1/5 vs 0/5 swings the estimate by 20 points). 10× is the common invertebrate-WGBS standard and materially reduces false DM calls; given mean ~38.5M deduplicated alignments/sample, 10× retains ample CpGs. Report the CpG-retention trade-off (sites surviving at 1×/5×/10×) so the choice is transparent.
- Also cap extreme coverage (PCR/duplication artefacts) with `filterByCoverage(..., hi.perc = 98)` as already done in `code/11`.
- **Global structure:** on the united methylation matrix (sites covered in all/most samples), compute and plot:
  - pairwise sample **correlation** (`getCorrelation`),
  - hierarchical **clustering** (`clusterSamples`, ward/correlation),
  - **PCA** (`PCASamples`).
- Use these to assess whether samples cluster by **site** or **PAH group** vs. by technical/coverage artefacts, and to formally evaluate **93M** as an outlier. Flag any individual that is a clustering/coverage outlier and decide retention for Goal 2.

### 1.6 Goal 1 deliverables
- Per-sample global methylation + coverage table; bimodality plot; CpG class proportions; feature-distribution barplot; CpG O/E distribution; correlation/clustering/PCA figures; documented 93M decision.

## Goal 2 — High vs low PAH differential methylation & biomarker discovery

Objective: identify CpGs (and regions) whose methylation differs between **low-PAH (HC, AP, BS = treatment 0)** and **high-PAH (EB, SA, SC = treatment 1)** mussels, annotate them, and frame robust candidates as biomarkers.

### 2.1 Object construction & filtering (23-sample set)
- `methRead` the 23 `cov` files with the treatment vector `c(0,0,0,0, 0,0,0,0, 0,0,0,0, 1,1,1,1, 1,1,1,1, 1,1,1)` (HC,AP,BS = 0; EB,SA,SC = 1), as in `code/11-methylkit-klone.Rmd`.
- Filter: `filterByCoverage(lo.count = 10, hi.perc = 98)` then `unite(min.per.group = 5L, destrand = FALSE)` (cov files are already CpG-merged, so `destrand = FALSE` is correct — destranding only applies to per-strand reports). `min.per.group = 5L` requires a CpG to be covered in ≥5 individuals per group, balancing power against missingness.

### 2.2 Per-CpG differential methylation (DML)
- `calculateDiffMeth()` (logistic regression per CpG).
- **Overdispersion:** field bivalves are overdispersed (biological + technical variance > binomial). Use `overdispersion = "MN"` with `test = "Chisq"` — this is a concrete, recommended deviation from the plain logistic-regression default used in the manuscript pipeline, and it directly mitigates inflated significance from pseudoreplication (see 2.5).
- **Significance/effect thresholds:** retain the manuscript's primary definition — **q < 0.01 and |methylation difference| ≥ 55%** — but also export the q<0.01/≥25% and ≥50% sets to characterize the full effect-size spectrum and to avoid over-reliance on a single hard cutoff.
- **Multiple testing:** methylKit defaults to **SLIM**; report **BH/FDR** as well for transparency and comparability with other studies (BH is more conventional and conservative-interpretable). Use BH-q as the headline.
- Split into hyper-/hypo-methylated (relative to low-PAH baseline) with `getMethylDiff(type = "hyper"/"hypo")`.

### 2.3 Differentially methylated regions (DMR) — recommended addition
- **Add a tiled-window DMR analysis in addition to single-CpG DMLs.** Justification: (i) bivalve methylation is organized at the gene-body/regional scale, so regional signal is biologically natural; (ii) DMRs aggregate evidence across neighboring CpGs, increasing robustness to single-CpG noise and to the small per-group n; (iii) regional calls are more portable to targeted biomarker assays.
- Implement with methylKit `tileMethylCounts(win.size = 1000, step.size = 1000, cov.bases = 5)` → `unite` → `calculateDiffMeth` (same overdispersion settings), then `getMethylDiff(difference = 25, qvalue = 0.01)` (region-level effect thresholds are conventionally lower than the 55% single-CpG cutoff).
- Optionally cross-check key regions with `dmrseq` (smoothing + region-level permutation inference), which provides an independent, model-based DMR significance estimate.

### 2.4 Annotation of DMLs/DMRs
- Convert DML/DMR coordinates to BED/bedGraph (the `code/12-IGV-prep.Rmd` → `code/14-intersectbed.Rmd` route already does this: `myDiff_1055p.tab` → `myDiff1055p.bedgraph`).
- **`bedtools intersect`** (BEDTools v2.30.0) DML/DMR coordinates against the `GCF_036588685.1` `genomic.gff` (and the per-feature GFFs from `12-IGV-prep`) to recover overlapping `gene`/`mRNA`/`exon`/`CDS` features and **gene IDs (LOC…)**. Use `bedtools closest` for loci that fall in intergenic space to report the nearest gene and distance.
- Record genomic context (gene body vs exon vs intron vs promoter vs intergenic) for each DML/DMR, mirroring the Goal-1 feature classes.
- Map LOC IDs to NCBI gene descriptions; assign functional categories as in the manuscript (membrane dynamics/vesicle trafficking, cell signaling, stress/apoptosis, cell adhesion, uncharacterized).

### 2.5 Pseudoreplication & design structure (critical)
The design is **4 individuals nested within 3 sites within each PAH group** — individuals are not independent of site, and "PAH exposure" is confounded with "site." Treating 12 vs 11 individuals as independent replicates (the default methylKit setup) is **pseudoreplication** and inflates significance. Recommended handling, in order of preference:

1. **Overdispersion correction** (`overdispersion = "MN"`, `test = "Chisq"`) — already specified in 2.2; partially absorbs the extra-binomial variance from within-/between-site clustering.
2. **Site as covariate:** supply a `covariates` data frame (site factor) to `calculateDiffMeth` so the group test adjusts for site-level differences, isolating the PAH-associated component.
3. **Site-as-replicate sensitivity analysis:** collapse to **site-level methylation means (n = 3 low vs 3 high)** and test at the site level. This is the statistically honest unit of replication for "site/exposure" effects; it will be low-powered but any locus significant here is far more credible as an exposure biomarker.
4. **Consistency filtering:** require a candidate DML/DMR to show the same direction of effect in **all 3 sites** within a group (not driven by a single site/individual). This per-site concordance check is the single most important guard against site-confounded false positives and feeds directly into the biomarker criteria (2.7).

**Honest framing:** with 3 sites/group this study cannot fully separate PAH exposure from site identity; the limitation must be stated explicitly (see Statistical considerations). Candidate loci are *associated with* the high-PAH sites, not proven PAH-causal.

### 2.6 Enrichment analyses (optional but recommended)
- **Gene-feature enrichment:** test whether DMLs/DMRs are over-represented in gene bodies/exons vs the genome-wide CpG background (Goal 1.3 distribution as the expected null; Fisher/χ²). Expect gene-body enrichment.
- **GO / functional enrichment:** map DML/DMR genes to GO terms (via NCBI annotation / orthology to a better-annotated bivalve such as *C. gigas* / *M. galloprovincialis* using eggNOG-mapper or InterProScan if NCBI GO is sparse) and test enrichment against all annotated genes carrying covered CpGs (the proper background — only genes that *could* be detected).
- **Targeted xenobiotic check:** explicitly query whether any DML/DMR maps to **CYP450, GST, sulfotransferase, ABC-transporter** genes. The manuscript found *none* in canonical biotransformation genes; confirm/extend this and report it as a substantive result (consistent with stable methylation of constitutively expressed detox genes in the mosaic-methylation framework).

### 2.7 Biomarker framing
A candidate locus is promoted to **biomarker candidate** only if it meets:

| Criterion | Operational definition |
|-----------|------------------------|
| **Effect size** | \|methylation difference\| ≥ 55% (DML) / ≥ 25% (DMR), q < 0.01 |
| **Consistency** | same direction in all 3 sites of the group; not driven by 1 individual; survives site-as-covariate model |
| **Coverage robustness** | ≥10× in ≥5 individuals/group; effect persists in 5× sensitivity pass |
| **Biological plausibility** | maps to an annotated gene with a coherent functional role (or a reproducible region) |
| **Assayability** | locus/region tractable for targeted bisulfite sequencing or MS-PCR |

- **Validation outlook:** shortlist the most robust loci/regions for an independent, cost-effective **targeted assay** (targeted bisulfite amplicon sequencing or methylation-sensitive PCR) deployable across the broader Nearshore Monitoring sample set — the manuscript's stated end goal.

### 2.8 Goal 2 deliverables
- DML table (with q, methylation difference, direction, per-site concordance); DMR table; annotated gene table with functional categories; feature/GO enrichment results; xenobiotic-gene assessment; ranked biomarker-candidate shortlist; IGV session(s) for the top loci (extending `code/13-IGV-klone.Rmd`).

## Proposed analysis steps

New scripts follow the repo convention `NN.NN-topic.Rmd` with a matching `output/NN.NN-topic/` directory (`code/README.md`, `output/README.md`). They continue from the highest existing number (`14-intersectbed.Rmd`).

| # | New script | Goal | Purpose | Output dir |
|---|-----------|------|---------|------------|
| ☐ 15 | `15-methylation-landscape.Rmd` | 1 | Load 24 `cov` files; per-sample global %CpG methylation, coverage stats, bimodality plots, CpG class proportions; correlation/clustering/PCA; **93M outlier decision** | `output/15-methylation-landscape/` |
| ☐ 16 | `16-genomic-feature-distribution.Rmd` | 1 | Build intron/promoter/intergenic BEDs from `genomic.gff`; `bedtools intersect` covered CpGs × features; feature × methylation-level summary (gene-body enrichment) | `output/16-genomic-feature-distribution/` |
| ☐ 17 | `17-cpg-oe-mosaic.Rmd` | 1 | Per-gene/CDS CpG O/E from reference FASTA; bimodal O/E distribution; O/E vs observed methylation | `output/17-cpg-oe-mosaic/` |
| ☐ 18 | `18-dml-dmr-methylkit.Rmd` | 2 | 23-sample filter (10×, hi.perc 98); DML with overdispersion + site covariate; tiled DMRs; thresholds q<0.01 / ≥55% (DML), ≥25% (DMR); BH + SLIM; per-site concordance | `output/18-dml-dmr-methylkit/` |
| ☐ 19 | `19-dml-dmr-annotation.Rmd` | 2 | DML/DMR → bedGraph; `bedtools intersect`/`closest` vs `genomic.gff`; recover LOC gene IDs + features + NCBI descriptions | `output/19-dml-dmr-annotation/` |
| ☐ 20 | `20-functional-enrichment.Rmd` | 2 | Feature & GO enrichment vs covered-CpG/gene background; explicit CYP450/GST/xenobiotic check | `output/20-functional-enrichment/` |
| ☐ 21 | `21-biomarker-candidates.Rmd` | 2 | Apply biomarker criteria; rank candidates; site-as-replicate sensitivity; shortlist for targeted assay; IGV sessions for top loci | `output/21-biomarker-candidates/` |

Existing `code/10`–`code/14` (methylKit DML at q<0.01/≥55% and BEDTools intersect) remain the reference implementation; the new scripts above generalize and harden them (overdispersion, site covariate, DMRs, landscape/QC, enrichment, biomarker ranking).

## Statistical considerations & limitations

- **Pseudoreplication / site–exposure confound (primary limitation):** PAH "treatment" is confounded with site identity (3 sites/group, 4 individuals/site). Individual-level methylKit tests overstate significance; mitigate with overdispersion correction, site as covariate, per-site concordance filtering, and a site-level (n=3 vs 3) sensitivity analysis. Report associations as exposure-*correlated*, not PAH-causal.
- **Unbalanced 23-sample design:** dropping 93M leaves SC with 3 individuals (others 4). Acceptable with `min.per.group` filtering but note the slight imbalance.
- **Coverage sensitivity:** 10× primary vs 5× sensitivity; report how DML/DMR counts depend on the threshold so conclusions are not artifacts of filtering.
- **Multiple testing:** report BH-FDR alongside methylKit's SLIM; the ≥55% effect-size filter is the main guard against trivial-but-significant calls.
- **Overdispersion is partial:** Methods-of-moments (`"MN"`) correction is not a full random-effects model; `dmrseq`/`PQLseq`-style mixed models (site as random effect) are the more rigorous alternative and are recommended for confirmation of top candidates.
- **Annotation completeness:** `GCF_036588685.1` is a recent assembly; many loci are uncharacterized `LOC…` genes. GO coverage may be sparse — supplement with orthology-based annotation, and treat functional categories as hypotheses.
- **Tissue & single time point:** mantle, one snapshot. Methylation differences reflect mantle cellular state, not necessarily digestive-gland biotransformation; generalization to other tissues/seasons is untested.
- **No matched expression data:** methylation–transcription coupling cannot be tested here; gene-body methylation changes are interpreted within the invertebrate framework but not validated against expression.

## Reproducibility notes

- **Tool versions (pin in each script):** FastQC v0.12.1, MultiQC v1.12, fastp v1.3.2, Bismark v0.24–0.25.1 (+ Bowtie2 v2.5.5), SAMtools v1.x, methylKit (R ≥ 4.4.2), BEDTools v2.30.0. Record `sessionInfo()` in every `.Rmd`.
- **Naming convention:** `NN.NN-topic.Rmd` in `code/` with a 1:1 matching `output/NN.NN-topic/` directory (per `code/README.md` / `output/README.md`).
- **Inputs:** standardize on the merged-CpG `cov` files (`*.CpG_report.merged_CpG_evidence.cov`) as the canonical methylKit input; keep large intermediates (BAMs, cov files) off-repo on gannet/owl as currently done and pull via `wget` in-script.
- **Determinism:** set seeds where stochastic (e.g., `dmrseq` permutations), save key R objects (`saveRDS` of `methylRawList`/`methylDiff`) so downstream steps are resumable, as already done for `myobj_23.rds`.
