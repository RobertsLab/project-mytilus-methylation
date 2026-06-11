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
