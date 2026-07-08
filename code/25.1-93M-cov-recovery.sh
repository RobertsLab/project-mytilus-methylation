#!/usr/bin/env bash
# 25.1 - Full-pipeline recovery of sample 93M's merged-CpG coverage file.
#
# WHY: 93M was never fully aligned. The gannet mirror holds only 93M's trimmed
# FASTQs; the repo's output/04-bismark-pipeline/93M* is a 10,000-read TEST run
# (see "-u 10000" in 05-bismark-cov.Rmd), which yields a ~17k-CpG cov ~1000x
# smaller than the real ~795 MB files. This script reproduces the production
# pipeline (scripts 22->25) for 93M only, matching every other sample's params:
#   bismark_genome_preparation
#   bismark (-score_min L,0,-0.6 --non_directional, paired)
#   deduplicate_bismark (--paired)
#   bismark_methylation_extractor
#   coverage2cytosine --merge_CpG --zero_based
#
# Idempotent: each step is skipped if its output already exists, so the job is
# resumable after interruption.
set -euo pipefail

BISMARK=/home/shared/Bismark-0.24.0
BOWTIE2=/home/shared/bowtie2-2.4.4-linux-x86_64
SAMTOOLS=/home/shared/samtools-1.12/samtools

ROOT="$(cd "$(dirname "$0")/.." && pwd)"          # repo root
GENOME_DIR="$ROOT/output/22-genome-prep"
ALIGN_DIR="$ROOT/output/23-bismark-align"
DEDUP_DIR="$ROOT/output/24-dedup-sort"
EXTRACT_DIR="$ROOT/output/25-meth-extract"
WORK="$ROOT/output/25.1-93M-cov-recovery"
READS="$WORK/reads"

GANNET="https://gannet.fish.washington.edu/seashell/bu-github/project-mytilus-methylation/data"
SAMPLE=93M
SCORE_MIN="L,0,-0.6"

mkdir -p "$GENOME_DIR" "$ALIGN_DIR" "$DEDUP_DIR" "$EXTRACT_DIR" "$READS"

log(){ echo "[$(date '+%F %T')] $*"; }

# ---------------------------------------------------------------------------
# 1. Trimmed reads (4.5 GB each from gannet)
# ---------------------------------------------------------------------------
R1="$READS/${SAMPLE}_R1.fastp-trim.fq.gz"
R2="$READS/${SAMPLE}_R2.fastp-trim.fq.gz"
for pair in R1 R2; do
  f="$READS/${SAMPLE}_${pair}.fastp-trim.fq.gz"
  if [ ! -s "$f" ]; then
    log "Downloading ${SAMPLE}_${pair} ..."
    wget -q -O "$f" "$GANNET/${SAMPLE}_${pair}.fastp-trim.fq.gz"
  fi
  # integrity: gzip must be valid and end cleanly
  gzip -t "$f"
done
log "Reads ready: $(du -h "$R1" "$R2" | awk '{print $1}' | tr '\n' ' ')"

# ---------------------------------------------------------------------------
# 2. Bismark genome preparation (bisulfite-converted Bowtie2 index)
# ---------------------------------------------------------------------------
# A COMPLETE bowtie2 index needs the reverse files too (.rev.1.bt2/.rev.2.bt2);
# a dir that only has .1-.4.bt2 is a half-built index (build was interrupted).
CT_REV="$GENOME_DIR/Bisulfite_Genome/CT_conversion/BS_CT.rev.1.bt2"
GA_REV="$GENOME_DIR/Bisulfite_Genome/GA_conversion/BS_GA.rev.1.bt2"
if [ ! -s "$CT_REV" ] || [ ! -s "$GA_REV" ]; then
  log "Bisulfite index missing/incomplete - (re)building from scratch ..."
  rm -rf "$GENOME_DIR/Bisulfite_Genome"
  "$BISMARK/bismark_genome_preparation" \
    --path_to_aligner "$BOWTIE2" --verbose "$GENOME_DIR/"
else
  log "Complete Bisulfite_Genome present - skipping prep."
fi

# ---------------------------------------------------------------------------
# 3. Bismark alignment (paired, non-directional) - the heavy step
# ---------------------------------------------------------------------------
ALN_BAM="$ALIGN_DIR/${SAMPLE}_pe.bam"
if [ ! -s "$ALN_BAM" ]; then
  log "Aligning ${SAMPLE} with Bismark/Bowtie2 ..."
  # NB: Bismark 0.24 forbids --basename together with --multicore. For a
  # non-directional run Bismark spawns 4 concurrent bowtie2 instances, so
  # -p 12 already saturates the 48-core box without --multicore.
  "$BISMARK/bismark" \
    --path_to_bowtie2 "$BOWTIE2" \
    -genome "$GENOME_DIR/" \
    -p 12 \
    -score_min "$SCORE_MIN" \
    --non_directional \
    -1 "$R1" -2 "$R2" \
    -o "$ALIGN_DIR" \
    --basename "${SAMPLE}" \
    > "$ALIGN_DIR/${SAMPLE}_stdout.log" 2> "$ALIGN_DIR/${SAMPLE}_stderr.log"
fi
# fail loudly if the aligner did not actually produce the BAM
[ -s "$ALN_BAM" ] || { log "ERROR: alignment produced no BAM; see ${SAMPLE}_stderr.log"; exit 1; }

# ---------------------------------------------------------------------------
# 4. Deduplicate (paired)
# ---------------------------------------------------------------------------
DEDUP_BAM="$DEDUP_DIR/${SAMPLE}_pe.deduplicated.bam"
if [ ! -s "$DEDUP_BAM" ]; then
  log "Deduplicating ${SAMPLE} ..."
  "$BISMARK/deduplicate_bismark" --bam --paired \
    --output_dir "$DEDUP_DIR" "$ALN_BAM"
else
  log "Deduplicated BAM exists - skipping."
fi

# ---------------------------------------------------------------------------
# 5. Methylation extraction -> bedGraph/.cov
# ---------------------------------------------------------------------------
COVGZ="$EXTRACT_DIR/${SAMPLE}_pe.deduplicated.bismark.cov.gz"
if [ ! -s "$COVGZ" ]; then
  log "Extracting methylation calls ..."
  "$BISMARK/bismark_methylation_extractor" \
    --bedGraph --counts --comprehensive --merge_non_CpG \
    --multicore 24 --buffer_size 50% \
    --output "$EXTRACT_DIR" "$DEDUP_BAM"
else
  log "bismark.cov.gz exists - skipping extraction."
fi

# ---------------------------------------------------------------------------
# 6. coverage2cytosine --merge_CpG  -> canonical methylKit input
# ---------------------------------------------------------------------------
MERGED="$EXTRACT_DIR/${SAMPLE}.CpG_report.merged_CpG_evidence.cov"
if [ ! -s "$MERGED" ]; then
  log "Merging CpG strands (coverage2cytosine) ..."
  "$BISMARK/coverage2cytosine" \
    --genome_folder "$GENOME_DIR/" \
    -o "$EXTRACT_DIR/${SAMPLE}" \
    --merge_CpG --zero_based \
    "$COVGZ"
  # drop the multi-GB genome-wide intermediate; keep only the merged cov
  rm -f "$EXTRACT_DIR/${SAMPLE}.CpG_report.txt"
fi

# ---------------------------------------------------------------------------
# 7. Report
# ---------------------------------------------------------------------------
log "DONE. Output: $MERGED"
ls -lh "$MERGED"
echo "records: $(wc -l < "$MERGED")"
awk '{m+=$5;u+=$6} END{printf "global weighted %%CpG meth: %.2f%%  (%d CpG)\n",100*m/(m+u),NR}' "$MERGED"
