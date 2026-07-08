#!/usr/bin/env bash
# 25.2 - Run the Goal 1 landscape analysis end-to-end (scripts 26 -> 29).
#
# Stages inputs from gannet/NCBI, then knits each Rmd in order. The source Rmds
# keep their eval=FALSE default; this driver knits an eval=TRUE TEMP COPY of each
# so the originals are never modified (safe against interruption). Real products
# (.rds/.csv/.png) are written to the matching output/ dirs by the chunks; the
# knitted <script>.md is a run log.
set -euo pipefail
cd "$(dirname "$0")"                    # -> code/

RMDS=(26-methylkit-object 27-methylation-landscape \
      28-genomic-feature-distribution 29-cpg-oe-mosaic)
COV_DIR=../output/25-meth-extract
GENOME_DIR=../output/22-genome-prep
GANNET="https://gannet.fish.washington.edu/seashell/bu-github/project-mytilus-methylation/data"
NCBI_GFF="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/036/588/685/GCF_036588685.1_PNRI_Mtr1.1.1.hap1/GCF_036588685.1_PNRI_Mtr1.1.1.hap1_genomic.gff.gz"

log(){ echo "[$(date '+%F %T')] $*"; }

# ---------------------------------------------------------------------------
# 1. Stage inputs up front so any fetch failure surfaces before the long run
# ---------------------------------------------------------------------------
mkdir -p "$COV_DIR" "$GENOME_DIR"

for s in 105M 106M 107M 109M 239M 240M 241M 242M 269M 270M 271M 272M \
         69M 70M 71M 72M 78M 79M 80M 81M 92M 94M 95M; do
  f="$COV_DIR/${s}.CpG_report.merged_CpG_evidence.cov"
  if [ ! -s "$f" ]; then log "downloading cov ${s}"; wget -q -O "$f" "$GANNET/${s}.CpG_report.merged_CpG_evidence.cov"; fi
done
# 93M is local (25.1); confirm it is present
[ -s "$COV_DIR/93M.CpG_report.merged_CpG_evidence.cov" ] || { log "ERROR: 93M cov missing"; exit 1; }

if [ ! -s "$GENOME_DIR/genomic.gff" ]; then
  log "downloading annotation GFF"
  wget -q -O "$GENOME_DIR/genomic.gff.gz" "$NCBI_GFF" && gunzip -f "$GENOME_DIR/genomic.gff.gz"
fi

n_cov=$(ls -1 "$COV_DIR"/*.CpG_report.merged_CpG_evidence.cov | wc -l)
log "cov files staged: $n_cov (expect 24); gff lines: $(wc -l < "$GENOME_DIR/genomic.gff")"
[ "$n_cov" -eq 24 ] || { log "ERROR: expected 24 cov files, found $n_cov"; exit 1; }

# ---------------------------------------------------------------------------
# 2. Knit each script in order (eval-on temp copy; original untouched)
# ---------------------------------------------------------------------------
for r in "${RMDS[@]}"; do
  tmp="_run_${r}.Rmd"
  sed 's/  eval = FALSE,/  eval = TRUE,/' "${r}.Rmd" > "$tmp"
  log "knitting ${r} ..."
  if Rscript -e "knitr::knit('${tmp}', output='${r}.md')"; then
    rm -f "$tmp"
    log "  -> ${r}.md done"
  else
    rm -f "$tmp"
    log "ERROR knitting ${r}; see console output above"
    exit 1
  fi
done

log "Goal 1 landscape analysis complete."
echo "Outputs:"
ls -1 ../output/26-methylkit-object ../output/27-methylation-landscape \
      ../output/28-genomic-feature-distribution ../output/29-cpg-oe-mosaic 2>/dev/null
