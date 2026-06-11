#!/usr/bin/env Rscript
# Audit of DML -> gene assignments for paper.md Table 2.
# Re-derives the mapping from scratch by intersecting all 16 methylKit DMLs
# against the NCBI GFF gene features (GCF_036588685.1).

suppressMessages({library(data.table)})

out_dir <- "../output/35-dml-gene-model-figure"

## ---- 1. Load all 16 DMLs from the methylKit headline set --------------------
tab <- fread("../output/11-methylkit-klone/myDiff_1055p.tab")
setnames(tab, c("rowid","chr","start","end","strand","pvalue","qvalue","meth.diff"))
# methylKit start is the 1-based CpG position
tab[, pos := start]
tab[, direction := ifelse(`meth.diff` >= 0, "Hyper", "Hypo")]

## ---- 2. Load gene features --------------------------------------------------
g <- fread(file.path(out_dir, "genes_all.gff"), sep = "\t", header = FALSE,
           quote = "", col.names = c("chr","src","type","start","end",
                                      "score","strand","frame","attr"))
get_attr <- function(attr, key) {
  out <- rep(NA_character_, length(attr))
  m   <- regexpr(paste0(key, "=[^;]+"), attr)
  hit <- m > 0
  out[hit] <- sub(paste0(key, "="), "", regmatches(attr, m))
  out
}
g[, gene := get_attr(attr, "gene")]
g[, description := get_attr(attr, "description")]
genes <- g[, .(chr, start, end, strand, gene, description)]

## ---- 3. For each DML: containing gene + nearest up/down genes ---------------
audit <- rbindlist(lapply(seq_len(nrow(tab)), function(i) {
  d <- tab[i]
  gc <- genes[chr == d$chr]
  inside <- gc[d$pos >= start & d$pos <= end]

  if (nrow(inside) > 0) {
    # if multiple (overlapping genes), take the one with smallest span
    inside[, span := end - start]
    setorder(inside, span)
    hit <- inside[1]
    data.table(
      chr = d$chr, pos = d$pos, qvalue = d$qvalue, meth.diff = d$`meth.diff`,
      direction = d$direction,
      context = "genic",
      containing_gene = hit$gene, containing_desc = hit$description,
      gene_span = paste0(hit$start, "-", hit$end), gene_strand = hit$strand,
      dist_to_gene = 0,
      up_gene = NA_character_, up_desc = NA_character_, up_dist = NA_integer_,
      down_gene = NA_character_, down_desc = NA_character_, down_dist = NA_integer_
    )
  } else {
    # intergenic: nearest gene whose END is < pos (upstream by coordinate)
    up   <- gc[end < d$pos][order(-end)][1]
    down <- gc[start > d$pos][order(start)][1]
    up_dist   <- if (nrow(up)   && !is.na(up$end))     d$pos - up$end       else NA_integer_
    down_dist <- if (nrow(down) && !is.na(down$start)) down$start - d$pos   else NA_integer_
    nearest_dist <- min(c(up_dist, down_dist), na.rm = TRUE)
    data.table(
      chr = d$chr, pos = d$pos, qvalue = d$qvalue, meth.diff = d$`meth.diff`,
      direction = d$direction,
      context = "intergenic",
      containing_gene = NA_character_, containing_desc = NA_character_,
      gene_span = NA_character_, gene_strand = NA_character_,
      dist_to_gene = nearest_dist,
      up_gene = up$gene, up_desc = up$description, up_dist = up_dist,
      down_gene = down$gene, down_desc = down$description, down_dist = down_dist
    )
  }
}))

setorder(audit, chr, pos)

## ---- 4. Mark which loci the paper's Table 2 used, and its assignment --------
paper_q <- c(5.77e-14, 1.72e-21, 1.16e-11, 2.92e-13, 4.17e-14,
             1.10e-13, 7.20e-13, 8.98e-18, 6.02e-19, 1.06e-12,
             1.32e-13, 2.52e-16)
paper_assign <- c("LOC134725187","LOC134725187","LOC134690221","LOC134714536",
                  "LOC134695637","LOC134707074","LOC134720671","LOC134711256",
                  "LOC134711256","LOC134685297","LOC134712818","LOC134724475")
paper_status <- c("Both","Both","Hyper","Hyper","Hypo","Hyper","Hypo",
                  "Both","Both","Hypo","Hyper","Hypo")
pmap <- data.table(qvalue = paper_q, paper_gene = paper_assign,
                   paper_status = paper_status)

# match by nearest qvalue (Table 2 q-values are rounded)
audit[, paper_gene := NA_character_]
audit[, paper_status := NA_character_]
for (j in seq_len(nrow(pmap))) {
  k <- which.min(abs(audit$qvalue - pmap$qvalue[j]))
  audit[k, paper_gene := pmap$paper_gene[j]]
  audit[k, paper_status := pmap$paper_status[j]]
}
audit[, in_paper_table2 := !is.na(paper_gene)]
audit[, paper_assignment_correct := in_paper_table2 &
        !is.na(containing_gene) & containing_gene == paper_gene]

## ---- 4b. Feature-level context (CDS / exon-UTR / intron / intergenic) -------
ov <- fread(file.path(out_dir, "dml_overlap_features.tsv"), header = FALSE,
            sep = "\t", quote = "",
            col.names = c("locus","type","start","end","strand","gene","desc","prod"))
feat_ctx <- ov[, .(
  has_exon = any(type == "exon"),
  has_cds  = any(type == "CDS"),
  has_gene = any(type == "gene")
), by = locus]
feat_ctx[, feature_context := fifelse(has_cds, "CDS (coding exon)",
                              fifelse(has_exon, "exon (UTR)",
                              fifelse(has_gene, "intron", "intergenic")))]
audit[, locus := paste0(chr, ":", pos)]
audit <- merge(audit, feat_ctx[, .(locus, feature_context)],
               by = "locus", all.x = TRUE, sort = FALSE)
audit[context == "intergenic", feature_context := "intergenic"]
setorder(audit, chr, pos)

fwrite(audit, file.path(out_dir, "dml_gene_assignment_audit.csv"))

## ---- 5. Console report ------------------------------------------------------
cat("\n================ FULL 16-DML AUDIT ================\n")
print(audit[, .(chr, pos, qvalue, meth.diff, direction, context,
                containing_gene, dist_to_gene, in_paper_table2,
                paper_gene, paper_assignment_correct)])

cat("\n---- Intergenic loci (nearest up/down genes) ----\n")
print(audit[context == "intergenic",
            .(chr, pos, meth.diff, up_gene, up_dist, up_desc,
              down_gene, down_dist, down_desc)])

genic <- audit[context == "genic"]
cat("\n================ CORRECTED SUMMARY (genic only) ================\n")
cat("Genic DMLs:", nrow(genic), " across unique genes:",
    uniqueN(genic$containing_gene), "\n")
cat("Hyper:", genic[direction=="Hyper", .N],
    " Hypo:", genic[direction=="Hypo", .N], "\n")
gene_n <- genic[, .(n_dml = .N,
                    dirs = paste(sort(unique(direction)), collapse="/")),
                by = .(containing_gene, containing_desc)]
cat("Genes with >1 DML (true bidirectional candidates):\n")
print(gene_n[n_dml > 1])
cat("\nAll genic gene -> DML count:\n")
print(gene_n[order(-n_dml)])

cat("\n================ INTERGENIC = the 4 excluded ================\n")
print(audit[context=="intergenic", .(chr, pos, qvalue, meth.diff,
                                      in_paper_table2)])
