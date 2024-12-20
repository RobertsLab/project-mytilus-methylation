library(methylKit)


analysisFiles <- list(
  "/home/shared/16TB_HDD_01/sr320/github/project-mytilus-methylation/output/09-meth-quant/105M_pe.deduplicated.sorted.bam",
  "/home/shared/16TB_HDD_01/sr320/github/project-mytilus-methylation/output/09-meth-quant/106M_pe.deduplicated.sorted.bam",
  "/home/shared/16TB_HDD_01/sr320/github/project-mytilus-methylation/output/09-meth-quant/107M_pe.deduplicated.sorted.bam",
  "/home/shared/16TB_HDD_01/sr320/github/project-mytilus-methylation/output/09-meth-quant/109M_pe.deduplicated.sorted.bam",
  "/home/shared/16TB_HDD_01/sr320/github/project-mytilus-methylation/output/09-meth-quant/239M_pe.deduplicated.sorted.bam",
  "/home/shared/16TB_HDD_01/sr320/github/project-mytilus-methylation/output/09-meth-quant/240M_pe.deduplicated.sorted.bam",
  "/home/shared/16TB_HDD_01/sr320/github/project-mytilus-methylation/output/09-meth-quant/241M_pe.deduplicated.sorted.bam",
  "/home/shared/16TB_HDD_01/sr320/github/project-mytilus-methylation/output/09-meth-quant/242M_pe.deduplicated.sorted.bam",
  "/home/shared/16TB_HDD_01/sr320/github/project-mytilus-methylation/output/09-meth-quant/269M_pe.deduplicated.sorted.bam",
  "/home/shared/16TB_HDD_01/sr320/github/project-mytilus-methylation/output/09-meth-quant/270M_pe.deduplicated.sorted.bam",
  "/home/shared/16TB_HDD_01/sr320/github/project-mytilus-methylation/output/09-meth-quant/271M_pe.deduplicated.sorted.bam",
  "/home/shared/16TB_HDD_01/sr320/github/project-mytilus-methylation/output/09-meth-quant/272M_pe.deduplicated.sorted.bam",
  "/home/shared/16TB_HDD_01/sr320/github/project-mytilus-methylation/output/09-meth-quant/69M_pe.deduplicated.sorted.bam",
  "/home/shared/16TB_HDD_01/sr320/github/project-mytilus-methylation/output/09-meth-quant/70M_pe.deduplicated.sorted.bam",
  "/home/shared/16TB_HDD_01/sr320/github/project-mytilus-methylation/output/09-meth-quant/71M_pe.deduplicated.sorted.bam",
  "/home/shared/16TB_HDD_01/sr320/github/project-mytilus-methylation/output/09-meth-quant/72M_pe.deduplicated.sorted.bam",
  "/home/shared/16TB_HDD_01/sr320/github/project-mytilus-methylation/output/09-meth-quant/78M_pe.deduplicated.sorted.bam",
  "/home/shared/16TB_HDD_01/sr320/github/project-mytilus-methylation/output/09-meth-quant/79M_pe.deduplicated.sorted.bam",
  "/home/shared/16TB_HDD_01/sr320/github/project-mytilus-methylation/output/09-meth-quant/80M_pe.deduplicated.sorted.bam",
  "/home/shared/16TB_HDD_01/sr320/github/project-mytilus-methylation/output/09-meth-quant/81M_pe.deduplicated.sorted.bam",
  "/home/shared/16TB_HDD_01/sr320/github/project-mytilus-methylation/output/09-meth-quant/92M_pe.deduplicated.sorted.bam",
  "/home/shared/16TB_HDD_01/sr320/github/project-mytilus-methylation/output/09-meth-quant/94M_pe.deduplicated.sorted.bam",
  "/home/shared/16TB_HDD_01/sr320/github/project-mytilus-methylation/output/09-meth-quant/95M_pe.deduplicated.sorted.bam"
)



myobj_23 = processBismarkAln(location = analysisFiles, sample.id = list("105M", "106M", "107M", "109M", "239M", "240M", "241M", "242M", "269M", "270M", "271M", "272M", "69M", "70M", "71M", "72M", "78M", "79M", "80M", "81M", "92M", "94M", "95M"), assembly = "Mt", read.context="CpG", mincov=2, treatment = c("HC", "HC", "HC", "HC", "AP", "AP", "AP", "AP", "BS", "BS", "BS", "BS", "EB", "EB", "EB", "EB", "SA", "SA", "SA", "SA", "SC", "SC", "SC"))