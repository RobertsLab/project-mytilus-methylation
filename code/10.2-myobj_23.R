library(methylKit)


analysisFiles2 <- list(
  "~/github/project-mytilus-methylation/data/70M_pe.deduplicated.sorted.bam",
  "~/github/project-mytilus-methylation/data/71M_pe.deduplicated.sorted.bam",
  "~/github/project-mytilus-methylation/data/72M_pe.deduplicated.sorted.bam",
  "~/github/project-mytilus-methylation/data/78M_pe.deduplicated.sorted.bam",
  "~/github/project-mytilus-methylation/data/79M_pe.deduplicated.sorted.bam",
  "~/github/project-mytilus-methylation/data/80M_pe.deduplicated.sorted.bam",
  "~/github/project-mytilus-methylation/data/81M_pe.deduplicated.sorted.bam"
) #Put all .bam files into a list for analysis.



myobj_07 = processBismarkAln(location = analysisFiles2, 
                             sample.id = list("70M", "71M", "72M", "78M", "79M", "80M", "81M"),
                             assembly = "Mt",
                             read.context="CpG",
                             mincov=2,
                             treatment = c("EB", "EB", "EB", "SA", "SA", "SA", "SA")
                             )