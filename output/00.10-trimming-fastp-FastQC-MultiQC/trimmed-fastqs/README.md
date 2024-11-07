`project-mytilus-methylation/output/00.10-trimming-fastp-FastQC-MultiQC/trimmed-fastqs`

Output from [00.10-trimming-fastp-FastQC-MultiQC.Rmd](../../code/00.10-trimming-fastp-FastQC-MultiQC.Rmd).

Due to size, trimmed FastQ files are not included in this repo. They can be downloaded from here:

- [https://gannet.fish.washington.edu/Atumefaciens/project-mytilus-methylation/output/00.10-trimming-fastp-FastQC-MultiQC/trimmed-fastqs/](https://gannet.fish.washington.edu/Atumefaciens/project-mytilus-methylation/output/00.10-trimming-fastp-FastQC-MultiQC/trimmed-fastqs/).

  - `*fastp-trim.fq.gz`
  - `*fastp-trim.fq.gz.md5`

---

- `*fastqc.html`: Individual FastQC HTML files.

- `*.fastp-trim.report.*`: Report files (in HTML and JSON) produced by fastp.

- [`fastp.stderr`](./fastp.stderr): Standard error output by fastp.

- [`multiqc_report.html`](./multiqc_report.html): MultiQC report in HTML format.

- [`multiqc_data/`](./multiqc_data/): MutliQC data directory. Contains anciallary files used by MultiQC to generate the HTML report file.
