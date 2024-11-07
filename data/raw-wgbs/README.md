`project-mytilus-methylation/data/raw-wgbs`

Due to size, raw sequencing data will not be stored in this directory.

Data can be downloaded from:

- [https://owl.fish.washington.edu/nightingales/M_trossulus/](https://owl.fish.washington.edu/nightingales/M_trossulus/)


Download all files and corresponding checksums with the following:

1.

```bash
wget \
--recursive \
--level=1 \
--no-directories \
--accept "[0-9]M*.fastq.gz,[0-9]M*.fastq.gz.md5sum" \
https://owl.fish.washington.edu/nightingales/M_trossulus/
```

2.

```bash
wget https://owl.fish.washington.edu/nightingales/M_trossulus/105M_1.fastq.gz
wget https://owl.fish.washington.edu/nightingales/M_trossulus/105M_1.fastq.gz.md5sum
wget https://owl.fish.washington.edu/nightingales/M_trossulus/105M_2.fastq.gz
wget https://owl.fish.washington.edu/nightingales/M_trossulus/105M_2.fastq.gz.md5sum
wget https://owl.fish.washington.edu/nightingales/M_trossulus/106M_1.fastq.gz
wget https://owl.fish.washington.edu/nightingales/M_trossulus/106M_1.fastq.gz.md5sum
wget https://owl.fish.washington.edu/nightingales/M_trossulus/106M_2.fastq.gz
wget https://owl.fish.washington.edu/nightingales/M_trossulus/106M_2.fastq.gz.md5sum
wget https://owl.fish.washington.edu/nightingales/M_trossulus/107M_1.fastq.gz
wget https://owl.fish.washington.edu/nightingales/M_trossulus/107M_1.fastq.gz.md5sum
wget https://owl.fish.washington.edu/nightingales/M_trossulus/107M_2.fastq.gz
wget https://owl.fish.washington.edu/nightingales/M_trossulus/107M_2.fastq.gz.md5sum
wget https://owl.fish.washington.edu/nightingales/M_trossulus/109M_1.fastq.gz
wget https://owl.fish.washington.edu/nightingales/M_trossulus/109M_1.fastq.gz.md5sum
wget https://owl.fish.washington.edu/nightingales/M_trossulus/109M_2.fastq.gz
wget https://owl.fish.washington.edu/nightingales/M_trossulus/109M_2.fastq.gz.md5sum
wget https://owl.fish.washington.edu/nightingales/M_trossulus/239M_1.fastq.gz
wget https://owl.fish.washington.edu/nightingales/M_trossulus/239M_1.fastq.gz.md5sum
wget https://owl.fish.washington.edu/nightingales/M_trossulus/239M_2.fastq.gz
wget https://owl.fish.washington.edu/nightingales/M_trossulus/239M_2.fastq.gz.md5sum
wget https://owl.fish.washington.edu/nightingales/M_trossulus/240M_1.fastq.gz
wget https://owl.fish.washington.edu/nightingales/M_trossulus/240M_1.fastq.gz.md5sum
wget https://owl.fish.washington.edu/nightingales/M_trossulus/240M_2.fastq.gz
wget https://owl.fish.washington.edu/nightingales/M_trossulus/240M_2.fastq.gz.md5sum
wget https://owl.fish.washington.edu/nightingales/M_trossulus/241M_1.fastq.gz
wget https://owl.fish.washington.edu/nightingales/M_trossulus/241M_1.fastq.gz.md5sum
wget https://owl.fish.washington.edu/nightingales/M_trossulus/241M_2.fastq.gz
wget https://owl.fish.washington.edu/nightingales/M_trossulus/241M_2.fastq.gz.md5sum
wget https://owl.fish.washington.edu/nightingales/M_trossulus/242M_1.fastq.gz
wget https://owl.fish.washington.edu/nightingales/M_trossulus/242M_1.fastq.gz.md5sum
wget https://owl.fish.washington.edu/nightingales/M_trossulus/242M_2.fastq.gz
wget https://owl.fish.washington.edu/nightingales/M_trossulus/242M_2.fastq.gz.md5sum
wget https://owl.fish.washington.edu/nightingales/M_trossulus/269M_1.fastq.gz
wget https://owl.fish.washington.edu/nightingales/M_trossulus/269M_1.fastq.gz.md5sum
wget https://owl.fish.washington.edu/nightingales/M_trossulus/269M_2.fastq.gz
wget https://owl.fish.washington.edu/nightingales/M_trossulus/269M_2.fastq.gz.md5sum
wget https://owl.fish.washington.edu/nightingales/M_trossulus/270M_1.fastq.gz
wget https://owl.fish.washington.edu/nightingales/M_trossulus/270M_1.fastq.gz.md5sum
wget https://owl.fish.washington.edu/nightingales/M_trossulus/270M_2.fastq.gz
wget https://owl.fish.washington.edu/nightingales/M_trossulus/270M_2.fastq.gz.md5sum
wget https://owl.fish.washington.edu/nightingales/M_trossulus/271M_1.fastq.gz
wget https://owl.fish.washington.edu/nightingales/M_trossulus/271M_1.fastq.gz.md5sum
wget https://owl.fish.washington.edu/nightingales/M_trossulus/271M_2.fastq.gz
wget https://owl.fish.washington.edu/nightingales/M_trossulus/271M_2.fastq.gz.md5sum
wget https://owl.fish.washington.edu/nightingales/M_trossulus/272M_1.fastq.gz
wget https://owl.fish.washington.edu/nightingales/M_trossulus/272M_1.fastq.gz.md5sum
wget https://owl.fish.washington.edu/nightingales/M_trossulus/272M_2.fastq.gz
wget https://owl.fish.washington.edu/nightingales/M_trossulus/272M_2.fastq.gz.md5sum
wget https://owl.fish.washington.edu/nightingales/M_trossulus/69M_1.fastq.gz
wget https://owl.fish.washington.edu/nightingales/M_trossulus/69M_1.fastq.gz.md5sum
wget https://owl.fish.washington.edu/nightingales/M_trossulus/69M_2.fastq.gz
wget https://owl.fish.washington.edu/nightingales/M_trossulus/69M_2.fastq.gz.md5sum
wget https://owl.fish.washington.edu/nightingales/M_trossulus/70M_1.fastq.gz
wget https://owl.fish.washington.edu/nightingales/M_trossulus/70M_1.fastq.gz.md5sum
wget https://owl.fish.washington.edu/nightingales/M_trossulus/70M_2.fastq.gz
wget https://owl.fish.washington.edu/nightingales/M_trossulus/70M_2.fastq.gz.md5sum
wget https://owl.fish.washington.edu/nightingales/M_trossulus/71M_1.fastq.gz
wget https://owl.fish.washington.edu/nightingales/M_trossulus/71M_1.fastq.gz.md5sum
wget https://owl.fish.washington.edu/nightingales/M_trossulus/71M_2.fastq.gz
wget https://owl.fish.washington.edu/nightingales/M_trossulus/71M_2.fastq.gz.md5sum
wget https://owl.fish.washington.edu/nightingales/M_trossulus/72M_1.fastq.gz
wget https://owl.fish.washington.edu/nightingales/M_trossulus/72M_1.fastq.gz.md5sum
wget https://owl.fish.washington.edu/nightingales/M_trossulus/72M_2.fastq.gz
wget https://owl.fish.washington.edu/nightingales/M_trossulus/72M_2.fastq.gz.md5sum
wget https://owl.fish.washington.edu/nightingales/M_trossulus/78M_1.fastq.gz
wget https://owl.fish.washington.edu/nightingales/M_trossulus/78M_1.fastq.gz.md5sum
wget https://owl.fish.washington.edu/nightingales/M_trossulus/78M_2.fastq.gz
wget https://owl.fish.washington.edu/nightingales/M_trossulus/78M_2.fastq.gz.md5sum
wget https://owl.fish.washington.edu/nightingales/M_trossulus/79M_1.fastq.gz
wget https://owl.fish.washington.edu/nightingales/M_trossulus/79M_1.fastq.gz.md5sum
wget https://owl.fish.washington.edu/nightingales/M_trossulus/79M_2.fastq.gz
wget https://owl.fish.washington.edu/nightingales/M_trossulus/79M_2.fastq.gz.md5sum
wget https://owl.fish.washington.edu/nightingales/M_trossulus/80M_1.fastq.gz
wget https://owl.fish.washington.edu/nightingales/M_trossulus/80M_1.fastq.gz.md5sum
wget https://owl.fish.washington.edu/nightingales/M_trossulus/80M_2.fastq.gz
wget https://owl.fish.washington.edu/nightingales/M_trossulus/80M_2.fastq.gz.md5sum
wget https://owl.fish.washington.edu/nightingales/M_trossulus/81M_1.fastq.gz
wget https://owl.fish.washington.edu/nightingales/M_trossulus/81M_1.fastq.gz.md5sum
wget https://owl.fish.washington.edu/nightingales/M_trossulus/81M_2.fastq.gz
wget https://owl.fish.washington.edu/nightingales/M_trossulus/81M_2.fastq.gz.md5sum
wget https://owl.fish.washington.edu/nightingales/M_trossulus/92M_1.fastq.gz
wget https://owl.fish.washington.edu/nightingales/M_trossulus/92M_1.fastq.gz.md5sum
wget https://owl.fish.washington.edu/nightingales/M_trossulus/92M_2.fastq.gz
wget https://owl.fish.washington.edu/nightingales/M_trossulus/92M_2.fastq.gz.md5sum
wget https://owl.fish.washington.edu/nightingales/M_trossulus/93M_1.fastq.gz
wget https://owl.fish.washington.edu/nightingales/M_trossulus/93M_1.fastq.gz.md5sum
wget https://owl.fish.washington.edu/nightingales/M_trossulus/93M_2.fastq.gz
wget https://owl.fish.washington.edu/nightingales/M_trossulus/93M_2.fastq.gz.md5sum
wget https://owl.fish.washington.edu/nightingales/M_trossulus/94M_1.fastq.gz
wget https://owl.fish.washington.edu/nightingales/M_trossulus/94M_1.fastq.gz.md5sum
wget https://owl.fish.washington.edu/nightingales/M_trossulus/94M_2.fastq.gz
wget https://owl.fish.washington.edu/nightingales/M_trossulus/94M_2.fastq.gz.md5sum
wget https://owl.fish.washington.edu/nightingales/M_trossulus/95M_1.fastq.gz
wget https://owl.fish.washington.edu/nightingales/M_trossulus/95M_1.fastq.gz.md5sum
wget https://owl.fish.washington.edu/nightingales/M_trossulus/95M_2.fastq.gz
wget https://owl.fish.washington.edu/nightingales/M_trossulus/95M_2.fastq.gz.md5sum
```