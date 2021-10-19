# Identifying kiwi Z-chromosome scaffolds

This document describes the identification of which scaffolds in the kiwi reference genome putatively belong to the Z-chromosome. These scaffolds will be excluded from downstream analyses.

## Download reference genomes from NCBI and prepare

For the kiwi reference genome for this project, we will use *Apteryx rowi* GCA_003343035.1, available from:
https://www.ncbi.nlm.nih.gov/assembly/GCF_003343035.1/

To identify putative Z-chromosomes, compare to a genome for which the Z-chromosome has already been identified. Use ostrich (*Struthio camelus australis*) GCA_000698965.1, available from:
https://www.ncbi.nlm.nih.gov/assembly/GCF_000698965.1/

Index the reference genomes using [SAMtools](http://www.htslib.org/).

```
samtools faidx kiwi_ref_genome.fna
samtools faidx ostrich_ref_genome.fna
```

## Align 
