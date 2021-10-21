# Observed heterozygosity (Ho)

Observed heterozygosity (Ho) will be calculated in [ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD). The the [ANGSD wiki](http://www.popgen.dk/angsd/index.php/Heterozygosity) suggests calculating Ho from all positions in the genome with minimal filtering, i.e., no filters for specific sites or minimum numbers of individuals. So, we will not use any pre-existing sites filter. However, we will provide a regions file to restrict our analysis to autosomes. The regions file (specified with -rf) is a list of autosomal scaffolds, one per line, followed by ":" to indicate to include all sites on that scaffold, for example:

```
scaffold_1:
scaffold_2:
etc.
```

The identity of autosomal vs. Z-chromosome scaffolds was [previously](https://github.com/jordanbemmels/kiwi-holocene/blob/main/01_Identify_Zchr_scaffolds.md) determined.

Here we will demonstrate the calculation of Ho for a single individual (KW36) from the aHaast population.

Calculate the folded site frequency spectrum (SFS) for a single diploid individual (single input BAM file specified with -i).

```
angsd -i australis__Haast__KW36__RA0997_sorted.bam -GL 1 -P 1 -minMapQ 20 -minQ 20 -anc kiwi_ref_genome.fna -rf kiwi_autosomal_scaffolds.txt -doSaf 1 -fold 1 -out australis__Haast__KW36__RA0997
misc/realSFS australis__Haast__KW36__RA0997.saf.idx > australis__Haast__KW36__RA0997.sfs
```

For a single individual, a folded SFS has only two entries, corresponding to homozygotes and heterozygotes. The estimate of Ho is thus the second number in the output file australis__Haast__KW36__RA0997.sfs divided by the total number of sites (the sum of all the numbers, i.e., the sum of homozygotes and heterozygotes):

```
# R script
cur_SFS <- read.table("australis__Haast__KW36__RA0997.sfs", header = F, sep = "", stringsAsFactors = F);
n_sites <- sum(cur_SFS);
Ho <- cur_SFS[ , 2] / sum(cur_SFS);
```
