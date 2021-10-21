# Genetic differentiation (Fst)

Population pairwise genetic differentiation (Fst) will be calculated in [ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD).

Use the unthinned SNPs whitelist [previously](https://github.com/jordanbemmels/kiwi-holocene/blob/main/03_Create_SNP_whitelists.md) generated from Whitelist 01, excluding singletons ("maf0136").

```
sites_bfv2_3x_maf0136_55ind_noZ.txt
```

ANGSD requires several steps to calculate Fst, here an example is between the two populations aHaast and aNorthFiordland.

Calculate a SAF (sample allele frequency) file for each population individually. Note that we must provide an ancestral state (```-anc kiwi_ref_genome.fna```); we can use the reference genome as it does not matter for Fst which is ancestral/derived as long as the same polarization is used across all populations.

```
FILTERS="-minMapQ 20 -minQ 20 -sites sites_bfv2_3x_maf0136_55ind_noZ.txt"
TODO="-dosaf 1"
angsd -b BAM__australis__Haast.txt -anc kiwi_ref_genome.fna -GL 1 -P 1 $FILTERS $TODO -out australis__Haast_noZ
angsd -b BAM__australis__NorthFiordland.txt -anc kiwi_ref_genome.fna -GL 1 -P 1 $FILTERS $TODO -out australis__NorthFiordland_noZ
```

Calculate the 2-dimensional SFS (site frequency spectrum) prior between the two populations.

```
misc/realSFS australis__Haast_noZ.saf.idx australis__NorthFiordland_noZ.saf.idx > aHaast_aNorthFiordland.ml
```

Prepare files to calculate Fst for windowed analysis, etc. Note that we use ```-whichFst 1``` (see ANGSD [manual](http://www.popgen.dk/angsd/index.php/Fst) for small sample sizes).

```
misc/realSFS fst index australis__Haast_noZ.saf.idx australis__NorthFiordland_noZ.saf.idx -whichFst 1 -sfs aHaast_aNorthFiordland.ml -fstout aHaast_aNorthFiordland.indexed
```

Calculate the global Fst. Note that we are interested in the outputted *weighted* Fst (not unweighted).

```
misc/realSFS fst stats aHaast_aNorthFiordland.indexed
```

Calculate Fst in 50-kbp sliding windows with 25-kbp step size. Use ```-type 2``` so that the window always starts on the first base pair of the scaffold (instead of starting where data begins), so that the windows are defined identically across multiple comparisons between other populations.

```
misc/realSFS fst stats2 aHaast_aNorthFiordland.indexed.fst.idx -win 50000 -step 25000 -type 2 > aHaast_aNorthFiordland.win50k_step25k.txt
```
