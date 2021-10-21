# Nucleotide diversity (pi) and Tajima's D

Nucleotide diversity (pi) and Tajima's D are calculated together using the thetas commands in [ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD).

The calculation of pi is done using an old method from older versions of ANGSD; this older method is now considered obsolete. However, the new replacement method (currently described on ANGSD Wiki) appears still buggy with folded data (see [https://github.com/ANGSD/angsd/issues/336](https://github.com/ANGSD/angsd/issues/336)) to the best of my abilities and comprehension. Thus, I have used the original older method from ANGSD.

An example is provided here for the population aHaast (australis__Haast).

## Nucleotide diversity (pi)

Generate a global estimate of the site frequency spectrum (SFS), using ```-fold 1``` to specify that the SFS is folded (i.e., it is unknown which allele is ancestral vs. derived). As input, use the sites file generated [previously]() from Whitelist 02 (```sites_bfv3_forDxy_maf00_noZ.txt```) for dXY.

```
FILTERS="-minMapQ 20 -minQ 20 -sites sites_bfv3_forDxy_maf00_noZ.txt -anc kiwi_ref_genome.fna"
TODO="-doSaf 1 -fold 1"
angsd -b BAM__australis__Haast.txt -GL 1 -P 1 $FILTERS $TODO -out BAM__australis__Haast
misc/realSFS BAM__australis__Haast.saf.idx > BAM__australis__Haast.sfs
```

Now, calculate thetas by re-running ANGSD with ```-doThetas 1```. Use ```-fold 1``` for the folded SFS. We supply priors with ```-pest``` using the SFS file produced in the previous step.

```
FILTERS="-minMapQ 20 -minQ 20 -sites sites_bfv3_forDxy_maf00_noZ.txt -anc kiwi_ref_genome.fna"
TODO="-doThetas 1 -doSaf 1 -fold 1"
angsd -b BAM__australis__Haast.txt -GL 1 -P 1 -pest BAM__australis__Haast.sfs $FILTERS $TODO -out australis__Haast
```

Print the estimates for parameters like pi, Tajima's D, other theta statistics. The outfile is ```australis__Haast.thetas.idx.pestPG```.

```
misc/thetaStat do_stat australis__Haast.thetas.idx
```

Adjust the estimate of pi relative to the total number of sites (variant + invariant) that met the same filtering criteria as used for the variant sites. The total number of sites was [previously](https://github.com/jordanbemmels/kiwi-holocene/blob/main/03_Create_SNP_whitelists.md) calculated to be 803,652,424 for Whitelist 02. However, some variable sites were excluded from the calculation of theta stats. The number of variable sites in our input sites file ```sites_bfv3_forDxy_maf00_noZ.txt``` was 26,759,111. In contrast, the number of sites in the ```australis__Haast.thetas.idx.pestPG``` file generated here (sum of the last column) is 26,385,298, so a total of 373,813 fewer sites were processed than expected. These excluded sites can be subtracted from the number of sites total when adjusting global pi. We can perform these calculations as follows:

```
# R script
wc_allSites_file <- 803652424; # total number of sites in genome that meet filtering criteria;
wc_variableSites_file <- 26759111 # total number of variable sites (SNPs) in the input sites file;
thetas <- read.table("australis__Haast.thetas.idx.pestPG", header = F, stringsAsFactors = F);
pi_raw <- sum(thetas[, 5]); # this corresponds to column tP;
variable_nSites_true <- sum(thetas[, 14]); # this corresponds to column nSites;
excluded_nSites <- wc_variableSites_file - variable_nSites_true;
total_nSites_adj <- wc_allSites_file - excluded_nSites;
pi_adj <- pi_raw / total_nSites_adj; # this is the adjusted global pi estimate;
```

Also print pi per 50-kbp window, with 25-kbp step size, to match the same designation as used previously in Fst.

```
misc/thetaStat do_stat australis__Haast.thetas.idx -win 50000 -step 25000 -type 2 -outnames aHaast.thetas_win50k_step25k
```

Adjust the windowed pi relative to the total number of sites passing filters, which was calculated [previously](https://github.com/jordanbemmels/kiwi-holocene/blob/main/07_dXY.md) for dXY and is available in the file ```FST_DIRECTORY/aHaast_aNorthFiordland.win50k_step25k.countSites.txt```. Use the script [adjustSlidingWindowPi_50kbp_step25kbp_git.R](https://github.com/jordanbemmels/kiwi-holocene/blob/main/adjustSlidingWindowPi_50kbp_step25kbp_git.R). The two command-line arguments are the input file and output file names, resectively.

```
Rscript adjustSlidingWindowPi_50kbp_step25kbp_git.R aHaast.thetas_win50k_step25k.pestPG aHaast.thetas_win50k_step25k_adj.pestPG
```

## Tajima's D

Tajima's D is cannot be calculated globally in ANGSD, so instead we will report mean Tajima's D across non-overlapping 10-kbp windows, following [Colella et al. (2017)](https://doi.org/10.1093/jhered/esab009).

Colella JP, Tigano A, Dudchenko O, Omer AD, Khan R, Bochkov ID, Aiden EL, MacManes MD. 2021 Limited evidence for parallel evolution among desert-adapted Peromyscus deer mice. J. Hered. 112, 286–302.

Summarize theta statistics in non-overlapping 10-kbp windows.

```
misc/thetaStat do_stat australis__Haast.thetas.idx -win 10000 -step 10000 -type 2 -outnames aHaast.thetas_win10k
```

The mean Tajima's D is the mean of the "Tajima" column in the output file ```aHaast.thetas_win10k```.
