# Divergence dating with SNAPPER

A quick time-calibrated phylogeny with divergence times can be estimated from SNPs from our whole-genome data using [SNAPPER](https://github.com/ForBioPhylogenomics/tutorials/blob/main/divergence_time_estimation_with_snp_data/README.md), a package add-on within [Beast2](http://www.beast2.org/). Note that robust divergence times in kiwi have previously been estimated with SNPs from genotyping-by-sequencing data plus mtDNA ([Weir et al. (2016)])(https://doi.org/10.1073/pnas.1603795113, using a coalescent approach and a much more detailed model allowing different population sizes and gene flow between lineages.

Weir JT, Haddrath O, Robertson HA, Colbourne RM, Baker AJ. 2016 Explosive ice age diversification of kiwi. Proc. Natl. Acad. Sci. U. S. A. 113, E5580–E5587.

We will use one individual to represent each lineage, for a total of 11 individuals.

## Identify and call SNPs

SNPs will be called in [ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD).

As input, begin with our list of filtered SNPs [previously identified](https://github.com/jordanbemmels/kiwi-holocene/blob/main/03_Create_SNP_whitelists.md) in Whitelist 01, with no minimum minor allele frequency (MAF): ```sites_bfv2_3x_maf00_55ind_10kbp_noZ.txt```. Since we are using 11 individuals, we will need to re-identify which sites are variable across the 11 focal individuals (```-SNP_pval 0.01 -rmTriallelic 0.01```). Require data for all individuals (```-minInd 11```), and only call SNPs with a posterior probability ≥0.99 (```-postCutoff 0.99```). Output as a VCF file (```-doBCF 1```).

```
FILTERS="-minMapQ 20 -minQ 20 -postCutoff 0.99 -minInd 11 -SNP_pval 0.01 -rmTriallelic 0.01 -sites sites_bfv2_3x_maf00_55ind_10kbp_noZ.txt -rf regions_bfv2_3x_maf00_55ind_10kbp_noZ.txt"
TODO="-doGeno 5 -dopost 1 -doMajorMinor 3 -doMAF 1 -doBCF 1 -doCounts 1"
angsd -b BAM_FILES_LIST_11ind.txt -GL 1 -P 1 $FILTERS $TODO -out kiwi11ind_maf00
```

The output file is ```kiwi11ind_maf00.bcf```.

## Prep and run SNAPPER

For dating, SNAPPER requires a constraints file. We will constrain the timing of the basal split between extant kiwi lineages, which was previously estimated (Weir et al. 2016) at 5.96 Ma, with 95% CI from 2.9-10.8 Ma. Following the [SNAPPER tutorial](https://github.com/ForBioPhylogenomics/tutorials/blob/main/divergence_time_estimation_with_snp_data/README.md), we will approximate the prior as a lognormal distribution centered at 5.96, with standard deviation of 0.33. Specify this in a new text file ```snapper_kiwiCrownConstraints.txt``` with a single line that reads as follows ("aHaast,aNorthFiordland,etc." are the names of the 11 lineages in the analysis):

```
lognormal(0,5.96,0.33)	crown	aHaast,aNorthFiordland,aSouthFiordland,aStewartIsland,hHaastii,mCoromandel,mEastern,mNorthland,mTaranaki,oKapiti,rOkarito
```




