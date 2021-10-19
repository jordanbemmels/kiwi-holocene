# Check relatedness of individuals

We will use [ngsRelateV2](https://github.com/ANGSD/NgsRelate) to check the pairwise relatedness of individuals. The purpose is to make sure there are no closely related individuals, which would need to be excluded from further analyses. In the empirical data with 55 kiwi, we did find 3 separate pairs of individuals that were closely related; with one individual excluded per pair, this resulted in retaining 52 individuals for downstream analyses.

As input, we need ANGSD VCF files (.bcf, extension has "b" not "v"). To generate these, we will calculate genotype likelihoods for each population individually and create output .bcf files (-dobcf 1). We can almost use the SNPs from [Whitelist 03](https://github.com/jordanbemmels/kiwi-holocene/blob/main/03_Create_SNP_whitelists.md); **however**, ngsRelateV2 should be run with unlinked sites, otherwise relatedness coefficients can be overestimated.

So, first create list of unlinked sites (minimum distance 10kbp) for each of the individual populations. An example is shown for aHaast (further details are explained [previously](https://github.com/jordanbemmels/kiwi-holocene/blob/main/03_Create_SNP_whitelists.md)): 

```
python thinSNPs.py sites_australis__Haast_bfv2_3x_maf00_noZ.txt sites_australis__Haast_bfv2_3x_maf00_10kbp_noZ.txt 10
```

We are ready to print the .bcf file using ANGSD. We require data for all 5 individuals, and a minimum depth of 8x per individual to include the site (optional).

```
FILTERS="-minQ 20 -minMapQ 20 -minInd 5 -setMinDepthInd 8"
TODO="-doCounts 1 -doMajorMinor 3 -doMaf 1 -doPost 1 -doGeno 1 -dobcf 1 -fold 1 --ignore-RG 0"
angsd -b BAM__australis__Haast.txt -sites sites_australis__Haast_bfv2_3x_maf00_10kbp_noZ.txt -GL 1 -P 1 $FILTERS $TODO -out aHaast_10kbp
```

Now, run ngsRelateV2 using the generated .bcf file. Note that ngsRelateV2 uses the genotype likelihoods in the .bcf file and not the genotype calls, which can be confirmed by inspecting the output of the program to terminal where is explains which part of the .bcf file is used in the calculation.

```
ngsRelate -h aHaast_10kbp.bcf -O aHaast_10kbp.res
```

Inspect the aHaast_10kbp.res output. Pairs of individuals with a KING coefficient ~0.25 or higher are parent-offspring or full-sibs and should be removed. According to [Manichaikul et al. (2010)](https://doi.org/10.1093/bioinformatics/btq559), the expected range around ~0.25 is [0.177, 0.354], so exclude one individual from any pairs with KING coefficient with KING coefficient ≥0.177.

Manichaikul A, Mychaleckyj JC, Rich SS, Daly K, Sale M, Chen WM. 2010 Robust relationship inference in genome-wide association studies. Bioinformatics 26, 2867–2873.

In our data there were three individuals to remove - they are completely excluded from the study in all subsequent analyses. The total number of individuals has now decreased from 55 to 52.
