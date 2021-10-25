# Detecting selective sweeps with RAiSD

Use [RAiSD]() to detect selective sweeps. The overarching logic here is to calculate the RAiSD *mu*-statistic on empirical data, then use simulated data based on the demographic history of each population (as recovered in PopSizeABC) to set a population-specific threshold for distinguishing what *mu* value is considered a selective sweep.

## Prepare empirical input files

Input whitelists of lineage-specific variable sites (no minimum minor allele frequency) have [previously]() been described in Whitelist 03. We will use ANGSD to create VCF files of called SNPs using these sites, requiring data for all 5 individuals per population, and a minimum depth of 8x per individual per site to call the genotype. An example is showng for aHaast (also written as australis__Haast).

```
FILTERS="-minQ 20 -minMapQ 20 -minInd 5 -setMinDepthInd 8"
TODO="-doCounts 1 -doMajorMinor 3 -doMaf 1 -doBCF 1 -doPost 1 -doGeno 3"
angsd -b BAM__australis__Haast.txt -sites sites_australis__Haast_bfv2_3x_maf00_noZ.txt -GL 1 -P 1 $FILTERS $TODO -out aHaast
```

The output file is ```aHaast.bcf```.

## Run RAiSD on empirical datasets

There are not many settings that are relevant to us to modify, so RAiSD was run using default parameters.

```
RAiSD -n aHaast_default -I aHaast.bcf -O -D
```

There are several output files, but the most relevant one with mu-statistics at each calculated position in the genome is ```RAiSD_Report.aHaast_default.tabular.txt```. Optionally, we can summarize the mu-statistic to report the maximum mu per 50-kbp window (step size 25-kbp), which was done for kiwi and saved as ```RAiSD_Report.aHaast_default.win_50kbp_step25kbp.txt```.

## Simulate neutral datasets for comparison

We will use [msprime](https://github.com/tskit-dev/msprime/) to simulate neutral genomic datasets for comparison, using the demographic scenarios [previously](https://github.com/jordanbemmels/kiwi-holocene/blob/main/12_PopSizeABC_demography.md) estimated with PopSizeABC for each population. Then, we will re-run RAiSD on the neutral datasets. This will give us an expectation of the mu-statistic values we might expect due to neutral dynamics (i.e., population bottlenecks) alone, which can mimic the effects of selective sweeps.

For each population, we will use the 500 retained simulations (regression-ajusted values, not raw values) from the ABC step of PopSizeABC. This was done previously in the [PopSizeABC/abc_kiwi_git.R](https://github.com/jordanbemmels/kiwi-holocene/blob/main/PopSizeABC/abc_kiwi_git.R) script on Line 213:

```
write.table(retained_aHaast, "abc_output/retainedSims_adj/retainedSims_aHaast.txt", sep = "\t", quote = F, row.names = F, col.names = T);
```



