# Detecting selective sweeps with RAiSD

Use [RAiSD](https://github.com/alachins/raisd) to detect selective sweeps. The overarching logic here is to calculate the RAiSD *mu*-statistic on empirical data, then use simulated data based on the demographic history of each population (as recovered in PopSizeABC) to set a population-specific threshold for distinguishing what *mu* value is considered a selective sweep.

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

For each population, we will use the 500 retained simulations (regression-adjusted values, not raw values) from the ABC step of PopSizeABC. A file with the population sizes for those 500 simulations was created previously in the [PopSizeABC/abc_kiwi_git.R](https://github.com/jordanbemmels/kiwi-holocene/blob/main/PopSizeABC/abc_kiwi_git.R) script on Line 213:

```
write.table(retained_aHaast, "abc_output/retainedSims_adj/retainedSims_aHaast.txt", sep = "\t", quote = F, row.names = F, col.names = T);
```

After we have this file, we are ready to simulate the neutral simulations. Do this using the [provided](https://github.com/jordanbemmels/kiwi-holocene/blob/main/simulateNeutralSNPs_10M_500retained_20reps_git.py) script ```simulateNeutralSNPs_10M_500retained_20reps_git.py``` (here, the command-line parameter ```1``` indicates to perform the simulation for the first population, which is aHaast):

```
python2 simulateNeutralSNPs_10M_500retained_20reps_git.py 1
```

Note that the above script cannot be run in parallel for multiple populations because temporary files with the same name may be overwritten, so please run different populations sequentially. To use the script, you would need to adjust the relative paths to input and output files given your particular file system.

The key lines in the python script showing the msprime simulation, and then printing to a VCF outfile, are as follows:

```
tree_sequence = msprime.simulate(sample_size=n, Ne=N, length=L, recombination_rate=myrec, mutation_rate=mymut, demographic_events=mydemo)
# [...]
tree_sequence.write_vcf(vcf_file, ploidy=2) 
```

These lines also indicate the relevant parameters, which are defined under the section labelled ```##### Set up input data that will be the same across all populations #####```:
- ```sample_size``` is haploid sample size, set to 10
- ```Ne``` is the initial population size, determined from the current PopSizeABC demography
- ```length``` is the length of the chromosome to be simulated, set to 10 Mbp
- ```recombination_rate``` and ```mutation_rate``` are fixed at the values used in previous simulations
- ```demographic_events``` is set to match the population size changes defined by the current PopSizeABC demography

The output of this script for each population would be 500 VCF files of simulated neutral data, each of which contains data for 20 independently simulated chromosomes, e.g.:

```
aHaast_neutral_10M_500retained_20reps_scenario0.vcf
aHaast_neutral_10M_500retained_20reps_scenario1.vcf
[...]
```

## Run RAiSD on simulated datasets

Run RAiSD on each of the 500 simulated neutral datasets using the same commands as above for the empirical datasets. The output will be 500 different sets of RAiSD output files, including a *RAiSD_Report* file for each of the 500 demographic scenarios. Each *RAiSD_Report* file contains 20 different chromosomes.

To set the mu threshold for each population as in the manuscript, inspect the output for each of the 500 demographic scenarios, and retain the maximum mu value for each of the 20 simulated chromosomes. Discard the 5% of simulations (25 out of 500) resulting in the highest mean value of maximum mu (i.e., mean of the 20 simulated chromosomes). This is to account for demographic uncertainty, as these very extreme scenarios are unlikely to represent the true scenario. After discarding those, that leaves 475 scenarios x 20 chromosomes = 9,500 simulated chromosomes remaining. The threshold is the 99th percentile of the maximum mu value observed across these 9,500 chromosomes. This value is chosen because it corresponds to ~1 false positive per 100 chromosomes. Each simulated chromosome is 10 Mbp, so that is ~1 false positive per 1,000 Mbp. The reference genome is ~1,144 Mbp, so the final expectation is 1144 Mbp x 1 outlier/1000 Mbp = 1.144 outliers, or 1.144 false positives due to demographic stochasiticty per population.

