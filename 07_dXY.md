# Genetic divergence (dXY)

Population pairwise genetic divergence (dXY) will be calculated using [ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD) and [ngsPopGen](https://github.com/mfumagalli/ngsPopGen). Specifically, dXY can be calculated from ANGSD output files using the script [calcDxy.R](https://github.com/mfumagalli/ngsPopGen/blob/master/scripts/calcDxy.R), distributed with ngsPopGen.

An example is shown here for comparing the aHaast and aNorthFiordland populations.

## Prepare files in ANGSD

Step 1 as indicated in the ```calcDxy.R``` script (prior to running the script) is to identify sites. This was [previously](https://github.com/jordanbemmels/kiwi-holocene/blob/main/03_Create_SNP_whitelists.md) accomplished for the generation of SNP Whitelist 02. Note that we are using only variable sites (across all individuals), then adjusting our estimates relative to the total number of sites (variant plus invariant).

The sites file to use as input from Whitelist 02 is:

```
sites_bfv3_forDxy_maf00_noZ.txt
```

We still need to do Step 2 as indicated in ```calcDxy.R``` prior to running that script. This step is to calculate a .maf file (allele frequencies) for each population, being careful to include exactly the same sites with exactly the same polarization for all populations. Removing the ```-SNP_pval``` flag ensures that sites which are monomorphic in particular populations (but variable across all individuals from multiple populations) are still included in the output. Using our sites file in combination with ```-doMajorMinor``` 3 will use a pre-defined major-minor designation that we can re-use for all populations, to ensure identical polarization.

```
FILTERS="-minMapQ 20 -minQ 20 -sites sites_bfv3_forDxy_maf00_noZ.txt"
TODO="-doMajorMinor 3 -doMaf 1"
angsd -b BAM__australis__Haast.txt -GL 1 -P 1 $FILTERS $TODO -out aHaast_noZ
angsd -b BAM__australis__NorthFiordland.txt -GL 1 -P 1 $FILTERS $TODO -out aNorthFiordland_noZ
```

## Run the calcDxy.R script

We are ready to run calcDxy.R script.

The output of this command will be the **total** dxy (sum across all sites), rather than the per-site dXY averaged across all sites (what is normally meant by "dXY"). 

```
Rscript calcDxy.R -p aHaast_noZ.mafs -q aNorthFiordland_noZ.mafs > aHaast_aNorthFiordland_noZ_terminalOutput.txt
```

The ```aHaast_aNorthFiordland_terminalOutput.txt``` file saved from the terminal output will report the unadjusted global dXY (sum across all sites). The ```Dxy_persite.txt``` output file has the dXY per site (i.e., dXY at every individual site). Neither of these are exactly what we want yet.

To keep track, we should rename the file ```Dxy_persite.txt``` to something more specific, e.g., ```aHaast_aNorthFiordland_Dxy_persite.txt```.

## Sum dXY per 50-kbp window

Sum the dXY per site into 50-kbp windows, using the provided script [slidingWindow_Dxy_kiwi55_threadable_git.R](https://github.com/jordanbemmels/kiwi-holocene/blob/main/slidingWindow_Dxy_kiwi55_threadable_git.R). For this script, we will re-use the 50-kbp window definitions from the Fst calculation above as a template (```FST_DIRECTORY/aHaast_aNorthFiordland.win50k_step25k.txt```, where "FST_DIRECTORY" is changed to wherever you have saved this file), so that our windows exactly match those from the Fst calculations. This is specified on Line 43. Note that the window size is assumed to be 50-kbp, but can be changed on Line 48.

```
Rscript slidingWindow_Dxy_kiwi55_threadable_git.R aHaast_aNorthFiordland_Dxy_persite.txt aHaast_aNorthFiordland_Dxy_win50k_step25k.txt
```

## Adjust dXY relative to the total sites (variant + invariant)

We are ready to adjust the total divergence to convert it into an average per-site estimate of dXY. Essentially, we divide the total divergence by the total number of sites (variant plus invariant sites). This needs to be done for both the global estimate and the 50-kbp windowed estimates.

To adjust global estimate, we know that the total number of variant + invariant sites (passing the same filtering criteria as used for the variant sites only) was [previously](https://github.com/jordanbemmels/kiwi-holocene/blob/main/03_Create_SNP_whitelists.md) calculated when generating Whitelist 02 and was 803,652,424. As an example, the total divergence from the aHaast_aNorthFiordland_terminalOutput.txt file above is 1717128.07028095. **Thus, the
global dXY = 1717128.07028095 / 803652424 = 0.00214 for this particular population pair.**

To adjust the windowed estimates, first we need to calculate the number of total sites (variant + invariant) per 50-kbp window. We already have [previously](https://github.com/jordanbemmels/kiwi-holocene/blob/main/03_Create_SNP_whitelists.md) calculated the number of sites per 5-kbp window in Whitelist 02 (which was done for modularity so that we could use any window size we wanted that is a multiple of 5kbp). To combine the 5-kbp windows into 50-kbp windows, use the script [countSitesPerWindow_git.R](https://github.com/jordanbemmels/kiwi-holocene/blob/main/countSitesPerWindow_git.R). Again, we will match to the the Fst windows as a template (```FST_DIRECTORY/aHaast_aNorthFiordland.win50k_step25k.txt```). Note that this only needs to be done once (not repeated for every possible population pair) assuming that the same sites are used for all population pairs. The output file is ```FST_DIRECTORY/aHaast_aNorthFiordland.win50k_step25k.countSites.txt```.

```
Rscript countSitesPerWindow_git.R
```

Now, we are ready to finally adjust the windowed estimates, using the script [adjustSlidingWindowDxy_countSites_50kbp_step25kbp_git.R](https://github.com/jordanbemmels/kiwi-holocene/blob/main/adjustSlidingWindowDxy_countSites_50kbp_step25kbp_git.R), and the total number of sites per window just calculated above in ```FST_DIRECTORY/aHaast_aNorthFiordland.win50k_step25k.countSites.txt``` (specified within the script).

```
Rscript adjustSlidingWindowDxy_countSites_50kbp_step25kbp_git.R aHaast_aNorthFiordland_Dxy_win50k_step25k.txt aHaast_aNorthFiordland_Dxy_win50k_step25k_countSites_adj.txt
```
