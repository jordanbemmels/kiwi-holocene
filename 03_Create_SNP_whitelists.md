# Create whitelists of fitlered SNPs

These whitelists of filtered SNPs will be used in [ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD) either directly (for analyses performed in ANGSD) or to print files in other formats to use for downstream applications. The SNP whitelists (01, 02, and 03) and the specific analyses they are used for are described in more details the supplementary methods of the publication.

# WHITELIST 01

Whitelist 01 provides input sets of SNPs for PCAngsd, Fst, and SNAPPER analyses (see publication for details).

## Initial run to find all SNPs, further filtering will be required

Here, 44 of 55 sequenced individuals (80%) are required, but subsequent analyses suggested that three individuals turned out to be closely related to others and were removed. This left an actual total of 52 individuals used in the project, but at this point we do not know which individuals are closely related so we must use all 55 individuals. The parameters ```-minInd 44``` and ```-setMinDepthInd 3``` are very relaxed; this will find a large number of SNPs initially. The input file ```BAM_FILES_LIST.txt``` should contain a list of sorted, indexed bam files for all the 55 individuals.

```
FILTERS="-minQ 20 -minMapQ 20 -minInd 44 -setMinDepthInd 3 -SNP_pval 0.01 -rmTriallelic 0.01"
TODO="-doCounts 1 -doDepth 1 -maxDepth 50000 -doMajorMinor 1 -doMAF 1 -doGeno 11 -doPost 2"
angsd -b BAM_FILES_LIST.txt -GL 1 -P 4 $FILTERS $TODO -out filter01
```

Next, we want to exclude SNPs with depth greater than the 99th percentile of all SNPs, as these may be gene duplications or other problems. To determine what depth corresponds to the 99th percentile, simply open and inspect the output file ```filter01.depthGlobal``` containing the depth distribution. In our real kiwi data the max depth to include was 926.

Also, exclude any SNPs with >75% heterozygosity, which also could be errors in alignment or assembly. To do this, we use the script 
[HetMajorityProb.py](https://github.com/z0on/2bRAD_denovo/blob/master/HetMajorityProb.py) by Nate S. Pope. However, the default for this script is to calculate the probability that more than 50% of calls for a given site (the majority) are heterozygotes. To update the script to calculate the probability that more than 75% of the calls at a site are heterozygotes, we need to change line 25 in the script (current line number as of 2021/10/15):

**OLD LINE:** ```utail_prob = pois_binom.pval(len(pr_heteroz)/2+1)```<br>
**NEW LINE:** ```utail_prob = pois_binom.pval(int(len(pr_heteroz)*3/4))```

After updating, we are ready to run the heterozoygosity estimator using the .geno output file from ANGSD above.

```
gunzip --keep filter01.geno.gz
python HetMajorityProb.py < filter01.geno > hetExcessProbs.txt
```

Convert the output into an ANGSD-format whitelist of sites that did not show >75% heterozygosity.

```
sed '' hetExcessProbs.txt | awk '{OFS="\t"; if ($6 < 0.05){print $1,$2}}' > sitesHetFilter.txt
angsd sites index sitesHetFilter.txt
```

## Re-run ANGSD with heterozygosity and depth filters applied

Re-run ANGSD, using the whitelist of sites (```-sites sitesHetFilter.txt```) that did not show excess heterozygosity and the maximum depth (```-setMaxDepth 926```) calculated above. The *SNP_pval*, *minInd*, *rmTriallelic*, and *setMinDepthInd* filters do not need to be repeated as the sitesHetFilter.txt includes only sites that previously passed all those filters.

```
FILTERS="-minMapQ 20 -minQ 20 -setMaxDepth 926 -sites sitesHetFilter.txt"
TODO="-doCounts 1 -doDepth 1 -maxDepth 1000 -doMajorMinor 1 -doMAF 1 -doGeno 3 -doPost 2"
angsd -b BAM_FILES_LIST.txt -GL 1 -P 4 $FILTERS $TODO -out filter02
```

Now, print a whitelist of sites that passed this second round of filtering. As input, use the output .maf file from the previous step.

```
gunzip --keep filter02.mafs.gz
sed '' filter02.mafs | awk '{OFS="\t"; if ($1 !="chromo"){print $1,$2,$3,$4}}' > sitesFilter02.txt
angsd sites index sitesFilter02.txt
```

## Additional filtering and subsetting

We will want to further filter the basic list of sites (```sitesFilter02.txt```) to prepare whitelists for specific analyses.

In all cases, we would like to print only sites with data at all 55 individuals.

```
sed '' filter02.mafs | awk '{OFS="\t"; if ($1 !="chromo" && $6 == "55"){print $1,$2,$3,$4}}' > sitesFilter02_55ind.txt
mv sitesFilter02_55ind.txt sites_bfv2_3x_maf00_55ind.txt # rename
angsd sites index sites_bfv2_3x_maf00_55ind.txt
```

We can also create a different filter also requiring complete data for all 55 individuals but additionally identifying only non-singletons (minimum minor allele frequency = 1.5 / (55\*2) alleles = 0.0136; we require 1.5 instead of 2 allele copies to account for ANGSD's probabilistic framework).

```
sed '' filter02.mafs | awk '{OFS="\t"; if ($1 !="chromo" && $5 >= 0.0136 && $6 == "55"){print $1,$2,$3,$4}}' > sitesFilter02_maf0136_55ind.txt
angsd sites index sites_bfv2_3x_maf0136_55ind.txt
```

In all cases, we need to remove sites putatively on the Z-chromosome, using the list of Z-chromosome scaffolds ([kiwi_zChr_scaffolds.txt](https://github.com/jordanbemmels/kiwi-holocene/blob/main/kiwi_zChr_scaffolds.txt)) generated in [01_Identify_Zchr_scaffolds.md](https://github.com/jordanbemmels/kiwi-holocene/blob/main/01_Identify_Zchr_scaffolds.md). Do this by using the small script [remove_kiwiZchr_fromSites_git.R](https://github.com/jordanbemmels/kiwi-holocene/blob/main/remove_kiwiZchr_fromSites_git.R), which will remove any lines from the whitelist that correspond to sites on putative Z-chromosome scaffolds.

```
Rscript remove_kiwiZchr_fromSites_git.R sites_bfv2_3x_maf00_55ind.txt
Rscript remove_kiwiZchr_fromSites_git.R sites_bfv2_3x_maf0136_55ind.txt
```

The outfiles are named *sites_bfv2_3x_maf00_55ind_noZ.txt* and *sites_bfv2_3x_maf0136_55ind_noZ.txt*, respectively.

Finally, for some analyses in Whitelist 01, we want to use SNPs thinned to a minimum distance of 10 kbp. We can use the small script [thinSNPs_git.py](https://github.com/jordanbemmels/kiwi-holocene/blob/main/thinSNPs_git.py), which will process read the (sorted) sites line by line and print only sites that are on a new scaffold or are at least 10kbp away from the previously printed site. If a different distance in kbp is desired, change the command-line argument to something other than 10. The argument format is: infile outfile minDistanceRequired.

```
python thinSNPs.py sites_bfv2_3x_maf00_55ind_noZ.txt sites_bfv2_3x_maf00_55ind_10kbp_noZ.txt 10
python thinSNPs.py sites_bfv2_3x_maf0136_55ind_noZ.txt sites_bfv2_3x_maf0136_55ind_10kbp_noZ.txt 10
```

Now we have available whitelists from Whitelist 01 to use for specific analyses:<br>
```sites_bfv2_3x_maf00_55ind_noZ.txt```: Fst<br>
```sites_bfv2_3x_maf0136_55ind_10kbp_noZ.txt```: PCAngsd<br>
```sites_bfv2_3x_maf00_55ind_10kbp_noZ.txt```: SNAPPER

# WHITELIST 02

Whitelist 02 describes creating input sets of SNPs for pi (nucleotide diversity) and dXY. For these calculations, filtering is different than above because we need two list of both VARIABLE and TOTAL sites. The two lists must be generated using exactly the same commands, but initial trials repeating whitelist 01 for TOTAL sites ran into difficulties, with the lists not being directly comparable due to different ways in which applying filters sequentially across multiple steps sometimes resulted in different behaviour for variable and non-variable sites. Instead, Whitelist 02 standardizes the commands and ensures that the filters for both sets used are exactly the same and applied in a straightforward way, and all at once, to get rid of the unwanted behaviour.

## Variable sites

We require data for all 55 individuals right from the beginning, implement a maximum depth directly, and importantly require no minimum MAF.

```
FILTERS="-minQ 20 -minMapQ 20 -minInd 55 -setMinDepthInd 3 -setMaxDepth 926 -rmTriallelic 0.01 -SNP_pval 0.01"
TODO="-doCounts 1 -dumpCounts 1 -doMajorMinor 1 -doMAF 1"
angsd -b BAM_FILES_LIST.txt -GL 1 -P 1 $FILTERS $TODO -out filter02_noHetFilter
```

Then, repeat the same steps as above to identify and remove sites with excess heterozygosity, and remove the Z-chromosome. The final outfile is:

```sites_bfv3_forDxy_maf00_noZ.txt```

## Total sites

Technically, we are more interested in the NUMBER of sites that pass the same filters as for the variable sites, rather than the actual identities of these sites themselves, but we will first specifically find their identities then count the number of sites in that list. We will count the number of sites globally (all autosomes), and in smaller windows so we can additionally perform windowed analyses on these statistics. We will use the number of sites identified here in downstream analyses to adjust the pi and dXY estimates, rather than calculate those statistics directly from a whitelist of total sites (very inefficient in ANGSD).

We repeat above command exactly (for the variable sites) **EXCEPT** do not include a SNP filter, so that all sites will be printed including invariant sites (i.e., not just SNPs). And of course, the outfile name is different.

```
FILTERS="-minQ 20 -minMapQ 20 -minInd 55 -setMinDepthInd 3 -setMaxDepth 926 -rmTriallelic 0.01"
TODO="-doCounts 1 -dumpCounts 1 -doMajorMinor 1 -doMAF 1"
angsd -b BAM_FILES_LIST.txt -GL 1 -P 1 $FILTERS $TODO -out allSites_noHetFilter
```

Again, follow the same steps above from Whitelist 01 to remove the sites with excess heterozygosity. Note that you do not need to re-determine a new set of excess heterozygous sites (relative to that created immediately above in Whitelist 02 for the variable sites), because none of the additional new sites added here will have excess heterozygosity, as the additional sites are all invariant. In other words, you can re-use the excess heterozygosity file from *Whitelist 02/Variable sites*. Depending on your machine, you may need to split up the sites file on multiple cores for efficiency/memory, then recombine the output into a single file. The final file with excess heterozygosity sites removed that we need for further processing is called:

```allSites_filtered.pos```

Note that the Z-chromosome was not yet removed from this file (although it perhaps could have been).

Finally, count the total number of sites in all 5-kbp windows of the genome. We use 5-kbp (quite small) so that we can employ modularity in downstream analyses if desired, i.e., we can quickly combine the windows to get the number of sites in windows that are any multiple of 5,000 bp, without having to tediously re-process the enormous file of all sites every single time. The summary of number of sites per 5-kbp is more wieldy. Use the script [countPer5kbpWindow_startPos1_git.py](https://github.com/jordanbemmels/kiwi-holocene/blob/main/countPer5kbpWindow_startPos1_git.py); the command-line arguments are an ANGSD-formatted .pos file, plus the output file name.

```
python countPer5kbpWindow_startPos1_git.py allSites_filtered.pos allSites_filtered_5kbpWindow_startPos1.txt
```

Finally, count the number of total sites in the genome that have now passed the equivalent set of filters used to identify variable sites for Whitelist 02. To do so, simply open the previous outfile ```allSites_filtered_5kbpWindow_startPos1.txt``` (e.g., in *R*) and print the number of sites that are or are not on kiwi Z-chromosome scaffolds (as identified in ```kiwi_zChr_scaffolds.txt```).

Total sites including Zchr: 844430544<br>
Sites on Zchr: 40778120<br>
Total sites excluding Zchr: 803652424<br> # the total number of autosomal sites we can use to adjust global estimates of pi and dXY
Percent sites on Zchr: 4.8<br>

# WHITELIST 03

Whitelist 03 describes creating input sets of SNPs that are **SPECIFIC TO EACH INDIVIDUAL LINEAGE**, for use with PopSizeABC and RAiSD. We want only SNPs that are truly variable in each individual lineage, to get the max number of SNPs and so as not to include invariant sites in population-specific analyes. Here, an example is provided for the population aHaast (also may be written out in full as australis__Haast), but these steps would need to be repeated for all of the populations in the project.

We will run ANGSD individually using our file from Whitelist 01 that has all SNPs across all 55 individuals: ```sites_bfv2_3x_maf00_55ind.txt```. Note that although we provide a sites file that already has major-minor allele designations across all 55 kiwi individuals, we use ```-doMajorMinor 1``` to re-identify major and minor allele designations for each population individually. We also require data for all 5 individuals in this population. The key step here is the ```-SNP_pval 0.01``` flag, which will re-identify which of our previously-filtered SNPs are still variable within each population.

```
FILTERS="-minMapQ 20 -minQ 20 -minInd 5 -SNP_pval 0.01 -sites sites_bfv2_3x_maf00_55ind.txt"
TODO="-doMajorMinor 1 -doMAF 1"
angsd -b BAM__australis__Haast.txt -GL 1 -P 1 $FILTERS $TODO -out australis__Haast
```

Create an indexed sites file for this population.

```
gunzip --keep australis__Haast.mafs.gz
sed '' australis__Haast.mafs | awk '{OFS="\t"; if ($1 !="chromo"){print $1,$2,$3,$4}}' > sites_australis__Haast_bfv2_3x_maf00.txt
angsd sites index sites_australis__Haast_bfv2_3x_maf00.txt
```

Follow the example from above in the Whitelist 01 section to remove the Z-chromosome, and save the outfile as:

```
sites_australis__Haast_bfv2_3x_maf00_noZ.txt
```

Repeat the above for the remainin 10 lineages. These sites files are used as whitelists for PopSizeABC and RAiSD.

# FILTERS FOR OTHER ANALYSES

Other analyses do not use pre-defined SNP whitelists. SNP filtering (if any) for other analyses is described individually in the respective sections for each of the other analyses.
