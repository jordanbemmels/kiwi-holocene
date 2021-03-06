# Divergence dating with SNAPPER

A quick time-calibrated phylogeny with divergence times can be estimated from SNPs from our whole-genome data using [SNAPPER](https://github.com/ForBioPhylogenomics/tutorials/blob/main/divergence_time_estimation_with_snp_data/README.md), a package add-on within [Beast2](http://www.beast2.org/). Note that robust divergence times in kiwi have previously been estimated with SNPs from genotyping-by-sequencing data plus mtDNA ([Weir et al. 2016](https://doi.org/10.1073/pnas.1603795113)), using a coalescent approach and a much more detailed model allowing different population sizes and gene flow between lineages.

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

Next, create an input population map file ```snapper_pops_11ind.txt```. This file lists all the individuals and what lineage ("species") they belong to. Here, we only have one individual per lineage. The contents of this file are as follows:

```
species	individual
aHaast	australis__Haast__KW36__RA0997
aNorthFiordland	australis__NorthFiordland__KW41__R32961
aSouthFiordland	australis__SouthFiordland__KW50__RA0202
aStewartIsland	australis__StewartIsland__KW51__RA0891
hHaastii	haastii__haastii__KW26__CD830
mCoromandel	mantelli__Coromandel__KW09__R32852
mEastern	mantelli__Eastern__KW11__NIBKOpo1
mNorthland	mantelli__Northland__KW04__41
mTaranaki	mantelli__Taranaki__KW18__kiwiCarcass1
oKapiti	owenii__Kapiti__KW24__O20581
rOkarito	rowi__Okarito__KW32__R32934
```

Finally, create an input .xml file for SNAPPER, using the [Ruby](https://www.ruby-lang.org/en/) script ```snapp_prep.rb``` that is [distributed](https://github.com/mmatschiner/snapp_prep) with SNAPP / SNAPPER. We specify that we want a SNAPPER and not a SNAPP analysis with ```-a SNAPPER```; provide our input VCF (```-v```), constraints (```-c```), and population map (```-t```) files; use a random sample of 1,000 SNPs (```-m 1000```); set the MCMC chain length to 2,000,000 (```-l 2000000```); and specify the output file name (```-x```).

```
ruby snapp_prep.rb -a SNAPPER -v kiwi11ind_maf00.bcf -t snapper_pops_11ind.txt -c snapper_kiwiCrownConstraints.txt -m 1000 -l 2000000 -x snapper_11ind_m1K_l2M.xml
```

Now, we are ready to run SNAPPER. The settings are all already specified in the .xml file.

```
beast/bin/beast -threads 23 snapper_11ind_m1K_l2M.xml
```

Repeat the Beast2 command for 5 replicate runs (in different directories), and we will check each for convergence and combine the output of the 5 individual runs below.

## Post-processing of individual runs

Population size (theta) is not printed to the SNAPPER output by defalt. Add theta using the ```add_theta_to_log.rb``` script that is [distributed](https://github.com/mmatschiner/snapp_prep) with SNAPP / SNAPPER. We must provide the generation time with ```-g```. The generation time for kiwi was estimated to be 19.1357 years (see manuscript text).

```
ruby add_theta_to_log.rb -l snapper.log -t snapper.trees -g 19.1357 -o snapper_w_popsize.log
```

We also want to know how much support there is for the branching pattern. We can quantify the posterior probabilities of clades on the maximum clade credibility tree using the [TreeAnnotator](https://www.beast2.org/treeannotator/) add-on for Beast2. Set the burn-in to 10% (```-burnin 10```).

```
beast/bin/treeannotator -burnin 10 -heights mean snapper.trees snapper_maxCladeCredibility.tree
```

The output files can be inspected using different visualization programs:<br>
  ```snapper_w_popsize.log``` -> [Tracer](https://www.beast2.org/tracer-2/)<br>
  ```snapper.trees``` -> [DensiTree](https://www.cs.auckland.ac.nz/~remco/DensiTree/)<br>
  ```snapper_maxCladeCredibility.tree``` -> [FigTree](http://tree.bio.ed.ac.uk/software/figtree/)

In particular, for each run we should check convergence with the .log file in Tracer, with the burn-in set to 10%. Visually check for convergence plus aim for effective sample size (ESS) values >200 for each run. If this is not achieved, try re-running with a longer MCMC chain.

## Combine runs for final result

If each of the individual runs have converged, we can combine them into a single final result, discarding the first 10% as burn-in. Use the provided script [combine_snapper_runs_git.R](https://github.com/jordanbemmels/kiwi-holocene/blob/main/combine_snapper_runs_git.R). Parameters can be changed within the *R* script if needed, but by default are set to the values corresponding to the analyses described above.

```
Rscript combine_snapper_runs_git.R
```

The output file are ```snapper_combined_11ind_m1K_l10M.log``` and ```snapper_combined_11ind_m1K_l10M.trees```.

Now, repeat the command from above to get the maximum clade credibility tree, this time for the combined dataset. However, change the burnin to zero (```-burnin 0```) because we have already removed the first 10% of each individual run prior to combining.

```
beast/bin/treeannotator -burnin 0 -heights mean snapper_combined_11ind_m1K_l10M.trees snapper_combined_11ind_m1K_l10M_maxCladeCredibility.tree
```

As above, visually check for convergence and aim for an ESS (viewed in Tracer) of >1000 on the combined runs.


 
