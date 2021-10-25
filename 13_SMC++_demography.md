# Demographic inference with SMC++

Another demographic inference (alternative to PopSizeABC) is done with [SMC++](https://github.com/popgenmethods/smcpp).

## Prepare input files

SMC++ requires info about *all* sites along a chromosome, so we do not use any pre-existing site filters. We will use all scaffolds at least 1 Mbp in length. For each scaffold, create a VCF file and a corresponding .bed mask file (to show which sites are masked due to having insufficient data for all individuals), using the [bamCaller.py](https://github.com/stschiff/msmc-tools/blob/master/bamCaller.py) script available from [msmc-tools](https://github.com/stschiff/msmc-tools). This script takes piped output from [bcftools](https://samtools.github.io/bcftools/). 

An example is shown for a single scaffold ("PTFB01001001.1"). Note that we are also requiring a minimum depth of 8x to retain the genotype call (+setGT command).

```
bcftools mpileup -q 20 -Q 20 -C 50 -r PTFB01001001.1 -a FORMAT/AD -d 10000 -f kiwi_ref_genome.fna -b BAM_FILES_LIST.txt | bcftools call -c -V indels | python bamCaller.py 708 maskIncluded_PTFB01001001.1.bed.gz | bcftools +setGT -- -t q -i 'FORMAT/AD[*:0] + FORMAT/AD[*:1] < 8' -n "./." | bcftools reheader -s BAM_FILES_LIST_55ind_renamed.txt | gzip -c > kiwi55_maskIncluded_PTFB01001001.1.vcf.gz
```

The output mask file ```kiwi55_maskIncluded_PTFB01001001.1.bed.gz``` is a list of sites that are *included* across all individuals, but SMC++ requires a mask of sites that are *excluded*. So, we will have to take the inverse of the maskIncluded file. This can be done using the [provided](https://github.com/jordanbemmels/kiwi-holocene/blob/main/convertMask_includedToExcluded.R) script ```convertMask_includedToExcluded_git.R```.

```
Rscript convertMask_includedToExcluded_git.R
```

The output ```maskExcluded_PTFB01000001.1.bed``` now has the sites that are *excluded*.

## Run SMC++

SMC++ first converts the VCF and .bed mask file into SMC++ file format. SMC++ also requires the designation of a distinguished individual (see SMC++ tutorial and publications). We can run the analysis with all possible combinations of distinguished individual. An example for aHaast for a single chromosome, with individual KW36 as the distinguished individual:

```
smc++ vcf2smc -d KW36 KW36 -m maskExcluded_PTFB01000001.1.bed.gz kiwi55_PTFB01000001.1.vcf.gz smc_files/aHaast/aHaast_PTFB01000001.1_di1.smc.gz PTFB01000001.1 aHaast:KW36,KW37,KW38,KW39,KW40
```

The relevant output file is ```smc_files/aHaast/aHaast_PTFB01000001.1_di1.smc.gz```. Repeat for all possible combinations of scaffold x population x distinguished individual.

Now, run SMC++, with all possible combinations of the .smc.gz output file (i.e., all scaffolds x distinguished individual) as input (```smc_files/aHaast/aHaast_*.smc.gz```). We use the ```cv``` command (and not ```estimate```) along with ```-folds 10``` to perform 10-fold cross validation. Increase the ```---knots``` parameter to potentially allow more sensitivity. Recombination and mutation rates used previously in other analyses are provided. We also begin at 5 generations ago (```timepoints 5 1568```), though the ending time of the analysis (1568) is not currently implemented (see [https://github.com/popgenmethods/smcpp/issues/177](https://github.com/popgenmethods/smcpp/issues/177)).

```
smc++ cv -r 2.1e-8 -o aHaast_cv_10folds_16knots_30ka --timepoints 5 1568 --folds 10 --knots 16 1.345868e-8 smc_files/aHaast/aHaast_*.smc.gz
```

To visualize the output (```aHaast_cv_10folds_16knots_30ka/model.final.json```), we supply a generation time (19.1357 years):

```
smc++ plot -g 19.1357 -c aHaast_cv_10folds_16knots_30ka/aHaast_cv_10folds_16knots_30ka.pdf aHaast_cv_10folds_16knots_30ka/model.final.json
```
