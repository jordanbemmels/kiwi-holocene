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
smc++ vcf2smc -d KW36 KW36 -m maskExcluded_PTFB01000001.1.bed.gz kiwi55_PTFB01000001.1.vcf.gz aHaast_PTFB01000001.1_di1.smc.gz PTFB01000001.1 aHaast:KW36,KW37,KW38,KW39,KW40
```

Repeat for all possible combinations of chromosome x population x distinguished individual.

