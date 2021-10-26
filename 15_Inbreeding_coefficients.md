# Inbreeding coefficients from runs of homozygosity

Inbreeding coefficients are calculated for each individual from runs of homozygosity (ROHs) using [ROHan](https://github.com/grenaud/ROHan).

## Calculate transition:transvertion ratio for kiwi

ROHan requires information on the transition/transvertion ratio. This can be estimated in [bcftools](https://samtools.github.io/bcftools/) using the VCF files (.bcf) [previously generated](https://github.com/jordanbemmels/kiwi-holocene/blob/main/14_Selective_sweeps.md) as input for RAiSD. An example is shown for aHaast:

```
bcftools stats aHaast.bcf | grep "TSTV"
```

Repeat for all 11 populations and take the average: 3.08.

## Prepare ROHan input files

ROHan uses raw BAM files with reads aligned to a reference genome, rather than called SNPs or files that have been filtered. We do not need to do any further preparation, as BAM files were [previously generated](https://github.com/jordanbemmels/kiwi-holocene/blob/main/02_Initial_read_processing.md) during initial read processing. Note that the Z-chromosome will be excluded below.

## Run ROHan

ROHan is run on a single individual at a time. The tolerance for a background heterozygosity rate in true ROHs is increased to ```--rohmu 5e-5```; the transition:transversion rate for kiwi was calculated above ```--tstv 3.08```; the ROH size in this example is at least 1.5 Mbp ```--size 1500000```; we provide a list of autosomes with ```-auto``` (determined [previously](https://github.com/jordanbemmels/kiwi-holocene/blob/main/01_Identify_Zchr_scaffolds.md)) to restrict the analysis to autosomes only.

An example is shown for a single individual, KW36:

```
rohan --rohmu 5e-5 --tstv 3.08 --size 1500000 --auto autosomes_aRowi_1500kbp.txt -t 1 -o aHaast__KW36/aHaast__KW36 kiwi_ref_genome.fna australis__Haast__KW36__RA0997_sorted.bam
```

The output file ```aHaast__KW36.summary.txt``` classifies the proportion of the genome into three categories: ROH, non-ROH, and unclassifiable. The inbreeding coefficient (F_ROH) is the proportion in ROH (exactly as reported; this value ignores unclassifiable segments so no adjustment is needed).

Repeat the analysis for all individuals, as well as for all individuals with ```--size 5000000```. The two sets of analyses are used to calculate F_ROH_1.5 and F_ROH_5.0, respectively.

