# Inbreeding coefficients from runs of homozygosity

Inbreeding coefficients are calculated for each individual from runs of homozygosity (ROHs) using [ROHan](https://github.com/grenaud/ROHan).

## Calculate transition:transvertion ratio for kiwi

ROHan requires information on the transition/transvertion ratio. This can be estimated in [bcftoos](https://samtools.github.io/bcftools/) using the VCF files (.bcf) [previously generated](https://github.com/jordanbemmels/kiwi-holocene/blob/main/14_Selective_sweeps.md) as input for RAiSD. An example is shown for aHaast:

```
bcftools stats aHaast.bcf | grep "TSTV"
```

Repeat for all 11 populations and take the average: 3.08.

## Prepare ROHan input files

ROHan uses raw BAM files with reads aligned to a reference genome, rather than called SNPs or files that have been filtered. We do not need to do any further preparation, as BAM files were [previously generated](https://github.com/jordanbemmels/kiwi-holocene/blob/main/02_Initial_read_processing.md) during initial read processing. Note that the Z-chromosome will be excluded below.

## Run ROHan

ROHan is run on a single individual at a time. The tolerance for a background heterozygosity rate in true ROHs is set with ```--rohmu 5e-5```; see the manuscript for a discussion. We set the transition:transversion rate ```--tstv 3.08```, a ROH size of ```--size 1500000```, and provide a list of autosomes with ```-auto``` (determined [previously](https://github.com/jordanbemmels/kiwi-holocene/blob/main/01_Identify_Zchr_scaffolds.md)) to restrict the analysis to autosomes.


```
rohan --rohmu 5e-5 --tstv 3.08 --size 1500000 --auto autosomes_aRowi_1500kbp.txt -t 1 -o aHaast__KW36/aHaast__KW36 kiwi_ref_genome.fna australis__Haast__KW36__RA0997_sorted.bam
```
