# Selection analyses with HyPhy

Code for this section was developed by manuscript author Else Mikkelsen.
Consolidation and minor edits by Jordan Bemmels.

We will use the *aBSREL* model within [HyPhy](https://stevenweaver.github.io/hyphy-site/) to test for selection in the most recent common ancestor of kiwi.

## Call kiwi genotypes

Call sites sites and index the output file using [BCFtools](https://samtools.github.io/bcftools/). Note that bam.list is a text file indicating the paths to the 11 input .bam files for each of the 11 individuals used in HyPhy analyses (see manuscript). The bam files were generated [previously](https://github.com/jordanbemmels/kiwi-holocene/blob/main/02_Initial_read_processing.md) during initial read processing.

```
bcftools mpileup -B -C 50 -a AD,DP,SP -Ou -f kiwi_ref_genome.fna --threads 40 -b bam.list | bcftools call --threads 20 -f GQ,GP -m -O z - > Kiwi_11genomes.allsites.vcf.gz
bcftools index Kiwi_11genomes.allsites.vcf.gz
```
