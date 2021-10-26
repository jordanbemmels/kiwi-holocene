# Selection analyses with HyPhy

Code for this section was developed by manuscript author Else Mikkelsen.
Consolidation and minor edits by Jordan Bemmels.

We will use the *aBSREL* model within [HyPhy](https://stevenweaver.github.io/hyphy-site/) to test for selection in the most recent common ancestor of kiwi.

## Call kiwi genotypes

Call sites sites and index the output file using [BCFtools](https://samtools.github.io/bcftools/). Note that ```bam.list``` is a text file indicating the paths to the 11 input .bam files for each of the 11 individuals used in HyPhy analyses (see manuscript). The bam files were generated [previously](https://github.com/jordanbemmels/kiwi-holocene/blob/main/02_Initial_read_processing.md) during initial read processing.

```
bcftools mpileup -B -C 50 -a AD,DP,SP -Ou -f kiwi_ref_genome.fna --threads 40 -b bam.list | bcftools call --threads 20 -f GQ,GP -m -O z - > Kiwi_11genomes.allsites.vcf.gz
bcftools index Kiwi_11genomes.allsites.vcf.gz
```

Filter the VCF file using [VCFtools](https://vcftools.github.io/man_latest.html) to a minimum quality of 20, a minimum mean depth of 5x, and a maximum mean depth of 21. Empirical mean depth is 14. The minimum mean depth of 5x is set lower than in other analyses, but is appropriate here because miscalling heterozygotes as homozygous due to low depth is not a major problem in this analysis.

```
vcftools --gzvcf Kiwi_11genomes.allsites.vcf.gz --remove-indels --minQ 20 --min-meanDP 5 --max-meanDP 21 --recode --stdout | gzip -c > Kiwi_11genomes.min5max21.vcf.gz 

#compress with bgzip since it appears to not have been compressed
mv Kiwi_11genomes.min5max21.vcf.gz Kiwi_11genomes.min5max21.vcf
bcftools view -I Kiwi_11genomes.min5max21.vcf -O z -o Kiwi_11genomes.min5max21.vcf.gz
bcftools index Kiwi_11genomes.min5max21.vcf.gz
```

## Extract genes

Coding sequence (CDS) sites are extracted using BCFtools. Importantly, the ```-M``` and ```--absent``` (```-a```) options should be set to avoid having missing data replaced with the reference allele (this requires bcftools version 1.11+).

```
#create a consensus genome sequence for each sample
cat samples | parallel --colsep " " 'bcftools consensus --fasta-ref kiwi_ref_genome.fna --sample {1} -M N -a N -H 1 Kiwi_11genomes.allsites.vcf.gz > genomes/{2}.fa'

mkdir filtered_genomes
#create a consensus genome sequence for each sample
cat samples | parallel --colsep " " 'time bcftools consensus --fasta-ref kiwi_ref_genome.fna --sample {1} -M N -a N -H 1 Kiwi_11genomes.min5max21.vcf.gz > filtered_genomes/{2}.fa'

#The samples were mapped to the genbank assembly but the annotation file is for the refseq assembly. The assemblies are identical, except that the names of the scaffolds are different. We just need to rename everything.
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/343/035/GCF_003343035.1_aptRow1/GCF_003343035.1_aptRow1_assembly_report.txt #lists the translation between Refseq and Genbank names
grep "RefSeq assembly and GenBank assemblies identical:" GCF_003343035.1_aptRow1_assembly_report.txt
# RefSeq assembly and GenBank assemblies identical: yes
grep -v "^#" GCF_003343035.1_aptRow1_assembly_report.txt | cut -f 5,7 > Genbank_Refseq_namingtable #needs to be a 2-column, tab-separated file, in which col1=oldname and col2=newname to be replaced
cp -r filtered_genomes temp_filtered_genomes #making a copy in case something goes wrong
conda activate maker
cat samples | parallel --colsep " " time /home/0_PROGRAMS/maker/bin/map_fasta_ids Genbank_Refseq_namingtable filtered_genomes/{2}.fa #took 11 seconds
#note if you had already indexed the fastas, you must now index them again
```

Obtain a .bed file listing the location of each gene in the genome, using [BEDOPS](https://bedops.readthedocs.io/en/latest/). Critically, this must list the coordinates for the coding sequences (CDS), without the introns or untranslated regions, since we will need to be able to translate them straight into protein sequence.

```
#get BED files for each protein
#loop through the list of genes. For each gene, look in the gff file for lines for the gene of interest, where column three indicates that that line pertains to the CDS. Then convert it to bed format.
mkdir bed
cat Aptrow.longest_RNA.txt | while read gene ; do grep "$gene" GCF_003343035.1/genomic.gff | awk '($3 == "CDS")' - | /home/0_PROGRAMS/bedops/bin/convert2bed -i gff - > bed/"$gene".bed ; done

#flip the order of the exons for genes on the minus strand, so that they are in the correct order relative to the gene, not the chromosome.
#first, get a file listing whether each gene is on the plus or minus strand
cat Aptrow.longest_RNA.txt | while read gene ; do head -n 1 bed/"$gene".bed | cut -f 6 | xargs -I {} echo "$gene" {} >> strandedness ; done
#reverse the order of the exons in the bed file for genes on the minus strand
cat strandedness | while read gene strandedness ; do if [[ "$strandedness" == "-" ]]; then sort --numeric-sort --reverse -k 2 bed/"$gene".bed > temp && mv temp bed/"$gene".bed ; fi ; done
```
