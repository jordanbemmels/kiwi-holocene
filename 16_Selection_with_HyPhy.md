# Selection analyses with HyPhy

## PART I: PREPARE CODING SEQUENCES IN KIWI AND OUTGROUPS

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

Coding sequence (CDS) sites are extracted for each individual using BCFtools and outputted as a .bed file. Importantly, the ```-M``` and ```--absent``` (```-a```) options should be set to avoid having missing data replaced with the reference allele (this requires bcftools version 1.11+).

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

Detour: there are multiple isoforms available for some genes. Ideally, we might pick the most abundantly expressed isoform of each gene. We do not have that information for our species and we will choose the longest isoform of each gene instead. We already made a list of the longest isoforms for each gene annotated for Apteryx (```Aptrow.longest_AA.txt```). Use this list. It lists the protein names, so we need to translate these into a list of RNA names that will be used for the CDS sequences.

```
cd /home/0_GENOMES5/kiwi/0_data
cat Aptrow.longest_AA.txt | while read protein ; do grep "$protein" GCF_003343035.1/genomic.gff | head -n 1 | sed 's/^.*;Parent=//g' | sed 's/;.*$//g' >> Aptrow.longest_RNA.txt ; done
```

Ready to continue.

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

## Extract sequences from genes

Now we can use bedtools to extract the CDS from each sample using the bed files.

```
mkdir genes
#bedtools will extract a sequence for each of the CDS exons listed in the file, producing separate fasta entries for each exon
#then, we can process that file right away: remove the headers, remove the newlines (concatenating the exons), add a newline to the end, then add a header. That gives a single fasta entry for the whole CDS
#then concatenate that sequence to a single fasta that will hold the gene sequence for all samples. Make sure the unique sample name makes it into the fasta header so you can tell which is which.

#first, start it out with just the reference genome
#cat Aptrow.longest_RNA.txt | while read gene ; do bedtools getfasta -fi kiwi_ref_genome.fna -bed bed/"$gene".bed | sed 's/>.*$//g' | tr -d '\n' | sed 's/$/\n/g' | sed "s/^/>$sample\n/g" >> genes/"$gene".fa ; done

#loop through each sample and each gene
cat samples | cut -f 2 -d " " | while read sample ; do time cat Aptrow.longest_RNA.txt | while read gene ; do bedtools getfasta -s -fi filtered_genomes/"$sample".fa -bed bed/"$gene".bed | sed 's/>.*$//g' | tr -d '\n' | sed 's/$/\n/g' | sed "s/^/>$sample\n/g" >> genes/"$gene".fa ; done ; done

#cat samples | cut -f 2 -d " " | while read sample ; do time cat Aptrow.longest_RNA.txt | parallel bedtools getfasta -fi filtered_genomes/"$sample".fa -bed bed/{1}.bed '|' sed 's/>.*$//g' '|' tr -d '\n' '|' sed 's/$/\n/g' '|' sed "s/^/>$sample\n/g" '>>' genes/{1}.fa ; done
```

Now we have fasta files containing (unaligned) gene sequences for each annotated kiwi gene.

## Transfer kiwi annotation to outgroups

Use [Liftoff](https://github.com/agshumate/Liftoff) to transfer the kiwi annotation to outgroups, in order to find orthologous gene sequences in the outgroups.

The outgroup reference genome info is saved in the filw ```palaeognathae_acessions.txt``` as follows (format: genome_accession_number genus species assembly_version_name):

```
GCA_013396415.1 Casuarius casuarius ASM1339641v1
GCA_003342915.1 Crypturellus cinnamomeus cryCin1
GCA_016128335.1 Dromaius novaehollandiae ZJU1.0
GCA_003342835.1 Rhea pennata rhePen1
```

Liftoff takes a .gff file and maps its features onto another genome. Before beginning, filter the kiwi .gff file to retain only the longest isoform of each gene. Use this file (```GCF_003343035.1/Longest_CDS.gff```) to locate their putative orthologs in the outgroup taxa.

```
mkdir -p liftoff
#First, run it once on just 1 sample to get the database file
cat palaeognathae_acessions.txt | head -n 1 | parallel --colsep " " liftoff -g GCF_003343035.1/Longest_CDS.gff -dir liftoff/{1}_{4}_liftoff -o liftoff/{1}_{4}_liftoff.gff -u liftoff/{1}_{4}.unmapped_features.txt ncbi_dataset/data/{1}/{1}_{4}_genomic.fna kiwi_ref_genome.fna

#now run for all the rest. use the newly-created *.gff_db instead of the gff in order to save time.
cat palaeognathae_acessions.txt | tail -n +1 | parallel --colsep " " liftoff -db GCF_003343035.1/Longest_CDS.gff_db -dir liftoff/{1}_{4}_liftoff -o liftoff/{1}_{4}_liftoff.gff -u liftoff/{1}_{4}.unmapped_features.txt ncbi_dataset/data/{1}/{1}_{4}_genomic.fna kiwi_ref_genome.fna
```

Now that we have the annotations, we need to extract the protein and CDS sequences for each gene. For this task, we can use [gffread](https://github.com/gpertea/gffread). Liftoff did not polish the codon structure of the CDS's, so there are frameshifts that preclude directly translating them into AA sequence. To fix this, we can use [Transdecoder](https://github.com/TransDecoder/TransDecoder/wiki) to identify the longest open reading frame in each CDS and extract the corresponding protein sequence, which is what we will use to assign orthology.

```
mkdir -p transcriptome

cat palaeognathae_acessions.txt | parallel --colsep " " gffread -x transcriptome/{1}_{4}.longest_CDS.fna -y transcriptome/{1}_{4}.longest_proteins.faa -g ncbi_dataset/data/{1}/{1}_{4}_genomic.fna liftoff/{1}_{4}_liftoff.gff

mkdir transcriptome/transdecoder
#get longest ORFs for each CDS (liftoff did not respect codon structure so the AA sequences are very frameshifted)
cat palaeognathae_acessions.txt | parallel --colsep " " mkdir transcriptome/transdecoder/{1}_{4} ';' /home/0_PROGRAMS/TransDecoder-TransDecoder-v5.5.0/TransDecoder.LongOrfs -t transcriptome/{1}_{4}.longest_CDS.fna -m 30 -S --output_dir transcriptome/transdecoder/{1}_{4}
export PERL5LIB=/home/0_PROGRAMS/TransDecoder-TransDecoder-v5.5.0/PerlLib
cat palaeognathae_acessions.txt | parallel --colsep " " perl /home/0_PROGRAMS/TransDecoder-TransDecoder-v5.5.0/get_longest_ORF_per_transcript.pl transcriptome/transdecoder/{1}_{4}/longest_orfs.pep '>' transcriptome/transdecoder/{1}_{4}.faa #dont use the script copy in the util folder, and make sure to set PERL5LIB to point to where transdecoder has its perl modules
```

## Assign orthology

Use [OrthoFinder](https://github.com/davidemms/OrthoFinder) to assigne orthologs between all taxa. We can use this info to select single-copy orthologs to analyze.

```
mkdir -p orthology

#combine all proteomes of interest into a single folder
cat palaeognathae_acessions.txt | parallel --colsep " " cp transcriptome/transdecoder/{1}_{4}.faa orthology
time /home/0_PROGRAMS/seqkit grep -r -n -f Aptrow.longest_AA.txt GCF_003343035.1/protein.faa > orthology/longest.Aptrow.AA.fa

#add names of taxa to the fasta headers in each proteome to keep track
cat palaeognathae_acessions.txt | parallel --colsep " " sed -i '"s/>/>{2} /g"' orthology/{1}_{4}.faa
#replace protein name with RNA name for kiwi
paste Aptrow.longest_AA.txt Aptrow.longest_RNA.txt > name_map
time /home/0_PROGRAMS/maker/bin/map_fasta_ids name_map orthology/longest.Aptrow.AA.fa

#run OrthoFinder
#set up environment for OrthoFinder
export PATH=$PATH:/home/0_PROGRAMS/fastme-2.1.5/binaries
export PATH=$PATH:/home/0_PROGRAMS/ #for Diamond
export PATH=$PATH:/home/0_PROGRAMS/OrthoFinder
export PATH=$PATH:/home/0_PROGRAMS/ncbi-blast-2.11.0+/bin/
/home/0_PROGRAMS/OrthoFinder/orthofinder -o ./orthofinder -f ./orthology -t 60 -S blast > ./orthofinder.log
```

