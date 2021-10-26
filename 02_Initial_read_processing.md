# Initial processing of raw sequencing reads

These filtering steps are performed prior to proceeding with analyses. A generic example is shown for a single individual with the identifier "SAMPLE". Due to paired-end sequencing, there are two input files, with suffixes "R1_001" and "R2_002".

## Remove Illumina adaptors and low-quality bases

This is performed using [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic). The file of Illumina adapter sequences *TruSeq3-PE-2.fa* is [distributed](https://github.com/timflutre/trimmomatic/blob/master/adapters/TruSeq3-PE-2.fa) with Trimmomatic.

```
java -jar trimmomatic-0.39.jar PE -threads 3 SAMPLE_R1_001.fasta.gz \
	SAMPLE_R2_001.fasta.gz \
	adaptorless_SAMPLE_R1.fastq.gz \
	adaptorless_SAMPLE_R1_unpaired.fastq.gz \
	adaptorless_SAMPLE_R2.fastq.gz \
	adaptorless_SAMPLE_R2_unpaired.fastq.gz \
	ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:7:5:true \
	TRAILING:3 LEADING:3 SLIDINGWINDOW:4:15 MINLEN:40 2>>SAMPLE_adaptrim_out.log
```

## Align reads to reference genome

This is done with [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml), using the ```--very-sensitive``` group of settings. The output is then piped through [SAMtools](http://www.htslib.org/) and saved as a .bam file.
                                      
```
bowtie2 -q --very-sensitive --no-unal --phred33 --time --threads 3 -x kiwi_ref_genome \
	-1 adaptorless_SAMPLE_R1.fastq.gz \
	-2 adaptorless_SAMPLE_R2.fastq.gz \
	-U adaptorless_SAMPLE_R1_unpaired.fastq.gz,adaptorless_SAMPLE_R2_unpaired.fastq.gz \
	2> SAMPLE__Bowtie2_log.txt | samtools view -S -b -@ 3 -h > SAMPLE_unsorted.bam;
```

Sort and index the bam file.
                                      
```
samtools sort -m 4G --threads=3 SAMPLE_unsorted.bam -o SAMPLE_sorted.bam
samtools index -b SAMPLE_sorted.bam
```

Calculate depth and breadth of coverage.
                                      
```
samtools depth -a SAMPLE_sorted.bam | awk '{c++;s+=$3}END{print s/c}' > depth_SAMPLE_sorted.txt
samtools depth -a SAMPLE_sorted.bam | awk '{c++; if($3>0) total+=1}END{print (total/c)*100}' > breadth_SAMPLE_sorted.txt
```
