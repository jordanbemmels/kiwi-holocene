# UCE phylogeny

A phylogeny of ultra-conserved elements (UCEs) will be generated in [Astral-III](https://github.com/smirarab/ASTRAL).

## Identify UCEs

We will identify UCEs in the kiwi reference genome using [phyluce](https://phyluce.readthedocs.io/en/latest/).

The reference genome should first be converted to 2bit format, using scripts downloaded from the [Kent Source Archive](http://hgdownload.soe.ucsc.edu/admin/exe/). To download the scripts:

```
rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/faToTwoBit ./ # downloads the faToTwoBit script for performing conversion
rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/twoBitInfo ./ # downloads the twoBitInfo script for viewing results
```

Now, perform the conversion and view the results. Note that we are re-naming the genome ```kiwi_ref_genome.fna``` to ```aRowi_ref.2bit``` to keep the original .fna and .2bit files clearly separate.

```
./faToTwoBit kiwi_ref_genome.fna aRowi_ref.2bit # perform conversion
./twoBitInfo aRowi_ref.2bit aRowi_ref_sizes.txt # view results
```

Download the vertebrate 5K probe set for UCE loci.

```
wget https://raw.githubusercontent.com/faircloth-lab/uce-probe-sets/master/uce-5k-probe-set/uce-5k-probes.fasta
```

Align the UCE probe set to the kiwi reference genome using phyluce.

```
phyluce_probe_run_multiple_lastzs_sqlite --db kiwi_UCE.sqlite --output kiwi_UCE-genome-lastz --scaffoldlist aRowi_ref --genome-base-path ./ --probefile uce-5k-probes.fasta --cores 1
```

The matches of the above command can be viewed in the outfile ```./kiwi_UCE-genome-lastz/uce-5k-probes.fasta_v_aRowi_ref.lastz.clean```. To understand what the output means and the format, inspect with [sqlite3](https://www.sqlite.org/index.html):

```
sqlite3 kiwi_UCE.sqlite
.mode columns
.headers on
.nullvalue .
.tables # prints which tables are in the database
SELECT * FROM aRowi_ref LIMIT 10; # prints the first 10 lines plus column names of aRowi_ref table
# result of this command, column names only: id  score  name1  strand1  zstart1  end1  length1  name2  strand2  zstart2  end2  length2  diff  cigar  identity  percent_identity  continuity  percent_continuity  coverage  percent_coverage
# the .lastz.clean file is the same, except is missing the "id" column
```

We need to convert the matches from the  output file ```uce-5k-probes.fasta_v_aRowi_ref.lastz.clean``` into a format ANGSD can use, such as an ANGSD regions file. When performing this conversion, note that phyluce defines site positions using a zero-origin half-open system, whereas ANGSD uses a one-indexed closed system, so we need to account for that in the conversion. We will also retain 1,000-bp flanking regions on either side of the UCE locus. If any separate UCEs ± 1000-bp flanking regions end up overlapping, then should be merged into a single larger locus. We should also remove any UCE positions on the Z-chromosome. Scaffolds corresponding to the Z-chromosome were identified [previously](https://github.com/jordanbemmels/kiwi-holocene/blob/main/01_Identify_Zchr_scaffolds.md). Perform the conversion taking these considerations into account using the script [convertLastzToAngsdRegions_git.R](https://github.com/jordanbemmels/kiwi-holocene/blob/main/convertLastzToAngsdRegions_git.R).

```
Rscript convertLastzToAngsdRegions.R
```

The input files are specified within the script. There are several output files produced corrseponding to UCE loci with different lengths of flanking regions. We are interested in those with 1000-bp flanking regions, so the relevant outfile is ```aRowi_ref_UCEs_plus_2x1000bp.txt```.

## Identify SNPs and call genotypes in ANGSD

We will use ANGSD to identify SNPs and call genotypes the UCE ± 1000-bp regions specified with ```-rf aRowi_ref_UCEs_plus_2x1000bp.txt```, corresponding to the file we just generated. We are also requiring sequencing data for all 52 individuals (```-minInd 52```), only calling biallelic SNPs (```-SNP_pval 0.01 -rmTriallelic 0.01```), and implementing a minimum minor allele frequency (MAF) of (2 / 2\*52) = 0.0192, which gets rid of singletons by requiring at least two allele copies. This minimum MAF is different than in some previous analyses with ANGSD where we required only 1.5 allele copies due to ANGSD's probabilistic genotype calling framework, but here we increased it to 2 allele copies since we are actually calling genotypes and have half an allele copy. To account for high uncertainty (especially for sites with low depth for particular individuals), we only call genotypes with posterior probaiblity greater than 0.99 (```-postCutoff 0.99```). We use ```-doBCF 1``` to output a VCF file (.bcf extension).

```
FILTERS="-minMapQ 20 -minQ 20 -postCutoff 0.99 -minInd 52 -minMAF 0.0192 -SNP_pval 0.01 -rmTriallelic 0.01 -rf aRowi_ref_UCEs_plus_2x1000bp.txt"
TODO="-doGeno 5 -dopost 1 -doMajorMinor 1 -doMAF 1 -doBCF 1 -doCounts 1"
angsd -b BAM_FILES_LIST_52ind.txt -GL 1 -P 1 $FILTERS $TODO -out UCEs_plus_2x1000bp_noOG_52ind_maf0192
```

Convert the .bcf file to phylip format using the script [vcf2phylip.py](https://github.com/edgardomortiz/vcf2phylip) by Edgardo M. Ortiz (DOI 10.5281/zenodo.1257057).

```
git clone https://github.com/edgardomortiz/vcf2phylip.git
python vcf2phylip/vcf2phylip.py -i UCEs_plus_2x1000bp_noOG_52ind_maf0192.bcf
```

The outfile is named ```UCEs_plus_2x1000bp_noOG_52ind_maf0192.bcf.phy```.

## 
