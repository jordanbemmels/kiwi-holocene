# Identifying kiwi Z-chromosome scaffolds

This document describes the identification of which scaffolds in the kiwi reference genome putatively belong to the Z-chromosome. These scaffolds will be excluded from downstream analyses.

## Download reference genomes from NCBI and prepare

For the kiwi reference genome for this project, we will use *Apteryx rowi* [**GCA_003343035.1**](https://www.ncbi.nlm.nih.gov/assembly/GCF_003343035.1/).

To identify putative Z-chromosomes, compare to a genome for which the Z-chromosome has already been identified. Use ostrich (*Struthio camelus australis*) [**GCA_000698965.1**](https://www.ncbi.nlm.nih.gov/assembly/GCF_000698965.1/).

Index the reference genomes using [SAMtools](http://www.htslib.org/).

```
samtools faidx kiwi_ref_genome.fna
samtools faidx ostrich_ref_genome.fna
```

## Align the kiwi to the ostrich reference genome

Use the promer tool in [MUMmer](http://mummer.sourceforge.net/). Kiwi is the query and ostrich is the reference.

```
promer ostrich_ref_genome.fna kiwi_ref_genome.fna -p ostrich_kiwi
```

Filter alignments using further tools from MUMmer. Retain the best possible alignment, print summary info.

```
delta-filter -q -i 50 -l 500 ostrich_kiwi.delta > ostrich_kiwi_q_i50_l500.delta
show-coords -l -q -T ostrich_kiwi_q_i50_l500.delta > ostrich_kiwi_q_i50_l500_coords.txt
```

## Identify which aligned kiwi scaffolds are putatively on the Z-chromosome

Identify any kiwi scaffolds aligning to known ostrich Z-chromosome scaffolds. A list of ostrich Z-chromosome scaffolds is available from Fig. 1 in [Zhang et al. (2015)](https://doi.org/10.1186/s13742-015-0062-9).  Exclude from this list the erroneously placed scaffolds from [Xu et al. (2019)](https://doi.org/10.1093/gbe/evz154), i.e., those which were determined to have been erroneously mapped to the Z-chromosome.

Zhang J, Li C, Zhou Q, Zhang G. 2015 Improving the ostrich genome assembly using optical mapping data. Gigascience 4, 4–6.<br>
Xu L, Wa Sin SY, Grayson P, Edwards S V., Sackton TB, Mank J. 2019 Evolutionary dynamics of sex chromosomes of Paleognathous birds. Genome Biol. Evol. 11, 2376–2390.

Save the list of ostrich Z-chromosome scaffolds as the following file:

[ostrich_Zscaffolds_withNames.txt](https://github.com/jordanbemmels/kiwi-holocene/blob/main/ostrich_Zscaffolds_withNames.txt)

Use the *R* script [determineWhichKiwiZ_git.R](https://github.com/jordanbemmels/kiwi-holocene/blob/main/determineWhichKiwiZ_git.R) to determine which kiwi scaffolds putatively map to ostrict Z-chromosome scaffolds. These will be those with ≥50% of the total alignment length or ≥500,000 bp matching. The MUMmer file plus the ostrich Z-chromosome scaffolds file are inputs. The output file is called:

[kiwi_zChr_scaffolds.txt](https://github.com/jordanbemmels/kiwi-holocene/blob/main/kiwi_zChr_scaffolds.txt)
