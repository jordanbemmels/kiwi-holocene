# Infer genetic clusters (STRUCTURE-like analysis) and PCA with PCAngsd

Use the sites whitelist ([previously](https://github.com/jordanbemmels/kiwi-holocene/blob/main/03_Create_SNP_whitelists.md) created in Whitelist 01), using only non-singleton SNPs thinned to 10-kbp:

```
sites_bfv2_3x_maf0136_55ind_10kbp_noZ.txt
```

The input for [PCAngsd](http://www.popgen.dk/software/index.php/PCAngsdv2) is a BEAGLE-formatted file. Use ANGSD to create the BEAGLE file (-doGlf 2). Note that because in our [previous analysis](https://github.com/jordanbemmels/kiwi-holocene/blob/main/03_Create_SNP_whitelists.md) we discovered three closely related individuals to remove, our list of input individuals now has only 52 individuals instead of all 55: ```BAM_FILES_LIST_52ind.txt```.

```
FILTERS="-minMapQ 20 -minQ 20 -sites sites_bfv2_3x_maf0136_55ind_10kbp_noZ.txt"
TODO="-doGlf 2 -doMajorMinor 3 -doMaf 1"
angsd -b BAM_FILES_LIST_52ind.txt -GL 1 -P 1 $FILTERS $TODO -out genolike_maf0136_10kbp_noZ_52ind
```

Run PCAngsd with the number of clusters determined automatically by the program (no -e parameter specified). Since we switched to 52 instead of 55 individuals without changing the input sites file, the minimum minor allele frequency (minMAF) that corresponds to non-singletons needs to be updated. The new cutoff for identifying non-singletons will be 1.5 / (2\*52) = 0.0144, set directly in PCAngsd with ```-minMAF 0.0144```. Otherwise, the PCAngsd default is 0.05.

```
python pcangsd.py -beagle genolike_maf0136_10kbp_noZ_52ind.beagle.gz -admix -o maf0136_min0144_52ind_eAuto -minMaf 0.0144 -threads 4
```

The outfile ```maf0136_min0144_52ind_eAuto.admix.9.Q``` contains the proportion membership of each individual in each of the clusters. The automatically selected number of clusters (K) was K = 9 clusters.

We can also manually set the number of clusters and re-run the analysis. This is done using the -e parameter, set to 1 less than the number of clusters desired. For example, to set the clusters manually to K = 5 then use -e 4, and to set to K = 11 then use -e 10

```
python pcangsd.py -beagle genolike_maf0136_10kbp_noZ_52ind.beagle.gz -admix -o maf0136_min0144_52ind_e4 -minMaf 0.0144 -e 4 -threads 4
python pcangsd.py -beagle genolike_maf0136_10kbp_noZ_52ind.beagle.gz -admix -o maf0136_min0144_52ind_e10 -minMaf 0.0144 -e 10 -threads 4
```

