# Create whitelists of fitlered SNPs

These whitelists of filtered SNPs will be used in [ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD) either directly (for analyses performed in ANGSD) or to print files in other formats to use for downstream applications. The SNP whitelists (01, 02, and 03) and the specific analyses they are used for are described in more details the supplementary methods of the publication.

# WHITELIST 01

Whitelist 01 provides input sets of SNPs for PCAngsd, Fst, and SNAPPER analyses (see publication for details).

## Initial run to find all SNPs, further filtering will be required

Here, 44 of 55 sequenced individuals (80%) are required, but subsequent analyses suggested that three individuals turned out to be closely related to others and were removed. This left an actual total of 52 individuals used in the project, but at this point we do not know which individuals are closely related so we must use all 55 individuals. The parameters ```-minInd 44``` and ```-setMinDepthInd 3``` are very relaxed; this will find a large number of SNPs initially. The input file ```BAM_FILES_LIST.txt``` should contain a list of sorted, indexed bam files for all the 55 individuals.

```
FILTERS="-minQ 20 -minMapQ 20 -minInd 44 -setMinDepthInd 3 -SNP_pval 0.01 -rmTriallelic 0.01"
TODO="-doCounts 1 -doDepth 1 -maxDepth 50000 -doMajorMinor 1 -doMAF 1 -doGeno 11 -doPost 2"
angsd -b BAM_FILES_LIST.txt -GL 1 -P 4 $FILTERS $TODO -out filter01
```

Next, we want to exclude SNPs with depth greater than the 99th percentile of all SNPs, as these may be gene duplications or other problems. To determine what depth corresponds to the 99th percentile, simply open and inspect the output file ```filter01.depthGlobal``` containing the depth distribution. In our real kiwi data the max depth to include was 926.

Also, exclude any SNPs with >75% heterozygosity, which also could be errors in alignment or assembly. To do this, we use the script 
[HetMajorityProb.py](https://github.com/z0on/2bRAD_denovo/blob/master/HetMajorityProb.py) by Nate S. Pope. However, the default for this script is to calculate the probability that more than 50% of calls for a given site (the majority) are heterozygotes. To update the script to calculate the probability that more than 75% of the calls at a site are heterozygotes, we need to change line 25 in the script (current line number as of 2021/10/15):

OLD LINE: ```utail_prob = pois_binom.pval(len(pr_heteroz)/2+1)```
NEW LINE: ```utail_prob = pois_binom.pval(int(len(pr_heteroz)*3/4))```


# ready to run on the .geno output file from ANGSD above

gunzip --keep filter01.geno.gz
python HetMajorityProb.py < filter01.geno > hetExcessProbs.txt

# print ANGSD-format whitelist of sites that did not show >75% heterozygosity

sed '' hetExcessProbs.txt | awk '{OFS="\t"; if ($6 < 0.05){print $1,$2}}' > sitesHetFilter.txt
angsd sites index sitesHetFilter.txt

#####

# re-run ANGSD, using the whitelist of sites (-sites sitesHetFilter.txt) that did not show excess heterozygosity and the maximum depth (-setMaxDepth 926) calculated above
# SNP_pval, minInd, rmTriallelic, setMinDepthInd filters do not need to be repeated as the sitesHetFilter.txt includes only sites that previously passed all those filters

FILTERS="-minMapQ 20 -minQ 20 -setMaxDepth 926 -sites sitesHetFilter.txt"
TODO="-doCounts 1 -doDepth 1 -maxDepth 1000 -doMajorMinor 1 -doMAF 1 -doGeno 3 -doPost 2"
angsd -b BAM_FILES_LIST.txt -GL 1 -P 4 $FILTERS $TODO -out filter02

# print whitelist of sites that passed this second round of filtering; use the output .maf file

gunzip --keep filter02.mafs.gz
sed '' filter02.mafs | awk '{OFS="\t"; if ($1 !="chromo"){print $1,$2,$3,$4}}' > sitesFilter02.txt
angsd sites index sitesFilter02.txt

#####

# apply additional filtering and subsetting for analyses

# print only sites with data at all 55 individuals

sed '' filter02.mafs | awk '{OFS="\t"; if ($1 !="chromo" && $6 == "55"){print $1,$2,$3,$4}}' > sitesFilter02_55ind.txt
mv sitesFilter02_55ind.txt sites_bfv2_3x_maf00_55ind.txt # rename
angsd sites index sites_bfv2_3x_maf00_55ind.txt

# also identify non-singletons (minimum minor allele frequency = 1.5 / (55*2) alleles = 0.0136; we require 1.5 instead of 2 allele copies to account for ANGSD's probabilistic framework)

sed '' filter02.mafs | awk '{OFS="\t"; if ($1 !="chromo" && $5 >= 0.0136 && $6 == "55"){print $1,$2,$3,$4}}' > sitesFilter02_maf0136_55ind.txt
angsd sites index sites_bfv2_3x_maf0136_55ind.txt

# remove sites putatively on Z-chromosome, using list of Z-chromosome scaffolds from above
# process e.g. in R, remove any lines corresponding to sites on scaffolds that are listed in kiwi_zChr_scaffolds.txt
# example for maf00 version:

//
# R script
zChrScaffolds <- read.table("kiwi_zChr_scaffolds.txt");
infile <- "sites_bfv2_3x_maf00_55ind.txt";
outfile <- "sites_bfv2_3x_maf00_55ind_noZ.txt";
#
con  <- file(infile, open = "r");
conOut <- file(outfile, open = "w");
count = 0;
while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
	if (count %% 100000 == 0) {
		print(paste0("working on line ", count));
	}
	temp_line <- (strsplit(oneLine, "\t"));
	if (temp_line[[1]][1] %in% zChrScaffolds) {
		count = count + 1;
	} else {
		writeLines(oneLine, conOut);
		count = count + 1;
	}
} 
close(con);
close(conOut);
//

# names of outfiles:

sites_bfv2_3x_maf00_55ind_noZ.txt
sites_bfv2_3x_maf0136_55ind_noZ.txt

# thin to a minimum distance of 10-kbp
# process line-by-line e.g. in python, print only sites that are on a new scaffold or are at least 10kbp away from the previously printed site
# example for maf00 version:

//
# python script
import sys
#
inputSites = sites_bfv2_3x_maf00_55ind_noZ.txt
outputFile = sites_bfv2_3x_maf00_55ind_10kbp_noZ.txt
kbp = 10 # minimum distance between SNPs in kbp for thinning
outfile = open(outputFile, 'w')
#
linecount = 0
with open(inputSites) as f:
	for line in f:
		linecount = linecount + 1
		if (linecount % 10000 == 0):
			print("working on line " + str(linecount))
		if (linecount > -1):
			# split into chrom and pos
			currentChrom, currentPos = line.split()[0], int(line.split()[1])
			#print (currentChrom)
			#print(currentPos)
			# 
			if (linecount == 1):
				# if it's the very first line, accept the SNP
				outfile.write(line)
				#print("condition1")
				previousChrom = currentChrom
				previousPos = currentPos
			elif (currentChrom != previousChrom):
				# otherwise, check if it's a new chromosome, and if it is, accept the SNP
				outfile.write(line)
				#print("condition2")
				previousChrom = currentChrom
				previousPos = currentPos
			elif (currentPos >= previousPos + kbp*1000):
				# otherwise, only write the SNP if it is at least the desired kbp away
				outfile.write(line)
				#print(previousPos + kbp*1000)
				#print(currentPos)
				#print("condition3")
				previousChrom = currentChrom
				previousPos = currentPos
#
outfile.close()
//

# names of outfiles:

sites_bfv2_3x_maf00_55ind_10kbp_noZ.txt
sites_bfv2_3x_maf0136_55ind_10kbp_noZ.txt

# input whitelists of SNPs from Group 01 for specific analyses:
	# Fst: sites_bfv2_3x_maf0136_55ind_noZ.txt
	# PCAngsd: sites_bfv2_3x_maf0136_55ind_10kbp_noZ.txt
	# SNAPPER: sites_bfv2_3x_maf00_55ind_10kbp_noZ.txt

#####
##### WHITELIST 02 AS DESCRIBED IN SUPPLEMENTARY METHODS ####
#####

# whitelist 02 describes creating input sets of SNPs for pi (nucleotide diversity) and dXY

# filtering is different than above because we need two list of both VARIABLE and TOTAL sites
# the two lists must be generated using exactly the same commands
# initial trials repeating whitelist 01 but for TOTAL sites ran into difficulties
# instead, this standardizes the commands and ensures that the filters for both sets used are exactly the same

### VARIABLE SITES

# require data for all 55 individuals right from the beginning, implement the maximum depth, no minimum MAF

FILTERS="-minQ 20 -minMapQ 20 -minInd 55 -setMinDepthInd 3 -setMaxDepth 926 -rmTriallelic 0.01 -SNP_pval 0.01"
TODO="-doCounts 1 -dumpCounts 1 -doMajorMinor 1 -doMAF 1"
angsd -b BAM_FILES_LIST.txt -GL 1 -P 1 $FILTERS $TODO -out filter02_noHetFilter

# follow the same steps above to remove the sites with excess heterozygosity, remove Z-chromosome, save final file as:

sites_bfv3_forDxy_maf00_noZ.txt

### TOTAL SITES

# here, we are more interested in the NUMBER of sites that pass the same filters as for the variable sites, rather than the actual identities of these sites themselves
# the logic will be to count the number of sites globally (all autosomes) and in smaller windows
# then, we can use the number of sites later (see pi and Dxy calculation sections) to adjust the pi and dXY estimates, rather than calculate directly from total sites (very inefficient in ANGSD)

# repeat above command exactly EXCEPT do not include a SNP filter, so that all sites will be printed (not only SNPs)

FILTERS="-minQ 20 -minMapQ 20 -minInd 55 -setMinDepthInd 3 -setMaxDepth 926 -rmTriallelic 0.01"
TODO="-doCounts 1 -dumpCounts 1 -doMajorMinor 1 -doMAF 1"
angsd -b BAM_FILES_LIST.txt -GL 1 -P 1 $FILTERS $TODO -out allSites_noHetFilter

# follow the same steps above to remove the sites with excess heterozygosity (may need to split file on multiple cores for efficiency / memory), save file as:
# note that the Z-chromosome was not yet removed (although it perhaps could have been here, as we have no need of the below info the Z-chromosome);

allSites_filtered.pos

# count the total number of sites initially per 5-kbp window

//
# python script
import sys
import math
posfile = allSites_filtered.pos
outfile = allSites_filtered_5kbpWindow_startPos1.txt
#
currentChr = "none"
currentStartPos = 1
currentCount = 0
with open(outfile, 'w') as o:
	o.write("chr\tstart\tend\tSites\n")
	with open(posfile) as f:
		next(f) # this skips the header line - assumes the file has text header
		for line in f:
			chr,pos = line.split("\t")[0:2]
			pos = int(pos)
			if (chr == currentChr):
				if (pos < currentStartPos + 5000):
					# if we are on the same chromosome and haven't stepped outside the window, add to the count
					currentCount = currentCount + 1
				else:
					# if we are on the same chromosome and HAVE stepped outside the window, print the line and start a new count
					o.write('\t'.join([str(currentChr), str(int(currentStartPos)), str(int(currentStartPos + 5000)), str(currentCount)]) + "\n")
					currentChr = chr
					currentStartPos = (math.floor(pos / 5000) * 5000) + 1
					currentCount = 1
			else:
				# if we are NOT on the same chromosome, print the line and start a new count
				o.write('\t'.join([str(currentChr), str(int(currentStartPos)), str(int(currentStartPos + 5000)), str(currentCount)]) + "\n")
				currentChr = chr
				currentStartPos = (math.floor(pos / 5000) * 5000) + 1
				currentCount = 1
		# need to print the very last line to record the final count for the last window of the last chromosome
		o.write('\t'.join([str(currentChr), str(int(currentStartPos)), str(int(currentStartPos + 5000)), str(currentCount)]) + "\n")
//

# use the outfile to count the total number of sites (that have passed filters) on the Z-chromosome and not on the Z-chromosome
# e.g., open in R, sum the total number of sites on scaffolds matching and not matching kiwi_zChr_scaffolds.txt

# total sites including Zchr: 844430544
# sites on Zchr: 40778120
# total sites excluding Zchr: 803652424
# percent sites on Zchr: 4.8

#####
##### WHITELIST 03 AS DESCRIBED IN SUPPLEMENTARY METHODS ####
#####

# whitelist 03 describes creating input sets of SNPs that are SPECIFIC TO EACH INDIVIDUAL LINEAGE, for use with PopSizeABC and RAiSD
# we want only SNPs that are truly variable in each individual lineage, to get the max number of SNPs and so as not to include invariant sites in population-specific analyes

# begin with the whitelist from 01 that has all SNPs across all 55 individuals: sites_bfv2_3x_maf00_55ind.txt
# note that although we provide a sites file that already has major-minor allele designations across all 55 kiwi individuals, we use -doMajorMinor 1 to re-identify major and minor allele designations for each population individually
# require data for all 5 individuals in this population
# call SNPs for each population individual; here an example for the population aHaast (also written out in full as australis__Haast)

FILTERS="-minMapQ 20 -minQ 20 -minInd 5 -SNP_pval 0.01 -sites sites_bfv2_3x_maf00_55ind.txt"
TODO="-doMajorMinor 1 -doMAF 1"
angsd -b BAM__australis__Haast.txt -GL 1 -P 1 $FILTERS $TODO -out australis__Haast

# create sites file specific to this population

gunzip --keep australis__Haast.mafs.gz
sed '' australis__Haast.mafs | awk '{OFS="\t"; if ($1 !="chromo"){print $1,$2,$3,$4}}' > sites_australis__Haast_bfv2_3x_maf00.txt
angsd sites index sites_australis__Haast_bfv2_3x_maf00.txt

# follow example from above in Whitelist 01 section to remove Z-chromosome, save outfile as:

sites_australis__Haast_bfv2_3x_maf00_noZ.txt

# repeat Whitelist 03 steps for remaining lineages
# these sites files are used as whitelists for PopSizeABC and RAiSD

#####
##### FILTERS FOR OTHER ANALYSES ####
#####

# other analyses do not use pre-defined whitelists
# SNP filtering (if any) for other analyses is described individually in the respective sections for each of the other analyses
