# Code is to convert a .lastz.clean file from phyluce output into a regions file for use with ANGSD. The LASTZ file contains the positions of the UCE elements in the kiwi reference genome (A. rowi) based on alignment to UCE probes. We will try several combinations of various amount of flanking region, including 0bp on either side 250bp on either side, 2x500bp, 2x750bp, and 2x1000bp.;
# We will also remove any UCEs on the Z-chromosome;

lastz <- read.table("uce-5k-probes.fasta_v_aRowi_ref.lastz.clean", sep = "\t", header = F, stringsAsFactors = F);

# add the header (which was not contained in the original file);

colnames(lastz) <- c("score", "name1", "strand1", "zstart1", "end1", "length1", "name2", "strand2", "zstart2", "end2", "length2", "diff", "cigar", "identity", "percent_identity", "continuity", "percent_continuity", "coverage", "percent_coverage");

# inspect some summary info about the UCE probe alignment;

nrow(lastz); # 5,249 UCE probes were identified;

hist(lastz$length1);
hist(lastz$length2);
hist(lastz$percent_identity);
hist(lastz$percent_continuity);
hist(lastz$percent_coverage);

#########################
#########################

# remove any duplicate probes that aligned multiple places in the kiwi reference genome - only retain the best alignment (highest score) for each individual probe, as it should not match more than one location;
# in one case we have two matches with exatly the same score, so we will actually just take the first alignment that is equal to the maximum;

lastz[duplicated(lastz$name2), ] # there will be 5 duplicate alignments to remove;

duplicated_probes <- lastz$name2[duplicated(lastz$name2)];

lastz$badDuplicate <- 0;
for (i in 1:length(duplicated_probes)) {
	temp_lastz <- lastz[lastz$name2 == duplicated_probes[i], ];
	#print("##########")
	#print(temp_lastz);
	# find the row name corresponding to the (first match of the) best alignment score;
	bestRow <- rownames(temp_lastz[temp_lastz$score == max(temp_lastz$score), ])[1];
	# all the duplicates that are not the best row will be changed to bad duplicates (1);
	badRows <- rownames(temp_lastz)[rownames(temp_lastz) != bestRow];
	lastz$badDuplicate[rownames(lastz) %in% badRows] <- 1;
}

lastz <- lastz[lastz$badDuplicate != 1, ];
lastz <- lastz[ , colnames(lastz) != "badDuplicate"];

# note, however, that some probes will be tiled - from Faircloth et al. (2012) (DOI:10.1093/sysbio/sys004):
	# "If UCEs were >180 bp, we tiled 120 bp probes across target regions at 2× density (i.e., probes overlapped by 60 bp). If UCEs were <180 bp total length, we selected a single probe from the center of the UCE."
# for example, compare these two probe matches:
	# 2355 PTFB01001320.1     4228     4348
	# 2354 PTFB01001320.1     4288     4408
		# here we can see that the matches are overlapping and offset by exactly 60bp, as expected - they are not really distinct locations, but rather, it was two separate probes (probes-locus:170,probes-probe:1 and probes-locus:170,probes-probe:2) that were tiled over the same locus;
		# we will have to account for this when printing our sets of loci for the ANGSD regions file - do not want to repeat partially-overlapping regions on multiple lines;

#########################
#########################

##### create columns corresponding to the ANGSD regions files positions for each UCE, with different amounts of flanking region;

# note that the positions in the LASTZ file are origin-zero half-open, whereas ANGSD uses an origin-one closed system;
	# this means that positions 0 to 120 in LASTZ are actually the 1st to 121th positions, but we exclude the final 121th position;
	# the ANGSD equivalent of LASTZ 0 to 120 would be positions 1 to 120 (inclusive);
	# see Kiwi55_10_UCE_phylogeny.docx for more info;

lastz$region_2x0bp <- paste0(lastz$name1, ":", lastz$zstart1 + 1, "-", lastz$end1);
lastz$region_2x250bp <- paste0(lastz$name1, ":", lastz$zstart1 + 1 - 250, "-", lastz$end1 + 250);
lastz$region_2x500bp <- paste0(lastz$name1, ":", lastz$zstart1 + 1 - 500, "-", lastz$end1 + 500);
lastz$region_2x750bp <- paste0(lastz$name1, ":", lastz$zstart1 + 1 - 750, "-", lastz$end1 + 750);
lastz$region_2x1000bp <- paste0(lastz$name1, ":", lastz$zstart1 + 1 - 1000, "-", lastz$end1 + 1000);

# note that if there is flanking region, we must not request a position below a value of 1 (the first position);
lastz$region_2x250bp[lastz$zstart1 + 1 - 250 < 1] <- paste0(lastz$name1[lastz$zstart1 + 1 - 250 < 1], ":", 1, "-", lastz$end1[lastz$zstart1 + 1 - 250 < 1] + 250);
lastz$region_2x500bp[lastz$zstart1 + 1 - 500 < 1] <- paste0(lastz$name1[lastz$zstart1 + 1 - 500 < 1], ":", 1, "-", lastz$end1[lastz$zstart1 + 1 - 500 < 1] + 500);
lastz$region_2x750bp[lastz$zstart1 + 1 - 750 < 1] <- paste0(lastz$name1[lastz$zstart1 + 1 - 750 < 1], ":", 1, "-", lastz$end1[lastz$zstart1 + 1 - 750 < 1] + 750);
lastz$region_2x1000bp[lastz$zstart1 + 1 - 1000 < 1] <- paste0(lastz$name1[lastz$zstart1 + 1 - 1000 < 1], ":", 1, "-", lastz$end1[lastz$zstart1 + 1 - 1000 < 1] + 1000);

lastz[lastz$zstart1 + 1 - 250 < 1, ];
lastz[lastz$zstart1 + 1 - 500 < 1, ];
lastz[lastz$zstart1 + 1 - 750 < 1, ]; # there were originally two matches, but both where duplicates and removed above!;
lastz[lastz$zstart1 + 1 - 1000 < 1, ]; # so, in a sense, this became irrelevant here!;

#########################
#########################

##### remove any regions on the Z-chromosome;

zChr <- read.table("kiwi_zChr_scaffolds.txt", stringsAsFactors = F);

lastz <- lastz[!(lastz$name1 %in% zChr$V1), ];

#########################
#########################

##### for each set of UCE regions (0bp flanking, 250bp flanking, etc.), create a sorted regions dataframe where all regions are non-overlapping;
	# see above: we expect some regions to be overlapping as they are from tiled probes over the same locus, but even some additional regions may overlap if there were two distinct loci that happened to be located very close to one another and the flanking regions start overlapping as well;

# strategy: (1) sort, (2) identify overlapping regions, (3) merge overlapping regions into a new dataframe;

# write a generic function

mergeUCEs <- function(uce) {
	
	uce$chr <- NA;
	uce$startPos <- NA;
	uce$endPos <- NA;
	
	for (i in 1:nrow(uce)) {
		
		uce$chr[i] <- strsplit(uce[i , 1], "[:-]")[[1]][1];
		uce$startPos[i] <- as.numeric(strsplit(uce[i , 1], "[:-]")[[1]][2]);
		uce$endPos[i] <- as.numeric(strsplit(uce[i , 1], "[:-]")[[1]][3]);
		
	}
	
	##### sort the uce dataframe;
	
	uce <- uce[order(uce$chr, uce$startPos), ];
	
	##### proceed through each of the overlapping regions and collect all the lines in a row that overlap;
	
	temp_uce <- uce[1, ]; # initialize a temp_uce dataframe with the first line;
	out_uce <- data.frame(angsd = character(0), chr = character(0), startPos = numeric(0), endPos = numeric(0)); # initialize an output dataframe;
	
	for (i in 2:nrow(uce)) {
		
		if ((uce$chr[i] != unique(temp_uce$chr)) | (uce$chr[i] == unique(temp_uce$chr) & uce$startPos[i] > temp_uce$endPos[nrow(temp_uce)])) {
			
			# if the new line is NOT on the same chromosome OR it is but is non-overlapping, we can merge and add the previous info to our output df (without incorporating the new line in any way);
			
			new_output_line <- data.frame(angsd = NA, chr = unique(temp_uce$chr), startPos = min(temp_uce$startPos), endPos = max(temp_uce$endPos));
			new_output_line$angsd <- paste0(new_output_line$chr, ":", new_output_line$startPos, "-", new_output_line$endPos);
			out_uce <- rbind(out_uce, new_output_line);
			
			# now, start a new temp_uce with the new line that was NOT an overlapping match;
			temp_uce <- uce[i, ];
			
		} else {
			
			# if the new line IS on the same chromosome AND it is overlapping, we will not merge anything just yet but add it to our growing temp dataframe for future processing;
			
			temp_uce <- rbind(temp_uce, uce[i, ]);
			
		}
		
		if (i == nrow(uce)) {
			
			# if it's the very last row, we need to add whatever we have got (after we have performed the above to see if our temp_uce will be just the last line or will be the last line plus whatever came before and hasn't yet been processed);
			new_output_line <- data.frame(angsd = NA, chr = unique(temp_uce$chr), startPos = min(temp_uce$startPos), endPos = max(temp_uce$endPos));
			new_output_line$angsd <- paste0(new_output_line$chr, ":", new_output_line$startPos, "-", new_output_line$endPos);
			out_uce <- rbind(out_uce, new_output_line);

		}
		
	}
	
	return(out_uce);
	
}

df_2x0bp <- mergeUCEs(as.data.frame(lastz$region_2x0bp, stringsAsFactors = F));
df_2x250bp <- mergeUCEs(as.data.frame(lastz$region_2x250bp, stringsAsFactors = F));
df_2x500bp <- mergeUCEs(as.data.frame(lastz$region_2x500bp, stringsAsFactors = F));
df_2x750bp <- mergeUCEs(as.data.frame(lastz$region_2x750bp, stringsAsFactors = F));
df_2x1000bp <- mergeUCEs(as.data.frame(lastz$region_2x1000bp, stringsAsFactors = F));

nrow(lastz);
nrow(df_2x0bp);
nrow(df_2x250bp);
nrow(df_2x500bp);
nrow(df_2x750bp);
nrow(df_2x1000bp);

###############
###############

# write the output in angsd regions-file format;

write.table(df_2x0bp$angsd, "aRowi_ref_UCEs_plus_2x0bp.txt", sep = "\t", row.names = F, col.names = F, quote = F);
write.table(df_2x250bp$angsd, "aRowi_ref_UCEs_plus_2x250bp.txt", sep = "\t", row.names = F, col.names = F, quote = F);
write.table(df_2x500bp$angsd, "aRowi_ref_UCEs_plus_2x500bp.txt", sep = "\t", row.names = F, col.names = F, quote = F);
write.table(df_2x750bp$angsd, "aRowi_ref_UCEs_plus_2x750bp.txt", sep = "\t", row.names = F, col.names = F, quote = F);
write.table(df_2x1000bp$angsd, "aRowi_ref_UCEs_plus_2x1000bp.txt", sep = "\t", row.names = F, col.names = F, quote = F);
