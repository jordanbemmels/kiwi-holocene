##### Code is to determine which kiwi scaffolds are on the Z-chromosome;

mummer <- read.table("ostrich_kiwi_q_i50_l500_coords.txt", header = F, skip = 4, sep = "\t", stringsAsFactors = F);
colnames(mummer) <- c("S1", "E1", "S2", "E2", "LEN_1", "LEN_2", "IDY", "SIM", "STP", "LEN_R", "LEN_Q", "FRM_1", "FRM_2", "TAGS_1", "TAGS_2");

ostrichZ <- read.table("ostrich_Zscaffolds_withNames.txt", sep = "\t", header = T, stringsAsFactors = F);

#########

# create a dataframe with a single line for each kiwi scaffold, and for each scaffold, count the total length of alignments (that have already been filtered for quality) that map to an ostrich Z-chromosome scaffold vs. to an ostrich putatively autosomal scaffold;

kiwi <- data.frame(scf = unique(mummer$TAGS_2), len_alignOstrichZ = 0, len_alignOstrichAuto = 0);
kiwi$scf <- kiwi$scf[order(kiwi$scf)];

for (i in 1:nrow(kiwi)) {
	print(paste0("working on kiwi scaffold ", i, " of ", nrow(kiwi)));
	#
	# reduce mummer to matches for the current kiwi scaffold only, and divide into matches tentatively corresponding to ostrich Zchr and to putative ostrich autozomes;
	mummer_red <- mummer[mummer$TAGS_2 == kiwi$scf[i], ];
	mummer_Z <- mummer_red[mummer_red$TAGS_1 %in% ostrichZ$chr, ]; # tentativley matches on Z chromosome, we have to check the position still (especially for scaffold9);
	mummer_auto <- mummer_red[!(mummer_red$TAGS_1 %in% ostrichZ$chr), ];
	#
	# for all the Z-matches, check if the midpoint of the ostrich chromosome falls within the min and max pos on the ostrich scaffolds that are truly on the ostrich Z (in reality, scaffold9 is the only scaffold where this could be an issue, because it was likely misassembled and a small portion is not on the Z-chromosome);
	if (nrow(mummer_Z) > 0) {
		mummer_Z$MIDPOS_1 <- (mummer_Z$S1 + mummer_Z$E1) / 2;
		mummer_Z$onZ <- NA;
		for (j in 1:nrow(mummer_Z)) {
			if (mummer_Z$MIDPOS_1[j] > ostrichZ$minPos[ostrichZ$chr == mummer_Z$TAGS_1[j]] & mummer_Z$MIDPOS_1[j] <= ostrichZ$maxPos[ostrichZ$chr == mummer_Z$TAGS_1[j]]) {
				mummer_Z$onZ <- 1;
			} else {
				mummer_Z$onZ <- 0;
			}
		}
		kiwi$len_alignOstrichZ[i] <- kiwi$len_alignOstrichZ[i] + sum(mummer_Z$LEN_2[mummer_Z$onZ == 1]); # add the length of the kiwi genome that is truly on the Z chromosome;
		kiwi$len_alignOstrichAuto[i] <- kiwi$len_alignOstrichAuto[i] + sum(mummer_Z$LEN_2[mummer_Z$onZ == 0]); # add the length of the kiwi genome that was on a chromosome that partly matches the Z, but it turned out this particular segment of the chromosome is not actually on the Z chromosome;		
	} else {
		kiwi$len_alignOstrichZ[i] <- kiwi$len_alignOstrichZ[i] + 0;
	}
	#
	# add the rest of the length of the kiwi genome that is on a putative autosome;
	if (nrow(mummer_auto) > 0) {
		kiwi$len_alignOstrichAuto[i] <- kiwi$len_alignOstrichAuto[i] + sum(mummer_auto$LEN_2);
	} else {
		kiwi$len_alignOstrichAuto[i] <- kiwi$len_alignOstrichAuto[i] + 0;
	}
}

# decide which scaffolds are on the Z-chromosome;
# FINAL CRITERIA: A KIWI SCAFFOLD IS CONSIERED PUTATIELY Z-CHROMOSOME IF ≥50% OF TOTAL ALIGNEMENT LENGTH CORRESPONDS TO THE OSTRICH Z-CHROMOSOME, OR IF TOTAL ALIGNMENT LENGTH TO OSTRICH Z-CHROMOSOME IS ≥500,000bp;

for (i in 1:nrow(kiwi)) {
	kiwi$percent_alignOstrichZ[i] <- kiwi$len_alignOstrichZ[i] / (sum(kiwi$len_alignOstrichZ[i], kiwi$len_alignOstrichAuto[i])) * 100;
}

kiwi_zChr <- kiwi[kiwi$percent_alignOstrichZ >= 50 | kiwi$len_alignOstrichZ >= 500000, ];
sum(kiwi_zChr$len_alignOstrichZ, kiwi_zChr$len_alignOstrichAuto); # kiwi Z-chr inferred to be 54.8 million bp;

##########

# output;

write.table(kiwi, "kiwi_ostrichZchr_matches.txt", sep = "\t", row.names = F, quote = F);
write.table(kiwi_zChr$scf, "kiwi_zChr_scaffolds.txt", sep = "\t", row.names = F, col.names = F, quote = F);
