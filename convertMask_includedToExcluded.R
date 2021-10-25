##### Code is to convert the mask files from a list of included sites (default printed by bamCaller.py) to a list of excluded sites (required by SMC++). Both input and output are in BED format. #####

indir <- "out_mask_includedSites";
outdir <- "out_mask_excludedSites";

scaffold_files <- list.files(indir, pattern = "*.bed.gz");
scaffolds <- substr(scaffold_files, 14, 27);

scaffold_lengths <- read.table("kiwi_ref_genome_assembly_report.txt", header = F, comment.char = "#", stringsAsFactors = F); # reference genome;

if (sum(!(scaffolds %in% scaffold_lengths$V5))) {
	stop("ERROR: not all scaffolds present in the reference genome.");
}

#

for (i in 1:length(scaffolds)) {

	print(i);

	mask <- read.table(gzfile(paste0(indir, "/", scaffold_files[i])));

	cur_length <- scaffold_lengths$V10[scaffold_lengths$V5 == scaffolds[i]];
	if (max(mask$V3) > cur_length) {
		stop("ERROR: maximum position for scaffold ", scaffolds[i], " is ", max(mask$V3), " whereas scaffold is only ", cur_length, " bp long. Is the correct reference genome specified?");
	}

	# depending whether the very first position (index 0) is included, we would need to process the first line differently;
	if (mask$V2[1] == 0) {
		cur_end_positions_exclusion <- mask$V2[2:nrow(mask)];
		cur_start_positions_exclusion <- mask$V3;
		print(paste0("Unusual scaffold ", scaffolds[i], ": reads begin at very first bp of scaffold"));
	} else {
		cur_end_positions_exclusion <- mask$V2;
		cur_start_positions_exclusion <- c(0, mask$V3);
	}
	# finalize the last line of the exclusion mask - will depend whether the very last position is equal to the total length of the chromosome or not;
	if (tail(cur_start_positions_exclusion, 1) == cur_length) {
		# in this case, we already reached the end of the scaffold and called base positions go right up to the very final base, so we do not complete the final line to indicate any final segment to exclude;
		cur_start_position_exclusion <- cur_start_positions_exclusion[length(cur_start_positions_exclusion) - 1];
		print(paste0("Unusual scaffold ", scaffolds[i], ": reads span to very final bp of scaffold"));
	} else {
		# in this case [expect to be the majority of cases], we need to complete the final line to indicate that up until the very last bp of the scaffold needs to be excluded;
		cur_end_positions_exclusion <- c(cur_end_positions_exclusion, cur_length);
	}

	# confirm that the lengths of the start and end positions match;
	if (length(cur_start_positions_exclusion) != length(cur_end_positions_exclusion)) {
		stop("ERROR: bug in code.");
	}

	# create the unmasked df to output as the final list of regions to EXCLUDE from SMC++ analysis;
	unmask <- data.frame(chrom = rep(scaffolds[i], length(cur_start_positions_exclusion)), chromStart = cur_start_positions_exclusion, chromEnd = cur_end_positions_exclusion);

	# in order to ensure that the output is written as 100000 instead of 1e+05, which can cause errors downstream, convert to character before outputting;
	unmask[ , 2] <- format(unmask[ , 2], scientific = FALSE);
	unmask[ , 3] <- format(unmask[ , 3], scientific = FALSE);

	write.table(unmask, paste0(outdir, "/maskExcluded_", scaffolds[i], ".bed"), sep = "\t", row.names = F, col.names = F, quote = F);


}




