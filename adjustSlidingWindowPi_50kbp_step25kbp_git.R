# Code is to to adjust the pi estimate relative to the number of total sites (variant + invariant) per 50-kbp sliding window;

# usage: Rscript adjustSlidingWindowPi_50kbp_step25kbp_git.R input_filename output_filename
args <- commandArgs(trailingOnly = T);
infile = args[1];
outfile = args[2];

# Template file with the number of sites total per 50-kbp window. This was calculated for aHaast and aNorthFiordland, but the population identity does not matter and the same template can be used for any pair of populations, provided the same number of sites is used;
template <- read.table("FST_DIRECTORY/aHaast_aNorthFiordland.win50k_step25k.countSites.txt", sep = "\t", header = T, stringsAsFactors = F);

# adjust pi;
df <- read.table(infile, sep = "\t", header = T, comment.char = "", stringsAsFactors = F);	
df$totalSites <- NA;
for (j in 1:nrow(df)) {
	if (length(template$totalSites[template$scaffold == df$Chr[j] & template$midPos == df$WinCenter[j]]) > 0) {
		# the if statement is needed because if the number of total sites is unknown, then we cannot update the value;
		df$totalSites[j] <- template$totalSites[template$scaffold == df$Chr[j] & template$midPos == df$WinCenter[j]];
	}
}
df$pi <- df$tP / df$totalSites;
df <- df[ , 2:ncol(df)]; # remove the first column, which is formatting info we don't need to keep;
write.table(df, outfile, row.names = F, quote = F, sep = "\t");
