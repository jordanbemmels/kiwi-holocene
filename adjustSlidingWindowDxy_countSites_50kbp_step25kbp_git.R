# code is to adjust the estimates of dXY (sum of all individual variant sites) relative to the total number of sites (variant + invariant) in each sliding window;

# usage: Rscript adjustSlidingWindowDxy_countSites_50kbp_step25kbp_git.R input_filename output_filename
args <- commandArgs(trailingOnly = T);
infile = args[1];
outfile = args[2];

# load the template containing the number of total sites (variant + invariant) per 50-kbp window;
# note that this template was generated for aHaast and aNorthFiordland populations, but can be re-used for any combination of population pairs, as long as the exact same set of sites are being used;
template <- read.table("FST_DIRECTORY/aHaast_aNorthFiordland.win50k_step25k.countSites.txt", sep = "\t", header = T, stringsAsFactors = F);

# adjust the windowed estimates;
cur_df <- read.table(paste0("slidingWindow_win50k_step25k/", infile), sep = "\t", header = T, stringsAsFactors = F);
cur_df$totalSites <- NA;
cur_df$Dxy <- NA;	
for (j in 1:nrow(cur_df)) {		
	cur_df$totalSites[j] <- template$totalSites[template$scaffold == cur_df$chr[j] & template$midPos == cur_df$midPos[j]];
	cur_df$Dxy[j] <- cur_df$Dxy_uncorrected[j] / cur_df$totalSites[j];		
}

write.table(cur_df, outfile, row.names = F, quote = F, sep = "\t");
