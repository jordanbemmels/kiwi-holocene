# code is to take a Dxy_persite.txt file and calculate mean Dxy per sliding window, to match the sliding windows in a template file (e.g., from Fst sliding window analyses in ANGSD);
		
##############################
##############################
##############################

#####  WRITE A FUNCTION #####		

slidingWindowDxy <- function(Dxy_per_SNP, outfile) {
	
	# get a list of all chromosomes to process;	
	allChr <- unique(win$chr);		
	
	# make sure to set to NA;
	win$Dxy_uncorrected <- NA;	
	win$Dxy_variableSites <- NA;	

	for (i in 1:length(allChr)) {
		print(paste0("chromosome", i, " of ", length(allChr)));
		currentChr <- Dxy_per_SNP[Dxy_per_SNP$chromo == allChr[i], ];
		currentWin <- win[win$chr == allChr[i], ];
		for (j in 1:nrow(currentWin)) {
			# note the ">=" for the lower bound of the window and the "<" for the upper bound of the window;
			currentWin$Dxy_uncorrected[j] <- sum(currentChr$dxy[(currentChr$position >= currentWin$midPos[j] - winSize/2) & (currentChr$position < currentWin$midPos[j] + winSize/2)]);
			currentWin$Dxy_variableSites[j] <- length(currentChr$dxy[(currentChr$position >= currentWin$midPos[j] - winSize/2) & (currentChr$position < currentWin$midPos[j] + winSize/2)]);
		}
		if (i == 1) {
			write.table(currentWin, outfile, sep = "\t", append = F, col.names = T, row.names = F, quote = F);
		} else {
			write.table(currentWin, outfile, sep = "\t", append = T, col.names = F, row.names = F, quote = F);
		}
	}
	
}

##############################
##############################
##############################

#####  REAL THING #####	

# template for sliding windows;
win <- read.table("FST_DIRECTORY/aHaast_aNorthFiordland.win50k_step25k.txt", header = F, sep = "\t", stringsAsFactors = F, skip = 1); # skip = 1 to skip the annoying first column;
colnames(win) <- c("region", "chr", "midPos", "Nsites", "Fst");
win$Fst_Nsites <- win$Nsites;
win <- win[ , !(colnames(win) %in% c("region", "Fst", "Nsites"))];

winSize <- 50000; # user must adjust if different

#####

# takes a few minutes to load each file;

args <- commandArgs(trailingOnly = T);
infile = args[1];
outfile = args[2];

print(infile)
print(outfile)

temp_persite_dat <- read.table(infile, header = T, sep = "\t", stringsAsFactors = F);
slidingWindowDxy(temp_persite_dat, outfile);
