sites <- read.table("allSites_filtered_5kbpWindow_startPos1.txt", sep = "\t", header = T, stringsAsFactors = F);

# LOAD FILE;

pairwise <- read.table("FST_DIRECTORY/aHaast_aNorthFiordland.win50k_step25k.txt", sep = "\t", skip = 1, stringsAsFactors = F);
pairwise <- pairwise[ , 2:ncol(pairwise)]; # remove the row names, which are not useful and get in the way;
colnames(pairwise) <- c("scaffold", "midPos", "Nsites", "Fst");

# CALCULATE TOTAL NUMBER OF SITES PER WINDOW - 50-kbp WINDOWS WITH 25-kbp STEP SIZE;

pairwise$totalSites <- NA;
for (k in 1:nrow(pairwise)) {
	pairwise$totalSites[k] <- sum(sites$Sites[sites$chr == pairwise$scaffold[k] & sites$start >= (pairwise$midPos[k] - 25000) & sites$end <= (pairwise$midPos[k] + 25000)]);
}

write.table(pairwise, "FST_DIRECTORY/aHaast_aNorthFiordland.win50k_step25k.countSites.txt", sep = "\t", row.names = F, col.names = T, quote = F);
