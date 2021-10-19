##### Code is to remove scaffolds on the kiwi Z-chromosome from any standard-formatted ANGSD sites file. Since thinning is on a scaffold basis, this code can be equally applied to either a thinned or un-thinned sites file; there is no need to remove sites first then re-thin. #####

# usage: Rscript remove_kiwiZChr_fromSites.R infile
# the outfile will automatically be named the same as the infile except instead of ending in .txt it will now be _noZ.txt

# no need to specify the Z-chromosome file in command line, the same standardized file will always be loaded:
zChrScaffolds <- read.table("kiwi_zChr_scaffolds.txt");
zChrScaffolds <- zChrScaffolds$V1;

# identify the input sites files;

args <- commandArgs(trailingOnly = T);

infile <- args[1];

# identify the name for the output file;

outfile <- paste0(substr(infile, 1, nchar(infile) - 4), "_noZ.txt");

# read the sites file line by line and process;

con  <- file(infile, open = "r");
conOut <- file(outfile, open = "w");

count = 0;
while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
	if (count %% 100000 == 0) {
		print(paste0("working on line ", count));
	}
	#
	temp_line <- (strsplit(oneLine, "\t"));
	if (temp_line[[1]][1] %in% zChrScaffolds) {
		# do not print the line;
		count = count + 1;
	} else {
		# print the line if the scaffolds is not on the Z-chromosome;
		writeLines(oneLine, conOut);
		count = count + 1;
	}
} 

close(con);
close(conOut);

