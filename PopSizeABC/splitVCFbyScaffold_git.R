# code is to split a single whole-genome vcf file for a kiwi population into separate vcf files for each scaffold;
# note that a .vcf and a .vcf.gz are BOTH required (also the tabix .tbi file for the .gz);

args <- commandArgs(trailingOnly = T);

infile <- args[1];
outdir <- args[2];

#print(infile)
#print(outdir)

# get a list of all the scaffolds - but remove any with <2 SNPs;

vcf <- read.table(infile, sep = "\t", header = F, comment.char = "#", stringsAsFactors = F);

scaffolds <- unique(vcf$V1);

n_snps_per_scaffolds <- data.frame(scaffold = scaffolds, n_snps = NA);
for (i in 1:length(scaffolds)) {
	n_snps_per_scaffolds$n_snps[i] <- sum(vcf$V1 == scaffolds[i]);
}

scaffolds2 <- scaffolds[n_snps_per_scaffolds$n_snps >= 2];

# split the vcf into each of scaffolds2 using bcftools;

for (i in 1:length(scaffolds2)) {

	print(paste0("working on scaffold ", i, " of ", length(scaffolds2)));
	splitCommand <- paste0("bcftools view -r ", scaffolds2[i], " ", infile, ".gz > ", outdir, "/", outdir, "_", scaffolds2[i], ".vcf");
	system(splitCommand);

}
