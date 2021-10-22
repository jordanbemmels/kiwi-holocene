# Demographic history as estimated in PopSizeABC

We will use [PopSizeABC](https://forge-dga.jouy.inra.fr/projects/popsizeabc) to infer changes in effective population size (Ne) in each population, focusing on the Holocene.

## Optional modifications to PopSizeABC code

PopSizeABC is written in python code. In the kiwi analysis, I made several minor and mostly optional modifications to the code. However, PopSizeABC is not my intellectual property and I do not have permissions to publically share PopSizeABC code here. Instead, for reproducibility of the kiwi project I will describe the changes that any user could also make to the PopSizeABC code if desired.

###### Modifications to gwas/summary_stat.py

Add hashes (```#```) at the beginning of the lines to comment out the entire if-statement beginning with ```if len(r2_list) < 2:``` (lines 160-182, i.e., where the original code indicates to ```# try a more exhaustive screening of SNP pairs```). This statement is removed because for a given distance bin, if there are fewer than 2 SNP pairs available, it would calculate the linkage disequilibrium (LD) across *every* possible SNP pair that meets the distance criteria, rather than ensuring that the SNP pairs are physically non-overlapping (the normal behaviour). Thus, the SNP pairs may be non-overlapping and non-independent. I prefer not to do this. If these lines are removed, then the LD mean and standard deviation for that distance bin will not be updated and will remain as -1. This is fine as PopSizeABC will recognize that something was wrong with these simulations (very rare), and they will fail the good_sims test and be discarded. 

