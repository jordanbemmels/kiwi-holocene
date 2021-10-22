#### Code is to combine the .log and .trees files of individual SNAPPER runs, while discarding an amount of burn-in (e.g., 10%) set by the user

##### THIS SECTION - USER SETS THESE OPTIONS #####

basedir <- "./snapper_11ind_m1K_l2M"; # name of the overall parent directory that contains results for each of the five independent runs in subdirectories;

subdirs <- c("run_01", "run_02", "run_03", "run_04", "run_05"); # names of the five subdirectories for the individual runs;

out_prefix <- "snapper_combined_11ind_m1K_l10M"; # desired prefix for output files
n_sims <- 2000000; # total number of sims per file - important to set correctly in order to retain the correct number of chain samples;
burnin <- 10; # burn-in percentage, can be changed if desired;
sampling_freq <- 1000; # how frequently the chain was sampled;
lines_to_skip_in_trees_file <- 32; # user needs to set - equal to the number of header lines in the .trees input files;

#### OPTIONS HAVE BEEN SET, CONTINUE #####

##### FIRST COMBINE THE LOG FILES #####

for (i in 1:length(subdirs)) {
	
	cur_log <- read.table(paste0(basedir, "/", subdirs[i], "/snapper_w_popsize.log"), sep = "", header = T, stringsAsFactors = F);
	
	if (i == 1) {		
		combined_log <- cur_log[cur_log$Sample > n_sims * burnin / 100, ];
	} else {
		combined_log <- rbind(combined_log, cur_log[cur_log$Sample > n_sims * burnin / 100, ]);
	}
	
}

combined_log$Sample <- (1:nrow(combined_log)) * 1000; # rename the Sample to be sequential numbers;

write.table(combined_log, paste0(out_prefix, ".log"), sep = "\t", quote = F, row.names = F, col.names = T);

##### NOW COMBINE THE TREE SAMPLES #####

# first just print the desired lines;

start_line <- lines_to_skip_in_trees_file;
end_line <- (n_sims / sampling_freq) + lines_to_skip_in_trees_file + 2;

for (i in 1:length(subdirs)) {
	
	intermediate_tree_command <- paste0("awk 'NR>", start_line, "&&NR<", end_line, "' '", basedir, "/", subdirs[i], "/snapper.trees' > intermediate_", i, ".txt");
	system(intermediate_tree_command);

}

# now combine them;

intermediate_combine_command <- paste0("cat ", paste0("intermediate_", 1:length(subdirs), ".txt", collapse = " "), " > intermediate_A.txt");

system(intermediate_combine_command);

# load and remove the burn-in from each of the runs;

all_trees <- read.table("intermediate_A.txt", header = F, sep = "", stringsAsFactors = F);

states <- all_trees$V2;
states <- as.numeric(substr(states, 7, nchar(states)));

all_trees <- all_trees[states > n_sims * burnin / 100, ];

# rename the sample to be sequential numbers;

for (i in 1:nrow(all_trees)) {
		all_trees$V2[i] <- paste0("STATE_", format(i * 1000, scientific = F)); 
}

# write output in NEXUS format;

initialize_NEXUS_cmd <- paste0("awk 'NR<", start_line + 1, "' '", basedir, "/", subdirs[1], "/snapper.trees' > ", out_prefix, ".trees"); # print header lines;
system(initialize_NEXUS_cmd);

write.table(all_trees, paste0(out_prefix, ".trees"), sep = " ", quote = F, row.names = F, col.names = F, append = T);

write("End;", paste0(out_prefix, ".trees"), append = T);
