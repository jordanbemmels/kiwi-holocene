require(abc);

##### load observed summary statistics;

# the purpose of [ , 1] is to convert to a single vector format rather than a data frame;
obs_aHaast <- read.table("aHaast_n10_mac1_macld1.stat")[ , 1];

#########################
#########################
#########################

##### load the simulations - summary statics and paramters;
# simulations from independent batches should be combined into a single file before loading here;

sstat <- read.table("combined__L4_n10_s100_mac1_macld1.stat", sep = " ", header = F, stringsAsFactors = F);
params <- read.table("combined__L4_n10_s100.params", sep = " ", header = F, stringsAsFactors = F);

# check that the number of rows is correct and equal between models, and that there is no missing data;

nrow(sstat);
nrow(params);

# check if there are any sstat rows with any NA values anywhere;
sum(!(complete.cases(sstat)));

##### systemtically check every row – there may be some additional cases where insufficient sstats were generated;
apply(sstat, 2, function(x) {sum(x == -1)});

# from the above, in my empirical data, there were 2 problematic examples where columns V45 and V20 do not have proper sstats - for all other columns the sstats are completed because there are no NAs above and no -1 (missing data) values;
# find these rows;
sstat[sstat$V20 == -1, ];
sstat[sstat$V45 == -1, ];
# they are identical (makes sense, as V45 is the number of SNP pairs in the most ancient / smallest distance bin and V20 is the LD between those pairs, so of course if V45 is missing then V20 will also be missing);
# IMPORTANT: the column identities mentioned here may not match those when running example scripts from GitHub - users should adjust the column names as necessary!;

# remove these two rows from BOTH the params and sstat (very important to remove both consistently to maintain order of df);
bad_rows <- which(sstat$V20 == -1);
sstat <- sstat[-1*bad_rows, ];
params <- params[-1*bad_rows, ];

nrow(sstat);
nrow(params);

##### SELECT EXACTLY 500,000 SIMULATIONS ###

# we would like to select exactly 500,000 sims, out of the good / completed ones;
# this assumes there were slightly more than 500,000 good sims completed to choose from;

sstat <- sstat[1:500000, ];
params <- params[1:500000, ];

nrow(sstat);
nrow(params);

# save these for future reference so that we will always use exactly the same set of 500,000 simulations for all analyses we might do, and not need to resort;

# only perform once;
#write.table(sstat, "combined__L4_n10_s100_mac1_macld1_filtered500K.stat", sep = "\t", row.names = F, col.names = F, quote = F);
#write.table(params, "combined__L4_n10_s100_filtered500K.params", sep = "\t", row.names = F, col.names = F, quote = F);

#########################
#########################
#########################

##### Reload the filtered 500,000 datasets so as not to have to repeat filtering steps above;

sstat <- read.table("combined__L4_n10_s100_mac1_macld1_filtered500K.stat", sep = "\t", header = F, stringsAsFactors = F);
params <- read.table("combined__L4_n10_s100_filtered500K.params", sep = "\t", header = F, stringsAsFactors = F);

head(sstat);
head(params);

nrow(sstat);
nrow(params);

apply(sstat, 2, function(x) {sum(x == -1)});
apply(sstat, 2, function(x) {sum(is.na(x))});

#########################
#########################
#########################

##### Select the summary stats #####

# this is done by formulae in the abc estim files with PopSizeABC, but we can specify it manually here;
	# the values can be double-checked by running manually the abc_ex1_[...].R code on the lab server;
# the first 6 columns are the SFS (V1-V6);
	# however, the very first column of this is the value for zero, which is just noise in these simulations and observed summary stats as we did not actually intend to have any SNPs included with zero frequency, so we would ignore that one;
# the next 14 columns are the LD statistics (V7-V20);
# the next 11 columns are the IBS statsitics, which will not be used (V21-V31);
# the final 14 columns are the number of SNP pairs in each distance bin, which was just for the sake of info and is not a summary stat (V32-V45);
  # IMPORTANT: THIS INFO IS NOT PRINTED IN EXAMPLE SCRIPTS PROVIDED ON GIT
# 6 + 14 + 11 equals 31 statistics total - we are interested in only the first 6 + 14 = 20, excluding the very first one, thus 2:20;

head(sstat);

ind_stat <- 2:20; # IMPORTANT: USER MUST INSPECT AND SET COLUMN INDICES MANUALLY - MAY BE DIFFERENT THAN SHOWN HERE;

#########################
#########################
#########################

##### Select the model parameters #####

# all parameters are variable exept for V1 and V2 (the mutation and recombination rates, respectively);
# V1 = mutation rate;
# V2 = recombination rate;
# V3-V16 = population sizes in different time bins, with V3 being the current population size at the present;
# V17 = ancestral population size;

# note that there are 15 population sizes in the params, but only 14 LD bins in the sstat: this is because for the most recent / largest-distance time bin, we could not simulate chromosomes long enough in a meaningful time, so there is no category for LD corresponding to this most recent time bin;

head(params);

ind_param <- c(3:17); # IMPORTANT: USER MUST INSPECT AND SET COLUMN INDICES MANUALLY - MAY BE DIFFERENT THAN SHOWN HERE;

# note that these are all population sizes, and all need to be log10-transformed prior to use;

#########################
#########################
#########################

##### Parameter inference #####

##### prepare parameters and summary statistics;

# note that we need to log10-transform all of the population size parameters (all of our selected parameters are population size in ind_param) (as in PopSizeABC);

params_reduced <- params[ , ind_param];
params_reduced <- log10(params_reduced);

##### perform the parameter inference;

res.nnet_aHaast <- abc(obs_aHaast[ind_stat], params_reduced, sstat[ , ind_stat], tol = 0.001, method = "neuralnet", numnet = 100);

# commented out to protect from overwriting files;
#save(res.nnet_aHaast, file = "abc_output/res.nnet_aHaast.R");

load("abc_output/res.nnet_aHaast.R");

summary(res.nnet_aHaast);

plot(res.nnet_aHaast, param = params_reduced);

summary(res.nnet_aHaast, print = F)[3,]; # median;
summary(res.nnet_aHaast, print = F, intvl = 0.9)[2,]; # 5% quantile (NOT the 2.5%);
summary(res.nnet_aHaast, print = F, intvl = 0.9)[6,]; # 95% quantile (NOT the 97.5%);

##############################
##############################
##############################

##### Cross-validation of ability to estimate parameters #####

# see ?cv4abc() and note that this function can be used to select the appropriate tolerance rate for parameter inference;

# we should use nnet as that is what is used in the main model;
# res.nnet_aHaast <- abc(obs_aHaast[ind_stat], params_reduced, sstat[ , ind_stat], tol = 0.001, method = "neuralnet", numnet = 100);

# extremely slow;
# note that it is not possible to actually change numnet - I tried many times trying to get "numnet = 100" to match the number of nnet regressions used for the empirical data, but it appears that this actually has no effect;
	# using numnet = 10 is certainly less time-consuming and using nnet = 10 instead of 100 could only cause us to overestimate our uncertainty not to underestimate it;
# use 200 replicates;
cv.res.nnet <- cv4abc(param = params_reduced, sumstat = sstat[ , ind_stat], tols = c(0.001), method = "neuralnet", statistic = "median",  nval = 200);
# timing, ~40-60min;

# overwrite-protected;
#save(cv.res.nnet, file = "abc_output/cv.res.nnet.R");

load("abc_output/cv.res.nnet.R");

# prediction error;

summary(cv.res.nnet);

# optional - rename the parameter names from V3-V17 to something less nonsensical;

cv.res.nnet$names$parameter.names <- paste0(round(10**years), "-", c(round(10**years[2:length(years)]) - 1, ""), " YBP");

##############################
##############################
##############################

##### Print a text file of the 500 (tol = 0.001) retained simulations, to use in simulations of the "neutral" expectation for RAiSD when detecting selective sweeps. Use the ADJUSTED rather than the unadjusted values, because the unadjusted values are heavily influenced by the (non-uniform) prior around the prior limits. In addition, CONSTRAIN the adjusted values to be within the prior limits, so that we never extrapolate beyond the prior;

# write a function to constrain them for a given set of retained simulations;

constrain_retained <- function(retainedSims) {
	
	for (i in 1:nrow(retainedSims)) {
		for (j in 1:ncol(retainedSims)) {
			if (retainedSims[i, j] < log10(20)) {
				retainedSims[i, j] <- log10(20);
			} else if (retainedSims[i, j] > log10(2000000)) {
				retainedSims[i, j] <- log10(2000000);
			}
		}
	}
	
	# for the final column (V17) corresponding to the most ancient time period, the max value of the prior is different;
	for (k in 1:nrow(retainedSims)) {
		if (retainedSims[k, ncol(retainedSims)] > log10(50000)){
			retainedSims[k, ncol(retainedSims)] <- log10(50000);
		}		
	}
	
	return(retainedSims);
	
}

retained_aHaast <- constrain_retained(res.nnet_aHaast$adj.values);

write.table(retained_aHaast, "abc_output/retainedSims_adj/retainedSims_aHaast.txt", sep = "\t", quote = F, row.names = F, col.names = T);

# sanity check - plot the median values;
plot(apply(retained_aHaast, 2, median), ylim = c(log10(20), log10(2000000)));
