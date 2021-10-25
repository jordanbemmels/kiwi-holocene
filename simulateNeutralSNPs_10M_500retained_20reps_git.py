#!/usr/bin/python
import sys
import os
import numpy as np
import copy
import msprime
import itertools

##########
##########
##########

pop_id = int(sys.argv[1]) # pop_id ranges from 1-11 for each of the 11 populations in alphabetical order
#n_rep = int(sys.argv[2]) # the number of independent chromosomes to simulate # UPDATE: for 500retained, we will always be simulationg 1000 chromosomes (2x500 different parameter combinations)
reps_per_param_set = 20 # 20reps version will simulate 20 replicates per each of the 500 retained demographic scenarios, and has the chromosome naming convention modified to be 1_1, 1_2, 1_3, ... 1_20, 2_1, 2_2, 2_3 ... etc., so that we can keep track for the scenario and replicate (scenario_replicate);

##### Set up input data that will be the same across all populations, v8_0 #####

n=10 # haploid sample size
L=10000000 # size of each segment, in bp
r_fixed = 2.1*(10**(-8)) # recombination rate used in v8_0
mmu=1.345868*10**(-8) # per generation per bp mutation rate used in v8_0

# time windows
nb_times=15 # number of time windows
Tmax=622 # the oldest time window will start at Tmax
a=0.09 # the length of time windows increases when time increases, at a speed that is determined by this coefficient  
# computation of time windows based on the above parameters
times=-np.ones(shape=nb_times,dtype='float')
for i in range(nb_times):
    times[i]=(np.exp(np.log(1+a*Tmax)*i/(nb_times-1))-1)/a
print "Population size changes at the following times (in generations):"
print times
print ""

# complete a params object in the same format as in PopSizeABC
params_generic=-np.ones(shape=[1,nb_times+2],dtype='float')
params_generic[:,0]=mmu
params_generic[:,1]=r_fixed

##########
##########
##########

##### Add the population size changes that are specific to the particular population #####

##### For 500retained.py version, we load the desired parameters of the 500 retained simulations from the PopSizeABC nnet regressions (adjusted values; constrained to be within the prior) #####

retained_filenames = np.array(["PopSizeABC_v8_0_retainedSims_adj/retainedSims_aHaast.txt",
	"PopSizeABC_v8_0_retainedSims_adj/retainedSims_aNorthFiordland.txt",
	"PopSizeABC_v8_0_retainedSims_adj/retainedSims_aSouthFiordland.txt",
	"PopSizeABC_v8_0_retainedSims_adj/retainedSims_aStewartIsland.txt",
	"PopSizeABC_v8_0_retainedSims_adj/retainedSims_hHaastii.txt",
	"PopSizeABC_v8_0_retainedSims_adj/retainedSims_mCoromandel.txt",
	"PopSizeABC_v8_0_retainedSims_adj/retainedSims_mEastern.txt",
	"PopSizeABC_v8_0_retainedSims_adj/retainedSims_mNorthland.txt",
	"PopSizeABC_v8_0_retainedSims_adj/retainedSims_mTaranaki.txt",
	"PopSizeABC_v8_0_retainedSims_adj/retainedSims_oKapiti.txt",
	"PopSizeABC_v8_0_retainedSims_adj/retainedSims_rOkarito.txt"])

# use [pop_id - 1] and the retained_filenames array to determine which of the above files to load;
popSizes = np.loadtxt(retained_filenames[pop_id - 1], delimiter="\t", skiprows = 1)

##########
##########
##########

##### Define a function to perform the simulation #####

def simul_neutral_SNPs(nb_seg,L,n,times,params):
    '''
    Simulates one sample with msprime and returns a tskit.TreeSequence object.
    JBB: function is based on simul_stats_one_rep_macld() from PopSizeABC popgen_abc.py.
    '''
    # creates demographic events
    mymut=params[0]
    myrec=params[1]
    N=params[2]
    mydemo=[]
    for i in range(1,len(times)):
	mydemo.append(msprime.PopulationParametersChange(time=times[i],initial_size=params[i+2]))
    # debugging of the demography - note that the initial population size shows up as "1", but this should not be a problem, because we specify it below when simulating tree_sequence with "Ne=N", and it would show up when printing the demography debugger because this only shows CHANGES from the initial size, not the initial size itself - I double-checked and manually changing Ne=N in msprime.simulate() to different values has a huge effect on the number of SNPs simulated, as expected for the most recent time bin, suggesting that setting the initial size with Ne=N below is indeed working as expected;
    #dd = msprime.DemographyDebugger(demographic_events=mydemo)
    #dd.print_history()
    # simulate each segment
	# CURRENTLY ONLY SIMULATING A SINGLE SEGMENT
    for i in range(nb_seg):
	tree_sequence = msprime.simulate(sample_size=n, Ne=N, length=L, recombination_rate=myrec, mutation_rate=mymut, demographic_events=mydemo)
	p=tree_sequence.num_mutations
	print "Number of SNP mutations:",p
    return tree_sequence

##########
##########
##########

##### Perform the simulations #####

# remove any existing .vcf fils in the temporary holding directory, so that there is no possibility of somehow using old ones from a previous run
os.system("rm sim_chr/temp_holding_dir/*.vcf")

for i in range(len(popSizes)):

	##### for the 500retained.py versions, we have not yet defined the parameters, and they will be different for each replicate, so we will define them here;
	
	# need to use copy.deepcopy() to create an entirely new object, instead of simply assigning one to equal the other (otherwise the two variables refer back to the same object and changing one also changes the other)
	params_currentRep=copy.deepcopy(params_generic)
	params_currentRep[:,2:(nb_times+2)]=10**popSizes[i]

	##### simulate a single chromosome, as many times as needed per parameter set #####
	for k in range(reps_per_param_set):
		tree_sequence=simul_neutral_SNPs(nb_seg=1,L=L,n=n,times=times,params=params_currentRep[0,:])
		with open("sim_chr/temp_holding_dir/rep"+str(i)+"_"+str(k)+".vcf", "w") as vcf_file:
			tree_sequence.write_vcf(vcf_file, ploidy=2) 
		##### it is not possible to concat two vcfs properly with identical chromosome names (assuming they are not truly supposed to be identical), and it is not possible to output the vcf using tree_sequence.write_vcf() with different chromsome names right from the beginning, so we will have to use awk to modify the chromosome names in the vcf files before concat-ing them - note that there is an extra "_1" appended to all chromosomes names, i.e., 1_1_1, 1_2_1, 1_3_1, ..., 1_20_1, 2_1_1, 2_1_1, ..., 2_20_1, 3_1_1, etc.;
		# pseudocode for the awk command
		# if it's the header line defining the chromosome name:
			# replace the chromosome name with a unique prefix corresponding to the replicate number
		# else if it's any other header line:
			# print as normal with no changes
		# else if it's a non-header line:
			# append the same unique prefix to match the renamed chromosome in the header line above
		os.system("awk '{if($0 == \"##contig=<ID=1,length="+str(L)+">\") print \"##contig=<ID=chr"+str(i)+"_"+str(k)+"_1,length="+str(L)+">\"; else if($0 !~ /^#/) print \"chr"+str(i)+"_"+str(k)+"_\"$0; else print $0}' sim_chr/temp_holding_dir/rep"+str(i)+"_"+str(k)+".vcf > sim_chr/temp_holding_dir/chr_rep"+str(i)+"_"+str(k)+".vcf")

##### All simulations are completed now #####

##### concatenate the vcf files into a single vcf

names_allpops = ["aHaast", "aNorthFiordland", "aSouthFiordland", "aStewartIsland", "hHaastii", "mCoromandel", "mEastern", "mNorthland", "mTaranaki", "oKapiti", "rOkarito"]

# this did not work for the 20reps version, when calling os.system it returned error code 32512;
#concat_fileList = ' '.join(["sim_chr/temp_holding_dir/chr_rep" + str(sub1) + "_" + str(sub2) +".vcf" for sub1,sub2 in itertools.product(range(len(popSizes)), range(reps_per_param_set))]);
#print concat_fileList
#os.system("/home/0_PROGRAMS/bcftools-1.10.2/bcftools concat -o sim_chr/"+names_allpops[pop_id - 1]+"_neutral_10M_500retained_20reps.vcf "+concat_fileList)

# instead, concatenate all the replicates for a single retained demographic scenario (we will have 500 total), then later we will need to modify the pipeline to run RAiSD on EACH of the 500 files, rather than a single RAiSD run across one gigantic concatenated dataset;
os.system("mkdir "+"sim_chr/neutral_10M_500retained_20reps/");
os.system("mkdir "+"sim_chr/neutral_10M_500retained_20reps/"+names_allpops[pop_id - 1]);
for i in range(len(popSizes)):
	concat_fileList = ' '.join(["sim_chr/temp_holding_dir/chr_rep" + str(i) + "_" + str(rep_index) +".vcf" for rep_index in range(reps_per_param_set)]);
	os.system("bcftools concat -o sim_chr/neutral_10M_500retained_20reps/"+names_allpops[pop_id - 1]+"/"+names_allpops[pop_id - 1]+"_neutral_10M_500retained_20reps_scenario"+str(i)+".vcf "+concat_fileList);

