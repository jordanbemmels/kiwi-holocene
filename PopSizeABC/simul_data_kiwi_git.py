# This script simulates genomic samples, and the corresponding AFS, LD and IBS summary statistics, as described in the manuscript.

# Although this is an input file for PopSizeABC simulations, some of the lines of code did not need modifications and are exactly the same as those from PopSizeABC example scripts distributed with the program. Unavoidably, I (Jordan Bemmels) needed to retain many lines of codes which I did not develop myself, in order that the example input file can be run successfully. The PopSizeABC code is all freely available at https://forge-dga.jouy.inra.fr/projects/popsizeabc. Full credit for unmodified lines of codes and for the general setup of this input file belongs to the authors of PopSizeABC (but any errors or inconsistencies are my own):

# Boitard S, Rodríguez W, Jay F, Mona S, Austerlitz F. 2016 Inferring population size history from large samples of genome-wide molecular data - an Approximate Bayesian Computation approach. PLoS Genet. 12, 1–36. (doi:10.1371/journal.pgen.1005877)

################################################
############### dependencies ###################
################################################

#!/usr/bin/python

import sys
import time
import datetime
import numpy as np
import popgen_abc
import tarfile
import os

start_time=time.time()

# JBB edits to add parameters on the command line
L=int(sys.argv[1])*1000000 # length of chromosome segments, in Mb
nb_seg=int(sys.argv[2]) # number of independent chromosome segments to simulate
nb_rep=int(sys.argv[3]) # number of replicates
batch=sys.argv[4] # batch so we can run in parallel with different outfile names

################################################
############### paramaters #####################
################################################

# general parameters
outfile_name='batch'+batch+'_L'+str(L/1000000) # root of the output files # JBB add batch number for unique file name
#nb_rep=100 # number of simulated datasets # JBB original way of defining number of replicates not from command line
#nb_seg=100 # number of independent segments in each dataset
#L=4000000 # size of each segment, in bp.
n=10 # haploid sample size
mac=1 # minor allele count threshold for AFS and IBS statistics computation
mac_ld=1 # minor allele count threshold for LD statistics computation
save_msp=False # if this parameter is set to True, snp positions and haplotypes corresponding to the same dataset will be stored in a compressed tar file.
	      # this allows to keep and potentially re-use the exact genomic samples, rather than just the summary statistics,
	      # but this require high memory ressources (approx 1 Mo per simulated dataset, on average, with current parameter values). 

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

# prior distributions
#r_min=0.7271*10**(-8) # lower bound for the per generation per bp recombination rate
#r_max=7.271*10**(-8) # upper bound for the per generation per bp recombination rate
r_fixed = 2.1*(10**(-8)) # fix recombination rate using the rate for rhea
mmu=1.345868*10**(-8) # per generation per bp mutation rate
Nmin=1.3010 # lower bound for the population size in each time window (in log10 scale)
Nmax=6.3010 # upper bound for the population size in each time window (in log10 scale)
Nmin_anc=1.3010 # use different bounds for ancestral population min and max
Nmax_anc=4.6990 # use has different bounds for ancestral population min and max

# LD statistics parameters 
per_err=5 # the length of each interval, as a percentage of the target distance
r=r_fixed # the per generation per bp recombination rate (an approximation is sufficient)
# creation of the bins of physical distance for which the average LD will be computed, based on the time windows defined above
interval_list=[]
for i in range(nb_times-1):
    t=(times[i+1]+times[i])/2
    d=1/(2*r*t)
    if d <= L:
        interval_list.append([d-per_err*d/100,d+per_err*d/100])
t=Tmax+times[nb_times-1]-times[nb_times-2]
d=1/(2*r*t) # JBB edits this line to use the same formula for d as a few lines immediately above (original code used an approximation of r instead of the actual value given; this is for calculating the distance bin for the oldest time window / shortest distance which should be calculated in the same way as the other windows)
interval_list.append([d-per_err*d/100,d+per_err*d/100])
print "Average LD will be computed for the following distance bins (in bp) :"
print interval_list
print ""

# IBS statistics parameters
prob_list=[0.0001,0.001,0.01,0.1,0.25,0.5,0.75,0.9,0.99,0.999,0.9999] # quantiles used to summarize the distribution of IBS segment length
size_list=[1] # number of diploid individuals that are used to define IBS segments (several values can be concatenated)

################################################
############### the programm ###################
################################################

# create the matrices where results (parameters and statistics) are stored
nb_dist=len(interval_list)
nb_prob=len(prob_list)
nb_size=len(size_list)
params=-np.ones(shape=[nb_rep,nb_times+2],dtype='float')
results=-np.ones(shape=[nb_rep,nb_dist+nb_size*nb_prob+n/2+1],dtype='float')  # missing statistic values will be set to -1.
									     # this mostly concerns LD statistics, because it is sometimes impossible to find SNP pairs in a given distance bin.
print('Total number of statistics : '+str(nb_dist+nb_size*nb_prob+n/2+1)+'\n')

# simulate parameters from the prior
mu=mmu*np.ones([nb_rep])
params[:,0]=mu
#r=10**(np.random.uniform(low=np.log10(r_min),high=np.log10(r_max),size=nb_rep)) # JBB edits: originally it sampled on a regular scale [r_min, r_max], but I want it to sample on a log-10 scale, so I added the log transformation [log10()] and back-transformation [10**] here 
r=r_fixed # fix recombination rate
params[:,1]=r
# EDIT FOR KIWI PROJECT use the original PopSizeABC code, except that we have a different Nmax_anc than the default N_max for other time windows - we will allow this to be anything and to not depend on the multiplication factor;
for i in range(nb_times):
    if i==0:
        pop_size=10**(np.random.uniform(low=Nmin,high=Nmax,size=nb_rep))
    elif i<(nb_times-1): # JBB edits to consider all but the N_anc bin;
        multiplication_factor=10**(np.random.uniform(-2,2,size=nb_rep)) # JBB edits: default was to sample in [-1,1] which restricts to a 10-fold increase/decrease between time-windows, but changed to [-2,2] to allow a 100-fold increase/decrease
        pop_size=params[:,1+i]*multiplication_factor
        pop_size=np.maximum(10**Nmin,pop_size)
        pop_size=np.minimum(pop_size,10**Nmax)
    else:
        # JBB adds this else-statement to have a different formula for the N_anc bin;
        pop_size=10**(np.random.uniform(low=Nmin_anc,high=Nmax_anc,size=nb_rep))
    params[:,2+i]=pop_size

# simulate summary statistics
outfile_name_2 = outfile_name + "_n" + str(n) + "_s" + str(nb_seg)
print 'Started the simualtions'
for i in range(nb_rep):
    elapsed_time=time.time()-start_time
    print 'Simulating replicate', i+1,", current time :", time.ctime(),", elapsed time :", datetime.timedelta(seconds=elapsed_time)
    #print "Population sizes :"
    #print params[i,:]
    sys.stdout.flush()
    try:
    	results[i,:]=popgen_abc.simul_stats_one_rep_macld(outfile_name_2,i,nb_seg,L,n,times,params[i,:],interval_list,size_list,prob_list,L,mac=mac,mac_ld=mac_ld,save_msp=save_msp)
    except:
        print 'Problem with replicate', i+1
        pass
    # JBB ADDS - print the results as each new simulation replicate is completed, rather than only at the end after all are completed, to help protect against server crashes and lost results;
    np.savetxt(outfile_name_2+"_mac"+str(mac)+"_macld"+str(mac_ld)+'.stat',results[0:nb_rep,:],fmt='%.3e')
    np.savetxt(outfile_name_2+'.params',params[0:nb_rep,:],fmt='%.3e')
    print "Printed the results"

# print the results
np.savetxt(outfile_name_2+"_mac"+str(mac)+"_macld"+str(mac_ld)+'.stat',results[0:nb_rep,:],fmt='%.3e')
np.savetxt(outfile_name_2+'.params',params[0:nb_rep,:],fmt='%.3e')
print "Printed the results"

# tar ms files (if they exist)
if save_msp==True:   
    mytar=tarfile.open(outfile_name_2+'.ms.tar.bz2','w:bz2')
    for i in range(nb_rep):
    	for j in range(nb_seg):
	    mytar.add(outfile_name_2+'_rep'+str(i+1)+'_seg'+str(j+1)+'.msp')
	    os.remove(outfile_name_2+'_rep'+str(i+1)+'_seg'+str(j+1)+'.msp')
    mytar.close()
    print "Created tar for msp files"

# the end
elapsed_time=time.time()-start_time
print "Finished job, ","current time :", time.ctime(), "elapsed time :", datetime.timedelta(seconds=elapsed_time)
