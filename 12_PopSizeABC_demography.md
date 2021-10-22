# Demographic history as estimated in PopSizeABC

We will use [PopSizeABC](https://forge-dga.jouy.inra.fr/projects/popsizeabc) to infer changes in effective population size (Ne) in each population, focusing on the Holocene.

## Modifications to PopSizeABC code

PopSizeABC is written in python code. In the kiwi analysis, I made several modifications to the code. However, PopSizeABC is not my intellectual property and I do not have permissions to publically share the core PopSizeABC code here. Instead, for reproducibility of the kiwi project I will describe the changes that any user could also make to the PopSizeABC code if desired.

#### Modifications to gwas/summary_stat.py

Add hashes (```#```) at the beginning of the lines to comment out the entire if-statement beginning with ```if len(r2_list) < 2:``` (lines 160-182, i.e., where the original code indicates to ```# try a more exhaustive screening of SNP pairs```). This statement is removed because for a given distance bin, if there are fewer than 2 SNP pairs available, it would calculate the linkage disequilibrium (LD) across *every* possible SNP pair that meets the distance criteria, rather than ensuring that the SNP pairs are physically non-overlapping (the normal behaviour). Thus, the SNP pairs may be non-overlapping and non-independent. I prefer not to do this. If these lines are removed, then the LD mean and standard deviation for that distance bin will not be updated and will remain as -1. This is fine as PopSizeABC will recognize that something was wrong with these simulations (very rare), and they will fail the good_sims test and be discarded.

#### Modifications to stat_from_vcf.py

The ```parameters``` section needs to be updated to be relevant to kiwi. This will be described further below in a subsequent section, when we are describing how to develop and set up the models for kiwi.

However, there is one generic change to the calculations performed in the parameters section that will be described here: the calculation of *d* for the final distance bin (oldest time window / shortest physical linkage distance) is changed to ```d=1/(2*r*t)```, to correspond with the same formula used for all the other distance bins other than the final one. This is instead of the default code ```d=10**8/(2*t)``` which fixes the recombination rate at 1x10e-8 rather than using the actual pre-defined recombination rate *r* for the final distance bin. This change is also implemented separately in the ```simul_data.py``` script for consistency with the simulations. To make this change:

```
### OLD LINE 58
d=10**8/(2*t)

### NEW LINE
d=1/(2*r*t)
```

In addition, under the ```the programm``` section (sic), there are two changes.

First is to process *all* scaffolds (not only the first two scaffolds):

```
### OLD LINES 90 [...] 93
for chro in range(1,2):
    [...]
    infile_vcf='../cattle_data/Chr'+str(chro)+'.vcf.gz'

### NEW LINES
all_chromo=os.listdir('kiwi_data/splitVCF/'+pop) # this allows processing all VCF files in the specified directory with this path, so that all scaffolds can be included
for chro in range(len(all_chromo)):
    [...]
    infile_vcf='kiwi_data/splitVCF/'+pop+'/'+all_chromo[chro]
```

Second is to ensure that we read the kiwi pedigree, not the default cattle pedigree.

```
### OLD LINES
    IO.parsePedFile_nogeno('../cattle_data/indiv.ped',mydata)  

### NEW LINES
    IO.parsePedFile_nogeno('kiwi_data/kiwi_indiv.ped',mydata) # JBB there is only a single pedigree file for all kiwis across all populations
```















The parameters section looks as follows:

```
#pop='populationName' # JBB edits - original way from example code of setting pop
pop=sys.argv[1] # JBB edits
list_ani=IO.read_list('kiwi_data/list_indiv_'+pop+'.txt') # list of diploid animals used for computing the summary statistics
n=len(list_ani)*2 # haploid sample size
mac=1 # minor allele count threshold for AFS and IBS statistics computation
mac_ld=1 # minor allele count threshold for LD statistics computation
L=4000000 # size of each segment, in bp.

# time windows
nb_times=15 # number of time window
Tmax=622 # the oldest time window will start at Tmax
a=0.09 # the length of time windows increases when time increases, at a speed that is determined by this coefficient  
# computation of time windows based on the above parameters
times=-np.ones(shape=nb_times,dtype='float')
for i in range(nb_times):
    times[i]=(np.exp(np.log(1+a*Tmax)*i/(nb_times-1))-1)/a
print "Population size changes at the following times (in generations):"
print times
print ""

# LD statistics parameters 
per_err=5 # the length of each interval, as a percentage of the target distance
r=2.1*(10**(-8)) # the per generation per bp recombination rate (an approximation is sufficient value)
# creation of the bins of physical distance for which the average LD will be computed, based on the time windows defined above.
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
size_list=[1] # number of diploid individuals that are used for to define IBS segments (several values can be concatenated)
```
