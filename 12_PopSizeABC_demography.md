# Demographic history as estimated in PopSizeABC

We will use [PopSizeABC](https://forge-dga.jouy.inra.fr/projects/popsizeabc) to infer changes in effective population size (Ne) in each population, focusing on the Holocene.

## Modifications to PopSizeABC code

PopSizeABC is written in python code. In the kiwi analysis, I made several modifications to the code. However, I do not have permissions to publicly share the core PopSizeABC code here (only example input files are shared below). Instead, for reproducibility of the kiwi project I will describe the changes that any user could also make to the default PopSizeABC code if desired.

#### Modifications to gwas/summary_stat.py

Add hashes (```#```) at the beginning of the lines to comment out the entire if-statement beginning with ```if len(r2_list) < 2:``` (lines 160-182, i.e., where the original code indicates to ```# try a more exhaustive screening of SNP pairs```). This statement is removed because for a given distance bin, if there are fewer than 2 SNP pairs available, it would calculate the linkage disequilibrium (LD) across *every* possible SNP pair that meets the distance criteria, rather than ensuring that the SNP pairs are physically non-overlapping (the normal behaviour). Thus, the SNP pairs may be non-overlapping and non-independent. I prefer not to do this. If these lines are removed, then the LD mean and standard deviation for that distance bin will not be updated and will remain as -1. This is fine as PopSizeABC will recognize that something was wrong with these simulations (very rare), and they will fail the *good_sims* test during ABC and be discarded.

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
### OLD LINE 95
    IO.parsePedFile_nogeno('../cattle_data/indiv.ped',mydata)  

### NEW LINES
    IO.parsePedFile_nogeno('kiwi_data/kiwi_indiv.ped',mydata) # JBB there is only a single pedigree file for all kiwis across all populations
```

## Call SNPs for each population and prepare population-specific input files

An example is shown here for the popultion aHaast (also sometimes written australis__Haast).

SNPs will be called for each population using [ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD). The creation of SNP whitelists for each population has been [previously](https://github.com/jordanbemmels/kiwi-holocene/blob/main/03_Create_SNP_whitelists.md) described under Whitelist 03. For example, for aHaast the filtered popultion-specific SNP file is ```sites_australis__Haast_bfv2_3x_maf00_noZ.txt1```. We will print a VCF file (```-dobcf 1```; output has .bcf file extension). Require data for all five individuals (```-minInd 5```), and a minimum depth of 8x to call a SNP for each individual (```-setMinDepthInd 8```).

```
FILTERS="-minQ 20 -minMapQ 20 -minInd 5 -setMinDepthInd 8"
TODO="-doCounts 1 -doMajorMinor 3 -doMaf 1 -doPost 1 -doGeno 1 -dobcf 1 -fold 1 --ignore-RG 0"
angsd -b BAM__australis__Haast.txt -sites sites_australis__Haast_bfv2_3x_maf00_noZ.txt -GL 1 -P 1 $FILTERS $TODO -out aHaast
```

The output file is ```aHaast.bcf```.

However, PopSizeABC takes a separate .bcf file for each scaffold, so we need to split the overall VCF into a separate VCF files for each scaffold. First, bgzip and index the VCF (commands are from [SAMtools](http://www.htslib.org/)).

```
mv aHaast.bcf aHaast.vcf # renme as VCF
bgzip -c aHaast.vcf > aHaast.vcf.gz
tabix -p vcf aHaast.vcf.gz
```

To perform the split by scaffold, use the provided script ```splitVCFbyScaffold_git.R```. The script requires the installation of [bcftools](https://samtools.github.io/bcftools/). The two command-line arguments are the input file (a corresponding .vcf, .vcf.gz, and .tbi file are all required) and output directory names, respectively. For the aHaast population:

```
Rscript splitVCFbyScaffold.R aHaast.vcf aHaast
```

Separate VCF files will be written to the ```./aHaast``` directory, e.g., ```./aHaast/aHaast_PTFB01000001.1.vcf``` for the first scaffold.

The VCF files are now ready. However, for each population, in addition to calling SNPs, we also need to create a file containing a list of individuals, a pedigree file (.ped). First create the list of individuals. For aHaast, the file ```list_indiv_aHaast.txt``` has the following five lines, corresponding to the names of the samples as appearing in the VCF file:

```
australis__Haast__KW36__RA0997_sorted
australis__Haast__KW37__TFCK_sorted
australis__Haast__KW38__R34159_sorted
australis__Haast__KW39__RA0921_sorted
australis__Haast__KW40__RA0201_sorted
```

A single pedigree file can be re-used for all populations, if it contains data for all 52 individuals. The population is given in the first column, second column is the individual, columns 3 / 4 / 5 are mother / father / sex (all 0 if unknown), column 6 is phenotype (-999) for missing. An example pedigree file ```kiwi_indiv.ped``` (only the first five lines for the aHaast individuals are shown):

```
aHaast australis__Haast__KW36__RA0997_sorted 0 0 0 -999
aHaast australis__Haast__KW37__TFCK_sorted 0 0 0 -999
aHaast australis__Haast__KW38__R34159_sorted 0 0 0 -999
aHaast australis__Haast__KW39__RA0921_sorted 0 0 0 -999
aHaast australis__Haast__KW40__RA0201_sorted 0 0 0 -999
```

## Set up stat_from_vcf.py file to calculate summary statistics

Calcultion of summary statistics is done using the ```stat_from_vcf.py``` file. Note that this script contains an extensive parameters section followed by the actual code (and thus I cannot provide a full example as the code is not mine to distribute). However, a few changes to the code section have already been described above. Here, I will provide text indicating how the ```parameters``` section was specified for the kiwi project, organized by the different subsections. The same parameters can be used for any kiwi population.

#### General parameters

```
#pop='populationName' # originally the population name is defined within the file
pop=sys.argv[1] # instead, make the script generic to take the population as a command argument
list_ani=IO.read_list('kiwi_data/list_indiv_'+pop+'.txt') # list of diploid animals used for computing the summary statistics (see for example the creation of the file list_indiv_aHaast.txt described above)
n=len(list_ani)*2 # haploid sample size
mac=1 # minor allele count threshold for AFS and IBS statistics computation
mac_ld=1 # minor allele count threshold for LD statistics computation
L=4000000 # size of each segment, in bp.
```

Note that we have re-defined the script to take the population name as a command-line argument. By setting ```mac``` and ```mac_ld``` to ```1```, we do not implement any minimum minor allele frequency. The size of each segment, in basepairs (```L=4000000```) is relevant in the simultions (described below), and is used here to set the largest physical distance will be for which LD statistics are calculated. It is set to 4,000,000 bp here because we have very few scaffolds larger than this size, so we would not be able to obtain an adequate number of non-physically-overlapping SNP pairs to have a robust estimate of LD at larger distances in the genome in the empirical data for kiwi.

#### Time windows

```
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
```

This section controls the number of time windows for which population sizes will be estimated, and their ages and lengths. We choose 15 windows, with the oldest time window starting at 622 generations (~11.9 ka for kiwi, assuming a generation time of 19.1357 years; see manuscript), with a rate parameter of a=0.09 for adjusting the relative size of windows. See [Boitard et al. 2016](https://doi.org/10.1371/journal.pgen.1005877) for more details. These values were set in order to give a reasonable number of windows with breakpoints roughly corresponding to important historical events of relevance to kiwi.

Boitard S, Rodríguez W, Jay F, Mona S, Austerlitz F. 2016 Inferring population size history from large samples of genome-wide molecular data - an Approximate Bayesian Computation approach. PLoS Genet. 12, 1–36.

#### LD statistics parameters

```
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
d=1/(2*r*t) # edited this line to use the same formula for d as a few lines immediately above (original code used an approximation of r instead of the actual value given; this is for calculating the distance bin for the oldest time window / shortest distance which should be calculated in the same way as the other windows)
interval_list.append([d-per_err*d/100,d+per_err*d/100])
print "Average LD will be computed for the following distance bins (in bp) :"
print interval_list
print ""
```

Here, we allow up to 5% error relative to the target distance in order to find SNPs for each physical distance bin. We provide the kiwi recombination rate of 2.1 x 10e-8 (see manuscript text). Note also the edit to the code on line ```d=1/(2*r*t)```, previously described above.

#### IBS statistics parameters

No changes needed and in fact we do not use IBS summary statistics (see Boitard et al. 2016) so this section is irrelevant.

## Calculate summary statistics on empirical data

Once the ```stat_from_vcf.py``` file is prepared, we are ready to run it on each of the populations. Here is an example for aHaast.

```
python2 comp_stat1/stat_from_vcf.py aHaast > terminalOutput_aHaast.txt
```

The outfile with summary statistics will have the extension ```.stat```. If you want to change the outfile to something specific, you can modify the ```# print the result``` section of ```stat_from_vcf.py```.

## Generate simulated data







