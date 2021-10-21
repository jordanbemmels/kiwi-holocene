# UCE phylogeny

A phylogeny of ultra-conserved elements (UCEs) will be generated in [Astral-III](https://github.com/smirarab/ASTRAL).

## Identify UCEs

We will identify UCEs in the kiwi reference genome using [phyluce](https://phyluce.readthedocs.io/en/latest/).

The reference genome should first be converted to 2bit format, using scripts downloaded from the [Kent Source Archive](http://hgdownload.soe.ucsc.edu/admin/exe/). To download the scripts:

```
rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/faToTwoBit ./ # downloads the faToTwoBit script for performing conversion
rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/twoBitInfo ./ # downloads the twoBitInfo script for viewing results
```

Now, perform the conversion and view the results. Note that we are re-naming the genome ```kiwi_ref_genome.fna``` to ```aRowi_ref.2bit``` to keep the original .fna and .2bit files clearly separate.

```
./faToTwoBit kiwi_ref_genome.fna aRowi_ref.2bit # perform conversion
./twoBitInfo aRowi_ref.2bit aRowi_ref_sizes.txt # view results
```

Download the vertebrate 5K probe set for UCE loci.

```
wget https://raw.githubusercontent.com/faircloth-lab/uce-probe-sets/master/uce-5k-probe-set/uce-5k-probes.fasta
```

Align the UCE probe set to the kiwi reference genome using phyluce.

```
phyluce_probe_run_multiple_lastzs_sqlite --db kiwi_UCE.sqlite --output kiwi_UCE-genome-lastz --scaffoldlist aRowi_ref --genome-base-path ./ --probefile uce-5k-probes.fasta --cores 1
```

The matches of the above command can be viewed in the outfile ```./kiwi_UCE-genome-lastz/uce-5k-probes.fasta_v_aRowi_ref.lastz.clean```. To understand what the output means and the format, inspect with [sqlite3](https://www.sqlite.org/index.html):

```
sqlite3 kiwi_UCE.sqlite
.mode columns
.headers on
.nullvalue .
.tables # prints which tables are in the database
SELECT * FROM aRowi_ref LIMIT 10; # prints the first 10 lines plus column names of aRowi_ref table
	# result of this command, column names only: id  score  name1  strand1  zstart1  end1  length1  name2  strand2  zstart2  end2  length2  diff  cigar  identity  percent_identity  continuity  percent_continuity  coverage  percent_coverage
	# the .lastz.clean file is the same, except is missing the "id" column
```

We need to convert the matches from the  output file ```uce-5k-probes.fasta_v_aRowi_ref.lastz.clean``` into a format ANGSD can use, such as an ANGSD regions file. When performing this conversion, note that phyluce defines site positions using a zero-origin half-open system, whereas ANGSD uses a one-indexed closed system, so we need to account for that in the conversion. We will also retain 1,000-kbp flanking regions on either side of the UCE locus. If any separate UCEs ± 1000-kbp flanking regions end up overlapping, then should be merged into a single larger locus. Perform the conversion taking these considerations into account using the script [convertLastzToAngsdRegions_git.R]().



