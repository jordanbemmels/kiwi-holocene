# before beginning, identify ultra-conserved elements (UCEs) in the kiwi reference genome, using phyluce
# the reference genome should first be converted to 2bit format, using scripts downloaded from the Kent Source Archive
	# http://hgdownload.soe.ucsc.edu/admin/exe/

rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/faToTwoBit ./ # downloads the faToTwoBit script for performing conversion
rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/twoBitInfo ./ # downloads the twoBitInfo script for viewing results
./faToTwoBit kiwi_ref_genome.fna aRowi_ref.2bit # perform conversion
./twoBitInfo aRowi_ref.2bit aRowi_ref_sizes.txt # view results

# download the vertebrate 5K probe set for UCE loci

wget https://raw.githubusercontent.com/faircloth-lab/uce-probe-sets/master/uce-5k-probe-set/uce-5k-probes.fasta

# align UCE probe set to the kiwi reference genome using phyluce
# the matches can be viewed in the outfile: ./kiwi_UCE-genome-lastz/uce-5k-probes.fasta_v_aRowi_ref.lastz.clean

time phyluce_probe_run_multiple_lastzs_sqlite --db kiwi_UCE.sqlite --output kiwi_UCE-genome-lastz --scaffoldlist aRowi_ref --genome-base-path ./ --probefile uce-5k-probes.fasta --cores 1

# to understand what the output means and the format, inspect with sqlite3

sqlite3 kiwi_UCE.sqlite
.mode columns
.headers on
.nullvalue .
.tables # prints which tables are in the database
SELECT * FROM aRowi_ref LIMIT 10; # prints the first 10 lines plus column names of aRowi_ref table
	# result of this command, column names only: id  score  name1  strand1  zstart1  end1  length1  name2  strand2  zstart2  end2  length2  diff  cigar  identity  percent_identity  continuity  percent_continuity  coverage  percent_coverage
	# the .lastz.clean file is the same, except is missing the "id" column

