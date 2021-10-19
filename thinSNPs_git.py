"""
This code does one thing:
1) thin SNPs to ensure a minimum distance (kbp) between sites, to avoid linkage disequilibrium.

Input is an ANGSD sites file.

Usage: python selectFscSNPs.py inSitesFile outFile requiredKbpDistance

"""

import sys

inputSites = sys.argv[1]
outputFile = sys.argv[2]
kbp = int(sys.argv[3])

### set up an outfile

outfile = open(outputFile, 'w')

### read the infile line by line

linecount = 0
with open(inputSites) as f:
	for line in f:
		linecount = linecount + 1
		if (linecount % 10000 == 0):
			print("working on line " + str(linecount))
		if (linecount > -1):
			# split into chrom and pos
			currentChrom, currentPos = line.split()[0], int(line.split()[1])
			#print (currentChrom)
			#print(currentPos)
			# 
			if (linecount == 1):
				# if it's the very first line, accept the SNP
				outfile.write(line)
				#print("condition1")
				previousChrom = currentChrom
				previousPos = currentPos
			elif (currentChrom != previousChrom):
				# otherwise, check if it's a new chromosome, and if it is, accept the SNP
				outfile.write(line)
				#print("condition2")
				previousChrom = currentChrom
				previousPos = currentPos
			elif (currentPos >= previousPos + kbp*1000):
				# otherwise, only write the SNP if it is at least the desired kbp away
				outfile.write(line)
				#print(previousPos + kbp*1000)
				#print(currentPos)
				#print("condition3")
				previousChrom = currentChrom
				previousPos = currentPos
				
		"""
		currentChrom = 
		if (linecount % 10000 == 0):
			print("working on line " + str(linecount))
		"""


### close the outfile

outfile.close()
