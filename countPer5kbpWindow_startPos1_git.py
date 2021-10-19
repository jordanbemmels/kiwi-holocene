"""
Code is to take a .pos sites from ANGSD and count the number of sites per 5kbp non-overlapping window per chromosome.
The purpose of this is to reduce the .pos file (enormous) by a factor of ~5,000, so that it becomes more wieldy.

Logic is to read the enormous file line by line, keep track of lines in windows, and print the count to the outfile when ready to move on to the new window. 

Usage: python countPer5kbpWindow_startPos1.py input_pos_file output_summary_file
"""

import sys
import math

posfile = sys.argv[1]
outfile = sys.argv[2]

#####

currentChr = "none"
currentStartPos = 1
currentCount = 0

# start by opening the outfile to which to write
with open(outfile, 'w') as o:
	o.write("chr\tstart\tend\tSites\n")
	# now open the file to read and process
	with open(posfile) as f:
		next(f) # this skips the header line - assumes the file has text header
		for line in f:
			chr,pos = line.split("\t")[0:2]
			pos = int(pos) # convert string to integer
			if (chr == currentChr):
				if (pos < currentStartPos + 5000):
					# if we are on the same chromosome and haven't stepped outside the window, add to the count
					currentCount = currentCount + 1
				else:
					# if we are on the same chromosome and HAVE stepped outside the window, print the line and start a new count
					o.write('\t'.join([str(currentChr), str(int(currentStartPos)), str(int(currentStartPos + 5000)), str(currentCount)]) + "\n")
					currentChr = chr
					currentStartPos = (math.floor(pos / 5000) * 5000) + 1
					currentCount = 1
			else:
				# if we are NOT on the same chromosome, print the line and start a new count
				o.write('\t'.join([str(currentChr), str(int(currentStartPos)), str(int(currentStartPos + 5000)), str(currentCount)]) + "\n")
				currentChr = chr
				currentStartPos = (math.floor(pos / 5000) * 5000) + 1
				currentCount = 1
		# need to print the very last line to record the final count for the last window of the last chromosome
		o.write('\t'.join([str(currentChr), str(int(currentStartPos)), str(int(currentStartPos + 5000)), str(currentCount)]) + "\n")
