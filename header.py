import argparse
from Bio import SeqIO
from collections import deque
import collections
import os
import gzip
import sys

def main():

##Use the get_args function
	INPUT, PREFIX, OUTPUT = get_args()

#Count the number of lines for later use.
	NUMBEROFLINES = sum(1 for LINE in gzip.open(INPUT))
#	NUMBEROFENTRIES = NUMBEROFLINES/4
	
#Set starting points for indexing the files
	ODD = 0
	EVEN = 4
	
#Initialize the header counters for the headers. First read will be 000000000000000. Second read will be 000000000000001. Etc.
	HEADER1 = 0
	HEADER2 = 0

#Create queues of index numbers to select lines to read and modify. One queue for each paired read.
	HEADER1INDEX = list(range(ODD, NUMBEROFLINES, 8))
	HEADER1QUEUE = deque(HEADER1INDEX)
	HEADER2INDEX = list(range(EVEN, NUMBEROFLINES, 8))
	HEADER2QUEUE = deque(HEADER2INDEX)

#Open the output file for writing.
	OUTFILE = open(OUTPUT, 'w+')

#Open the input file.
	with gzip.open(INPUT, 'rt') as FILE:
#Label each line in the input file (enumerate) and iterate through each line (for ...)
		for NUMBEREDLINE, LINE in enumerate(FILE):
#If the queues are still populated...
			if HEADER1QUEUE and HEADER2QUEUE:
#And if the line number we're on matches the number at the first spot in the _first_ header queue...
				if NUMBEREDLINE == HEADER1QUEUE[0]:
#And if the line starts with '@'
					if LINE[0] == '@':
#Change the header to this using information from the header counters and the prefix input by the user.
						OUTFILE.write('@' + '{:0>15}'.format(HEADER1) + '\t' + PREFIX + ':i:1\n')
#Iterate forward the header counters
						HEADER1 += 1
#Pop off the first number in the queue.
						HEADER1QUEUE.popleft()
#If the line doesn't start with a '@', something is wrong.
					else:
						print('Something is wrong with your file.')
						sys.exit(1)
#If the line number matches the number in the first spot in the _second_ header queue...
				elif NUMBEREDLINE == HEADER2QUEUE[0]:
#Do all the same stuff but with different numbers
					if LINE[0] == '@':
						OUTFILE.write('@' + '{:0>15}'.format(HEADER2) + '\t' + PREFIX + ':i:2\n')
						HEADER2 += 1
						HEADER2QUEUE.popleft()
#If the line is not a header, just write it to the output file.
				else: 
					OUTFILE.write(LINE)
#Finish writing all remaining lines even though HEADER1QUEUE is empty
			else:
				OUTFILE.write(LINE)
#Close your output file when you're done going through both queues
	OUTFILE.close()
	

##Get arguments function
def get_args():
	#What this script does
	parser = argparse.ArgumentParser(description="Replace headers for python class", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	#Argument of the input data file
	parser.add_argument('-i', '--input', help='Name of your input gzipped fastq file, include path if outside current directory', required=True)
	#Argument for prefix after the header ID
	parser.add_argument('-p', '--prefix', type=str, help='Indicator text before read ID (1 or 2)', default='OP')
	#Argument for the name of the output file
	parser.add_argument('-o', '--output', type=str, help='Name of the output file.', required=True)

	args = parser.parse_args()
	INPUT = args.input
	if args.prefix:
		PREFIX = args.prefix
	else:
		BASENAME = os.path.basename(INPUT).split(".")[0]
		PREFIX = BASENAME
	OUTPUT = args.output
	
	return INPUT, PREFIX, OUTPUT

if __name__ =="__main__":main()	