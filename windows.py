from Bio import SeqIO
from Bio.SeqUtils import GC
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.switch_backend('Agg')
import argparse
from pylab import savefig

def main():
##Use the get_args function
	GENOME, WINDOWSIZE, STEP = get_args()
	print('File is = ' + GENOME + '.')
	print('Window size is ' + str(WINDOWSIZE) + '.')
	print('Step size is ' + str(STEP) + '.')

#initialize a list of GC content
	GCLIST = []
	SCAFFCOUNT = 0
	WINDOWCOUNT = 0
#populate the list
#for every record in your genome
	for SEQRECORD in SeqIO.parse(GENOME, 'fasta'):
#extract the sequence
		SEQUENCE = SEQRECORD.seq
#check to see if it's over 30000 br
		if len(SEQUENCE) > 30000:
			SCAFFCOUNT = SCAFFCOUNT + 1
#if so, use the windows function to get each window
			for SUBSEQ in windows(SEQUENCE, WINDOWSIZE, STEP):
				WINDOWCOUNT = WINDOWCOUNT + 1
#calculate the GC content for that window
				WINDOWGC = GC(SUBSEQ)
#add that number to the growing list
				GCLIST.append(WINDOWGC)

#build histogram
	plt.hist(GCLIST)
	plt.savefig("GChist.png")

#output basic stats
	print('Scaffolds over 30kb = ' + str(SCAFFCOUNT) + '.')
	print('Window count = ' + str(WINDOWCOUNT) + '.')

#create function that will extract each window from the input file
def windows(SEQUENCE, WINSIZE, STEPSIZE):
	SEQLEN = len(SEQUENCE)
	NUMBEROFWINDOWS = int(SEQLEN/STEPSIZE)
	for WINDOWNUMBER in range(0, NUMBEROFWINDOWS * STEPSIZE, STEPSIZE):
		yield SEQUENCE[WINDOWNUMBER:WINDOWNUMBER+WINSIZE]
	
def get_args():
	parser = argparse.ArgumentParser(description="Will generate a histogram of GC content in a genome.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-g', '--genome', type=str, help='genome file in .fasta format.', required=True)
	parser.add_argument('-w', '--windowsize', type=int, help='Window size in bp.', required=True)
	parser.add_argument('-s', '--step', type=int, help='Step size for moving windows in bp.', required=True)
	args = parser.parse_args()
	GENOME = args.genome
	WINDOWSIZE = args.windowsize
	STEP = args.step
	
	return GENOME, WINDOWSIZE, STEP

if __name__ =="__main__":main()
        

