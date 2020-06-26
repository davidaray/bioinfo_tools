import os
import argparse
from Bio import SeqIO
import re
import csv

####################################
## Usage: newname_fasta.py -p <filename> -f <fasta file> -o <fasta file>
##
## NOTE: make sure there are no trailing end of line characters at the bottom of your pull_list
## You will just end up copying the original fasta to a new file if there are.

def get_args():
	#What this script does
	parser = argparse.ArgumentParser(description="Pull set of sequences from a larger fasta file using a list of headers.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	#Argument of the input .clstr file
	parser.add_argument('-p', '--pull_file', type=str, help='Name of your file with the lists of old and new names.', required=True)
	#Argument for the taxon list.
	parser.add_argument('-f', '--fasta', type=str, help='Fasta file that contains entries to be renamed.', required=True)
	parser.add_argument('-o', '--outfasta', type=str, help='Name of the output file with renamed entries.', required=True) 
	parser.add_argument('-b', '--blast_friendly', type=str, help='Make blast-friendly by replacing # and /. y or n.', default = 'n')

	args = parser.parse_args()
	PULLLIST = args.pull_file
	FASTA = args.fasta
	OUTFILE = args.outfasta
	BLASTFRIENDLY = args.blast_friendly
	
	return PULLLIST, FASTA, OUTFILE, BLASTFRIENDLY

PULLLIST, FASTA, OUTFILE, BLASTFRIENDLY = get_args()

print('Pulling entries in ' + FASTA + ' to create ' + OUTFILE +'.')

#Create list to store sequences
RECORDLIST = []
#Open the fasta input file.
for record in SeqIO.parse(FASTA, 'fasta'):
	#Sanity check
	print('Checking ' + record.id)
	#Import the list of headers to search for.
	with open(PULLLIST, 'r') as LIST: 
		#For every line in the list.
		for LINE in LIST:
			LINE = LINE.rstrip() 
			if LINE in record.id:
				print('Match between ' + LINE + ' and ' + record.id)
				RECORDLIST.append(record)
#Write the list to a new file.
with open('tmp.fa', 'w') as TMPFILE:
	SeqIO.write(RECORDLIST, TMPFILE, 'fasta')

if BLASTFRIENDLY == 'n':
	os.rename('tmp.fa', OUTFILE)
else:
	print('Making headers blast-friendly.')
	RECORDLIST = []
	for record in SeqIO.parse('tmp.fa', 'fasta'):
		HEADER = record.id
		NEWHEADER = HEADER.replace('#', '__')
		NEWHEADER = NEWHEADER.replace('/', '___')
		record.id = NEWHEADER
		record.description = ''
		RECORDLIST.append(record)
	with open(OUTFILE, 'w') as OUT:
		SeqIO.write(RECORDLIST, OUT, 'fasta')

os.remove('tmp.fa')
	


			

