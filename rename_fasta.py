import os
import argparse
from Bio import SeqIO
import re
import csv

####################################
## Usage: newname_fasta.py -r <filename> -f <fasta file> -n <fasta file>
##
## Replace headers in a large fasta file using a text file that contains original names
## in the first column and new names in the second. Will find only the headers of interest 
## and leave any other headers as-is.

def get_args():
	#What this script does
	parser = argparse.ArgumentParser(description="Replace headers in a large fasta file using a text file that contains original names in the first column and new names in the second.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	#Argument of the input .clstr file
	parser.add_argument('-r', '--replace_file', type=str, help='Name of your file with the lists of old and new names.', required=True)
	#Argument for the taxon list.
	parser.add_argument('-f', '--fasta', type=str, help='Fasta file that contains entries to be renamed.', required=True)
	parser.add_argument('-o', '--outfasta', type=str, help='Name of the output file with renamed entries.', required=True) 

	args = parser.parse_args()
	REPLACEMENTS = args.replace_file
	FASTA = args.fasta
	RENAMED = args.outfasta
	
	return REPLACEMENTS, FASTA, RENAMED

REPLACEMENTS, FASTA, RENAMED = get_args()

print('Replacing names in ' + FASTA + ' to create ' + RENAMED +'.')

#Create list to store sequences
RECORDLIST = []
#Open the fasta output file for writing
with open(RENAMED, 'w') as OUT:
	#Open the fasta input file.
	for record in SeqIO.parse(FASTA, 'fasta'):
		#Import the list of old and new headers.
		with open(REPLACEMENTS, 'r') as LIST: 
			#Convert it to a list.
			THISLIST = csv.reader(LIST, delimiter='\t')
			#For every line in the list.
			for LINE in THISLIST:
				#Check if the old header is present in the current record.
				if LINE[0] in record.id:
					#If it is, replace with the new header.
					record.id = LINE[1]
					record.description = ''
			#Otherwise, just copy the record.
			RECORDLIST.append(record)
#Write the list to a new file.
SeqIO.write(RECORDLIST, RENAMED, 'fasta')

			

