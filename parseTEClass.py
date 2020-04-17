import os
import argparse
from Bio import SeqIO
import re

####################################
## Usage: python parseTEClass.py -f <Name of the library of analyzed TEs from TEClass.> 
##

def get_args():
	#What this script does
	parser = argparse.ArgumentParser(description="Compare header and result of TEClass analysis.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	#Argument for input library used for cd-hit-est
	parser.add_argument('-f', '--fasta', type=str, help='Name of the library of analyzed TEs from TEClass.', required=True)

	args = parser.parse_args()
	FASTA = args.fasta
	
	return FASTA

def parserecorddesc(ID):
	IDREPLACE = ID.replace('|', ',')
	IDREPLACE = IDREPLACE.replace(': ', ',')
	IDSPLIT = re.split(',', IDREPLACE)
	TECLASSRESULT = IDSPLIT[2]
	HEADER = IDSPLIT[0]
	SPLITHEADER = re.split('#', HEADER)
	SPLITHEADER = SPLITHEADER[1].replace('/', ',')
	SPLITHEADER = re.split(',', SPLITHEADER)
	RMASKCALL = SPLITHEADER[0]
	
	return HEADER, RMASKCALL, TECLASSRESULT
	

FASTA = get_args()
	
PREFIX = os.path.splitext(FASTA)[0]
#PREFIX = PREFIX[0]
with open(PREFIX + '_teClass_comparison.txt', 'w') as OUT:
	with open(PREFIX + '_tocheck.fa', 'w') as OUTFASTA:
		with open(PREFIX + '_valid.fa', 'w') as VALIDFASTA:
			for record in SeqIO.parse(FASTA, 'fasta'):
				HEADER, RMASKCALL, TECLASSRESULT = parserecorddesc(record.description)
				if RMASKCALL == TECLASSRESULT:
					MATCH = 'exact_match'
					SeqIO.write(record, VALIDFASTA, 'fasta')
				else:
					MATCH = 'no_match'
					SeqIO.write(record, OUTFASTA, 'fasta')
				OUT.write(HEADER + '\t' + RMASKCALL + '\t' + TECLASSRESULT + '\t' + MATCH + '\n')
			

			
