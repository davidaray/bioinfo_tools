import os
import argparse
from Bio import SeqIO
import re

####################################
## Usage: python parseTEs.py -f <fasta file with renamed TEs from 200 mammals project> 
##
## Script will also output a file with name equivalencies, "new name <tab> old name".

def get_args():
	#What this script does
	parser = argparse.ArgumentParser(description="Replace headers from TE curation efforts using results from 95% cd-hit-est cluster files and a list of taxon names.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	#Argument for input library used for cd-hit-est
	parser.add_argument('-f', '--fasta', type=str, help='Name of the library of renamed TEs from 200 mammals.', required=True)

	args = parser.parse_args()
	FASTA = args.fasta
	
	return FASTA

def main():

#This doesn't
	FASTA = get_args()
	
	
	with open('DNA.fa', 'w') as DNA, \
	open('LTR.fa', 'w') as LTR, \
	open('RC.fa', 'w' ) as RC, \
	open('Unknown.fa', 'w') as UNKNOWN, \
	open('bigLINEs.fa', 'w') as BIGLINES, \
	open('medLINEs.fa', 'w') as MEDLINES, \
	open('smallLINEs.fa', 'w') as SMALLLINES, \
	open('bigSINEs.fa', 'w') as BIGSINES, \
	open('smallSINEs.fa', 'w') as SMALLSINES, \
	open('Other.fa', 'w') as OTHER:
		for record in SeqIO.parse(FASTA, 'fasta'):
			if 'DNA' in record.id:
				SeqIO.write(record, DNA, 'fasta')
			elif 'LTR' in record.id:
				SeqIO.write(record, LTR, 'fasta')
			elif 'RC' in record.id:
				SeqIO.write(record, RC, 'fasta')		
			elif 'Unknown' in record.id:
				SeqIO.write(record, UNKNOWN, 'fasta')		
			elif ('LINE' in record.id) and (len(record.seq) >= 3000):
				SeqIO.write(record, BIGLINES, 'fasta')		
			elif ('LINE' in record.id) and (len(record.seq) >= 700) and (len(record.seq) < 3000):
				SeqIO.write(record, MEDLINES, 'fasta')
			elif ('LINE' in record.id) and (len(record.seq) < 700):
				SeqIO.write(record, SMALLLINES, 'fasta')		
			elif ('SINE' in record.id) and (len(record.seq) >= 700):
				SeqIO.write(record, BIGSINES, 'fasta')		
			elif ('SINE' in record.id) and (len(record.seq) < 700):
				SeqIO.write(record, SMALLSINES, 'fasta')	
			else:
				SeqIO.write(record, OTHER, 'fasta')

if __name__ =="__main__":main()	
			
