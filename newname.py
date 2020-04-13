import os
import argparse
from Bio import SeqIO
import re

####################################
## Usage: python newname.py -c <cluster file from cd-hit-est> -l <list of taxon names to search for>
## -f <fasta file submitted to cd-hit-est> -r <name of the file to hold renamed fasta entries>
##
## Script will also output a file with name equivalencies, "new name <tab> old name".

def get_args():
	#What this script does
	parser = argparse.ArgumentParser(description="Replace headers from TE curation efforts using results from 95% cd-hit-est cluster files and a list of taxon names.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	#Argument of the input .clstr file
	parser.add_argument('-c', '--cluster_file', type=str, help='Name of your .clstr file from cd-hit-est.', required=True)
	#Argument for the taxon list.
	parser.add_argument('-l', '--list', type=str, help='List file of taxon names. One name per line..', required=True)
	#Argument for input library used for cd-hit-est
	parser.add_argument('-f', '--fasta', type=str, help='Name of the library used for cd-hit-est.', required=True)
	parser.add_argument('-r', '--renamed', type=str, help='Name of the output file with renamed .fa entries.', required=True) 

	args = parser.parse_args()
	CLUSTER = args.cluster_file
	LIST = args.list
	FASTA = args.fasta
	RENAMED = args.renamed
	
	return CLUSTER, LIST, FASTA, RENAMED

#Function to read in a line from a cluster file, modify to new header format and get the taxon name
def linetoheader(LINE):
	if '_family-' in LINE:
		LINE = LINE.replace('...', ',')
		LINE = re.split(', >', LINE)
		LINE = LINE[1]
		LINE = re.split(', ', LINE)
		PREHEADER = LINE[0]
		HEADER = PREHEADER.replace('_family-', '.')
		HEADER = HEADER.replace('___', '/')
		HEADER = HEADER.replace('__', '#')
		if '-1_family' in PREHEADER:
		    SPLIT = re.split('-1_family', PREHEADER)
		    TAXON = SPLIT[0]
		elif '-2_family' in PREHEADER:
		    SPLIT = re.split('-2_family', PREHEADER)
		    TAXON = SPLIT[0]
		elif '-3_family' in PREHEADER:
		    SPLIT = re.split('-3_family', PREHEADER)
		    TAXON = SPLIT[0]
		elif '-4_family' in PREHEADER:
		    SPLIT = re.split('-4_family', PREHEADER)
		    TAXON = SPLIT[0]
		elif '-5_family' in PREHEADER:
		    SPLIT = re.split('-5_family', PREHEADER)
		    TAXON = SPLIT[0]
		elif '-6_family' in PREHEADER:
		    SPLIT = re.split('-6_family', PREHEADER)
		    TAXON = SPLIT[0]
		elif '-7_family' in PREHEADER:
		    SPLIT = re.split('-7_family', PREHEADER)
		    TAXON = SPLIT[0]
		return PREHEADER, HEADER, TAXON
	#else:
		#pass

def main():

#This doesn't
	CLUSTER, LIST, FASTA, RENAMED = get_args()
	print('Creating best_hits.txt')
	#Create a best_hits.txt file with only the lines that have the taxa in our LIST
	#Open best_hits.txt for writing
	with open('best_hits.txt', 'w') as BEST:
		#open cluster file
		with open(CLUSTER, 'rt') as CLUSTERFILE:
			#open list file for reading
			with open(LIST, 'r') as F:
				#read in list.txt as list and strip \n from ends.
				THISLIST = F.read().splitlines()
				#for line in CLUSTER
				for LINE in CLUSTERFILE:
					LINE = re.sub("\*", 'longest', LINE)
					#Find lines of interest
					if '_family-' in LINE:
						#get TAXON value
						PREHEADER, HEADER, TAXON = linetoheader(LINE)
						#If TAXON is in the list
						if TAXON in THISLIST:
							#write to a new file for use in next block
							BEST.write(LINE)
	
	#Create a table and a new fasta file with the hew headers
	#open the best_hits.txt file for reading
	print('Opening best_hits.txt, name_equivalents, and renamed.fa')
	with open('best_hits.txt', 'rt') as TARGETS:
		#open the output name equivalencies file
		with open('name_equivalents.txt', 'w') as NAMES:
			#open the new fasta for writing
			with open(RENAMED, 'w') as OUT:
				#for every line in the best_hits.txt file
				print('Staring loop to find lines of interest')
				for LINE in TARGETS:
					#get the header information and taxon
					PREHEADER, HEADER, TAXON = linetoheader(LINE)
					#Write new header and original header to the name equivalency file
					NAMES.write(HEADER + '\t' + PREHEADER + '\n')
					#for every record in the original fasta submitted to cd-hit-est
					for record in SeqIO.parse(FASTA, 'fasta'):
						#Check if it's one we're interested in
						if (PREHEADER in record.id) and ('longest' in LINE):
							print('Found ' + PREHEADER)
							#change the header and save as the new file
							record.id = HEADER
							record.description = ''
							SeqIO.write(record, OUT, 'fasta')

if __name__ =="__main__":main()	
			

