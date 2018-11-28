import argparse
import gzip
import subprocess
import re
from Bio.Blast.Applications import NcbiblastnCommandline
import os

####MAIN function
def main():

##Use the get_args function
	GENOME, LIBRARY, BLASTOUTFILE, SEQBUFFER, SEQNUMBER = get_args()
	print('Input file = ' + GENOME)
	PREFIX = re.split("[.]", GENOME)[0]

##Determine file type using extension and create a tmp.file to use.
	if GENOME.lower().endswith('gz'):
		print('Using .gz file.')
		TMP = PREFIX + '.fa'
		with gzip.open(GENOME, 'rt') as IN:
			with open(TMP, 'w') as OUT:
				for LINE in IN:
					OUT.write(LINE)

#	moduleload = 'module load intel rmblast'
#	os.system(moduleload)
					
##Attempt 1
#	makedb = 'makeblastdb -in ' + TMP + '-dbtype nucl'
#	os.system(makedb)
#	blastn = 'blastn –query ' + LIBRARY + ' -db ' + TMP + ' -out ' + BLASTOUTFILE + ' -outfmt 6'
#	os.system(blastn)

##Attempt 2
	blastn_cline = NcbiblastnCommandline(query=LIBRARY, out=BLASTOUTFILE, db=TMP, outfmt=6)
	
##Attempt 3	
#	subprocess.check_call('makeblastdb -in {} -dbtype nucl'.format(TMP), shell=True)
#	subprocess.check_call('blastn –query {} -db {} -out {} -outfmt 6'.format(LIBRARY, TMP, BLASTOUTFILE), shell=True)

##Call subprocess to make the blast database and run blast
#	subprocess.check_call('perl /lustre/work/daray/software/extractAlignTEs.pl --genome {} --blast {} --consTEs {} --seqBuffer {} --seqNum {} --align'.format (TMP, BLASTOUTFILE, LIBRARY, SEQBUFFER, SEQNUMBER), shell=True)

##Delete tmp file.
#	print('Removing tmp file.')
#	os.remove(TMP)		

##Get arguments function
def get_args():
	parser = argparse.ArgumentParser(description="Assignment 4", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-g', '--genome', type=str, help='Name of the .fa or .fa.gz file to be analyzed', required=True)
	parser.add_argument('-l', '--library', type=str, help='Library of putative TEs', required=True)
	parser.add_argument('-b', '--blast', type=str, help='Name of blast file to generate and use')
	parser.add_argument('-sb', '--seqbuffer', type=int, help='Number of flanking bases to extract, default = 500 bp', default=500)
	parser.add_argument('-sn', '--seqnumber', type=int, help='Number of hits from each putative TE to analyze')


	args = parser.parse_args()
	GENOME = args.genome
	LIBRARY = args.library
	BLASTOUTFILE = args.blast
	SEQBUFFER = args.seqbuffer
	SEQNUMBER = args.seqnumber

	return GENOME, LIBRARY, BLASTOUTFILE, SEQBUFFER, SEQNUMBER

if __name__ =="__main__":main()
		
