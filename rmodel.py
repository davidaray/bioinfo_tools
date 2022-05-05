import os
import subprocess
import glob
import logging
import time
import datetime
import gzip
import shutil
from Bio import SeqIO
import argparse

LOGGER = logging.getLogger(__name__)

# Function 1: Set up input arguments
def get_args():
	parser = argparse.ArgumentParser(description="Will take in a genome (genomeabbreviation.fa.gz), run RepeatModeler and, subsequently, RepeatAfterMe to generate a pre-curation directory consisting of potential curated consensus sequences as well as image files for evaluation.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-g', '--genomepath', type=str, help='Path to a gzippped, fasta formatted genome assembly to be queried.', required=True)
	parser.add_argument('-p', '--processors', type=int, help='Number of processors to use for multi-processor enabled steps.', required = True)
	parser.add_argument('-b', '--batches', type=int, help='Number of batches to run during RepeatMasker step.', required = True)
	parser.add_argument('-w', '--workpath', type=str, help='Directory path where work is to be done.', required = True)
#	parser.add_argument('-rb', '--rightbuffer', type=int, help='Right beffer size. The number of bp of flanking sequence for each hit to be extracted along with the hit. Optional, Default = 1000', default = 1000)
#	parser.add_argument('-n', '--hitnumber', type=int, help='The number of hits to be exracted. Optional. Default = 50.', default = 50)
#	parser.add_argument('-a', '--align', type=str, help='Align the output fasta file, y or n?. Default is y.', default = 'y')
#	parser.add_argument('-t', '--trimal', type=str, help='Use trimal to remove low-aligning regions, y or n? Trimal can sometimes encounter an error that prevents it from working, this results in an empty file in downstream analyses. Default is y.', default = 'y')
#	parser.add_argument('-e', '--emboss', type=str, help='Generate a trimal/emboss consensus, y or n. Optional. Default=y.', default = 'y')
#	parser.add_argument('-m', '--maxiters2', type=str, help='Limit muscle iterations to 2? y or n. Optional. Default=n.', default = 'n')
	parser.add_argument('-log', '--log_level', default='INFO')
#	parser.add_argument('-beds', '--keep_beds', type=str, help='Keep the bed files used for extractions. y or n. Optional. Default=n.', default='n')

	args = parser.parse_args()
	GENOMEPATH = args.genomepath
	PROCESSORS = args.processors
	BATCHES = args.batches
	WORKPATH = args.workpath
#	RBUFFER = args.rightbuffer
#	HITNUM = args.hitnumber
#	ALIGN = args.align
#	TRIMAL = args.trimal
#	EMBOSS = args.emboss
#	MAXITERS = args.maxiters2
	LOG = args.log_level
#	BEDS = args.keep_beds

	return WORKPATH, GENOMEPATH, PROCESSORS, BATCHES, LOG

# Function 2: Directory establishment
def DIRS(DIR):
	if os.path.exists(DIR):
        	shutil.rmtree(DIR)
	os.mkdir(DIR)

# Function 3: Get genome from source to working directory 
def GETGENOME(GENOMEPATH, WORKPATH):
	shutil.copy2(GENOMEPATH, WORKPATH + 'assemblies_dir')

# Function 4: Run RepeatModeler on target genome 
def RMODEL(GENOMEPREFIX, PROCESSORS, WORKPATH, REPEATMODELERPATH):
	os.chdir(WORKPATH + 'assemblies_dir')
	with gzip.open(GENOMEPREFIX + '.fa.gz', 'rb') as GENOMEIN:  
		with open(GENOMEPREFIX + '.fa', 'wb') as GENOMEOUT:
			shutil.copyfileobj(GENOMEIN, GENOMEOUT)
			os.chdir(WORKPATH + 'rmodeler_dir') 
			subprocess.run(REPEATMODELERPATH + 'BuildDatabase -name {} {}'.format(GENOMEPREFIX, WORKPATH + 'assemblies_dir/' + GENOMEPREFIX + '.fa'), shell=True)
			subprocess.run(REPEATMODELERPATH + 'RepeatModeler -rscout_dir /lustre/work/daray/software/RepeatScout-1.0.11/ -database {} -pa {} > {}.RMrun.out'.format (GENOMEPREFIX, PROCESSORS, GENOMEPREFIX), shell=True)

# Function 5: Reformat RepeatModeler output 
def REFORMATRM(GENOMEPREFIX, WORKPATH):
#    os.chdir(WORKPATH)
########Reinstate these lines when not testing
	os.chdir(WORKPATH + 'rmodeler_dir')
	with open(GENOMEPREFIX + '.RMrun.out') as SEARCH:
		for LINE in SEARCH:
			LINE = LINE.rstrip()
			if 'Using output directory' in LINE:
				RMOUTDIR = LINE.split()[4] 
	os.chdir(RMOUTDIR)
	ORIGINAL = GENOMEPREFIX + '-families.fa'
	EDITED = GENOMEPREFIX + '-families_mod.fa'
	with open(ORIGINAL) as ORIGINALFILE, open(EDITED, 'w') as EDITEDFILE:
		records = SeqIO.parse(ORIGINALFILE, 'fasta')
		for record in records:
			HEADER = record.id
			NEWHEADER = HEADER.split()[0]
			NEWHEADER = NEWHEADER.replace('#', '__')
			NEWHEADER = NEWHEADER.replace('/', '___')
			NEWHEADER = NEWHEADER.replace('rnd', GENOMEPREFIX + '-rnd')
			record.id = NEWHEADER
			record.description = NEWHEADER 
			SeqIO.write(record, EDITEDFILE, 'fasta')

#### End of subsidiary functions


####MAIN function
def main():	
##Get input arguments
	WORKPATH, GENOMEPATH, PROCESSORS, BATCHES, LOG = get_args()
    
	print('It is working')


# Setup logging and script timing
	handlers = [logging.FileHandler('extract_align.log'), logging.StreamHandler()]
	logging.basicConfig(format='', handlers = handlers)
	logging.getLogger().setLevel(getattr(logging, LOG.upper()))

	start_time = time.time()
    
# Establish paths
	SOFTWARE = '/lustre/work/daray/software/'
	REPEATMODELERPATH = '/lustre/work/daray/software/RepeatModeler/'
	CLUSTERRUNPATH = '/home/daray/gitrepositories/bioinfo_tools/slurm_clusterrun.py'
################Directories for testing only
#    WORKPATH = '/lustre/scratch/daray/dVir_TEs/'
#    GENOMEPATH = '/lustre/scratch/daray/dVir_TEs/testdir/dVir_160.fa.gz'

# Set up directories
	DIRS(WORKPATH + 'assemblies_dir')
	DIRS(WORKPATH + 'rmodeler_dir')
	DIRS(WORKPATH + 'rmasker_dir')

# Get genome assembly to the assemblies directory
	GETGENOME(GENOMEPATH, WORKPATH)
    
# Determine genome abbreviation
	GENOMEPREFIX = os.path.basename(GENOMEPATH)    
	GENOMEPREFIX = os.path.splitext(GENOMEPREFIX)[0]
	GENOMEPREFIX = os.path.splitext(GENOMEPREFIX)[0]

# Print information for user
	LOGGER.info('#\n# TE_curate.py\n#')

	LOGGER.info('Genome file: ' + GENOMEPATH)
	LOGGER.info('Number of processors: ' + str(PROCESSORS))
	LOGGER.info('Number of batches for RepeatMasker: ' + str(BATCHES))
	LOGGER.info('Working directory: ' + WORKPATH)
	LOGGER.info('Genome prefix to use: ' + GENOMEPREFIX)
#	LOGGER.info('Number of hits evaluated: ' + str(HITNUM))
#	LOGGER.info('Muscle alignment = ' + ALIGN)
#	LOGGER.info('Trimal processing = ' + TRIMAL)
#	LOGGER.info('Emboss consensus = ' + EMBOSS)
#	LOGGER.info('Keep bed files = ' + BEDS)
	LOGGER.info('Log level: ' + LOG)

    
# Run RepeatModeler
	RMODEL(GENOMEPREFIX, PROCESSORS, WORKPATH, REPEATMODELERPATH)

# Reformat RepeatModeler output
	REFORMATRM(GENOMEPREFIX, WORKPATH)
    
# Run RepeatMasker using the RepeatModeler library
#	REPEATMASK(GENOMEPREFIX, BATCHES)

if __name__ =="__main__":main()	