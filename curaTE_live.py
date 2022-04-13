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
import pysed

LOGGER = logging.getLogger(__name__)

# Function 1: Set up input arguments
def get_args():
	parser = argparse.ArgumentParser(description="Will take in a genome (genomeabbreviation.fa.gz), run RepeatModeler and, subsequently, RepeatAfterMe to generate a pre-curation directory consisting of potential curated consensus sequences as well as image files for evaluation.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-g', '--genomepath', type=str, help='Path to a gzippped, fasta formatted genome assembly to be queried.', required=True)
	parser.add_argument('-p', '--processors', type=int, help='Number of processors to use for multi-processor enabled steps.', required = True)
	parser.add_argument('-b', '--batches', type=int, help='Number of batches to run during RepeatMasker step.', required = True)
	parser.add_argument('-w', '--workpath', type=str, help='Directory path where work is to be done. Must have the trailing "/".', required = True)
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

# Function 4: Run RepeatModeler on target genome and copy output to rmodeler directory
def RMODEL(GENOMEPREFIX, PROCESSORS, WORKPATH, REPEATMODELERPATH):
	os.chdir(WORKPATH + 'assemblies_dir')
	with gzip.open(GENOMEPREFIX + '.fa.gz', 'rb') as GENOMEIN:  
		with open(GENOMEPREFIX + '.fa', 'wb') as GENOMEOUT:
			shutil.copyfileobj(GENOMEIN, GENOMEOUT)
			os.chdir(WORKPATH + 'rmodeler_dir') 
			subprocess.run(REPEATMODELERPATH + 'BuildDatabase -name {} {}'.format(GENOMEPREFIX, WORKPATH + 'assemblies_dir/' + GENOMEPREFIX + '.fa'), shell=True)
			subprocess.run(REPEATMODELERPATH + 'RepeatModeler -database {} -pa {} > {}.RMrun.out'.format (GENOMEPREFIX, PROCESSORS, GENOMEPREFIX), shell=True)
	# Find name of RM directory with RepeatModeler run
	with open(GENOMEPREFIX + '.RMrun.out') as SEARCH:
		for LINE in SEARCH:
			LINE = LINE.rstrip()
			if 'Using output directory' in LINE:
				RMOUTDIR = LINE.split()[4] 
	LIBRARYPATH = RMOUTDIR + '/consensi.fa.classified'
	shutil.copy(LIBRARYPATH, WORKPATH + 'rmodeler_dir')
    
# Function 5: Reformat RepeatModeler output 
def REFORMATRM(GENOMEPREFIX, WORKPATH):
	os.chdir(WORKPATH + 'rmodeler_dir')
	# Find name of RM directory with RepeatModeler run
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
			if len(record.seq) >= 50:
				SeqIO.write(record, EDITED, 'fasta')
    
# Function 6: Build bridge file for submitting RepeatMasker doLift script to queue
def BUILDBRIDGE(GENOMEPREFIX, WORKPATH):
	os.chdir(WORKPATH + 'rmasker_dir')
	print("Creating bridge.sh file...\n")
	OUT = open('bridge.sh', 'w+')
	OUT.write( '#!/bin/bash\n')
	OUT.write( '#SBATCH --job-name=' + 'bridge_' + GENOMEPREFIX + '\n')
	OUT.write( '#SBATCH --output=%x.%j.out\n')
	OUT.write( '#SBATCH --error=%x.%j.err\n')
	OUT.write( '#SBATCH --partition=nocona\n')
	OUT.write( '#SBATCH --nodes=1\n')
	OUT.write( '#SBATCH --ntasks=1\n')
	OUT.write( '\n')
	OUT.write( 'cd ' + WORKPATH + 'rmasker_dir\n\n')
	OUT.write( '\n')
	OUT.write( 'sh qsub.sh\n')
	OUT.write( '\n')
	OUT.write( '#Creates a list of jobIDs to keep track of the RepeatMasker batches being run.\n')
	OUT.write( 'jobIDs=\"\"; for i in `squeue  | grep ' + GENOMEPREFIX + ' | awk \'{print $1}\'`; do jobIDs=$jobIDs,$i; done; jobIDs=\"${jobIDs:1}\"\n')
	OUT.write( '\n')
	OUT.write( '#Submits the doLift script to the queue but holds it until all jobs in the jobIDs list (the RepeatMasker batches) have cleared.\n')
	OUT.write( 'sbatch --dependency=afterok:$jobIDs doLift.sh\n')
	OUT.write( '\n')
	OUT.close()

# Function 7: Repeatmasker 
def REPEATMASK(GENOMEPREFIX, BATCHES, WORKPATH, CLUSTERRUNPATH):
	os.chdir(WORKPATH + 'rmodeler_dir')
	with open(GENOMEPREFIX + '.RMrun.out') as SEARCH:
		for LINE in SEARCH:
			LINE = LINE.rstrip()
			if 'Using output directory' in LINE:
				RMOUTDIR = LINE.split()[4] 
	os.chdir(WORKPATH + 'rmasker_dir')
	ASSEMBLYPATH = WORKPATH + 'assemblies_dir/' + GENOMEPREFIX + '.fa'
	LIBRARYPATH = RMOUTDIR + '/consensi.fa.classified'
	subprocess.run('python ' + CLUSTERRUNPATH + ' -i {} -b {} -lib {} -dir . -xsmall '.format (ASSEMBLYPATH, BATCHES, LIBRARYPATH), shell=True)
	subprocess.run('sbatch bridge.sh', shell=True)
	while not os.path.exists(WORKPATH + 'rmasker_dir/alldone.ok'):
		time.sleep(10)
		print('Waiting for dLift.sh to complete.')
	if os.path.isfile(WORKPATH + 'rmasker_dir/alldone.ok'):
		return
        
# Function 8: Generate Stockholm files
def STOCKHOLM(GENOMEPREFIX, WORKPATH, REPEATMODELERPATH):
	os.chdir(WORKPATH + 'ram_dir')
	ASSEMBLYPATH = WORKPATH + 'rmasker_dir/' + GENOMEPREFIX + '.2bit'
	ALIGNPATH = WORKPATH + 'rmasker_dir/' + GENOMEPREFIX + '.fa.align.gz'
	with gzip.open(ALIGNPATH, 'rb') as ALIGNIN:  
		with open(GENOMEPREFIX + '.fa.align', 'wb') as ALIGNOUT:
			shutil.copyfileobj(ALIGNIN, ALIGNOUT)
	subprocess.run(REPEATMODELERPATH + 'util/generateSeedAlignments.pl -taxon {} -assemblyFile {} {}.fa.align >generateSeeds.log'.format(GENOMEPREFIX, ASSEMBLYPATH, GENOMEPREFIX), shell=True)

# Function 9: Convert Stockholm files to fasta
def STOCK2FASTA(GENOMEPREFIX, WORKPATH, REPEATMODELERPATH):
	os.chdir(WORKPATH + 'ram_dir')
	for FILE in glob.glob('*.stk'):
		STKNAME = os.path.basename(FILE)
		TENAME = FILE.split('.')[0]
		FANAME = TENAME + '.fa'
		subprocess.run(REPEATMODELERPATH + 'util/Linup -i {} -fasta >{}-{}'.format(FILE, GENOMEPREFIX, FANAME), shell=True)
		os.rename(STKNAME, GENOMEPREFIX + '-' + STKNAME)

# Function 10: Prepare files and file structure for running RAM
def RAMPREP(WORKPATH):
	os.chdir(WORKPATH + 'ram_dir')
	DIRS('extensions')
	DIRS('extendlogs')
	for STKFILE in glob.glob('*.stk'):
		# Get file names and TE ID
		TENAME = STKFILE.split('.')[0]
		EXTENSIONDIR = WORKPATH + 'ram_dir/extensions/' + TENAME
		DIRS(EXTENSIONDIR)
#		STOCKHOLMDIR = WORKPATH + 'ram_dir/stockholm_files/'
#		EDITED = TENAME + '_mod.fa'
#		with open(FILE) as ORIGINALFILE, open(EDITED, 'w') as EDITEDFILE:
#			records = SeqIO.parse(ORIGINALFILE, 'fasta')
#			for record in records:
#				HEADER = record.id
#				NEWHEADER = HEADER.split('#')[0]
#				record.id = NEWHEADER
#				record.description = NEWHEADER 
#				SeqIO.write(record, EDITEDFILE, 'fasta')
		shutil.move(STKFILE, EXTENSIONDIR)

# Function 11: Run RAM
def EXTEND_STK (GENOMEPREFIX, WORKPATH, RMODUTILPATH):	
	for FOLDER in glob.glob(WORKPATH + 'ram_dir/extensions/*'):
		os.chdir(WORKPATH + 'ram_dir/extensions')
		ASSEMBLYPATH = WORKPATH + 'rmasker_dir/' + GENOMEPREFIX + '.2bit'
		CURRENT_TE = os.path.basename(FOLDER)
		print('Extending ' + CURRENT_TE)
		os.chdir(CURRENT_TE)
		FAMILYFILE = CURRENT_TE + '.stk'
		ALIGNMENTOUT = CURRENT_TE + '-extended.out'
		OUTPUTFILE = CURRENT_TE + '-extended.stk'
		subprocess.run('{}extend-stk.pl -assembly {} -input {} -msaout {} -output {} -d'.format(RMODUTILPATH, ASSEMBLYPATH, FAMILYFILE, ALIGNMENTOUT, OUTPUTFILE), shell=True)


# Function 12: Generate .png files for visualizations
def PNG (GENOMEPREFIX, WORKPATH, RMODUTILPATH):	
	for FOLDER in glob.glob(WORKPATH + 'ram_dir/extensions/*'):
		os.chdir(WORKPATH + 'ram_dir/extensions')
		ASSEMBLYPATH = WORKPATH + 'rmasker_dir/' + GENOMEPREFIX + '.2bit'
		CURRENT_TE = os.path.basename(FOLDER)
		print('Creating msa and img file: ' + CURRENT_TE)
		os.chdir(CURRENT_TE)
		FAMILYFILE = CURRENT_TE + '-extended.stk'
		ALIGNMENTOUT = CURRENT_TE + '-extended.out'
		OUTPUTFA = CURRENT_TE + '-extended.fa'
		subprocess.run('{}Linup -msa -genome {} -includeFlanking 50 {} >{}'.format(RMODUTILPATH, ASSEMBLYPATH, ALIGNMENTOUT, OUTPUTFA), shell=True)
		subprocess.run('{}visualizeAlignPNG.pl -genome {} -fasta {}'.format(RMODUTILPATH, ASSEMBLYPATH, FAMILYFILE), shell=True)





######Extra stuff. Ignore

		# Archive original files
		# Replace generic names with TE IDs in relevant files
#		MSAFILEIN = 'MSA-extended_with_rmod_cons.fa'
#		MSAFILEOUT = TENAME + '_MSA_extended.fa'
#		INFILE = open(MSAFILEIN, 'r')
#		OUTFILE = open(MSAFILEOUT, 'w')
#		for LINE in INFILE:
#			OUTFILE.write(LINE.replace('repam-newrep', TENAME))
#			OUTFILE.write(LINE.replace('CORECONS', 'CONSENSUS-' + TENAME))
#		INFILE.close()
#		OUTFILE.close()
#		INFILE = open(REPFILEIN, 'r')
#		OUTFILE = open(REPFILEOUT, 'w')
#		for LINE in INFILE:
#			OUTFILE.write(LINE.replace('repam-newrep', TENAME))
#		INFILE.close()
#		OUTFILE.close()
#		shutil.copyfileobj('img.png', TENAME + '.png')
#		ENTRIES = list(SeqIO.parse('repseq.unextended', 'fasta'))
#		COUNT = len(ENTRIES)
#		# Sort rejects, repeats with fewer than 10 hits
#		if COUNT < 10:
#			shutil.copyfileobj(TENAME + '.png', WORKPATH + 'ram_dir/images_and_alignments/rejects/' + TENAME + '.png')
#			shutil.copyfileobj(TENAME + '_MSA_extended.fa', WORKPATH + 'ram_dir/images_and_alignments/rejects/' + TENAME + '_MSA_extended.fa')
#		if COUNT > 9:
#			shutil.copyfileobj(TENAME + '_rep.fa', WORKPATH + 'ram_dir/images_and_alignments/rejects/' + TENAME + '_rep.fa')
#		# Sort possible segmental duplications, >10,000 bp consensus
#		if ((COUNT > 9) and (LENGTH > 13000)):
#			shutil.copyfileobj(TENAME + '.png', WORKPATH + 'ram_dir/images_and_alignments/possible_SD/' + TENAME + '.png')
#			shutil.copyfileobj(TENAME + '_MSA_extended.fa', WORKPATH + 'ram_dir/images_and_alignments/possible_SD/' + TENAME + '_MSA_extended.fa')
#			shutil.copyfileobj(TENAME + '_rep.fa', WORKPATH + 'ram_dir/images_and_alignments/possible_SD/' + TENAME + '_rep.fa')
#		if ((COUNT > 9) and (LENGTH <13000)):
#			shutil.copyfileobj(TENAME + '.png', WORKPATH + 'ram_dir/images_and_alignments/likely_TEs/' + TENAME + '.png')
#			shutil.copyfileobj(TENAME + '_MSA_extended.fa', WORKPATH + 'ram_dir/images_and_alignments/likely_TEs/' + TENAME + '_MSA_extended.fa')
#			shutil.copyfileobj(TENAME + '_rep.fa', WORKPATH + 'ram_dir/images_and_alignments/likely_TEs/' + TENAME + '_rep.fa')
#### End of subsidiary functions
#### End ignore stuff
    
####MAIN function
def main():	
##Get input arguments
	WORKPATH, GENOMEPATH, PROCESSORS, BATCHES, LOG = get_args()
    
	print('It is working')

# Establish paths
	SOFTWARE = '/lustre/work/daray/software/'
	REPEATMODELERPATH = '/lustre/work/daray/software/RepeatModeler/'
	CLUSTERRUNPATH = '/home/daray/gitrepositories/bioinfo_tools/slurm_clusterrun.py'
	RMODUTILPATH = '/lustre/work/daray/software/RepeatModeler/util/'
################Directories for testing only
#	WORKPATH = '/lustre/scratch/daray/dVir_TEs/'
#	GENOMEPATH = '/lustre/scratch/daray/dVir_TEs/testdir/dVir_160.fa.gz'
#	BATCHES = 50
#	PROCESSORS = 36

# Determine genome abbreviation
	GENOMEPREFIX = os.path.basename(GENOMEPATH)    
	GENOMEPREFIX = os.path.splitext(GENOMEPREFIX)[0]
	GENOMEPREFIX = os.path.splitext(GENOMEPREFIX)[0]
	
# Setup logging and script timing
	handlers = [logging.FileHandler('curaTE.log'), logging.StreamHandler()]
	logging.basicConfig(format='', handlers = handlers)
	logging.getLogger().setLevel(getattr(logging, LOG.upper()))

	start_time = time.time()

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

# Confirm and set up directories
	if not os.path.isfile(WORKPATH + 'assemblies_dir/' + GENOMEPREFIX + '.fa.gz'):
		DIRS(WORKPATH + 'assemblies_dir')
		GETGENOME(GENOMEPATH, WORKPATH)
    
	CLASSIFIEDFILE = WORKPATH + 'rmodeler_dir/consensi.fa.classified'
	if not os.path.isfile(CLASSIFIEDFILE):
		DIRS(WORKPATH + 'rmodeler_dir')
	
	ALIGNPATH = WORKPATH + 'rmasker_dir/' + GENOMEPREFIX + '.fa.align.gz'
	if not os.path.isfile(ALIGNPATH):
		DIRS(WORKPATH + 'rmasker_dir')

	RAMPATH = WORKPATH + 'ram_dir/extensions'
	if not os.path.isdir(RAMPATH):
		DIRS(WORKPATH + 'ram_dir')
    
# Run RepeatModeler if consensi.fa.classified does not exist
	CLASSIFIEDFILE = WORKPATH + 'rmodeler_dir/consensi.fa.classified'
	if not os.path.isfile(CLASSIFIEDFILE):
		RMODEL(GENOMEPREFIX, PROCESSORS, WORKPATH, REPEATMODELERPATH)

# Reformat RepeatModeler output
#	REFORMATRM(GENOMEPREFIX, WORKPATH)
    
# Run RepeatMasker using the RepeatModeler library if align.gz does not exist
	ALIGNPATH = WORKPATH + 'rmasker_dir/' + GENOMEPREFIX + '.fa.align.gz'
	if not os.path.isfile(ALIGNPATH):
		BUILDBRIDGE(GENOMEPREFIX, WORKPATH)
		REPEATMASK(GENOMEPREFIX, BATCHES, WORKPATH, CLUSTERRUNPATH)
    
# Get Stockholm files from RepeatMasker .align file
	STOCKHOLM(GENOMEPREFIX, WORKPATH, REPEATMODELERPATH)

# Convert all stockholm files to fasta
#	STOCK2FASTA(GENOMEPREFIX, WORKPATH, REPEATMODELERPATH)

# Run RAM on all fasta files
# First create the necessary directories
	os.chdir(WORKPATH + 'ram_dir')
	DIRS('images_and_alignments')
	DIRS('possible_SD')
	DIRS('likely_TEs')
	DIRS('rejects')

# Now create all the files to run RAM on
	RAMPREP(WORKPATH)

# Now run RAM
	EXTEND_STK(GENOMEPREFIX, WORKPATH, RMODUTILPATH)
    
# Generate visualizations
	PNG (GENOMEPREFIX, WORKPATH, RMODUTILPATH)

# Future home of AI evaluation of alignments

if __name__ =="__main__":main()	