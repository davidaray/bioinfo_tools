#############################################
# usage: slurm_clusterrun.py [-h] -i INPUT [-sp SPECIES] 
# [-b BATCH_COUNT] -dir GENOME_DIR [-od OUTDIR] 
# [-lib LIBRARY] [-xsmall] [-nolow] [-s {q,s}]
# [-p PROCESSORS] [-q {quanah,nocona}]
#
# Replaces the generateSGEclusterrun perl scripts from Robert Hubley.

import sys
from pathlib import Path
import os
import argparse
import itertools
import subprocess
from pyfaidx import Faidx
from pyfaidx import Fasta
import time
import re
import stat
import shutil
from shutil import copyfile
import errno
from Bio import SeqIO
from Bio import SeqRecord
from Bio import Seq
import Bio.SeqIO as IO

# Where RepeatMasker is stored
REPEATMASKER = "/lustre/work/daray/software/RepeatMasker"
# Where this script can find liftUp, twoBitInfo and twoBitToFa
BIN_DIR = "/lustre/work/daray/software"

# Define arguments
def get_args():
    #What this script does
	parser = argparse.ArgumentParser(description="Generate SGE cluster runs for RepeatMasker; built in RepeatMasker parameters are -xsmall [softmasks repetitive regions] -a [.align output file] -gff [generates a GFF format output] -pa [runs in parallel], please see RepeatMasker for details of these run options", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	required = parser.add_argument_group('required arguments')
    #Give input genome FASTA
	parser.add_argument('-i', '--input', type=str, help='genome file in FASTA format', required=True)
    #Argument of species name
	parser.add_argument('-sp', '--species', type=str, help='Species library from RepBase to be used.', required=False)
    # Desired batch number
	parser.add_argument('-b', '--batch_count', type=int, help='Batch count', default=50)
    # Input genome directory
	parser.add_argument('-dir', '--genome_dir', type=str, help='Path to genome FASTA', required=True)
    # Argument for output directory
	parser.add_argument('-od', '--outdir', type=str, help='Location of directory for the output subdirectory', default='.')
    # Which queue to use
	parser.add_argument('-q', '--queue', type=str, help='Select the queue to run RepeatMasker in [quanah|nocona].', choices=['quanah', 'nocona'], default='nocona')
    #Argument of RepeatMasker run parameter
	parser.add_argument('-lib', '--library', type=str, help='RepeatMasker run parameter custom library "-lib [filename]" option', required=False)
    #Argument of RepeatMasker run parameter
	parser.add_argument('-xsmall', action='store_true', help='Select a RepeatMasker masking option as lowercase bases [-xsmall], default is to mask as Ns')
    #Argument of RepeatMasker run parameter
	parser.add_argument('-nolow', action='store_true', help='RepeatMasker parameter flag "-nolow" option; does not mask low complexity DNA or simple repeats')
    #Argument of RepeatMasker run parameter
	parser.add_argument('-s', '--speed', type=str, help='RepeatMasker run parameter "-q" or "-s" option; q=quick search; 5-10% less sensitive, 3-4 times faster than default; s=slow search; 0-5% more sensitive, 2.5 times slower than default', choices=['q', 's'], required=False)
	#Argument for number of processors to use during repeatmasker runs
	parser.add_argument('-p', '--processors', type=int, help='Number of processors to use for each RepeatMasker run. Default = 10.', required=False, default=10)
    
	args = parser.parse_args()
	GENOME = args.input
	SPECIES = args.species
	BATCH_COUNT = args.batch_count
	GENOME_DIR = args.genome_dir
	OUTDIR = args.outdir
	QUEUE = args.queue
	LIBRARY = args.library
	XSMALL = args.xsmall
	NOLOW = args.nolow
	SPEED = args.speed
	PROC=args.processors

	return GENOME, SPECIES, BATCH_COUNT, GENOME_DIR, OUTDIR, LIBRARY, XSMALL, NOLOW, SPEED, PROC, QUEUE
    
REPEATMASKERPATH = '/lustre/work/daray/software/RepeatMasker'
SOFTWARE = '/lustre/work/daray/software'
ACTCONDA = '. ~/miniforge3/etc/profile.d/conda.sh'

# Function: build doLift.sh with UCSC liftUp and preserve '*' in .out files
def buildDoLift(GENOME_NAME, OUTDIR, QUEUE, REPEATMASKERPATH):
    import os

    PARTITION_DIR = os.getcwd()
    print("Creating doLift.sh file...\n")
    out_path = os.path.join(OUTDIR, 'doLift.sh')

    with open(out_path, 'w') as OUT:
        OUT.write('#!/bin/bash\n')
        OUT.write(f'#SBATCH --job-name={GENOME_NAME}.2Bit-doLift\n')
        OUT.write('#SBATCH --output=%x.%j.out\n')
        OUT.write('#SBATCH --error=%x.%j.err\n')
        OUT.write(f'#SBATCH --partition={QUEUE}\n')
        OUT.write('#SBATCH --nodes=1\n')
        OUT.write('#SBATCH --ntasks=12\n')
        OUT.write('#SBATCH --time=24:00:00\n\n')
        OUT.write('. ~/miniforge3/etc/profile.d/conda.sh\n')
        OUT.write('conda activate\n\n')
        OUT.write(f'cd {PARTITION_DIR}\n\n')
        OUT.write(f'REPEATMASKERPATH={REPEATMASKERPATH}\n')
        OUT.write('SOFTWARE=/lustre/work/daray/software\n')
        OUT.write('CURATIONDIR=/lustre/scratch/daray/bat1k_TE_analyses/curation_templates\n\n')
        
        OUT.write('# Lift each partition\n')
        OUT.write('for d0 in RMPart/???/\n')
        OUT.write('do\n')
        OUT.write('  d0=${d0::-1}\n')
        OUT.write('  bNum=$(basename $d0)\n')
        OUT.write('  $SOFTWARE/ucscTools/liftUp -type=.out stdout $d0/$bNum.lft error $d0/$bNum.fa.out > $d0/$bNum.fa.liftedOut\n')
        OUT.write('  if [ -f $d0/$bNum.fa.align ]; then\n')
        OUT.write('    $SOFTWARE/ucscTools/liftUp -type=.align stdout $d0/$bNum.lft error $d0/$bNum.fa.align > $d0/$bNum.fa.liftedAlign\n')
        OUT.write('  fi\n')
        OUT.write('done\n\n')

        OUT.write('# Combine all lifted outputs\n')
        OUT.write(f'$SOFTWARE/ucscTools/liftUp {GENOME_NAME}.fa.out /dev/null carry RMPart/???/*.liftedOut\n')
        OUT.write(f'mv {GENOME_NAME}.fa.out {GENOME_NAME}.fa.out.original\n')
        OUT.write(f'$SOFTWARE/ucscTools/liftUp {GENOME_NAME}.fa.align /dev/null carry RMPart/???/*.liftedAlign\n\n')
        
        OUT.write('# Re-attach * from original .out files\n')
        OUT.write(f'python $CURATIONDIR/replace_stars_v2.py -d RMPart/ -i {GENOME_NAME}.fa.out.original -o {GENOME_NAME}.fa.raw.out\n\n')

        OUT.write('# Process for overlaps\n')
        OUT.write(f'python $CURATIONDIR/process_repeatmasker_v7.py -i {GENOME_NAME}.fa.raw.out -p {GENOME_NAME} -ot both -ov lower_divergence -t $SLURM_NTASKS --progress-dir progress_files\n\n')

        OUT.write('# Build RepeatMasker summary\n')
        OUT.write(f'perl $REPEATMASKERPATH/util/buildSummary.pl -useAbsoluteGenomeSize -genome {GENOME_NAME}.2bit {GENOME_NAME}.out > {GENOME_NAME}.summary\n\n')

        OUT.write('# Generate RepeatLandscape\n')
        OUT.write(f'perl $REPEATMASKERPATH/util/calcDivergenceFromAlign.pl -s {GENOME_NAME}.divsum {GENOME_NAME}.fa.align\n')
        OUT.write(f'perl $SOFTWARE/createRepeatLandscape.mod.pl -div {GENOME_NAME}.divsum -twoBit {GENOME_NAME}.2bit > {GENOME_NAME}-landscape.html\n\n')

        OUT.write('# Zip up final RepeatMasker outputs and genome inputs\n')
        OUT.write('conda activate pigz\n')
        OUT.write(f'pigz -p $SLURM_NTASKS {GENOME_NAME}.summary\n')
        OUT.write(f'pigz -p $SLURM_NTASKS {GENOME_NAME}.divsum\n')
        OUT.write(f'pigz -p $SLURM_NTASKS {GENOME_NAME}.out\n')
        OUT.write(f'pigz -p $SLURM_NTASKS {GENOME_NAME}.fa.out.original\n')
        OUT.write(f'pigz -p $SLURM_NTASKS {GENOME_NAME}.fa.raw.out\n')
        OUT.write(f'pigz -p $SLURM_NTASKS {GENOME_NAME}.bed\n')
        OUT.write(f'pigz -p $SLURM_NTASKS {GENOME_NAME}.fa.align\n')
        OUT.write(f'pigz -p $SLURM_NTASKS {GENOME_NAME}.fa\n')
        OUT.write(f'pigz -p $SLURM_NTASKS {GENOME_NAME}.2bit\n\n')

        OUT.write('# Signal completion\n')
        OUT.write('touch alldone.ok\n')

    print(f"doLift.sh written to {out_path}")

#Function 2: Creating batch files for RepeatMasker runs.
def create_batch(RECORDS, CHUNK_SIZE):
    # Create an object out of each record and go through them iteratively
	RECORD_IT = iter(RECORDS)
    # The "next" object in the list of "records"...
    # Basically going through each contig one at a time
	#print('The next object in the list of records')
	RECORD = next(RECORD_IT)
    # Initialize base pair counting
	CURRENT_BASE = 0
    # Create a dictionary for batches and initialize the batch size
    # "Batch" is defined as the output fasta file that has a collection of "chunks"
	#print('Create a dictionary for batches and initialize the batch size')
	BATCH = []
	BATCH_SIZE = 0
    # While there are still records left in the list, keep creating new batches
	while RECORD:
		# Loop over records until the batch is full (i.e. reached the max chunk size), or there are no new records
		while BATCH_SIZE != CHUNK_SIZE and RECORD:
			# Define the end... which sums up to the chunk size
			END = CURRENT_BASE + CHUNK_SIZE - BATCH_SIZE
			# Define the output sequence, which is the current base (beginning base), seperated by a ":" and the end base of the contig
			SEQ = RECORD[CURRENT_BASE:END]
            # Define where to cut the contig off, which is the current base + the length of the output sequence defined above
			END_OF_SLICE = CURRENT_BASE + len(SEQ)
			# Create the fasta headers to match that of the original SGE script
            # <original_contig_name> ":" <beginning_base> "-" <end_base>
			FASTA_HEADER = RECORD.id + ":{}-{}".format(CURRENT_BASE, END_OF_SLICE)
            # Change the seq.id to the fasta header defined above.
			SEQ.id = SEQ.name = FASTA_HEADER
            # Set a blank description for the sequence.
            # For some reason this throws off Biopython if there is nothing present in the description object.
			SEQ.description = ''
            # Add the sequence to the current batch
			BATCH.append(SEQ)
            # This is where we start doing the math.
            # Add the lenth of the current sequence we are iterating through to the current base.
            # When doing this, we also need to keep track of the batch_size... we want to make everything as equal as possible.
			CURRENT_BASE += len(SEQ)
			BATCH_SIZE += len(SEQ)
            # When we have "added" all of the bases from the current sequence we are iterating through,
            # then we need to go and grab the next sequence in the list.
			if CURRENT_BASE >= len(RECORD):
				RECORD = next(RECORD_IT, None)
				CURRENT_BASE = 0
	# Once we have a batch with the correct size, yield the batch.
	# OR... we have run out of sequences in the genome, so stop.
		yield BATCH
		BATCH = []
		BATCH_SIZE = 0

####End of functions

###The main section of the script
#Get input arguments
GENOME, SPECIES, BATCH_COUNT, GENOME_DIR, OUTDIR, LIBRARY, XSMALL, NOLOW, SPEED, PROC, QUEUE = get_args()

#### Sanity checks======================== ####
#Announce inputs
print("The query genome is {}.".format(GENOME))
print("{} batches will be made.".format(str(BATCH_COUNT)))
print("The genome FASTA is located in '{}'.".format(GENOME_DIR))
print("The output directory is '{}'.".format(OUTDIR))
print("The job queue is {}.".format(QUEUE))
print("The number of processors to be used during each RepeatMasker run is " + str(PROC) + '\n')

#Check species and/or library provided and yell if not
if not SPECIES and not LIBRARY:
	sys.exit("Must supply value for option 'species' or 'lib'!")
if SPECIES and LIBRARY:
	sys.exit("Only supply a value for one option: 'species' or 'lib'! Not both!")
if not os.path.isdir(GENOME_DIR):
	sys.exit("The given genome directory, '{}', does not exist.".format(GENOME_DIR))
#Check if library is accessible and valid
if LIBRARY:
	try:
		if not os.path.getsize(LIBRARY) > 0:
			sys.exit("The library file, '{}', is empty.".format(LIBRARY))
	except OSError as e:
		sys.exit("The library file '{}' does not exist or is inaccessible.".format(LIBRARY))
#Check if library is valid
#try:
#	if not os.path.getsize(LIBRARY) > 0:
#		sys.exit("The library file, '{}', is empty.".format(LIBRARY))
#except OSError as e:
#	sys.exit("The library file '{}' does not exist or is inaccessible.".format(LIBRARY))

#Create output directory if it doesn't exist
if not os.path.isdir(OUTDIR):
	sys.exit("The output directory '{}' does not exist.".format(OUTDIR))

#Check for optional flags
FLAGS = [LIBRARY, XSMALL, NOLOW, SPEED]
if not FLAGS:
	print("Default RepeatMasker parameters used, no custom library, -inv -a -gff -pa options used.\n")
else:
	print("Custom parameters used:")
	if XSMALL:
		print("-xsmall flag used.")
	if NOLOW:
		print("-nolow flag used.")
	if LIBRARY:
		print("-lib flag used. Custom library is '{}'.".format(os.path.basename(LIBRARY)))
	if SPEED:
		print("-{} flag used. Search sensitivity has changed.".format(SPEED))
	print('\n')

#### Define files and directories======================== ####
#Get genome filename and path
GENOME_FASTA = os.path.join(GENOME_DIR, GENOME)
#Create path for genome partitions directory
PARTITION_DIR = os.path.join(GENOME_DIR, "RMPart")
#Check if partition directory exists, is empty. If not empty, empty it.
if not os.path.exists(PARTITION_DIR):
	try:
		os.mkdir(PARTITION_DIR)
	except OSError as e:
		if e.errno != errno.EEXIST:
			raise
	print("Made '{}' directory.\n".format(PARTITION_DIR))
else:
	if not os.listdir(PARTITION_DIR):
		print("'{}' is empty. Continuing.".format(PARTITION_DIR))
	else:
		print("'{}' is not empty. Removing contents and continuing.\n".format(PARTITION_DIR))
        ## To remove anything in the RMPart folder from a previous run (all symbolic links (not expected here) and files and subdirectories) without deleting the RMPart directory itself
		for FILE in os.listdir(PARTITION_DIR):
			FILE_PATH = os.path.join(PARTITION_DIR, FILE)
			try:
				shutil.rmtree(FILE_PATH)
			except OSError:
				os.remove(FILE_PATH)

#Set optional parameters for RepeatMasker run.				
LIB_OR_SPECIES = ""
if LIBRARY:
	LIB_OR_SPECIES = ' -lib ' + LIBRARY + ' '
else:
	LIB_OR_SPECIES = ' -species ' + SPECIES + ' '
ADD_PARAMS = str(LIB_OR_SPECIES)
if NOLOW:
	ADD_PARAMS = ADD_PARAMS + '-nolow '
if XSMALL:
	ADD_PARAMS = ADD_PARAMS + '-xsmall '
if SPEED:
	ADD_PARAMS = ADD_PARAMS + '-' + str(SPEED) + ' '
else: ADD_PARAMS = ADD_PARAMS + '-s '          
print('Setting ADD_PARAMS as: ' + ADD_PARAMS)

#Copy library into partition directory
if LIBRARY:
	shutil.copy(LIBRARY, PARTITION_DIR)
	LIB_FILE = os.path.basename(LIBRARY)
	LIBRARY = os.path.join(PARTITION_DIR, LIB_FILE)

#Prep for generating qsub.sh
THISDIR=os.getcwd()
QSUB_FILE_PATH = os.path.join(THISDIR, 'qsub.sh')
if os.path.isfile(QSUB_FILE_PATH):
		os.remove(QSUB_FILE_PATH)
	
#### Set some variables======================== ####
GENOME_NAME = os.path.basename(GENOME_FASTA).split(".fa")[0]
print('This is the genome_name: ' + GENOME_NAME)
SLOTS_PER_BATCH = PROC	
# Create a list of all the seq records inside of the genome
RECORDS = list(SeqIO.parse(GENOME, "fasta"))
# Sum of all lengths of all records diveded by "batch_number"
CHUNK_SIZE = sum(len(i) for i in RECORDS) // BATCH_COUNT + 1
print('CHUNK_SIZE = ' + str(CHUNK_SIZE))

####Build dictionary from genome index
RECORD_DICT = IO.to_dict(IO.parse(GENOME_FASTA, "fasta"))

#### Write out the batches as new fasta files============== #### 
for i, BATCH in enumerate(create_batch(RECORDS, CHUNK_SIZE)):
    #Name the file and keep track of the numbering.
	FILENAME = "{:03d}.fa".format(i)
    #Create the filepath to put the fasta file in the appropriate RMPart subdirectory
	NUM_DIR = "{:03d}".format(i)
	PATH2_BATCH_DIR = os.path.join(PARTITION_DIR, NUM_DIR)
	os.mkdir(PATH2_BATCH_DIR)
	FILEPATH = os.path.join(PARTITION_DIR, NUM_DIR, FILENAME)
    # Write all the batches' sequences and their appropriate headers to the output fasta file.
	SeqIO.write(BATCH, FILEPATH, "fasta")
    #Create the filepath to put the fasta file in the appropriate RMPart subdirectory
	NUM_DIR = "{:03d}".format(i)
    #Write the .lft file as well using the .fa file just created.
	LIFT_FILE = "{:03d}.lft".format(i)
	LIFT_FILEPATH = os.path.join(PARTITION_DIR, NUM_DIR, LIFT_FILE)
	OUTPUT = open(LIFT_FILEPATH,'w')
	with open(FILEPATH, 'r') as FIN:
		for LINE in FIN:
			if LINE.startswith('>'):
				HEADER = re.sub('>','',LINE)
				HEADER2 = re.sub('\n','',HEADER)
				CONTIG = HEADER2.split(":")[0]
				PART2_HEADER = HEADER2.split(":")[1]
				START = int(PART2_HEADER.split("-")[0])
				END = int(PART2_HEADER.split("-")[1])
				LENGTH = END-START
				ORIGINAL_CONTIG_LENGTH = len(RECORD_DICT[str(CONTIG)])
				OUTPUT.write(str(START) + '\t' + str(HEADER2) + '\t' + str(LENGTH) + '\t' + str(CONTIG) + '\t' + str(ORIGINAL_CONTIG_LENGTH) + '\n')
	OUTPUT.close()
    #write the batch-***.sh file as well using the same system
	SH_FILE = "batch-{:03d}.sh".format(i)
	SLURMBATCH_PATH = os.path.join(PARTITION_DIR, NUM_DIR)
	SLURMBATCH_FILE = os.path.join(PARTITION_DIR, NUM_DIR, SH_FILE)
	with open(SLURMBATCH_FILE, 'w') as BATCH_FILE, open(QSUB_FILE_PATH, 'a') as QSUB_WRITE:
		#Write batch file.
		BATCH_FILE.write( '#!/bin/bash\n')
		BATCH_FILE.write( '#SBATCH --job-name=' + GENOME_NAME + '.batch.' + NUM_DIR +'\n')
		BATCH_FILE.write( '#SBATCH --output=%x.%j.out\n')
		BATCH_FILE.write( '#SBATCH --error=%x.%j.err\n')
		BATCH_FILE.write( '#SBATCH --partition=' + QUEUE + '\n')
		BATCH_FILE.write( '#SBATCH --nodes=1\n')
		BATCH_FILE.write( '#SBATCH --ntasks=' + str(PROC) + '\n')
		BATCH_FILE.write( '\n')
		BATCH_FILE.write( '\n')
		BATCH_FILE.write( '\n')
		BATCH_FILE.write( ACTCONDA + '\n')
		BATCH_FILE.write( 'conda activate repeatmasker\n')
		BATCH_FILE.write( '\n')
		BATCH_FILE.write( "cd {}\n\n".format(SLURMBATCH_PATH))
		BATCH_FILE.write( "RepeatMasker{}-inv -a -gff -pa {} {}.fa >& run.log\n".format(ADD_PARAMS, str(PROC - 1), NUM_DIR))
		
		#Write to qsub.sh
		QSUB_WRITE.write('sbatch ' + SLURMBATCH_FILE + '\n')

#Generate a 2bit file for doLift.		
if os.path.isfile(GENOME_NAME + '.2bit'):
	print('2bit file already exists.\n')
else:
	print('Generating a 2bit file.\n')
	subprocess.run(SOFTWARE + '/faToTwoBit {} {} '.format(GENOME, GENOME_NAME + '.2bit'), shell=True)

buildDoLift(GENOME_NAME, OUTDIR, QUEUE, REPEATMASKERPATH)
