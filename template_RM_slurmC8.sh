#!/bin/bash 
#SBATCH --job-name=<NAME>_RM
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=nocona
#SBATCH --mem-per-cpu=100G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-user david.a.ray@ttu.edu


###########################3
## Usage:
## This script will run the Centos-8 version of RepeatMasker on nocona. 
## You will need to replace <NAME> with the genome assembly identifier as appropriate.
## DIR may need modifying.
## sed "s|DIRPATH|path to your directory|g" template_RM_slurmC8.sh ><new file name>
## You will need to include the path to your library.
## sed "s|LIBFILE|path to your library file|g" <new file name> ><new new file name>
## You will need to include the path to your genome assembly.
## sed "s|ASSEMBLYPATH|path to your assembly file|g" <new file name> ><new new file name>
## If running lots of genomes, here's a quick way to do it using a 'list.txt' file that has all of your taxon names in it.
## Sed the DIRPATH, LIBFILE, and ASSEMBLYPATH, then:
## cat list.txt | while read i; do sed "s/<NAME>/$i/g" <new new file name> >${i}_RM_C8.sh; done
## The final script to run ${GENOME}_rm2bed.sh should be prepared before running this
## Then, submit the scripts using a similar loop.

#This script uses a python program. You need to make sure python is loaded and ready to go by activating conda.
#Your conda installation must have biopython and pyfaidx installed.
. ~/conda/etc/profile.d/conda.sh
conda activate

#makes variables used in this script
GENOME=<NAME>
RUNTYPE=${GENOME}_RM
DIR=DIRPATH/$RUNTYPE
LIBRARYPATH=LIBFILE
GENOMEPATH=ASSEMBLYPATH/${GENOME}.fa

# If GenomeAbbreviation is greater than 10 characters, will cause trouble
# with queueing of doLift.sh.
# Shorten GENOME and to 10 characters and call it SUBNAME to resolve this problem.
# Note that the qstat queries in the bottom few lines look for SUBNAME instead of GENOME.
CHARLEN=${#GENOME} #Finds length of GenomeAbbreviation
if test $CHARLEN -gt 10; then SHORTNAME=${GENOME:0:10}; else SHORTNAME=$GENOME; fi

#Creates a working directory and then goes into it to make all the files
mkdir $DIR
cd $DIR

#Creates a link to the genome assembly and unzip it. Then create a .2bit version for RepeatMasker.
ln -s $GENOMEPATH/${GENOME}.fa
#gunzip -c $GENOME".fa.gz" > $GENOME".fa"
/lustre/work/daray/software/faToTwoBit ${GENOME}.fa ${GENOME}.2bit
#creates a symbolic link to the 200 mammals TE library file so it can be used.
#cp /lustre/scratch/aosmansk/250mammals_RM_v2/final_mammal_library.fa .

#Use sge_clusterrun.py to generate all batches needed to run RepeatMasker and the doLift.sh to complile the results.
python /home/daray/gitrepositories/bioinfo_tools/slurm_clusterrunC8.py \
        -i ${GENOME}.fa \
        -b 50 \
        -lib $LIBRARYPATH \
        -dir . \
        -xsmall 

#Submits the RepeatMasker jobs created by sge_clusterrun.py
sh qsub.sh

#Creates a list of jobIDs to keep track of the RepeatMasker batches being run.
jobIDs=""; for i in `squeue  | grep $SHORTNAME | awk '{print $1}'`; do jobIDs=$jobIDs,$i; done; jobIDs="${jobIDs:1}"

#Submits the doLift script to the queue but holds it until all jobs in the jobIDs list (the RepeatMasker batches) have cleared.
sbatch --dependency=afterok:$jobIDs doLift.sh

#Creates a list of job IDs to keep track of the doLift script.
jobIDs=""; for i in `squeue  | grep $SHORTNAME | awk '{print $1}'`; do jobIDs=$jobIDs,$i; done; jobIDs="${jobIDs:1}"

#Submits the rm2bed script to the queue but holds it until all jobs in the jobIDs list (the doLift job) have cleared.
sbatch --dependency=afterok:$jobIDs ../${GENOME}_rm2bed.sh
