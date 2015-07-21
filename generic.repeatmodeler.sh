#!/bin/bash
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N gator.RMod
#$ -o $JOB_NAME.o$JOB_ID
#$ -e $JOB_NAME.e$JOB_ID
#$ -q Chewie,R2D2
#$ -pe sm 3
#$ -P communitycluster
 

#Preparation:
#Create a working directory
mkdir <path to working directory>

#Migrate the the working directory
cd <path to working directory>

#Create a soft link to the genome you plan to analyze in the working directory make sure extension of the file is .fa. 
ln -s <path to genome file> .

#Create abbreviation for genome name
GENOME_NAME=$(basename <genome file> .masked.fa)

#This script assumes you have already masked the genome for known elements and have a file called XXX.masked.fa
/lustre/work/daray/software/RepeatModeler-1.0.8/BuildDatabase -name $GENOME_NAME -engine ncbi $GENOME_NAME.masked.fa
/lustre/work/daray/software/RepeatModeler-1.0.8/RepeatModeler -engine ncbi -pa 3 -database $GENOME_NAME
