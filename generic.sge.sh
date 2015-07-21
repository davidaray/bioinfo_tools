#!/bin/bash
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N SGE
#$ -o $JOB_NAME.o$JOB_ID
#$ -e $JOB_NAME.e$JOB_ID
#$ -q R2D2
#$ -pe sm 1
#$ -P communitycluster


#Preparation:
#Create a working directory
mkdir <path to working directory>

#Migrate the the working directory
cd <path to working directory>

#Create a soft link to the genome you plan to analyze in the working directory make sure extension of the file is .fa. 
ln -s <path to genome file> .

#Create abbreviation for genome name
GENOME_NAME=$(basename <genome file> .fa)

/lustre/work/daray/software/faToTwoBit $GENOME_NAME".fa" $GENOME_NAME".2bit"
		
#run the generate RM cluster script and set to run against RepBase
/lustre/work/daray/software/generateSGEClusterRun_Chewie.pl \
				-nolow	\
              	-twoBit $GENOME".2bit" \
				-batch_count 10 \
				-species   #use -lib option if you have a custom library, otherwise decide the appropriate species designation

#change the permissions and submit the qsub script
chmod u+x qsub.sh doLift.sh			
#./qsub.sh


