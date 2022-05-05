#!/bin/bash 
#SBATCH --job-name=rmodel_Ssci
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --ntasks=36


. ~/conda/etc/profile.d/conda.sh
conda activate

#Define variables
#PROCESSORS = number of processors you want to use
PROCESSORS=36
#BATCHES = number of batches to use for RepeatMasker run, if using
BATCHES=50

cd /path/to/working/directory

python rmodel.py -g /path/to/genome/assembly.$NAME.fa.gz -p $PROCESSORS -b $BATCHES -w /lustre/scratch//npaulat/arachnids/Ssci/
