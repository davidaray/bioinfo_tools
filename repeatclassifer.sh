#!/bin/bash 
#SBATCH --job-name=rclassify
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=60G

module --ignore-cache load gcc/10.1.0 r/4.0.2
. ~/conda/etc/profile.d/conda.sh
conda activate repeatmodeler

cd /lustre/scratch/daray/ixodes/iRic1/repeatclassifier

RepeatClassifier -consensi iRic1_extended_rep.fa


