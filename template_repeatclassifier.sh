#!/bin/bash 
#SBATCH --job-name=<NAME>_classify
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=60G

module --ignore-cache load gcc/10.1.0 r/4.0.2
. ~/conda/etc/profile.d/conda.sh
conda activate repeatmodeler

TAXON=<NAME>
WORKDIR=/path/to/working/directory
mkdir -p $WORKDIR/repeatclassifier

#Fix the names & create the file for classification
cd $WORKDIR/extensions/final_consensuses
for FILE in *.fa; do
	cp $FILE ${TAXON}-${FILE}
	sed -i "s/rnd/${TAXON}-rnd/g" ${TAXON}-$FILE
done

cat ${TAXON}*.fa >$WORKDIR/${TAXON}_extended_rep.fa


#Run RepeatClassifier
cd $WORKDIR
RepeatClassifier -consensi ${TAXON}_extended_rep.fa


