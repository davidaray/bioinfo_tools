#/bin/bash
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N whale_split
#$ -o $JOB_NAME.o$JOB_ID
#$ -e $JOB_NAME.e$JOB_ID
#$ -q Chewie
#$ -pe fill 5
#$ -P communitycluster

BASDIR=/lustre/scratch/daray/whales
mkdir $BASDIR/faSplit
WORKDIR=$BASDIR/faSplit
WHLGEN=$BASDIR/genomes
SOFT=/lustre/work/daray/software

cd $WORKDIR

##### Split the genomes into chunks of ~500,000,000 bp #####

perl $SOFT/file_split.fasta.pl $WHLGEN/BalAcu.fasta -l 500000000               
perl $SOFT/file_split.fasta.pl $WHLGEN/OrcOrc.fasta -l 500000000
perl $SOFT/file_split.fasta.pl $WHLGEN/PhyCat.fasta -l 500000000
perl $SOFT/file_split.fasta.pl $WHLGEN/TurTru.fasta -l 500000000
perl $SOFT/file_split.fasta.pl $WHLGEN/LipVex.fasta -l 500000000

##### then combine the first chunks into a single file #####

cat $WHLGEN/*.fasta.0 > whale_cat.fa     

##### Create files to check number of bases in .0 files afterward to ensure the sizes are around 500 Mb #####

for FIRSTCHUNK in $WHLGEN/*.fasta.0
	do ABBREV=$(basename $FIRSTCHUNK .fasta.0)
	wc -m $WHLGEN/$ABBREV".fasta.0" >$ABBREV"_count.txt"
done

##### Check number of bases in cat files afterward to ensure the sizes are correct #####

perl $SOFT/file_analyze.fasta.pl *.fa >cat_info.txt
wc -m whale_cat.fa >whale_cat_count.txt
