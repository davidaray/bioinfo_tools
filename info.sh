#Custom shell script to process Kevin's paired genome survey sequences from Peromyscus
#Will take the assembled (peared) files in .fastq format and create a table file with the following columns
#filename, #seqs, 'average', <value>, 'maximum', <value>, 'minimum', <value>

for FASTQ in /lustre/work/daray/pear/*_.fastq.assembled.fastq
	do ID=$(basename $FASTQ _.fastq.assembled.fastq)
	fastq_to_fasta -i $FASTQ -o $ID"_assembled.fa"
	perl /lustre/work/daray/software/file_analyze.fasta.pl $ID"_assembled.fa" >$ID"_assembled.info"
	sed ':a;N;$!ba;s/\n/\t/g' $ID"_assembled.info" |	cut -d" " -f1,3,6,10,12,16,18,22 >$ID"_assembled.info.tbl"	
done

cat *.tbl >info_all.tbl

