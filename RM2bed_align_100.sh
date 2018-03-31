#USAGE: RM_post.sh <filename>
#Will process a .align.gz file from RepeatMasker to generate multiple sorted output files for downstream processing

ABBREV=$(basename $1 .align)
	
#gunzip $ABBREV".align"

perl /lustre/work/daray/software/RepeatMasker/util/calcDivergenceFromAlign.Ray.pl -noCpG $ABBREV".align" >$ABBREV"_div.out" 

#perl /lustre/work/daray/software/RepeatMasker/util/calcDivergenceFromAlign.Ray.pl $ABBREV".align" >$ABBREV"_wCpG_div.out"

perl /lustre/work/daray/software/RM.addsize.div.pl $ABBREV"_div.out" | grep -v Satellite | grep -v Simple_repeat | sed 's/kimura=//g' | sed 's/\#/\t/g' | sed 's/\//\t/g' >$ABBREV"_div_size.out"

#perl /lustre/work/daray/software/RM.addsize.div.pl $ABBREV"_wCpG_div.out" | grep -v Satellite | grep -v Simple_repeat | sed 's/kimura=//g' | sed 's/\#/\t/g' | sed 's/\//\t/g' >$ABBREV"_wCpG_div_size.out"

awk '{print $5 "\t" $6 "\t" $7 "\t" $11 "\t" $19 "\t" $10 "\t" $8 "\t" $12 "\t" $13 "\t" $18}' $ABBREV"_div_size.out" | sed 's/\tC\t/\t-\t/g' | sort -g -k 6,6 >$ABBREV"_sort_size.bed"

#awk '{print $5 "\t" $6 "\t" $7 "\t" $11 "\t" $19 "\t" $10 "\t" $8 "\t" $12 "\t" $13}' $ABBREV"_wCpG_div_size.out" | sed 's/\tC\t/\t-\t/g' >$ABBREV"_wCpG_processed.out"

awk '{if($7>=100)print;}' $ABBREV"_sort_size.bed" >$ABBREV"_select_size.bed"  

#sort -g -k 9,9 $ABBREV"_processed.out" >$ABBREV"_sort_div.out"

#awk '{if($4>=100)print;}' $ABBREV"_sort_div.out" >$ABBREV"_select_div.out"  

#sort -g -k 6,6 $ABBREV"_processed.out" >$ABBREV"_sort_name.out"

#awk '{if($4>=100)print;}' $ABBREV"_sort_name.out" >$ABBREV"_select_name.out"  

#sort -g -k 7,7 $ABBREV"_processed.out" >$ABBREV"_sort_class.out"

#awk '{if($4>=100)print;}' $ABBREV"_sort_class.out" >$ABBREV"_select_class.out"  

#sort -g -k 8,8 $ABBREV"_processed.out" >$ABBREV"_sort_family.out"

#awk '{if($4>=100)print;}' $ABBREV"_sort_family.out" >$ABBREV"_select_family.out"  

#sort -g -k 4,4 $ABBREV"_wCpG_processed.out" >$ABBREV"_wCpG_sort_size.out"

#awk '{if($4>=100)print;}' $ABBREV"_wCpG_sort_size.out" >$ABBREV"_wCpG_select_size.out"  

#sort -g -k 9,9 $ABBREV"_wCpG_processed.out" >$ABBREV"_wCpG_sort_div.out"

#awk '{if($4>=100)print;}' $ABBREV"_wCpG_sort_div.out" >$ABBREV"_wCpG_select_div.out"  

#sort -g -k 6,6 $ABBREV"_wCpG_processed.out" >$ABBREV"_wCpG_sort_name.out"

#awk '{if($4>=100)print;}' $ABBREV"_wCpG_sort_name.out" >$ABBREV"_wCpG_select_name.out"  

#sort -g -k 7,7 $ABBREV"_wCpG_processed.out" >$ABBREV"_wCpG_sort_class.out"

#awk '{if($4>=100)print;}' $ABBREV"_wCpG_sort_class.out" >$ABBREV"_wCpG_select_class.out"  

#sort -g -k 8,8 $ABBREV"_wCpG_processed.out" >$ABBREV"_wCpG_sort_family.out"

#awk '{if($4>=100)print;}' $ABBREV"_wCpG_sort_family.out" >$ABBREV"_wCpG_select_family.out"  

#gzip $ABBREV".align" &

#gzip $ABBREV"_div.out" &

#gzip $ABBREV"_wCpG_div.out" &

