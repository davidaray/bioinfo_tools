#6 Feb 2015
use lib '/lustre/work/daray/software/perlmodules/share/perl5/Bio';
use Bio::SearchIO; 
use Bio::SeqIO;
use Getopt::Long;

#TO DO
# get around using bioperl
# include options for repeatMasker, BED, gFF, SAM
# packaged alignment
# clean up unaligned files
# fix usage statement
 



#sets all input values to undef
my $genome_in 	= undef;
my $blast_in 	= undef;
my $te_in 	= undef;
my $buffer_in 	= undef;
my $seq_num_in 	= undef; 
my $align 	= undef; #sets align default to unalign
my $help	= undef;

GetOptions('genome=s'	 => \$genome_in,
	   'blast=s'	 => \$blast_in,
           'consTEs=s'	 => \$te_in,
	   'seqBuffer=i' => \$buffer_in,
           'seqNum=i'	 => \$seq_num_in,
           'align'       => \$align, 
           'h'	         => \$help );


#if input values not appropriate print usage statement.
if(defined $help){
	usage();
}

if( ! defined $genome_in ){
	print "\nERROR: No genome file given\n";
	usage();
}

if( ! defined $blast_in ){
	print "\nERROR: No blast file given\n";
	usage();
}

if( ! defined $te_in ){
	print "\nERROR: No file of consensus TEs given\n";
	usage();
}

if( ! defined $buffer_in ){
	print "\nERROR: No sequence buffer given for extractions\n";
	usage();
}

if( ! defined $seq_num_in ){
	print "\nError: Number of sequences to extract are not indicated\n";
	usage();
}

#step 1 organizing blast output----------------------------------------------
 &organize_blast_hits($blast_in, $buffer_in, $seq_num_in);

#step two creating ouput files with the consensus ---------------------------
 &create_te_out_files($te_in);

#step three loading genome into memory------------------------------


my $inseq = Bio::SeqIO->new(-file   => "<$genome_in", -format =>"fasta");
my %genome;
while (my $seq = $inseq->next_seq) {
	my $id = $seq->display_id;
	$genome{$id} = $seq->seq;
	}




    open (PROCESSED_BLAST_FILE, "<processed_blast.txt") or die "(BLAST) processed_blast.txt cannot be opened.\n";
    while(<PROCESSED_BLAST_FILE>){
        ($query, $subject, $extract_start, $extract_length, $orient) = split;
      
        $seq_extract = substr($genome{$subject}, $extract_start, $extract_length);
			
		# sequences in reverse orientatation reverse complemented to fit 
		#    the original query sequence
		if($orient == -1){
			$seq_extract =revTrans($seq_extract);
				
		}
        
		$to_align_list{$query}++;
        open (TE_OUT, ">>$query.fas") or die "$id.fas cannot be opened to insert new sequence to file.\n";
        print TE_OUT ">$subject-$extract_start-$extract_length($orient)\n$seq_extract\n";
    }

    my %genome;

#step 4 align if desired-----------------------------------
if( $align eq "1"){
    foreach $te_unaligned_file (keys %to_align_list){
        $te_unaligned_file =~ s/.fas//gi;
        system("/lustre/work/daray/software/muscle/muscle -in $te_unaligned_file.fas -out $te_unaligned_file.muscle.fas");
    }
}

################################################################################
sub revTrans
{
	local($sequence);
	
	$sequence = $_[0];	
	$sequence = reverse $sequence;
	$sequence =~ tr/ACGTacgt/TGCAtgca/;
	return $sequence;
}


################################################################################
sub create_te_out_files
{
    local($TE_library);
    $TE_library         = $_[0];



    my $inseq = Bio::SeqIO->new(-file   => "<$TE_library", -format =>"fasta");
    while (my $seq = $inseq->next_seq) {
        my $id = $seq->display_id;
	    my $seq = $seq->seq;
        open (TE_OUT, ">$id.fas") or die "$id.fas cannot be opened.\n";
        print TE_OUT ">CONSENSUS-$id\n$seq\n";
	}
}

################################################################################
sub organize_blast_hits
{
    local($blast_file, $buffer, $seqs_to_extract);
    $blast_file         = $_[0];
    $buffer             = $_[1];
    $seqs_to_extract    = $_[2];

    #open blastfile
    open (BLASTFILE, "<$blast_file") or die "blast file $blast_file cannot be opened.\n";


    my $element=0;
    while(<BLASTFILE>) {
	    ($query,        $subject,       $percent_id, 
         $align_length, $mismatches,    $gap_openings, 
         $query_start,  $query_end,     $subj_start, 
         $subj_end,     $e_value,       $bit_score    ) = split;
        
        #set orientation to 1    
        $orient = 1;
        
        #if hit on rev strand reorient and set to -1
        if( $subj_end < $subj_start ){
            ($subj_end, $subj_start) = ($subj_start, $subj_end);
            $orient = -1;
        }

        #set extract start with buffer (set to zero if extend beyond contig)
        $extract_start = $subj_start - $buffer;
        if($extract_start < 0){$extract_start=0;}
        
        #set extract length by subtracting end (with buffer from buff start) 
        $extract_length = (($subj_end + $buffer) - $extract_start);
        

        #create array of arrays with the process blast input data
        push(@processed_blast, [$query, $subject, $extract_start, $extract_length, $orient, $e_value, $bit_score]);
        ##############print "$processed_blast[$element][1]\n";

        $element++;

    }
 
    #sort all processes blast hits by evalue (small to large) then
    # sort by bit score large to small.
    my @evalue_sorted           = sort { $a->[5] <=> $b->[5] } @processed_blast;
    my @bit_and_evalue_sorted   = sort { $b->[6] <=> $a->[6] } @evalue_sorted;



    #reset counter
    my $iter;
    open (PROCESSED_BLAST_FILE, ">processed_blast.txt") or die "(processed blastFile <$blast_file> cannot be opened.\n";
    #cycle through sorted blast output
    while($iter <= $#bit_and_evalue_sorted){
    
        #if hit is one of the best 40 (or first) from sorted blast output -> print out
        if(($hits{$bit_and_evalue_sorted[$iter][0]}+1) <= $seqs_to_extract){
            print PROCESSED_BLAST_FILE "$bit_and_evalue_sorted[$iter][0]\t$bit_and_evalue_sorted[$iter][1]\t$bit_and_evalue_sorted[$iter][2]\t$bit_and_evalue_sorted[$iter][3]\t$bit_and_evalue_sorted[$iter][4]\n";
        }

    #increment hits and iter
    $hits{$bit_and_evalue_sorted[$iter][0]}++;
    $iter++;
    }
    
    close PROCESSED_BLAST_FILE;
}



################################################################################
sub usage
{
    print "\n";
	print " #################################################################\n";
	print " #                                                               #\n";	
	print " #  extract.align.pl - a program to extract TEs from a genome    #\n";
	print " #                                                               #\n";
	print " #  Required values:                                             #\n";  
	print " #    --genome    input genome file (in fasta format)            #\n";
	print " #    --blast     input blast file (tab format [6])              #\n"; 
	print " #    --consTEs   file of consensus elements                     #\n";
	print " #    --seqBuffer additional 5'/3' bases in extracted sequence   #\n";
	print " #    --seqNum    number of sequences to extract                 #\n";
	print " #                                                               #\n";	
	print " #  Optional values:                                             #\n";
	print " #    --align     aligns the sequences with MUSCLE               #\n";
	print " #                                                               #\n";
	print " # Questions: neal.platt at gmail.com		  	     	 #\n";
	print " #################################################################\n";
	
	exit;
}	
