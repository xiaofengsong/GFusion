#!/usr/bin/perl
# This programm find fusion genes;
use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case pass_through);
use Bio::SeqIO::fastq;
use Bio::DB::Fasta;

#//-----------------------------------//
# USAGE
#					GFusion.pl 1.0.0
#	GFusion is a software package to detect fusion genes using RNA-Seq data.
my $USAGE =<<USAGE;
	
	[example]
	perl Gfusion.pl -o out_file -r 0 -p 6 -i /path/to/bowtie1_index -g /path/to/genes.gtf -1 test_1.fastq -2 test_2.fastq
	

	[options]
	-o		output directory (default: output);
	-1		pair-end_fastq_1
	-2		pair-end_fastq_2
	-f		single-end_fastq
	-r		<int> inner distance between mate pairs (default = 0);
	-p		<int> threads (default = 4);
	-n		<int> min split num (default = 1);
	-L		<int> the length of anchors(more than 20, default = 0.4*read_length);
	-d		<int> min distance between adjacent genes (default = 50000); 
	-i		bowtie_index
	-g		gene annotation file(GTF with known transcripts))
	--help		usage
	

USAGE
#contact: xfsong@nuaa.edu.cn
#//-----------------------------------//


	my ( $bowtie_index, $fastq_1, $fastq_2, $opt_gtf, $opt_help, $fastq);
	my $out_file = "output";
	my $split_num = 1;
	my $thread = 4;
	my $r=0;
	my $pos_length=15;
	my $intron_length = 50000;
	my $fastq_in;
	my $pos_num=1;
	my $pair_num=2;
	my %hash_S;
	GetOptions(
		"out_file|o=s"=>\$out_file,
		"1=s"=>\$fastq_1,
		"2=s"=>\$fastq_2,
		"r:i"=>\$r,
		"n:i"=>\$split_num,
		"p:i"=>\$thread,
		'i=s' => \$bowtie_index,
		'g=s' => \$opt_gtf,
		'help!' => \$opt_help,
		"L:i"=>\$pos_length,
		"d:i"=>\$intron_length,
		"f=s"=>\$fastq,
		
 	);
	if(@ARGV>0) {
		print "\nerror unknown option: ","@ARGV\n","\n"; 
		exit 0;
	} 
	
	
	my $fusion_out = $out_file."/fusion_out";
	my $ref = $fusion_out."/ref";
	my $ref_index = $ref."/index";
	
	if ($opt_help) { print "$USAGE\n"; exit 0; }	
	
	if (!$bowtie_index ) { print "\n"," warning: missing parameters -i!","\n","$USAGE\n"; exit 0; }

	
	if (!$opt_gtf ) { print "\n"," warning: missing parameters -g!","\n","$USAGE\n"; exit 0; }

	
	

	if($fastq_1) {
		$fastq_in = Bio::SeqIO->new(-format => 'fastq', -file => "$fastq_1");
	} elsif ($fastq) {
		$fastq_in = Bio::SeqIO->new(-format => 'fastq', -file => "$fastq");
	}	
	
		my $fastq_obj=$fastq_in->next_seq();
		my $read_length = $fastq_obj->length;
		$pos_length = int(0.4*$read_length) if $pos_length<=20;
		$pos_length = 20 if $pos_length<=20;
		my $segment_length = $read_length-2*$pos_length;
		if($segment_length <10 or $pos_length<0) { print "error about anchor length, we suggest the best anchor length is 0.4*read_length!","\n","$USAGE\n"; exit 0; }
	

		system("mkdir $out_file") unless -e $out_file;
		system("mkdir $fusion_out") unless -e $fusion_out;
		system("mkdir $ref") unless -e $ref;
		system("mkdir $ref_index") unless -e $ref_index;
		my $start_when = time();
		my $when;
		$when = localtime();
		print STDERR "[",$when,"]","\n";

	if($fastq_1 and $fastq_2 and !$fastq) {

		system("tophat -o $out_file --bowtie1 -p $thread -r $r -I100000 --no-coverage-search $bowtie_index $fastq_1 $fastq_2");
		system("samtools view -q 30 -F4 -f8 $out_file/accepted_hits.bam > $out_file/unmate.sam");
		system("samtools view $out_file/unmapped.bam  >$out_file/unmapped.sam");
		&mate("$out_file/unmate.sam","$out_file/unmapped.sam");   #>$out_file/un_mate.sam
		unlink "$out_file/unmate.sam" if -e "$out_file/unmate.sam";
		$when = localtime();
		print STDERR "[",$when,"]","\n";
		
		#split ends
		system("bowtie -p $thread $bowtie_index $out_file/un.fastq -S $fusion_out/un.sam");

		system("samtools view -F4 -b -S $fusion_out/un.sam > $fusion_out/un.bam");
		system("samtools sort -n $fusion_out/un.bam -o $fusion_out/un_S");
		system("samtools view $fusion_out/un_S.bam > $fusion_out/un_S.sam");
		system("samtools view -H $fusion_out/un_S.bam > $fusion_out/pos_header.sam");
		&both_filter;	

	
		&sortAll("$fusion_out/un_both.sam","$fusion_out/pos_S.sam");
		&locating("$fusion_out/pos_S.sam","$fusion_out/pos_locale.sam","$fusion_out/pos_header.sam");
	
		system("samtools view -b -S $fusion_out/pos_locale.sam > $fusion_out/pos_locale.bam");
		system("samtools sort -n $fusion_out/pos_locale.bam -o $fusion_out/pos_locale_S");
		system("samtools view $fusion_out/pos_locale_S.bam > $fusion_out/pos_locale_S.sam");

		&pos_filter;  #>choice pos_both.sam

		#choose mate
		&other_mate("$fusion_out/pos_both.sam","$fusion_out/pos_other.sam");
		&sortAll("$fusion_out/pos_other.sam","$fusion_out/pos_other_S.sam");	
		&locating("$fusion_out/pos_other_S.sam","$fusion_out/pos_other_locale.sam");

				
		&split_find("$fusion_out/pos_both.sam","$fusion_out/pos_other_locale.sam","$fusion_out/pos_gene.sam");
		
		#paired-end
		system("samtools view -H $out_file/accepted_hits.bam > $fusion_out/pair_header.sam");
		system("samtools view -q30 -F12 $out_file/accepted_hits.bam > $fusion_out/pair_map.sam");
		&pair_filter;  #paired-end filter
		&sortAll("$fusion_out/pair_filter.sam","$fusion_out/pair_S.sam");
		&locating("$fusion_out/pair_S.sam","$fusion_out/pair_locale.sam","$fusion_out/pair_header.sam");

		system("samtools view -b -S $fusion_out/pair_locale.sam > $fusion_out/pair_locale.bam");
		system("samtools sort -n $fusion_out/pair_locale.bam -o $fusion_out/pair_locale_S");
		system("samtools view $fusion_out/pair_locale_S.bam > $fusion_out/pair_locale_S.sam");



		&span_find;  #spanning filter
		&no_result;
		&ref_Index;
		
		&srr_fq($fastq_1,$fastq_2);

		
		system("tophat -o $fusion_out/final -p $thread --bowtie1 -r $r $ref_index/re $ref/er_1.fq $ref/er_2.fq");

		system("samtools view -F12 $fusion_out/final/accepted_hits.bam > $fusion_out/final/map.sam");

		&result("$out_file/temp","$fusion_out/final/map.sam","$out_file/fusion_gene_list.txt");
		
=pod		#total fastq
	#	system("tophat -o $fusion_out/total -p $thread --bowtie1 -r $r $ref_index/re $fastq_1 $fastq_2");

		system("samtools view -q30 -F12 $fusion_out/total/accepted_hits.bam > $fusion_out/total/map.sam");
		&result("$out_file/temp","$fusion_out/total/map.sam","$out_file/total_result");
=cut
		$when = localtime();
		print STDERR "[",$when,"]";
		print STDERR " Completed successfully!  The time elapsed: about ",0.01*int ( (time()-$start_when)/36 )," hours.\n";			


	} elsif(!$fastq_1 and !$fastq_1 and $fastq) {

		&toUn("$fastq","$out_file/un_1.fastq","$out_file/un_2.fastq");	
		
		system("tophat -o $fusion_out --bowtie1 -p $thread -r $segment_length --no-coverage-search $bowtie_index $out_file/un_1.fastq $out_file/un_2.fastq");
		
		system("samtools view -bh -q30 -F12 $fusion_out/accepted_hits.bam > $fusion_out/top_accQfi.bam");	
		system("samtools sort -n $fusion_out/top_accQfi.bam -o $fusion_out/un_S");
		system("samtools view $fusion_out/un_S.bam > $fusion_out/un_S.sam");
		system("samtools view -H $fusion_out/un_S.bam > $fusion_out/pos_header.sam");

		&both_filter;	


		&sortAll("$fusion_out/un_both.sam","$fusion_out/pos_S.sam");
		&locating("$fusion_out/pos_S.sam","$fusion_out/pos_locale.sam","$fusion_out/pos_header.sam");

		
		system("samtools view -b -S $fusion_out/pos_locale.sam > $fusion_out/pos_locale.bam");
		system("samtools sort -n $fusion_out/pos_locale.bam -o $fusion_out/pos_locale_S");
		system("samtools view $fusion_out/pos_locale_S.bam > $fusion_out/pos_locale_S.sam");

		&pos_filter;  #>choice pos_both.sam		
		&seed_filter("$fusion_out/pos_both.sam", "$fusion_out/pos_gene.sam");
		
		&seed_find;  #split reads
		&no_result;
		&ref_Index;
		&srr_fq($fastq);

		system("tophat -o $fusion_out/final -p $thread --bowtie1 -r $r $ref_index/re $ref/er_1.fq");

		system("samtools view -F4 $fusion_out/final/accepted_hits.bam > $fusion_out/final/map.sam");
		&single_result("$out_file/temp","$fusion_out/final/map.sam","$out_file/fusion_gene_list.txt");
		
		$when = localtime();
		print STDERR "[",$when,"]";
		print STDERR " Completed successfully!  The time elapsed: about ",0.01*int ( (time()-$start_when)/36 )," hours.\n";			
	
	} else {
		 print "\n"," warning: fastq file error!","\n","$USAGE\n"; exit 0;
	}


sub toUn {
	open(IN,"$_[0]") || die "Can't open $_[0]:$!\n";	
	open(UN_1,">$_[1]") || die "Can't open $_[1]:$!\n";
	open(UN_2,">$_[2]") || die "Can't open $_[2]:$!\n";
	my $long = $read_length-$pos_length;
	my ($line, @line, $srr_line, $qu_line);
	while($line =<IN>) {
		@line = split('\s+',$line); 
		$srr_line =<IN>;
		<IN>;
		$qu_line =<IN>;
		print UN_1 $line[0],"\n";			
		print UN_1 substr($srr_line, 0, $pos_length),"\n";
		print UN_1 "+","\n";
		print UN_1 substr($qu_line, 0, $pos_length),"\n";
		print UN_2 $line[0],"\n";	
		print UN_2 substr($srr_line, $long, $pos_length),"\n";
		print UN_2 "+","\n";
		print UN_2 substr($qu_line, $long, $pos_length),"\n";

	}
	#print "123","\n";
	close(IN);
	close(UN_1);	close(UN_2);
}





sub mate {
	open(UNMATE,"$_[0]") || die "Can't open $_[0]:$!\n";	
	open(UNMAP,"$_[1]") || die "Can't open $_[1]:$!\n";
	open(UN,">$out_file/un.fastq") || die "Can't open $out_file/un.fastq:$!\n";
	open(OUT_M,">$out_file/un_mate.sam") || die "Can't open $out_file/un_mate.sam:$!\n";
	my $long = $read_length-$pos_length;
	my $unmateLine;
	my @unmateArr;
	my %mateHash;
	my $unmapLine;
	my @unmapArr;
	while($unmateLine = <UNMATE>) {
		@unmateArr = split(/\s+/,$unmateLine);	
		$mateHash{$unmateArr[0]} = $unmateLine unless exists $mateHash{$unmateArr[0]};
		
	}
	while($unmapLine = <UNMAP>) {
		@unmapArr = split(/\s+/,$unmapLine);
		if(exists $mateHash{$unmapArr[0]}) {
			my $srr_line = $unmapArr[9];
			my $srr_line2 = $unmapArr[10];
			print UN "@".$unmapArr[0]."/1","\n";			
			print UN substr($srr_line, 0, $pos_length),"\n";
			print UN "+","\n";
			print UN substr($srr_line2, 0, $pos_length),"\n";
			print UN "@".$unmapArr[0]."/2","\n";	
			print UN substr($srr_line, $long, $pos_length),"\n";
			print UN "+","\n";
			print UN substr($srr_line2, $long, $pos_length),"\n";
			print OUT_M $mateHash{$unmapArr[0]};
						
			delete $mateHash{$unmapArr[0]};
		}
	
	}
		
	undef %mateHash; close(UNMATE); close(UNMAP); close(UN); close(OUT_M);
	


}


sub both_filter {
	open(POS,"$fusion_out/un_S.sam");
	open(BOTH,">$fusion_out/un_both.sam");
	my $temp = <POS>;
	my @temp;
	my $line;
	my @line;
	my $line_lo;

	while($line=<POS>) {
		@temp = split(/\s+/,$temp);
		$temp[0] =~ s/\/(\d)$//;
	
		@line= split(/\s+/,$line);
		$line[0] =~ s/\/(\d)$//;
		if($temp[0] eq $line[0]) {
			if(($temp[2] ne $line[2]) || (abs($line[3]-$temp[3]) >= $intron_length)){
			

				print BOTH $temp;
				print BOTH $line;
				$temp = <POS>;	
			
			} else {
								
				$temp = $line;		
			}	
			
		} else {	
			$temp = $line;
		}
	}
	close(POS);
	close(BOTH);
}	

sub pair_filter {
	open(PAIR,"$fusion_out/pair_map.sam");
	open(FILTER,">$fusion_out/pair_filter.sam");
	my $line;
	my @line;
	while($line=<PAIR>) {
		@line= split(/\s+/,$line);
		unless(($line[6] eq "=") and (abs($line[8]) <= $intron_length)) {  #$lone[8]
			print FILTER $line;
		}
	}
	close(PAIR);
	close(FILTER);
}	

sub sortAll {
	open(IN,"$_[0]") || die "Can't open $_[0]:$!\n";	
	open(OUT,">$_[1]") || die "Can't open $_[1]:$!\n";
	&decide("$_[0]","$fusion_out/sort1.sam","$fusion_out/sort2.sam","$fusion_out/sort3.sam");
	&sorting("$fusion_out/sort1.sam");
	&sorting("$fusion_out/sort2.sam");
	&sorting("$fusion_out/sort3.sam");
	unlink "$fusion_out/sort1.sam" if -e "$fusion_out/sort1.sam";
	unlink "$fusion_out/sort2.sam" if -e "$fusion_out/sort2.sam";
	unlink "$fusion_out/sort3.sam" if -e "$fusion_out/sort3.sam";

	close(IN);
	close(OUT);

}
sub decide
{
	open(IN,"$_[0]");
	open(OUT1,">$_[1]");
	open(OUT2,">$_[2]");
	open(OUT3,">$_[3]");
	my @arrays ;
	while(my $line = <IN>) { 
	@arrays = split(/\s/,$line);
	if('chr15' ge $arrays[2]) {
		print OUT1 $line;	
	} elsif ('chr3' ge $arrays[2]) {
		print OUT2 $line;	
	} else {
		print OUT3 $line;	
	}
}

		
	close(IN);
	close(OUT1);
	close(OUT2);
	close(OUT3);
	
}

sub sorting
{
	
	my $key;
	my $count = 0;
	my @arrays;
	open(IN,"$_[0]")|| die "Can't open $_[0]:$!\n";	
	while(my $line = <IN>) { 
		@arrays = split(/\s/,$line);
		$hash_S{$count}{'line'} = $line;
		$hash_S{$count}{'chr'} = $arrays[2];
		$hash_S{$count}{'local'} = $arrays[3];
		#print $hash{$count}{'line'};
		$count++;
	}

	foreach $key(sort sort_chome(keys(%hash_S))) {
		print OUT $hash_S{$key}{'line'};
		
	}

	sub sort_chome
	{
		$hash_S{$a}{'chr'} cmp $hash_S{$b}{'chr'} or $hash_S{$a}{'local'}  <=> $hash_S{$b}{'local'};
	}

	undef %hash_S;	
	close(IN);
	
	
}
sub locating {
	my $head;	
	open (GTF,"$opt_gtf")|| die "Can't open $_[1]:$!\n";
	open (IN, "$_[0]")|| die "Can't open $_[0]:$!\n";
	open (OUT,">$_[1]")|| die "Can't open $_[1]:$!\n";
	if (defined $_[2]) {
		open (HEAD,"$_[2]");
		while($head = <HEAD>) {
			print OUT $head;	
		}	
	}

	my $count;
	my $sam_line;
	my $gtf_line;
	my @sam_array;
	my @gtf_array;
	my $bit;
	my @flag;
	my $flag_length;
	my $strand;
	my $last_sam_chr = "chr1";
	my $index;
	my $seg;
	$count = 0;

	$gtf_line = <GTF>;
	@gtf_array = split(/\s/, $gtf_line);
		
		
	
	while ( $sam_line = <IN> )
	{	
		
		chomp $sam_line;
		@sam_array = split(/\s/, $sam_line);
		if ( $sam_array[2] eq "chrM" ) 	{ next; }
		
		LABEL:
		
		if ( $sam_array[2] eq $gtf_array[0] )
		{
			if ( $gtf_array[2] eq "exon" ) 
			{
				if ( ( $sam_array[3] >= $gtf_array[3] - 4 ) and ($sam_array[3] < $gtf_array[4] + 4)  )
				{
					$count += 1;
					#the strand;below is not important;
					$bit = reverse ( sprintf("%b",$sam_array[1])+0 ); #decimal to binary; 
					@flag = split ("", $bit);
					$flag_length = @flag;
					if ($flag_length<6) {
						$flag[6] = 0;
					}
					if ($flag_length<4) {
						$flag[4] = 0;
					}
					if ( $flag[4] == 1) { $strand = "-"; }# the strand is forward;
					else {$strand = "+";}
					if($flag[6] == 1) {$seg = 1;} 
					else {$seg = 2;}
					
					#the strand;above is not important;
					$gtf_array[11] =~ s/^\"|\";$//g;
					if ( $gtf_array[12] eq "transcript_id")
					{
						$index = 13;
					}
					else
					{
						$index = 15;
					}
					$gtf_array[$index] =~ s/^\"|\";$//g;
					print OUT $sam_line,"\t","XC:Z:",join(" ",$gtf_array[0], join("","coor:",$gtf_array[3], "*", $gtf_array[4]), join("","gene_id:",$gtf_array[11]), join("","read",$strand), join("","strand",$gtf_array[6]), join("","seg",$seg),$gtf_array[$index] ),"\n";	
					
				}
				else
				{
					if ( $sam_array[3] > $gtf_array[4] )
					{
						$gtf_line = <GTF>;
			
						if ( defined $gtf_line ) 
						{
							@gtf_array = split(/\s+/, $gtf_line);
							goto LABEL;
						}
						else { last; }
					}
				}
			}
			else
			{
				$gtf_line = <GTF>;
				if ( defined $gtf_line )
				{
					@gtf_array = split(/\s+/, $gtf_line);
					goto LABEL;
				}
				else { last; }
			}
		}
		else
		{
			if ( $sam_array[2] ne $last_sam_chr )
			{
				$gtf_line = <GTF>;
				if ( defined $gtf_line ) 
				{
					@gtf_array = split(/\s+/, $gtf_line);
					goto LABEL;
				}
				else {last;}
			}
		}
		$last_sam_chr = $sam_array[2];
	}
	close(IN);
	close(GTF);
	close(HEAD);
	close(OUT);
}


sub pos_filter {
	open(POS,"$fusion_out/pos_locale_S.sam");
	open(BOTH,">$fusion_out/pos_both.sam");
	my $temp = <POS>;
	my @temp;
	my $temp_flag;
	my $temp_lo;
	my $temp_gene;
	my $line;
	my @line;
	my $line_lo;
	my $line_gene;
	my $line_flag;
	my $strand_1 = -4;
	my $strand_2 = $segment_length+4;
	my $strand_3 = $pos_length-4;
	my $strand_4 = $read_length-$pos_length+4;
	my $coor_1;
	my $coor_2;

	while($line=<POS>) {
		@temp = split(/\s+/,$temp);
		if($temp[0] =~ s/\/(\d)$//) { $temp_lo = $1;}
		elsif($temp =~ /\sseg(\d)\s/) { $temp_lo = $1;}
		
		@line= split(/\s+/,$line);
		if($line[0] =~ s/\/(\d)$//) {$line_lo = $1;}
		elsif($line =~ /\sseg(\d)\s/) { $line_lo = $1;}
		if($temp[0] eq $line[0]) {
			
			$temp_gene = $1 if $temp =~ /\sgene_id:(\S+)\s/;
			$line_gene = $1 if $line =~ /\sgene_id:(\S+)\s/;
			if($temp_gene ne $line_gene) {
				
				$temp_flag = 0;
				$line_flag = 0;			
				$temp =~ /\scoor:(\d+)\*(\d+)\s/;
				$coor_1 = $1;
				$coor_2 = $2;
				if($temp_lo == 1) {
					
					if(($temp =~ /\sread-\s/) and (($temp[3]-$coor_1<=$strand_2)and($temp[3]-$coor_1>=$strand_1))) {
						$temp_flag = 1;	
						
					} elsif(($temp =~ /\sread\+\s/) and (($coor_2-$temp[3]<=$strand_4 )and($coor_2-$temp[3]>=$strand_3))) {
						$temp_flag = 1;	
					}
				} elsif($temp_lo == 2) {
					if(($line =~ /\sread-\s/) and (($coor_2-$temp[3]<=$strand_4 )and($coor_2-$temp[3]>=$strand_3))) {
						$temp_flag = 1;		
					} elsif(($temp =~ /\sread\+\s/) and (($temp[3]-$coor_1<=$strand_2)and($temp[3]-$coor_1>=$strand_1))) {
						$temp_flag = 1;	
					}
				}

				$line=~ /\scoor:(\d+)\*(\d+)\s/;
				$coor_1 = $1;
				$coor_2 = $2;
				if($line_lo == 1) {
					if(($line =~ /\sread-\s/) and (($line[3]-$coor_1<=$strand_2)and($line[3]-$coor_1>=$strand_1))) {
						$line_flag = 1;	
						
					} elsif(($line =~ /\sread\+\s/) and (($coor_2-$line[3]<=$strand_4 )and($coor_2-$line[3]>=$strand_3))) {
						$line_flag = 1;	
					}
				} elsif($line_lo == 2) {
					if(($line =~ /\sread-\s/) and (($coor_2-$line[3]<=$strand_4 )and($coor_2-$line[3]>=$strand_3))) {
						$line_flag = 1;		
					} elsif(($line =~ /\sread\+\s/) and (($line[3]-$coor_1<=$strand_2)and($line[3]-$coor_1>=$strand_1))) {
						$line_flag = 1;	
					}
				}
				if(($temp_flag == 1) and ($line_flag == 1)) {
					
					print BOTH $temp;
					print BOTH $line;
					$temp = <POS>;	
				}
			}	
		}		
			
		$temp = $line;
	}
	close(POS);
	close(BOTH);
}	

sub other_mate {
	open(UN,"$out_file/un_mate.sam") || die "Can't open $out_file/unmate_1.sam:$!\n";
	open(BOTH,"$_[0]") || die "Can't open BOTH:$!\n";
	open(OUT,">$_[1]") || die "Can't open SRR:$!\n";
	my %hash;
	my $count=0;
	my $num=0;
	#print OUT $open,"\n";
	while(my $line = <BOTH>) {
		my @line = split('\t',$line);	
		if ($line[0] =~ s/\/1$//) {
			$hash{$line[0]}=0;
			$num++;
		}
		
	}
	while(my $line = <UN>) {
		my @line = split('\t',$line);
		if(exists $hash{$line[0]}) {
			print OUT $line;
			$count++;
		}	

	}
	undef %hash; close(UN); close(BOTH); close(OUT);

}

sub seed_filter {
	my (%hash, @arrays, $line, $key, $r_2, $r_2_chr, $r_2_gene, $r_2_cor1, 
      		$r_2_cor2, $r_2_strand, $r_2_read, $r_1, $r_1_chr, $r_1_gene, $r_1_cor1, $r_1_cor2, $r_1_strand, $r_1_read);
	
	#print OUT $open,"\n";
	open (IN, "$_[0]")|| die "Can't open $_[0]:$!\n";
	open (OUT_GENE,">$_[1]")|| die "Can't open $_[1]:$!\n";

	while($line = <IN>) {
		@arrays = split(/\s+/,$line);
		$key = $arrays[0];
		if(exists $hash{$key}) {
			if ($line=~ /seg1/) {
				$r_2 = $hash{$key};
				$r_1 = $line;
			} else {
				$r_1 = $hash{$key};
				$r_2 = $line;
			}
			$r_2 =~ /(chr\w+)\scoor:(\d+)\*(\d+)\sgene_id:(\S+)\s(read(\+|-))\s(strand(\+|-))\s/i;
			$r_2_chr = $1;
			$r_2_cor1 = $2;
			$r_2_cor2 = $3;
			$r_2_gene = $4;
			$r_2_read = $5;
			$r_2_strand = $7;

			$r_1 =~ /(chr\w+)\scoor:(\d+)\*(\d+)\sgene_id:(\S+)\s(read(\+|-))\s(strand(\+|-))\s/i;
			$r_1_chr = $1;
			$r_1_cor1 = $2;
			$r_1_cor2 = $3;
			$r_1_gene = $4;
			$r_1_read = $5;
			$r_1_strand = $7;
			if ($r_1_gene ne $r_2_gene) {
				if($r_1_read eq 'read-' and $r_2_read eq 'read-' and $r_1_strand eq 'strand+' and $r_2_strand eq 'strand+') {
					print OUT_GENE 	join("\t",$r_2_gene, $r_2_chr, $r_2_cor2, $r_1_gene, $r_1_chr, $r_1_cor1, "ff", $key, "\n");		
				
				} elsif ($r_1_read eq 'read+' and $r_2_read eq 'read+' and $r_1_strand eq 'strand+' and $r_2_strand eq 'strand+') {
					print OUT_GENE 	join("\t",$r_1_gene, $r_1_chr, $r_1_cor2, $r_2_gene,$r_2_chr,$r_2_cor1,"ff",$key,"\n");		
				
				} elsif ($r_1_read eq 'read+' and $r_2_read eq 'read+' and $r_1_strand eq 'strand-' and $r_2_strand eq 'strand-') {
					print OUT_GENE 	join("\t",$r_2_gene,$r_2_chr,$r_2_cor1,$r_1_gene, $r_1_chr, $r_1_cor2, "rr",$key,"\n");		
				
				} elsif ($r_1_read eq 'read-' and $r_2_read eq 'read-' and $r_1_strand eq 'strand-' and $r_2_strand eq 'strand-') {
					print OUT_GENE 	join("\t",$r_1_gene, $r_1_chr, $r_1_cor1, $r_2_gene,$r_2_chr,$r_2_cor2,"rr",$key,"\n");		
				
				} elsif ($r_1_read eq 'read+' and $r_2_read eq 'read-' and $r_1_strand eq 'strand+' and $r_2_strand eq 'strand-') {
					print OUT_GENE 	join("\t",$r_1_gene, $r_1_chr, $r_1_cor2, $r_2_gene, $r_2_chr, $r_2_cor2,  "fr",$key,"\n");		
				
				} elsif ($r_1_read eq 'read+' and $r_2_read eq 'read-' and $r_1_strand eq 'strand-' and $r_2_strand eq 'strand+') {
					print OUT_GENE 	join("\t",$r_2_gene, $r_2_chr, $r_2_cor2, $r_1_gene, $r_1_chr, $r_1_cor2, "fr",$key,"\n");		
				
				} elsif ($r_1_read eq 'read-' and $r_2_read eq 'read+' and $r_1_strand eq 'strand+' and $r_2_strand eq 'strand-') {
					print OUT_GENE 	join("\t",$r_2_gene, $r_2_chr, $r_2_cor1, $r_1_gene, $r_1_chr, $r_1_cor1, "rf",$key,"\n");
					
				} elsif ($r_1_read eq 'read-' and $r_2_read eq 'read+'and $r_1_strand eq 'strand-' and $r_2_strand eq 'strand+') {
					print OUT_GENE 	join("\t",$r_1_gene, $r_1_chr, $r_1_cor1, $r_2_gene, $r_2_chr, $r_2_cor1, "rf",$key,"\n");
					
				}
			}

  
		} else {
			$hash{$key}=$line;			
		}
		
		
	}
	close(IN); close(OUT_GENE); undef %hash;
}









sub split_find {
	my (%hash, @arrays, $line, $key, $key2, $key1, $num);
	my ($mate_gene, $mate_strand, $r_2, $r_2_chr, $r_2_gene, $r_2_cor1, $r_2_cor2, $r_2_strand, 
		$r_2_read, $r_1, $r_1_chr, $r_1_gene, $r_1_cor1, $r_1_cor2, $r_1_strand, $r_1_read);
	
	#print OUT $open,"\n";
	open (IN, "$_[0]")|| die "Can't open $_[0]:$!\n";
	open (MATE,"$_[1]")|| die "Can't open $_[1]:$!\n";
	open (OUT_GENE,">$_[2]")|| die "Can't open $_[2]:$!\n";

	while($line = <IN>) {
		@arrays = split(/\s+/,$line);
		$hash{$arrays[0]}=$line;
	}
		
	while($line = <MATE>) {
		$num++;
		@arrays = split('\s',$line);
		$key = $arrays[0];
		$key1 = $key."/1";
		$key2 = $key."/2";
		$line =~ /\sgene_id:(\S+)\s(read(\+|-))\s/i;
		$mate_gene = $1;
		$mate_strand = $2;
		#print $mate_id,"	",$mate_strand,"\n";
		if(exists $hash{$key2} and exists $hash{$key1}) {	
			$r_2 = $hash{$key2};
			$r_2 =~ /(chr\w+)\scoor:(\d+)\*(\d+)\sgene_id:(\S+)\s(read(\+|-))\s(strand(\+|-))\s/i;
			$r_2_chr = $1;
			$r_2_cor1 = $2;
			$r_2_cor2 = $3;
			$r_2_gene = $4;
			$r_2_read = $5;
			$r_2_strand = $7;
			$r_1 = $hash{$key1};
			$r_1 =~ /(chr\w+)\scoor:(\d+)\*(\d+)\sgene_id:(\S+)\s(read(\+|-))\s(strand(\+|-))\s/i;
			$r_1_chr = $1;
			$r_1_cor1 = $2;
			$r_1_cor2 = $3;
			$r_1_gene = $4;
			$r_1_read = $5;
			$r_1_strand = $7;
			if ($r_2_gene eq $mate_gene and $r_1_gene ne $r_2_gene) {

				if($r_1_read eq 'read-' and $r_2_read eq 'read-' and $mate_strand eq 'read+' and $r_1_strand eq 'strand+' and $r_2_strand eq 'strand+') {
					print OUT_GENE 	join("\t",$r_2_gene, $r_2_chr, $r_2_cor2, $r_1_gene, $r_1_chr, $r_1_cor1, "ff", $key, "\n");		
				
				} elsif ($r_1_read eq 'read+' and $r_2_read eq 'read+' and $mate_strand eq 'read-'and $r_1_strand eq 'strand+' and $r_2_strand eq 'strand+') {
					print OUT_GENE 	join("\t",$r_1_gene, $r_1_chr, $r_1_cor2, $r_2_gene,$r_2_chr,$r_2_cor1,"ff",$key,"\n");		
				
				} elsif ($r_1_read eq 'read+' and $r_2_read eq 'read+' and $mate_strand eq 'read-'and $r_1_strand eq 'strand-' and $r_2_strand eq 'strand-') {
					print OUT_GENE 	join("\t",$r_2_gene,$r_2_chr,$r_2_cor1,$r_1_gene, $r_1_chr, $r_1_cor2, "rr",$key,"\n");		
				
				} elsif ($r_1_read eq 'read-' and $r_2_read eq 'read-' and $mate_strand eq 'read+'and $r_1_strand eq 'strand-' and $r_2_strand eq 'strand-') {
					print OUT_GENE 	join("\t",$r_1_gene, $r_1_chr, $r_1_cor1, $r_2_gene,$r_2_chr,$r_2_cor2,"rr",$key,"\n");		
				
				} elsif ($r_1_read eq 'read+' and $r_2_read eq 'read-' and $mate_strand eq 'read+'and $r_1_strand eq 'strand+' and $r_2_strand eq 'strand-') {
					print OUT_GENE 	join("\t",$r_1_gene, $r_1_chr, $r_1_cor2, $r_2_gene, $r_2_chr, $r_2_cor2,  "fr",$key,"\n");		
				
				} elsif ($r_1_read eq 'read+' and $r_2_read eq 'read-' and $mate_strand eq 'read+'and $r_1_strand eq 'strand-' and $r_2_strand eq 'strand+') {
					print OUT_GENE 	join("\t",$r_2_gene, $r_2_chr, $r_2_cor2, $r_1_gene, $r_1_chr, $r_1_cor2, "fr",$key,"\n");		
				
				} elsif ($r_1_read eq 'read-' and $r_2_read eq 'read+' and $mate_strand eq 'read-'and $r_1_strand eq 'strand+' and $r_2_strand eq 'strand-') {
					print OUT_GENE 	join("\t",$r_2_gene, $r_2_chr, $r_2_cor1, $r_1_gene, $r_1_chr, $r_1_cor1, "rf",$key,"\n");
					
				} elsif ($r_1_read eq 'read-' and $r_2_read eq 'read+' and $mate_strand eq 'read-'and $r_1_strand eq 'strand-' and $r_2_strand eq 'strand+') {
					print OUT_GENE 	join("\t",$r_1_gene, $r_1_chr, $r_1_cor1, $r_2_gene, $r_2_chr, $r_2_cor1, "rf",$key,"\n");
					
				}
			
			}

		} 
	}


	undef %hash; close(IN); close(MATE); close(OUT_GENE);
	
}

sub seed_find { 
	open(POS,"$fusion_out/pos_gene.sam");
	open(GTF,"$opt_gtf");
	open(TEMP,">$out_file/temp");
	open(SRR,">$fusion_out/srr");

	my ($pos, %split, %gene_H, $key, @pos, $gene, %srr, $srr, @temp, $line, @line, $line_lo, $line_gene,
		$temp_gene, $count, %gene_local, $gtf, @gtf);
	
	while($pos=<POS>) {
		chomp($pos);
		@pos = split(/\s+/,$pos);
		$pos =~ s/\s+($pos[7])//;
		print SRR $1,"\n";
		if(exists $split{$pos}) {   
			$split{$pos}++;		
		} else {
			$split{$pos} = 1; 
			$gene = $pos[0]."\t".$pos[3];		
			$gene_local{$pos[0]}{"start"} = 0 unless exists $gene_local{$pos[0]}{"start"};
			$gene_local{$pos[0]}{"stop"} = 0 unless exists $gene_local{$pos[0]}{"stop"};
			$gene_local{$pos[3]}{"start"} = 0 unless exists $gene_local{$pos[3]}{"start"};
			$gene_local{$pos[3]}{"stop"} = 0 unless exists $gene_local{$pos[3]}{"stop"};
			
		}		
	
	}

	#pair spanning reads SRR 
		
			
	
	while($gtf = <GTF>) {
		$gtf =~ /gene_name\s+"(\S+)"/;
		$gene = $1;		
		if(exists $gene_local{$gene}) {
			if($gene_local{$gene}{"start"} == 0 and $gene_local{$gene}{"stop"} == 0) {
				@gtf = split(/\s+/, $gtf);				
				$gene_local{$gene}{"start"} = $gtf[3];
				$gene_local{$gene}{"stop"} = $gtf[4];	
			} else {
				@gtf = split(/\s+/, $gtf);				
				if($gene_local{$gene}{"start"} > $gtf[3]) {
					$gene_local{$gene}{"start"} = $gtf[3];				
				}
				if($gene_local{$gene}{"stop"} < $gtf[4]) {
					$gene_local{$gene}{"stop"} = $gtf[4];	
				}				
			}		
		}

	
	}
		

	foreach $key (sort(keys %split)) {
		@pos = split('\s+',$key);
		unless ($gene_local{$pos[3]}{"start"} < $gene_local{$pos[0]}{"stop"} and $gene_local{$pos[3]}{"stop"} > $gene_local{$pos[0]}{"start"}) {
			$count++;
			print TEMP $count,"\t",$key,"\t",$gene_local{$pos[0]}{"start"},"\t",$gene_local{$pos[0]}{"stop"},"\t",$gene_local{$pos[3]}{"start"},"\t",$gene_local{$pos[3]}{"stop"},"\t",$split{$key},"\n";	
			
		
		}	
				
	
	}
	undef %srr;undef %split;undef %gene_H;undef %gene_local;close(POS);close(PAIR);	close(TEMP);close(OUT);close(SRR);	

}

sub no_result {
	open(TEMP,"$out_file/temp");
	my $temp =<TEMP>;
	unless(defined $temp) {
		$when = localtime();
		print STDERR "[",$when,"]","\n";
		print STDERR " Result: No Fusion Genes!  The time elapsed: about ",0.01*int ( (time()-$start_when)/36 )," hours.\n";
		exit 0;
	}
	
	
}




sub span_find { 
	open(POS,"$fusion_out/pos_gene.sam");
	open(PAIR,"$fusion_out/pair_locale_S.sam");
	open(GTF,"$opt_gtf");
	open(TEMP,">$out_file/temp");
	open(SRR,">$fusion_out/srr");

	my ($pos, %split, %gene_H, $key, @pos, $gene, %srr, $srr, @temp, $line, @line, $line_lo, $line_gene,
		$temp_gene, $count, %gene_local, $gtf, @gtf);
	my $temp = <PAIR>;

	while($pos=<POS>) {
		chomp($pos);
		@pos = split(/\s+/,$pos);
		$pos =~ s/\s+($pos[7])//;
		$srr = $1;
		if(exists $split{$pos}) {   
			$split{$pos}++;	
			$srr{$pos} = $srr{$pos}.$srr."\n";	
		} else {
			$split{$pos} = 1;
			$srr{$pos} = $srr."\n"; 
			$gene = $pos[0]."\t".$pos[3];		
			$gene_H{$gene} = 0 unless exists $gene_H{$gene};
			
		}		
	
	}

	
	while($line=<PAIR>) {	
		@temp = split(/\s+/,$temp);
		@line= split(/\s+/,$line);
		if(($temp[0] eq $line[0]) and ($temp=~ /NH:i:1/) and ($line=~ /NH:i:1/)) {
			$temp_gene = $1 if $temp =~ /\sgene_id:(\S+)\s/;
			$line_gene = $1 if $line =~ /\sgene_id:(\S+)\s/;
			$gene = $temp_gene."\t".$line_gene;
			if(exists $gene_H{$gene}){
				$gene_H{$gene}++;
				print SRR $temp[0]."\n";	 
				$gene_local{$temp_gene}{"start"} = 0 unless exists $gene_local{$temp_gene}{"start"};
				$gene_local{$temp_gene}{"stop"} = 0 unless exists $gene_local{$temp_gene}{"stop"};
				$gene_local{$line_gene}{"start"} = 0 unless exists $gene_local{$line_gene}{"start"};
				$gene_local{$line_gene}{"stop"} = 0 unless exists $gene_local{$line_gene}{"stop"};
				$temp = <PAIR>;	
			
			} else {
				$gene = $line_gene."\t".$temp_gene;
				if(exists $gene_H{$gene}){
					$gene_H{$gene}++;
					print SRR $temp[0]."\n";    
					$gene_local{$temp_gene}{"start"} = 0 unless exists $gene_local{$temp_gene}{"start"};
					$gene_local{$temp_gene}{"stop"} = 0 unless exists $gene_local{$temp_gene}{"stop"};
					$gene_local{$line_gene}{"start"} = 0 unless exists $gene_local{$line_gene}{"start"};
					$gene_local{$line_gene}{"stop"} = 0 unless exists $gene_local{$line_gene}{"stop"};
					$temp = <PAIR>;	
				} else	{
					$temp = $line;				
				}
				 
				
			}
		} else {	
			$temp = $line;
		}			
			
	}
	while($gtf = <GTF>) {
		$gtf =~ /gene_name\s+"(\S+)"/;
		$gene = $1;		
		if(exists $gene_local{$gene}) {
			if($gene_local{$gene}{"start"} == 0 and $gene_local{$gene}{"stop"} == 0) {
				@gtf = split(/\s+/, $gtf);				
				$gene_local{$gene}{"start"} = $gtf[3];
				$gene_local{$gene}{"stop"} = $gtf[4];	
			} else {
				@gtf = split(/\s+/, $gtf);				
				if($gene_local{$gene}{"start"} > $gtf[3]) {
					$gene_local{$gene}{"start"} = $gtf[3];				
				}
				if($gene_local{$gene}{"stop"} < $gtf[4]) {
					$gene_local{$gene}{"stop"} = $gtf[4];	
				}				
			}		
		}

	
	}
		

	foreach $key (sort(keys %split)) {
		@pos = split('\s+',$key);
		$gene = $pos[0]."\t".$pos[3];
		if($gene_H{$gene}>0) {	
			unless ($gene_local{$pos[3]}{"start"} < $gene_local{$pos[0]}{"stop"} and $gene_local{$pos[3]}{"stop"} > $gene_local{$pos[0]}{"start"}) {
				$count++;
				print TEMP $count,"\t",$key,"\t",$gene_local{$pos[0]}{"start"},"\t",$gene_local{$pos[0]}{"stop"},"\t",$gene_local{$pos[3]}{"start"},"\t",$gene_local{$pos[3]}{"stop"},"\t",$split{$key},"\t",$gene_H{$gene},"\n";	 
				print SRR $srr{$key};	
			
			}	
		}		
	
	}
	undef %srr;undef %split;undef %gene_H;undef %gene_local;close(POS);close(PAIR);	close(TEMP);close(OUT);close(SRR);	

}


sub ref_Index{
	open(IN,"$out_file/temp");
	open(OUT,">$ref_index/re.fa");
	#open (GTF,"$opt_gtf");
	#my $index = $bowtie_index.".fa";
	#print "\n",$index,"\n";
	my $dna_db = Bio::DB::Fasta->new($bowtie_index.".fa");				
	my ($line, @line, $r_1_gene, $r_1_chr, $r_1_cor, $r_1_start, $r_1_stop, $r_2_gene, $r_2_chr, $r_2_cor,
		 $r_2_start, $r_2_stop, $r_strand, $seq, $count, $boundary);
	
	while($line =<IN>) {
		@line = split(/\s+/,$line);
		$count = $line[0];
		$r_1_gene = $line[1];
		$r_1_chr = $line[2];
		$r_1_cor = $line[3];	
		$r_1_start = $line[8];	
		$r_1_stop = $line[9];	
		$r_2_gene = $line[4];
		$r_2_chr = $line[5];
		$r_2_cor = $line[6];
		$r_2_start = $line[10];	
		$r_2_stop = $line[11];	
		$r_strand = $line[7];
		if($r_strand eq "ff") {
			$seq = $dna_db->seq($r_1_chr, $r_1_start, $r_1_cor).$dna_db->seq($r_2_chr,$r_2_cor, $r_2_stop);
			$boundary = $r_1_cor-$r_1_start+1;
			print OUT ">",$count,"/",$r_1_gene,"/",$r_2_gene,"/",$boundary,"\n";
			while($seq=~s/(\S{50})//) {
				print OUT $1,"\n";
			}
			print OUT $seq,"\n";
		} elsif($r_strand eq "rr") {
			$seq = $dna_db->seq($r_2_chr, $r_2_start, $r_2_cor).$dna_db->seq($r_1_chr,$r_1_cor, $r_1_stop);
			$boundary = $r_2_cor-$r_2_start+1;
			print OUT ">",$count,"/",$r_1_gene,"/",$r_2_gene,"/",$boundary,"\n";
			while($seq=~s/(\S{50})//) {
				print OUT $1,"\n";
			}
			print OUT $seq,"\n";	
		} elsif($r_strand eq "fr") {
			$seq = $dna_db->seq($r_2_chr, $r_2_start, $r_2_cor);
			$seq =~ tr/ATGC/TACG/;
			$seq = reverse($seq);
			$seq = $dna_db->seq($r_1_chr, $r_1_start, $r_1_cor).$seq;
			$boundary = $r_1_cor-$r_1_start+1;
			print OUT ">",$count,"/",$r_1_gene,"/",$r_2_gene,"/",$boundary,"\n";
			while($seq=~s/(\S{50})//) {
				print OUT $1,"\n";
			}
			print OUT $seq,"\n";
		} elsif($r_strand eq "rf") {
			$seq = $dna_db->seq($r_2_chr, $r_2_cor, $r_2_stop);
			$seq =~ tr/ATGC/TACG/;
			$seq = reverse($seq);
			$seq = $seq.$dna_db->seq($r_1_chr,$r_1_cor, $r_1_stop);
			$boundary = $r_2_stop-$r_2_cor+1;
			print OUT ">",$count,"/",$r_1_gene,"/",$r_2_gene,"/",$boundary,"\n";
			while($seq=~s/(\S{50})//) {
				print OUT $1,"\n";
			}
			print OUT $seq,"\n";	
		} 



	}
	system("bowtie-build $ref_index/re.fa $ref_index/re");
	close(IN);close(OUT);
}

sub srr_fq {
	open(SRR,"$fusion_out/srr") || die "Can't open SRR:$!\n";
	open(IN_1,$_[0]) || die "Can't open $_[0]:$!\n";
	open(IN_2,$_[1]) if defined $_[1];
	open(OUT_1,">$ref/er_1.fq") || die "Can't open $ref/er_1.fq:$!\n";
	open(OUT_2,">$ref/er_2.fq") if defined $_[1];
	my $line;
	my %srr;
	my @line;
	while($line=<SRR>) {
		chomp($line);
		$srr{$line} = 1;
	}
	while($line=<IN_1>) {
		@line = split(/\s+/,$line);
		$line[0] =~ s/\@(\S+)/$1/;
		if(exists $srr{$line[0]}) {
			print OUT_1 $line;
			$line=<IN_1>;
			print OUT_1 $line;	
			$line=<IN_1>;
			print OUT_1 $line;	
			$line=<IN_1>;
			print OUT_1 $line;	
		} else {
			<IN_1>;
			<IN_1>;
			<IN_1>;
		}
	}
	if(defined $_[1]) {	
		while($line=<IN_2>) {
			@line = split(/\s+/,$line);
			$line[0] =~ s/\@(\S+)/$1/;
			if(exists $srr{$line[0]}) {
				print OUT_2 $line;
				$line=<IN_2>;
				print OUT_2 $line;	
				$line=<IN_2>;
				print OUT_2 $line;	
				$line=<IN_2>;
				print OUT_2 $line;	
			} else {
				<IN_2>;
				<IN_2>;
				<IN_2>;
			}
		}	

	}
	close(OUT_2) if defined $_[1];
	close(SRR);close(IN_1);close(IN_2);close(OUT_1);undef %srr;


}

sub result {
	open(TEMP,$_[0]) || die "Can't open $_[0]:$!\n";
	open(MAP,$_[1]) || die "Can't open $_[1]:$!\n";
	open(RESULT,">$_[2]") || die "Can't open $_[2]:$!\n";
	#open(OUT,">out") || die "Can't open $_[2]:$!\n";
	my ($line, @line, @temp, $temp, %id, $key, $count, $gene, $cor, $flag, %count, %gene, $xa, $log, $nh);

	
	while($line =<MAP>) {
		
		$xa = 0;
		
		#print $xa,"\n";
		@line = split(/\s+/,$line);
		$line[2] =~ /(\d+)\/(\S+)\/(\S+)\/(\d+)/;
		$count = $1;
		$gene = $2."\t".$3;
		$cor = $4;
		if(exists $id{$line[0]}) {
			
			#print $xa,"\n";
			@temp = split(/\s+/, $id{$line[0]});
			$temp = $temp[1]."\t".$temp[2];
			#print $temp[4],"\n";	
			
			$log = $xa+$temp[4];
			
			if (($gene ne $temp)||($log>2))  {
				
				$id{$line[0]}=join("\t", $1, $2, $3, 2,$log);					
			}elsif ($cor-$line[3]<$read_length and $cor-$line[3]>0) {
				$id{$line[0]}=join("\t", $1, $2, $3, 1,$log);	
				#print OUT $line;
			}
		}elsif ($cor-$line[3]<$read_length and $cor-$line[3]>0) {
			
			$id{$line[0]}=join("\t", $1, $2, $3, 1, $xa);
			#print OUT $line;
			
					
		} else {	
			
			$id{$line[0]}=join("\t", $1, $2, $3, 0, $xa);	
				
		}
		
		
	}
	
	foreach $key(keys(%id)) {
		@line = split(/\s+/,$id{$key});
		$count = $line[0];
		$gene = $line[1]."\t".$line[2];
		$flag = $line[3];
		
		if($flag eq '1') {
			
			$count{$count} = 0 unless exists $count{$count};
			$count{$count}++;	
		} elsif($flag eq '0' ) {
			
			$gene{$gene} = 0 unless exists $gene{$gene};
			$gene{$gene}++;	
			
		}		
	}
	
	while($line=<TEMP>) {
		@line = split(/\s+/,$line);
		$gene = $line[1]."\t".$line[4];
		if(defined $count{$line[0]} and $count{$line[0]}>= $split_num ) {
			print RESULT join("\t",$line[1],$line[2],$line[3],$line[4],$line[5],$line[6],$line[7],$count{$line[0]},$gene{$gene}),"\n" if (exists $count{$line[0]} and exists $gene{$gene});
		}
			
	}

	close(TEMP);close(MAP);close(RESULT);undef %count;undef %id;undef %gene;


}

sub single_result {
	open(TEMP,$_[0]) || die "Can't open $_[0]:$!\n";
	open(MAP,$_[1]) || die "Can't open $_[0]:$!\n";
	open(RESULT,">$_[2]") || die "Can't open $_[2]:$!\n";
	my ($line, @line, @temp, $temp, %id, $key, $count, $gene, $cor, $flag, %count, %gene);
	
	while($line =<MAP>) {
		@line = split(/\s+/,$line);
		$line[2] =~ /(\d+)\/(\S+)\/(\S+)\/(\d+)/;
		$count = $1;
		$gene = $2."\t".$3;
		$cor = $4;
		if ($cor-$line[3]<$read_length and $cor-$line[3]>0) {
			$id{$line[0]}=join("\t", $1, $2, $3, 1);	
		}
		else {		
			$id{$line[0]}=join("\t", $1, $2, $3, 0);		
		}
		
		
	}

	foreach $key(keys(%id)) {
		@line = split(/\s+/,$id{$key});
		$count = $line[0];
		$gene = $line[1]."\t".$line[2];
		$flag = $line[3];
		
		if($flag eq '1') {
			
			$count{$count} = 0 unless exists $count{$count};
			$count{$count}++;	
		} 	
	}
	
	while($line=<TEMP>) {
		@line = split(/\s+/,$line);
		print RESULT join("\t",$line[1],$line[2],$line[3],$line[4],$line[5],$line[6],$line[7],$count{$line[0]}),"\n" if (exists $count{$line[0]});
			
	}

	close(TEMP);close(MAP);close(RESULT);undef %count;undef %id;undef %gene;


}

