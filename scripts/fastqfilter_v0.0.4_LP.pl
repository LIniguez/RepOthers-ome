#!/usr/bin/perl
use Getopt::Long;
$fastqfil="";
$minbasequal=20;
$maxN=0.2;
$wsize=5;
$trimq=20;
$minread=100;
$smallsqs=0;
$loquant=0.25;
$midquant=0.5;
$minavgqual=28;
$lowqTres=15;
$midqTres=25;
$highqTres=30;
$encoding=33;
$noens=0;
$output= "fastq,fasta,summary";
$trm="0";

GetOptions(
	'fastqfil=s' => \$fastqfil,
	'encoding=i' => \$$encoding,
	'minread=s' => \$minread,
	'minbasequal=i' => \$minbasequal,
	'minavgqual=i' => \$minavgqual,
	'maxN=s' => \$maxN	,
	'output=s' => \$output,
	'loquant=f' => \$lodquant,
	'midquant=f' => \$midquant,
	'loqmin=i' => \$lowqTres,
	'midqmin=i' => \$midqTres,
	'hiqmin=i' => \$highqTres,
	'noens!' => \$noens,
	'trm!' => \$trm,
	'triwinsize=i'=>\$wsize,
	'trimin=i' =>\$trimq,
	'h|help!' =>\$help

)or die "Usage: $0 [options] --fastqfil <fastq_file>\n";

if(!$fastqfil && !$help){die "Usage: $0 [options] --fastqfil <fastq_file>\n";}
if($midquant>=1||$midquant<=0||$loquant>=1||$loquant<=0||$midquant<=$loquant){
	print "Incorrect arguments";
	$help=1;
}else{@quantiles=(0.1,0.25,0.5,0.75,0.9,1,$loquant,$midquant);}

@temp=split(",",$output);
$cont=0;
foreach $elem (@temp){
	if($elem eq "fastq" || $elem eq "fq"){
		$ofq=1;$cont++;
	}
	if($elem eq "fasta" || $elem eq "fa"){
		$ofa=1;$cont++;
	}
	if($elem eq "summary" || $elem eq "sum"){
		$summary=1;$cont++;
	}
}
if(!$cont){
	print "Incorrect arguments";
	$help=1;
}


if($help){
	print "\nUsage:\n\n$0 [options] --fastqfil <fastq_file>\n\n";
	print "Options:\n";
	print "  --fastqfil\tFILE\tfilename (with path if necessary) of the FASTQ file\n";
	print "  --encoding\tINT\tASCII phred quality encoding [33]\n";
	print "  --minread\tINT\tminimum read lenght to be processed for quality [100]\n";
	print "  --minbasequal\tINT\tbase phred quality of bases to be changed to N, unless noens = FALSE [20]\n";
	print "  --minavgqual\tINT\tminimum overall average phred quality to accept a read [28]\n";
	print "  --maxN\tSTR\tmaximum number of N's allowed in a read (if noens=FALSE), either as an integer or a proportion or each read lenght [0.2]\n";
	print "  --output\tSTR\ta string, comma separeted, with either or all of \"fastq\", \"fas\", \"summary\" (may be abbrviated \"fa\", \"fq\", \"sum\") defining output files [fq,fa,sum]\n";
	print "  --loquant\tFLT\tlower quantile in the read quality distribution [0.25]\n";
	print "  --midquant\tFLT\tlarger quantile in the read quality distribution than loquant [0.5]\n";
	print "  --loqmin\tINT\tminium base phred quality average on the loquant quantile [15]\n";
	print "  --midqmin\tINT\tminium base phred quality average on the quantile beteween loquant and midquant [25]\n";
	print "  --hiqmin\tINT\tminium base phred quality average on quantile > miquant [30]\n";
	print "  --triwinsize\tINT\twindow size used for triming sequences for its 3' end bases on an average sliding window method [5]\n";
	print "  --trimin\tINT\tminium average base phred quality on the window for trimming sequences [20]\n";
	print "  --notrm\t\tdo not trim sequences based on an average sliding window mehtod with triwinsize window size and a minimum phred average of trimin\n";
	print "  --noens\t\tdo not allow minbasequal N substitutions nor N's in the original file\n";
	exit(-1);
}

if($fastqfil=~/(.+)\.fq\.gz$/||$fastqfil=~/(.+)\.FQ\.gz$/||$fastqfil=~/(.+)\.fastq\.gz$/||$fastqfil=~/(.+)\.FASTQ\.gz$/){
	open(FA,"gunzip -c $fastqfil |")||die "can't open $fastqfil\n";
	$name_file=$1;
}elsif($fastqfil=~ /(.+)\.fq\.bz2$/||$fastqfil=~ /(.+)\.FQ\.bz2$/||$fastqfil=~ /(.+)\.fastq\.bz2$/||$fastqfil=~ /(.+)\.FASTQ\.bz2$/){
	open(FA,"bunzip2 -c $fastqfil |")||die "can't open $fastqfil\n";
	$name_file=$1;
}elsif($fastqfil=~/(.+)\.fq$/||$fastqfil=~/(.+)\.FQ$/||$fastqfil=~/(.+)\.fastq$/||$fastqfil=~/(.+)\.FASTQ$/){
	open(FA, "$fastqfil")||die "can't open $fastqfil\n";
	$name_file=$1;
}else{
	open(FA, "$fastqfil")||die "can't open $fastqfil\n";
	$name_file=$fastqfil;
}

if($ofq){open(OFAQ, ">$name_file"."_filtered_LP.fastq")||die "cant't create $name_file"."_filtered.fastq";}
if($ofa){open(OFAT, ">$name_file"."_filtered_LP.fa")||die "cant't create $name_file"."_filtered.fasta";}

$seqwN=0;
$fil_maxN=0;
$fil_lowquan=0;
$fil_midquan=0;
$fil_highquan=0;
$fil_avrg=0;
$seq_an=0;
$seq_pass=0;
if($summary){
	@quantiles_out=(0.1,0.25,0.5,0.75,0.9,1);
}

while(<FA>){
	$seq_an++;
	$id=$_;
	chomp $id;
	$seq=<FA>; 
	chomp $seq;
	<FA>;
	$qual=<FA>;
	chomp $qual;
	undef @vec_seq;
	undef @vec_qual;
	@vec_seq=split("",$seq);
	@vec_qual=split("",$qual);
	$qsum=0;
	$cumqsum[0]=$qsum;
	$fminbqual=0;
	$postrim=0;
	undef @rsum;
	$nopass=0;
	for($cont=0; $cont<=$#vec_qual; $cont++){
		$qint=(ord($vec_qual[$cont])-$encoding); 
		$qualint[$cont]=$qint;
		$qsum+=$qint;
		$cumqsum[$cont+1]=$qsum;
	} 
	for($cont=0;$cont<=($#vec_qual - $wsize +1); $cont++){
		$rsum[$cont]=($cumqsum[$cont + $wsize] -$cumqsum[$cont] )/$wsize; 
	}
	for($cont=$#rsum; $cont>=0;$cont--){ 
		if($rsum[$cont]>$trimq){
			$postrim=$cont+$wsize-1;
			last;
		}
	} 
	if($postrim!=$#vec_qual){$postrim=$postrim-$wsize;}
	for($cont=$#vec_qual;$cont>$postrim;$cont--){
		pop @vec_seq; pop @vec_qual; pop @qualint;
	}
	$seq=join("",@vec_seq);
	$qual=join("",@vec_qual);
	if($noens){
		if($seq =~/N/){ $nopass=1; $seqwN++;}
	}
	if($#vec_seq>=$minread && !$nopass){
		@perc=@quantiles;
		@qualint_sort=sort { $a <=> $b } @qualint;
		for($cont=0;$cont<=$#quantiles;$cont++){
			$perc[$cont]=$qualint_sort[int(($quantiles[$cont]*($#qualint_sort+1))+0.5)-1]
		}
		$midquanres=pop(@perc);
		$lowquanres=pop(@perc);
		
		$sumlow=0;
		$summed=0;
		$sumhigh=0;
		$contlow=0;
		$contmed=0;
		$conthigh=0;
		$sumtot=0;
		$fminbqual=0;
		if($maxN<1){$tempmaxN=int(($maxN*($#qualint_sort+1))+0.5);}else{$tempmaxN=$maxN;}
		for($cont=0; $cont<=$#qualint;$cont++){
			if($qualint[$cont]<$lowquanres){
				$sumlow+=$qualint[$cont];
				$contlow++;
			}elsif($qualint[$cont]<$midquanres){
				$summed+=$qualint[$cont];
				$contmed++;
			}else{	$sumhigh+=$qualint[$cont];
				$conthigh++;
			}
			$sumtot+=$qualint[$cont];
			if(($qualint[$cont]<$minbasequal && $fminbqual<$tempmaxN )||$vec_seq[$cont]eq"N"){
				$vec_seq[$cont]="N";
				$fminbqual++;
			}
		}
		if($fminbqual>=$tempmaxN){$fil_maxN++;$nopass=1;}
		if(!$nopass){
			if(!$contlow){$promlow=$lowquanres;}else{ $promlow=$sumlow/$contlow;}
			if(!$contmed){$prommed=$midquanres;}else{ $prommed=$summed/$contmed;}
			$promhigh=$sumhigh/$conthigh;
			$promtot=$sumtot/($#qualint+1);

			if($promlow<$lowqTres){$fil_lowquan++;$nopass=1;}
			if($prommed<$midqTres&&!$nopass){$fil_midquan++;$nopass=1;}
			if($promhigh<$highqTres&&!$nopass){$fil_highquan++;$nopass=1;}
			if($promtot<$minavgqual&&!$nopass){$fil_avrg++;$nopass=1;}
		}
			
	}else{$smallsqs++; $nopass=1;}
		
	if(!$nopass){
		$seq=join("",@vec_seq);
		if($ofq){print OFAQ "$id\n$seq\n+\n$qual\n";}
		if($ofa){print OFAT ">$id\n$seq\n";}
		if($summary){
		for($cont=0;$cont<=$#quantiles_out;$cont++){
			$quan[$cont][$seq_pass]=$perc[$cont];
		}
		$len[$seq_pass]=$postrim+1;
		}
		$seq_pass++;
	}
}
close(FA, OFAQ, OFAS);

if($summary){

	open(SUM, ">$name_file"."_summary_LP.txt")||die "cant't create $name_file"."_summary.txt";
	for($cont=0;$cont<=$#quantiles_out;$cont++){ 
		@temp=@{$quan[$cont]};
		@temp_sort=sort {$a <=> $b } @temp;
		for($cont2=0;$cont2<=$#quantiles_out;$cont2++){
			$nuerito=int(($quantiles_out[$cont]*($#temp_sort+1))+0.5); 
			$perc[$cont][$cont2]=$temp_sort[int(($quantiles_out[$cont2]*($#temp_sort+1))+0.5)-1];
		}
		@temp=@len;
		@temp_sort=sort {$a<=> $b} @temp;
		$perc_len[$cont]=$temp_sort[int(($quantiles_out[$cont]*($#temp_sort+1))+0.5)-1];
	}
	print SUM "Number of reads in file:\t\t\t\t$seq_an\n";
	print SUM "Quality quantile assessment (columns\: bases, rows: reads):\n";
	print SUM "\t";
	for($cont=0;$cont<=$#quantiles_out;$cont++){
		$gout=$quantiles_out[$cont] * 100;
		print SUM "\t$gout\%";
	}
	print SUM "\n";
	for($cont=0;$cont<=$#quantiles_out;$cont++){
		$gout=$quantiles_out[$cont] * 100;
		print SUM "\t$gout\%";
		for($cont2=0;$cont2<=$#quantiles_out;$cont2++){
			print SUM "\t$perc[$cont][$cont2]";
		}
		print SUM "\n";
	}
	print SUM "Read lengths distribution:\n";
	for($cont=0;$cont<=$#quantiles_out;$cont++){
		$gout=$quantiles_out[$cont] * 100;
		print SUM "\t$gout\%";
	}
	print SUM "\n";
	for($cont=0;$cont<=$#quantiles_out;$cont++){
		print SUM "\t$perc_len[$cont]";
	}
	print SUM "\n";
	$porlength=int(($smallsqs/$seq_an)*100 +0.5);
	print SUM "Number of reads < $minread bp:\t\t\t$smallsqs\t($porlength %)\n";
	$porlowquan=int(($fil_lowquan/$seq_an)*100 +0.5);
	$oloquant=int($loquant*100);
	print SUM "Reads filtered by quantile $oloquant % < $lowqTres:\t$fil_lowquan\t($porlowquan %)\n";
	$pormidquan=int(($fil_midquan/$seq_an)*100 +0.5);
	$omidquant=int($midquant*100);
	print SUM "Reads filtered by quantile $omidquant % < $midqTres:\t$fil_midquan\t($pormidquan %)\n";
	$porhighquan=int(($fil_highquan/$seq_an)*100 +0.5);
	print SUM "Reads filtered by quantile >$omidquant % < $highqTres:\t$fil_highquan\t($porhighquan %)\n";
	$poravrg=int(($fil_avrg/$seq_an)*100+0.5);
	print SUM "Reads filtered by a minimum average of $minavgqual:\t$fil_avrg\t($poravrg %)\n";
	$pormaxN=int(($fil_maxN/$seq_an)*100+0.5);
	print SUM "Reads with more than $maxN N base substitutions at <= $minbasequal threshold:\t$fil_maxN\t($pormaxN %)\n";
	$totfiltered=$smallsqs+$fil_lowquan+$fil_midquan+$fil_highquan+$fil_avrg+$fil_maxN;
	$portotfil=int(($totfiltered/$seq_an)*100+0.5);
	print SUM "Total filtered reads:\t$totfiltered\t($portotfil %)\n";
	$remreads=$seq_an - $totfiltered;
	print SUM "Total remaining reads:\t$remreads\n";

	close(SUM);
}




