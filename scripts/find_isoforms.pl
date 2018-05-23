#!/usr/bin/perl
open(INT, "$ARGV[0]")|| die "I can't open intersect file (gene/all possible transcripts)\n";
$thres=$ARGV[1];
$flagini=0;
while($l=<INT>){
	chomp $l; 
	@vec=split("\t", $l);
	$gene=$vec[0]."_".$vec[1]."_".$vec[2];
	$isofo=$vec[3]."_".$vec[4]."_".$vec[5];
	$size=$vec[7];
	$h{$gene}{$isofo}=$size;
	if($gene ne $passedgene && $flagini==1){ 
		@isofosbr = sort {$h{$passedgene}{$b} <=> $h{$passedgene}{$a} or $a cmp $b} keys (%{$h{$passedgene}});
		for($cont=0;$isofosbr[$cont+1]=~/.+/;$cont++){
			for($cont2=$cont+1;$isofosbr[$cont2]=~/.+/;$cont2++){
				if($h{$passedgene}{$isofosbr[$cont2]}){
					@aux1=split("_",$isofosbr[$cont]);
					@aux2=split("_",$isofosbr[$cont2]);
					$dif=abs($aux1[1]-$aux2[1])+abs($aux1[2]-$aux2[2]);
					if($dif<=$thres){$h{$passedgene}{$isofosbr[$cont2]}=0;}
				}
			}
		}
		foreach $isopas(keys %{$h{$passedgene}}){
			if($h{$passedgene}{$isopas}){
				$isopas =~ tr/_/\t/;
				print "$isopas\n";
			}
		}
	 	undef (%{$h{$passedgene}});
	}
	$passedgene=$gene;
	if($flagini==0){$flagini=1;}
}

@isofosbr = sort {$h{$passedgene}{$b} <=> $h{$passedgene}{$a} or $a cmp $b} keys (%{$h{$passedgene}});
for($cont=0;$isofosbr[$cont+1]=~/.+/;$cont++){
	for($cont2=$cont+1;$isofosbr[$cont2]=~/.+/;$cont2++){
		if($h{$passedgene}{$isofosbr[$cont2]}){
			@aux1=split("_",$isofosbr[$cont]);
			@aux2=split("_",$isofosbr[$cont2]);
			$dif=abs($aux1[1]-$aux2[1])+abs($aux1[2]-$aux2[2]);
			if($dif<=$thres){$h{$passedgene}{$isofosbr[$cont2]}=0;}
		}
	}
}
foreach $isopas(keys %{$h{$passedgene}}){
	if($h{$passedgene}{$isopas}){
		$isopas =~ tr/_/\t/;
		print "$isopas\n";
	}
}
undef (%{$h{$passedgene}});
