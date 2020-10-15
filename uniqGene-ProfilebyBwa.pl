use FindBin qw($Bin);

die "perl uniqGene-ProfilebySoap.pl gene.fa sam out\n" if(@ARGV!=3);
my ($line,@inf,%genelen,%abu,$sum);

open IN, "$ARGV[0]" or die "can not open $ARGV[0]\n";
$/=">";
<IN>;
while($line=<IN>){
	chomp $line;
	@inf=split /\n/,$line;
	my @id=split /\s/,$inf[0];
	my $seq="";
	for(my $i=1;$i<=$#inf;$i++){
		$seq.=$inf[$i];
	}
	$genelen{$id[0]}=length($seq);
	$abu{$id[0]}=0;
}
close IN;
$/="\n";

#sam results
if($ARGV[1]=~/\.bam$/){
	open IN,"$Bin/samtools view $ARGV[1]|" or die "can not open $ARGV[1]\n";
}
else{
	open IN,"$ARGV[1]" or die "can not open $ARGV[1]\n";
}
while($line=<IN>){
	chomp $line;
	@inf=split /\t/,$line;
	next if($line=~/^@/ || $inf[2] eq "*" ||$inf[5] eq "*");
	my ($m,$l)=(0,0);
	my @ele=split /[A-Z]/,$inf[5];
	for(my $i=0;$i<=$#ele;$i++){
		$l+=$ele[$i];
	}
	my $mn=$inf[5]=~s/M/M/g;
	if($mn==1){
		$inf[5]=~/(\d+)M/;
		$m=$1;
	}
	if($mn==2){
		$inf[5]=~/(\d+)M\d+\w(\d+)M/;
		$m+=$1+$2;
	}
	if($mn==3){
		$inf[5]=~/(\d+)M\d+\w(\d+)M\d+\w(\d+)M/;
		$m+=$1+$2+$3;
	}
#	print "$inf[5]\t$m\t$l\n";
	if($m>=50 && ($m/$l)>0.95){
		$sum++;
		$abu{$inf[2]}++;
	}
}
close IN;

open OA, ">$ARGV[2]" or die "can not open $ARGV[2]\n";
print OA "GeneID\tGeneLen\tReadsNum\tAbundance\tRelativeAbundance\tTotalAbundance\n";
my $temp=0;
foreach my $i (sort keys %genelen){
	$temp+=$abu{$i}/$genelen{$i}/$sum;
}
foreach my $i (sort keys %genelen){
	printf OA "$i\t%d\t%ld\t%.8e\t%.8e\t%.8e\n",$genelen{$i},$abu{$i},$abu{$i}/$genelen{$i}/$sum,$abu{$i}/$genelen{$i}/$sum/$temp,$temp;
}
close OA;

