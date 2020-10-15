die "perl Profile-Merge.pl sam.list out\n" if(@ARGV!=2);
my ($line,@inf,%ff,@sample);

open IN, "$ARGV[0]" or die "can not open $ARGV[0]\n";
open OA, ">$ARGV[1]" or die "can not open $ARGV[1]\n";
print OA "GeneID";
while($line=<IN>){
	chomp $line;
	@inf=split /\s+/,$line;
	print OA "\t$inf[0]";
	push @sample,$inf[0];
	open $ff{$inf[0]}, "$inf[1]" or die "can not open $inf[1]\n";
	my $temp=$ff{$inf[0]};<$temp>;
}
print OA "\n";
close IN;

open IN, "$inf[1]" or die "can not open $inf[1]\n";
<IN>;
while($line=<IN>){
	chomp $line;@inf=split /\t/,$line;
	print OA "$inf[0]";
	for(my $j=0;$j<=$#sample;$j++){
		my $temp=$ff{$sample[$j]};
		$line=<$temp>;chomp $line;@inf=split /\t/,$line;
#		print "$j $sample[$j] $inf[4] $line\n";
		if($inf[4]=~/e[+-]/){
			printf OA "\t%.6f",$inf[4]*1000000;
		}
		else{
			printf OA "\t$inf[4]";
		}
	}
	print OA "\n";
}
close OA;
close IN;
