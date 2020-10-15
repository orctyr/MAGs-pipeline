die ".pl input.fa output.fa name\n" if(@ARGV!=3);
open IN, "$ARGV[0]" or die "can not open file: $ARGV[0]\n";
open OA, ">$ARGV[1]" or die "can not open file: $ARGV[1]\n";

$/=">";<IN>;
my $n=0;
while($line=<IN>){
	chomp $line;
	my @ele=split /\n/,$line;
	$n++;
	my $seq="";
	for(my $i=1;$i<=$#ele;$i++){
		$seq.=$ele[$i];
	}
	print OA ">$ARGV[2]__c$n\n$seq\n";
}
close IN;
close OA;
