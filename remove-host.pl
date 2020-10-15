die "perl remove-host.pl host.temp name\n" if(@ARGV!=2);

my ($line,@inf,%id,%host);
open IN, "$ARGV[0]" or die "can not open: $ARGV[0]\n";
while($line=<IN>){
	chomp $line;
	@inf=split /\t/,$line;
	if(!defined $id{$inf[0]}){
		$id{$inf[0]}=1;
		my @ele=split /_/,$inf[1];
		$host{$ele[0]}++;
	}
}
close IN;

my $total=0;
open IA,"gunzip -dc $ARGV[1].clip.1.fq.gz|" or die "can not open $ARGV[1].clip.1.fq.gz\n";
open IB,"gunzip -dc $ARGV[1].clip.2.fq.gz|" or die "can not open $ARGV[1].clip.2.fq.gz\n";
open OA,"|gzip > $ARGV[1].clean.1.fq.gz" or die "can not open $ARGV[1].clean.1.fq.gz\n";
open OB,"|gzip > $ARGV[1].clean.2.fq.gz" or die "can not open $ARGV[1].clean.2.fq.gz\n";
while($line=<IA>){
	$total++;
	my @r1=();
	my @r2=();
	push @r1,$line;
	$line=<IB>;	push @r2,$line;
	for(my $i=1;$i<=3;$i++){
		$line=<IA>;push @r1,$line;
		$line=<IB>;push @r2,$line;
	}
	@inf=split /\s+/,$r1[0];
	$inf[0]=~s/^@//;
	if(!defined $id{$inf[0]}){
		print OA "$r1[0]$r1[1]$r1[2]$r1[3]";
		print OB "$r2[0]$r2[1]$r2[2]$r2[3]";
	}
}
close (IA,IB,OA,OB);

open OA, ">$ARGV[1].host.info" or die "can not open $ARGV[1].host.info\n";
printf OA "$ARGV[1]\tTotalInputReadsPair\t%ld\n",$total;
my $con=0;
foreach my $i (sort keys %host){
	print OA "$ARGV[1]\t$i\t$host{$i}\n";
	$con+=$host{$i};
}
printf OA "$ARGV[1]\tTotalCleanReadsPair\t%ld\n",$total-$con;
close OA;

