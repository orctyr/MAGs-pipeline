use File::Basename;
use Getopt::Long;
use FindBin qw($Bin);

sub usage{
	print STDERR <<USAGE;
	Version 1.0 2018-12-20 by TaoYe
	Plus Two Binning Methods (Concoct & MaxBin)
	Plus DAStools

	Options 
		-input  <s> : Required (Absolute Dir), Four Columns:
			ID Ass.fa Metabat2-OutDir DepthFile(Metabat2 input) 
		-thread <n> : Thread number, default:20
		-minlen <n> : Contig MinLength, default: 1500 
		-outdir <s> : Output Dir, default: pwd
USAGE
}

my ($input,$thread,$outdir,$minlen);
GetOptions(
	"input:s"=>\$input,
	"thread:n"=>\$thread,
	"minlen:n"=>\$minlen,
	"outdir:s"=>\$outdir,
);
$thread||=20;
$minlen||=1500;
$outdir||=`pwd`;chomp $outdir;

if(!defined($input)){
	usage;
	exit;
}

`mkdir -p $outdir/shell`;

my ($line,@inf);
open IA, "$input" or die "can not open file: $input\n";
while($line=<IA>){
	chomp $line;
	@inf=split /\t/,$line;
	open OA, ">$outdir/shell/S_$inf[0].sh" or die "can not open:S_$inf[0].sh\n";
	`mkdir -p $outdir/DAS_Bin_$inf[0]`;
	print OA "cd $outdir/DAS_Bin_$inf[0]\n";
	print OA "source activate metawrap\n";
	print OA "cut -f 1,3 $inf[3] \| sed -e \'1d\' > $inf[0].depth\n";
	print OA "run_MaxBin.pl  -contig $inf[1] -abund $inf[0].depth -out $inf[0].maxbin -min_contig_length $minlen -thread $thread \n";
	print OA "mkdir -p maxbin maxbin-tempfile concoct concoct-tempfile\n";
	print OA "mv $inf[0].maxbin\*.fasta maxbin\/\n";
	print OA "mv \*.log \*.marker \*.tar.gz \*.noclass \*.summary \*.tooshort maxbin-tempfile\/\n ";
	print OA "concoct --coverage_file $inf[0].depth --composition_file $inf[1] -l $minlen -b concoct-out\n";
	print OA "/lustre/sdb/taoye/miniconda3/envs/metawrap/bin/metawrap-scripts/split_concoct_bins.py concoct-out_clustering_gt$minlen.csv  $inf[1]  .\/ \n";
	print OA "ls bin\*.fa \> concoct.binlist\n";
	print OA "perl $Bin/bin-rename.pl concoct.binlist bin $inf[0].concoct concoct\n";
	print OA "mv concoct-out\*.csv concoct-out\*.txt concoct-tempfile\/\n";
	print OA "rm unbinned.fa\n";
	print OA "mkdir -p metabat2\n";
	print OA "cp $inf[2]/\*.fa .; ls \*.fa \> metabat2.binlist\n";
	print OA "perl $Bin/bin-rename.pl metabat2.binlist $inf[0] $inf[0].metabat2 metabat2\n";
	#print OA "/lustre/sdb/taoye/miniconda3/envs/metawrap/bin/metawrap-scripts/binning_refiner.py -1 concoct/ -2 maxbin/ -3 metabat2/ -o refined \n";
	print OA "source deactivate metawrap\n";
	print OA "source activate DAS_Tool\n";
	print OA "Fasta_to_Scaffolds2Bin.sh -e fasta -i concoct\/ \> concoct.scaffolds2bin.tsv\n";
	print OA "Fasta_to_Scaffolds2Bin.sh -e fasta -i metabat2\/ \> metabat2.scaffolds2bin.tsv\n";
	print OA "Fasta_to_Scaffolds2Bin.sh -e fasta -i maxbin\/ \> maxbin.scaffolds2bin.tsv\n";
	print OA "DAS_Tool -i concoct.scaffolds2bin.tsv,metabat2.scaffolds2bin.tsv,maxbin.scaffolds2bin.tsv -l concoct,metabat2,maxbin --search_engine diamond --write_bins 1 -t $thread --score_threshold 0 -c $inf[1] -o $inf[0] \n";
	print OA "source deactivate DAS_Tool\n";
	#print OA "cd $outdir/DAS_Bin_$inf[0]\n";
	#print OA "source activate metawrap\n";
	#print OA "checkm lineage_wf $inf[0]_DASTool_bins DAS_checkm/ -x fa -t $thread --pplacer_threads 4 --tab_table -f $inf[0].dastool.checkm.summary\n";
	print OA "checkm lineage_wf metabat2/ metabat2_checkm/ -x fa -t $thread --pplacer_threads 4 --tab_table -f $inf[0].metabat2.checkm.summary\n";
	print OA "checkm lineage_wf maxbin/ maxbin_checkm/ -x fasta -t $thread --pplacer_threads 4 --tab_table -f $inf[0].maxbin.checkm.summary\n";
	#print OA "perl $Bin/bin-summary-DAStool.pl $inf[0]\n";
	#close OA;
}
close IA;

