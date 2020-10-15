use File::Basename;
use Getopt::Long;
use FindBin qw($Bin);

my $java="/xiefei/tools/jdk1.8.0_181/bin/java";
my $trimmomatic="/xiefei/tools/bin/module/Trimmomatic-0.33/";
my $bwa="/xiefei/tools/bin/bwa";
my $samtools="/xiefei/tools/bin/samtools";
my $fqstat="/xiefei/tools/bin/module/FastqStat.jar";
my $megahit="/xiefei/tools/bin/module/software/megahit/megahit";
my $idba="/xiefei/tools/bin/idba/bin/";
my $minimus2="/xiefei/tools/bin/minimus2";
my $toAmos="/xiefei/tools/bin/toAmos";

sub usage{
	print STDERR <<USAGE;
	Version 1.0 2018-10-25 by YaoYe
	MetaGenome Pipeline. 
	Including: 
	1.QC; 
	2.Remove Host Genome
	3.Reads Taxnomy
	4.Assembly
	5.Assembly Merge by Minimus2
	6.Contig Correction
	7.Binning
	8.Gene Prediction

	Options 
		-fqlist <s> : Required Input
					  column1: ID  
					  column2: fastq1 file
					  column3: fastq2 file
		-outdir <s> : Output Directory
        -host   <s> : host genome sequences
		-runid  <s> : ID for every run, default: Test1
		-thread <n> : thread number, default: 10
		-set    <s> : group information for geneset(ID sample1 sample2...)
		-help       : show this help
USAGE
}

my ($fqlist,$outdir,$host,$runid,$thread,$set,$help);
GetOptions(
	"fqlist:s"=>\$fqlist,
	"outdir:s"=>\$outdir,
	"host:s"=>\$host,
	"runid:s"=>\$runid,
	"thread:n"=>\$thread,
	"set:s"=>\$set,
	"help"=>\$help,
);
$runid||="Test1";
if(!defined($fqlist)){
	usage;
	exit;
}
$outdir||=`pwd`;chomp $outdir;
my ($line,@inf);
`mkdir -p $outdir/shell`;
`mkdir -p $outdir/stat`;
`mkdir -p $outdir/01.QC`;
`mkdir -p $outdir/02.Tax`;
`mkdir -p $outdir/03.Assembly`;
`mkdir -p $outdir/04.Predict`;
`mkdir -p $outdir/04.PredictGeneSet`;
`mkdir -p $outdir/05.Anno`;
`mkdir -p $outdir/06.BIN`;

open IN, "$fqlist" or die "can not open file: $fqlist\n";
open OA, ">$outdir/01.QC/clean.fq.list" or die "can not open file: $outdir/01.QC/clean.fq.list\n";
open OB, ">$outdir/01.QC/clip.fq.list" or die "can not open file: $outdir/01.QC/clip.fq.list\n";
while($line=<IN>){
	chomp $line;
	my @inf=split /\t/,$line;
	open SH, ">$outdir/shell/S1.QC.$inf[0].$runid.sh" or die "can not open file: $outdir/shell/S1.QC.$inf[0].$runid.sh\n";
	open SH1, ">$outdir/shell/S2.idba.$inf[0].$runid.sh" or die "can not open file: $outdir/shell/S2.idba.$inf[0].$runid.sh\n";
	open SH2, ">$outdir/shell/S2.megahit.$inf[0].$runid.sh" or die "can not open file: $outdir/shell/S2.megahit.$inf[0].$runid.sh\n";
	open SH3, ">$outdir/shell/S3.merge.$inf[0].$runid.sh" or die "can not open file: $outdir/shell/S3.merge.$inf[0].$runid.sh\n";
	open SH4, ">$outdir/shell/S4.correct1.$inf[0].$runid.sh" or die "can not open file: $outdir/shell/S4.correct1.$inf[0].$runid.sh\n";
	open SH5, ">$outdir/shell/S4.correct2.$inf[0].$runid.sh" or die "can not open file: $outdir/shell/S4.correct2.$inf[0].$runid.sh\n";

	print SH "cd $outdir/01.QC/\n";
	print SH "$java -jar $trimmomatic/trimmomatic-0.33.jar PE -threads $thread -phred33 $inf[1] $inf[2] $outdir/01.QC/$inf[0].clip.1.fq.gz $outdir/01.QC/$inf[0].single.R1.fastq.gz $outdir/01.QC/$inf[0].clip.2.fq.gz $outdir/01.QC/$inf[0].single.R2.fastq.gz ILLUMINACLIP:$trimmomatic/adapters/TruSeq2-PE.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:75\n";
	print OA "$inf[0]\t$outdir/01.QC/$inf[0].clean.1.fq.gz\t$outdir/01.QC/$inf[0].clean.2.fq.gz\n";
	print OB "$inf[0]\t$outdir/01.QC/$inf[0].clip.1.fq.gz\t$outdir/01.QC/$inf[0].clip.2.fq.gz\n";
	print SH "echo \"========QC Finished, Removing Host Genome=========\"\n\n";
	#remove host
	print SH "$bwa mem -t $thread -M $host $outdir/01.QC/$inf[0].clip.1.fq.gz $outdir/01.QC/$inf[0].clip.2.fq.gz \| awk \'\$3\~\/chr\/\' \| cut -f 1\,3 \> $outdir/01.QC/$inf[0].host.temp \n";
	print SH "perl $Bin/remove-host.pl $outdir/01.QC/$inf[0].host.temp $inf[0]\n";

	#Megahit pipeline
	print SH2 "cd $outdir/03.Assembly\n";
	print SH2 "$megahit -t $thread -1 $outdir/01.QC/$inf[0].clean.1.fq.gz -2 $outdir/01.QC/$inf[0].clean.2.fq.gz --min-contig-len 500 -o $inf[0]-megahit\n";
	print SH2 "ln -s $inf[0]-megahit/final.contigs.fa $inf[0].megahit.contigs\n";
	print SH2 "perl $Bin/assemle-summary.pl $inf[0].megahit.contigs > $inf[0].megahit.contigs.summary\n";
	#IDBA pipeline
	print SH1 "cd $outdir/03.Assembly\n";
	print SH1 "gunzip -dc $outdir/01.QC/$inf[0].clean.1.fq.gz > $inf[0].clean.1.fq\n";
	print SH1 "gunzip -dc $outdir/01.QC/$inf[0].clean.2.fq.gz > $inf[0].clean.2.fq\n";
	print SH1 "$idba/fq2fa --merge $inf[0].clean.1.fq $inf[0].clean.2.fq $inf[0].clean.fa\n";
	print SH1 "$idba/idba_ud -l $inf[0].clean.fa --pre_correction --min_contig 500 --num_threads $thread -o $inf[0].idba --mink 90 --maxk 124 \n";
	print SH1 "rm $inf[0].clean.1.fq $inf[0].clean.2.fq\n";
	print SH1 "mv $inf[0].idba/contig.fa $inf[0].idba/contig.fa.temp\n";
	print SH1 "perl $Bin/fa_check.pl $inf[0].idba/contig.fa.temp $inf[0].idba/contig.fa\n";
	print SH1 "rm $inf[0].idba/contig.fa.temp\n";
	print SH1 "ln -s $inf[0].idba/contig.fa $inf[0].idba.contigs\n";
	print SH1 "perl $Bin/assemle-summary.pl $inf[0].idba.contigs > $inf[0].idba.contigs.summary\n";
	#merge by minimus2
	print SH3 "cd $outdir/03.Assembly\n";
	print SH3 "cat $inf[0].megahit.contigs $inf[0].idba.contigs > $inf[0].temp.merge.fa\n";
	print SH3 "$toAmos -s $inf[0].temp.merge.fa -o $inf[0].afg\n";
	print SH3 "$minimus2 $inf[0] -D CONSERR=0 -D OVERLAP=100 -D MINID=100 -D REFCOUNT=\`grep -c \"\>\" $inf[0].megahit.contigs\` \n";
	print SH3 "cat $inf[0].fasta $inf[0].singletons.seq > $inf[0].merge.contigs\n";
	print SH3 "rm -rf $inf[0].idba.contigs.filter $inf[0].temp.merge.fa $inf[0].frg $inf[0].bnk $inf[0].runAmos.log $inf[0].afg $inf[0].ref.seq $inf[0].qry.seq $inf[0].delta $inf[0].coords $inf[0].ovl $inf[0].OVL $inf[0].contig $inf[0].fasta $inf[0].singletons.seq $inf[0].singletons\n";
	print SH3 "perl $Bin/assemle-summary.pl $inf[0].merge.contigs> $inf[0].merge.contigs.summary\n";
	#correct by pilon
	print SH4 "cd $outdir/03.Assembly\n";
	print SH4 "sed -e \'s\/n\/A\/g\' $inf[0].merge.contigs > $inf[0].temp.fa \n";
	print SH4 "$Bin/bwa index $inf[0].temp.fa\n";
	print SH4 "$Bin/bwa aln $inf[0].temp.fa -n 0.02 -e 2 -t $thread -f $inf[0].clean.1.fq.gz.sai $outdir/01.QC/$inf[0].clean.1.fq.gz\n";
	print SH4 "$Bin/bwa aln $inf[0].temp.fa -n 0.02 -e 2 -t $thread -f $inf[0].clean.2.fq.gz.sai $outdir/01.QC/$inf[0].clean.2.fq.gz\n";
	print SH4 "$Bin/bwa sampe -a 1000 -N 5 $inf[0].temp.fa $inf[0].clean.1.fq.gz.sai $inf[0].clean.2.fq.gz.sai $outdir/01.QC/$inf[0].clean.1.fq.gz $outdir/01.QC/$inf[0].clean.2.fq.gz | $Bin/samtools view -bS -F 4 -T $inf[0].temp.fa - \> $inf[0].bam\n";
	print SH4 "$Bin/samtools sort $inf[0].bam $inf[0].sorted\n";
	print SH4 "$Bin/samtools index $inf[0].sorted.bam\n";
	print SH4 "rm $inf[0].clean.1.fq.gz.sai $inf[0].clean.2.fq.gz.sai $inf[0].bam\* \n";

	print SH5 "cd $outdir/03.Assembly\n";
	print SH5 "$Bin/samtools mpileup -guSDf $inf[0].temp.fa $inf[0].sorted.bam | $Bin/bcftools view -cvNg - > $inf[0].vcf\n";
	print SH5 "perl $Bin/seq_correct_vcf.pl -input $inf[0].temp.fa -vcf $inf[0].vcf -output $inf[0].correct.contigs\n";
	print SH5 "perl $Bin/contig-rename.pl $inf[0].correct.contigs $inf[0].final.contigs $inf[0]\n";
	print SH5 "perl $Bin/assemle-summary.pl $inf[0].final.contigs > $inf[0].final.contigs.summary\n";
	print SH5 "rm $inf[0].sorted.bam\* $inf[0].vcf $inf[0].temp.fa\* $inf[0].correct.contigs\n";

	#Binning-Checkm
	open SH6, ">$outdir/shell/S5.BIN.$inf[0].$runid.sh" or die "can not open file: $outdir/shell/S5.BIN.$inf[0].$runid.sh\n"; 
	print SH6 "cd $outdir/06.BIN/ \n";
	print SH6 "ln -s $outdir/03.Assembly/$inf[0].final.contigs\n";
	print SH6 "$Bin/bwa index $inf[0].final.contigs\n";
	print SH6 "$Bin/samtools faidx $inf[0].final.contigs\n";
	print SH6 "$Bin/bwa mem  -t $thread -M -R \'\@RG\\tID:$inf[0]\\tSM:$inf[0]\\tLB:$inf[0]\\tPL:Illumina\\tPI:500\' $inf[0].final.contigs $outdir/01.QC/$inf[0].clip.1.fq.gz $outdir/01.QC/$inf[0].clip.2.fq.gz \| $Bin/samtools view -bS -t $inf[0].final.contigs.fai - > $inf[0].bam\n";
	print SH6 "$Bin/samtools sort -m 6960000000 $inf[0].bam $inf[0].sort\n";
	print SH6 "$Bin/samtools index $inf[0].sort.bam\n";
	print SH6 "$Bin/metabat/jgi_summarize_bam_contig_depths --outputDepth $inf[0].depth $inf[0].sort.bam\n";
	print SH6 "rm $inf[0].bam\n";
	print SH6 "$Bin/metabat/metabat2 -t $thread --inFile $inf[0].final.contigs --abdFile $inf[0].depth  --outFile $outdir/06.BIN/$inf[0]/$inf[0] \n";
	print SH6 "ls $outdir/06.BIN/$inf[0]/$inf[0]\*fa > $inf[0].bins.list\n";
	print SH6 "source activate metawrap\ncheckm lineage_wf $outdir/06.BIN/$inf[0]/ $outdir/06.BIN/$inf[0]_checkm -x fa -t $thread --pplacer_threads 4 --tab_table -f $inf[0].checkm.summary\n";

	#gene predict
	open SH7, ">$outdir/shell/S6.predict.$inf[0].$runid.sh" or die "can not open file: $outdir/shell/S6.predict.$inf[0].$runid.sh\n";
	print SH7 "cd $outdir/04.Predict\n";
	print SH7 "$Bin/prodigal -a $inf[0].temp.orf.faa -i $outdir/03.Assembly/$inf[0].final.contigs -f gff -o $inf[0].gff -p meta -q -d $inf[0].temp.orf.ffn\n";
	print SH7 "perl $Bin/gene-filter.pl $inf[0].temp.orf.ffn $inf[0].temp.orf.faa $inf[0].orf.ffn $inf[0].orf.faa \n";
	print SH7 "rm $inf[0].temp.orf.faa $inf[0].temp.orf.ffn \n";

	#bwa mapping geneset
	open SH8, ">$outdir/shell/S8.Abund.$inf[0].$runid.sh" or die "S8.Abund.$inf[0].$runid.sh\n"; 
	print SH8 "cd $outdir/07.Profile \n";
	print SH8 "$Bin/bwa mem -t $thread -M -R \'\@RG\\tID:$inf[0]\\tSM:$inf[0]\\tLB:$inf[0]\\tPL:Illumina\\tPI:500\' $outdir/04.PredictGeneSet/$runid.geneSet.ffn $outdir/01.QC/$inf[0].clean.1.fq.gz $outdir/01.QC/$inf[0].clean.2.fq.gz  > $inf[0].sam\n";
	print SH8 "perl $Bin/uniqGene-ProfilebyBwa.pl $outdir/04.PredictGeneSet/$runid.geneSet.ffn $inf[0].sam $inf[0].sam.abd\n";
	print SH9 "$inf[0]\t$inf[0].sam.abd\n";
	
	close (SH,SH1,SH2,SH3,SH4,SH5,SH6,SH7,SH8,SH9);
}
close OA;
close OB;
close IN;

#merge sample cmd
open SH, ">$outdir/shell/S1.STAT.$runid.sh" or die "can not open file: $outdir/shell/S1.STAT.$runid.sh\n";
print SH "cd $outdir/stat/\n";
print SH "$java -jar $fqstat -i $outdir/$fqlist > $outdir/stat/raw.$runid.stat.xls\n";
print SH "$java -jar $fqstat -i $outdir/01.QC/clean.fq.list > $outdir/stat/clean.$runid.stat.xls\n";
print SH "$java -jar $fqstat -i $outdir/01.QC/clip.fq.list > $outdir/stat/clip.$runid.stat.xls\n";
close SH;

#Create Geneset using CDHIT
open IN,"$set" or die "Not DO GeneSet!\n";
open SH, ">$outdir/shell/S7.GeneSet.$runid.sh" or die "can not open $outdir/shell/S7.GeneSet.$runid.sh\n";
print SH "cd $outdir/04.PredictGeneSet\n";
print SH "\#rm \*all.ffn\n";
while($line=<IN>){
	chomp $line;
	@inf=split /\s+/,$line;
	for(my $i=1;$i<=$#inf;$i++){
		print SH "cat $outdir/04.Predict/$inf[$i].orf.ffn.gz >> $inf[0].all.ffn.gz\n";
	}
	print SH "\#$Bin/cd-hit-est -i $inf[0].all.ffn -o $inf[0].geneSet.ffn -n 9 -c 0.95 -G 0 -M 0 -d 0 -aS 0.9 -r 1 -T 80\n";
	print SH "\#rm $inf[0].all.ffn\n";
}
close SH;


