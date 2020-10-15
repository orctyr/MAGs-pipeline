use File::Basename;
use Getopt::Long;
use FindBin qw($Bin);

my $java="/lustre/sdb/zengl/tools/jdk1.8.0_181/bin/java";
my $trimmomatic="/lustre/sdb/zengl/bin/module/Trimmomatic-0.33/";
my $bwa="/lustre/sdb/taoye/bin/bwa";
my $samtools="/lustre/sdb/taoye/bin/samtools";
my $fqstat="/lustre/sdb/zengl/bin/module/FastqStat.jar";
my $megahit="/lustre/sdb/zengl/bin/module/software/megahit/megahit";

sub usage{
	print STDERR <<USAGE;
	Version 1.0 2019-01-02 by YaoYe
	MetaGenome Pipeline. 
	Including: 
	1.QC; 
	2.Remove Host Genome
	3.Reads Taxnomy
	4.Assembly:  Megahit (> 500bp)
	5.Binning
	6.Gene Prediction and Annotation
	7.Abundance

	Options 
		-fqlist <s> : Required Input
					  column1: ID  
					  column2: fastq1 file
					  column3: fastq2 file
		-outdir <s> : Output Directory
        -host   <s> : host genome sequences
		-runid  <s> : ID for every run, default: Test1
		-thread <n> : thread number, default: 10
		-help       : show this help
USAGE
}

my ($fqlist,$outdir,$host,$runid,$thread,$help);
GetOptions(
	"fqlist:s"=>\$fqlist,
	"outdir:s"=>\$outdir,
	"host:s"=>\$host,
	"runid:s"=>\$runid,
	"thread:n"=>\$thread,
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
`mkdir -p $outdir/07.Profile`;
`mkdir -p $outdir/08.MAGsPorfile`;
`mkdir -p $outdir/08.MAGsPorfile/BINs`;

open IN, "$fqlist" or die "can not open file: $fqlist\n";
open OA, ">$outdir/01.QC/clean.fq.list" or die "can not open file: $outdir/01.QC/clean.fq.list\n";
open OB, ">$outdir/01.QC/clip.fq.list" or die "can not open file: $outdir/01.QC/clip.fq.list\n";
#gene cd-hit shell
open CD, ">$outdir/shell/S7.GeneSet.$runid.sh" or die "can not open $outdir/shell/S7.GeneSet.$runid.sh\n";
print CD "cd $outdir/04.PredictGeneSet\n";
print CD "rm \*all.ffn\n";
open SH4, ">$outdir/07.Profile/bam.list" or die "$outdir/07.Profile/bam.list\n";
open SH5, ">$outdir/shell/S8-Merge.Abund.$runid.sh" or die "$outdir/shell/S8-Merge.Abund.$runid.sh\n";
open ANNO, ">$outdir/shell/S9.Anno.$runid.sh" or die "$outdir/shell/S9.Anno.$runid.sh\n";

my $cleanfastq="";
while($line=<IN>){
	chomp $line;
	my @inf=split /\t/,$line;
	open SH, ">$outdir/shell/S1.QC.$inf[0].$runid.sh" or die "can not open file: $outdir/shell/S1.QC.$inf[0].$runid.sh\n";
	open SH2, ">$outdir/shell/S2.megahit.$inf[0].$runid.sh" or die "can not open file: $outdir/shell/S2.megahit.$inf[0].$runid.sh\n";

	print SH "cd $outdir/01.QC/\n";
	print SH "$java -jar $trimmomatic/trimmomatic-0.33.jar PE -threads $thread -phred33 $inf[1] $inf[2] $outdir/01.QC/$inf[0].clip.1.fq.gz $outdir/01.QC/$inf[0].single.R1.fastq.gz $outdir/01.QC/$inf[0].clip.2.fq.gz $outdir/01.QC/$inf[0].single.R2.fastq.gz ILLUMINACLIP:$trimmomatic/adapters/merge.fa:2:30:10 SLIDINGWINDOW:4:15 MINLEN:75\n";
	print OA "$inf[0]\t$outdir/01.QC/$inf[0].clean.1.fq.gz\t$outdir/01.QC/$inf[0].clean.2.fq.gz\n";
	print OB "$inf[0]\t$outdir/01.QC/$inf[0].clip.1.fq.gz\t$outdir/01.QC/$inf[0].clip.2.fq.gz\n";
	print SH "echo \"========QC Finished, Removing Host Genome=========\"\n\n";
	#remove host
	print SH "$bwa mem -t $thread -M $host $outdir/01.QC/$inf[0].clip.1.fq.gz $outdir/01.QC/$inf[0].clip.2.fq.gz \| awk \'\$3\~\/chr\/\' \| cut -f 1\,3 \> $outdir/01.QC/$inf[0].host.temp \n";
	print SH "perl $Bin/remove-host.pl $outdir/01.QC/$inf[0].host.temp $inf[0]\n";
	print SH "cd $outdir/08.MAGsPorfile\n";
	print SH "ln -s $outdir/01.QC/$inf[0].clean.1.fq.gz $inf[0]_1.fastq.gz\n";
	print SH "ln -s $outdir/01.QC/$inf[0].clean.2.fq.gz $inf[0]_2.fastq.gz\n";
	$cleanfastq.=" $inf[0]_1.fastq.gz $inf[0]_2.fastq.gz ";
	#Megahit pipeline
	print SH2 "cd $outdir/03.Assembly\n";
	print SH2 "$megahit -t $thread -1 $outdir/01.QC/$inf[0].clean.1.fq.gz -2 $outdir/01.QC/$inf[0].clean.2.fq.gz --min-contig-len 500 -o $inf[0]-megahit\n";
	print SH2 "perl $Bin/contig-rename.pl $inf[0]-megahit/final.contigs.fa $inf[0].contigs $inf[0]\n";
	print SH2 "perl $Bin/assemle-summary.pl $inf[0].contigs > $inf[0].contigs.summary\n";

	#Binning-Checkm
	open SH6, ">$outdir/shell/S5.BIN.$inf[0].$runid.sh" or die "can not open file: $outdir/shell/S5.BIN.$inf[0].$runid.sh\n"; 
	print SH6 "cd $outdir/06.BIN/ \n";
	print SH6 "ln -s $outdir/03.Assembly/$inf[0].contigs \n";
	print SH6 "$Bin/bwa index $inf[0].contigs\n";
	print SH6 "$Bin/samtools faidx $inf[0].contigs\n";
	print SH6 "$Bin/bwa mem -t $thread -M -R \'\@RG\\tID:$inf[0]\\tSM:$inf[0]\\tLB:$inf[0]\\tPL:Illumina\\tPI:500\' $inf[0].contigs $outdir/01.QC/$inf[0].clip.1.fq.gz $outdir/01.QC/$inf[0].clip.2.fq.gz \| $Bin/samtools view -bS -t $inf[0].contigs.fai - > $inf[0].bam\n";
	print SH6 "$Bin/samtools sort -m 6960000000 $inf[0].bam $inf[0].sort\n";
	print SH6 "$Bin/samtools index $inf[0].sort.bam\n";
	print SH6 "$Bin/metabat/jgi_summarize_bam_contig_depths --outputDepth $inf[0].depth $inf[0].sort.bam\n";
	print SH6 "rm $inf[0].bam\n";
	print SH6 "$Bin/metabat/metabat2 --minContig 1500 -t $thread --inFile $inf[0].contigs --abdFile $inf[0].depth  --outFile $outdir/06.BIN/$inf[0]/$inf[0]\n";
	print SH6 "ls $outdir/06.BIN/$inf[0]/$inf[0]\*fa > $inf[0].bins.list\n";
	print SH6 "cp $outdir/06.BIN/$inf[0]/$inf[0]\*fa $outdir/08.MAGsPorfile/BINs \n";
	print SH6 "source activate metawrap\ncheckm lineage_wf $outdir/06.BIN/$inf[0]/ $outdir/06.BIN/$inf[0]_checkm -x fa -t $thread --pplacer_threads 4 --tab_table -f $inf[0].checkm.summary\n";

	#gene predict
	open SH7, ">$outdir/shell/S6.predict.$inf[0].$runid.sh" or die "can not open file: $outdir/shell/S6.predict.$inf[0].$runid.sh\n";
	print SH7 "cd $outdir/04.Predict\n";
	print SH7 "$Bin/prodigal -a $inf[0].temp.orf.faa -i $outdir/03.Assembly/$inf[0].contigs -f gff -o $inf[0].gff -p meta -q -d $inf[0].temp.orf.ffn\n";
	print SH7 "perl $Bin/gene-filter.pl $inf[0].temp.orf.ffn $inf[0].temp.orf.faa $inf[0].orf.ffn $inf[0].orf.faa \n";
	print SH7 "rm $inf[0].temp.orf.faa $inf[0].temp.orf.ffn \n";

	#gene cd-hit shell
	print CD "cat $outdir/04.Predict/$inf[0].orf.ffn >> $runid.all.ffn\n";
	
	#bwa mapping geneset
	open SH3, ">$outdir/shell/S8.Abund.$inf[0].$runid.sh" or die "S8.Abund.$inf[0].$runid.sh\n"; 
	print SH3 "cd $outdir/07.Profile \n";
	print SH3 "$Bin/bwa mem -t $thread -M -R \'\@RG\\tID:$inf[0]\\tSM:$inf[0]\\tLB:$inf[0]\\tPL:Illumina\\tPI:500\' $outdir/04.PredictGeneSet/$runid.geneSet.ffn $outdir/01.QC/$inf[0].clean.1.fq.gz $outdir/01.QC/$inf[0].clean.2.fq.gz  > $inf[0].sam\n";
	print SH3 "perl $Bin/uniqGene-ProfilebyBwa.pl $outdir/04.PredictGeneSet/$runid.geneSet.ffn $inf[0].sam $inf[0].sam.abd\n";
	print SH4 "$inf[0]\t$inf[0].sam.abd\n";
	close (SH,SH2,SH3,SH4,SH6,SH7);
}

#cd-hit
print CD "$Bin/cd-hit-est -i $runid.all.ffn -o $runid.geneSet.ffn -n 9 -c 0.95 -G 0 -M 0 -d 0 -aS 0.9 -r 1 -T 80\n";
print CD "$Bin/transeq -sequence $runid.geneSet.ffn -table 11 -trim -outseq $runid.geneSet.faa\n";
print CD "perl /lustre/sdb/zengl/bin/module/Meta/RefLeng.pl $runid.geneSet.faa\n";
print CD "$bwa index  $runid.geneSet.ffn\n";
print CD "$Bin/samtools faidx $runid.geneSet.ffn\n";
print CD "rm $runid.all.ffn\n";

#annotation
print ANNO "cd $outdir/05.Anno\n";
print ANNO "mkdir -p KEGG NR COG CAZy\n";
#NR
print ANNO "cd $outdir/05.Anno/NR/\n";
print ANNO "/lustre/sdb/zengl/bin/module/software/diamond blastp -d /mnt/data/nr/nr -q $outdir/04.PredictGeneSet/$runid.geneSet.faa -o uniqGeneSet.nr.m8 -f 6 --evalue 0.00001 -k 10 -t ./ -b 12\n";
print ANNO "perl /lustre/sdb/zengl/bin/module/best_m8.pl uniqGeneSet.nr.m8 >uniqGeneSet.nr.m8.fil\n";
print ANNO "perl /lustre/sdb/zengl/bin/module/nr_anno_overlap_func.pl uniqGeneSet.nr.m8.fil /lustre/sdb/Database/nr/gi.nr.func.xls >uniqGeneSet.nr.m8.fil.func.xls\n";
print ANNO "perl /lustre/sdb/zengl/bin/module/Meta/4-5.tax_byM8.pl uniqGeneSet.nr.m8 /lustre/sdb/zengl/bin/module/Meta/gi_2tax.prot.xls.gz nr.tax.xls\n";
print ANNO "perl /lustre/sdb/zengl/bin/module/Meta/4-6.tax_profile.pl -i nr.tax.xls -p $outdir/07.Profile/$runid.geneabundance.txt -o nr.tax.profile.xls\n";
print ANNO "perl /lustre/sdb/zengl/bin/module/Meta/4-7.taxLevel_profile.pl -i nr.tax.profile.xls -o summary\n";
print ANNO "cd summary\n";
print ANNO "perl /lustre/sdb/zengl/bin/module/Meta/plot-heatmap.pl -i phylum.xls -o heatmap_phylum.pdf -rtop 100 -ct 1 -rt 1 -clc 0.8 -rlc 0.7 -slas 2 -col darkblue-darkgreen-yellow-darkred -marble 3-0-0-5 -h 14 \n";
print ANNO "perl /lustre/sdb/zengl/bin/module/Meta/plot-heatmap.pl -i genus.xls  -o heatmap_genus.pdf  -rtop 100 -ct 1 -rt 1 -clc 0.8 -rlc 0.7 -slas 2 -col darkblue-darkgreen-yellow-darkred -marble 3-0-0-5 -h 14\n";
#KEGG
print ANNO "cd $outdir/05.Anno/KEGG/\n";
print ANNO "/lustre/sdb/zengl/bin/module/software/diamond blastp -d /mnt/data/seq_pep/ko -q $outdir/04.PredictGeneSet/$runid.geneSet.faa -o $runid.uniqGeneSet.faa.m8 -f 6 --evalue 0.00001 -k 10 -t ./ -b 12\n";
print ANNO "source activate py2.7\n";
print ANNO "perl /lustre/sdb/zengl/bin/module/Meta/4-8.kobas.pl -i $runid.uniqGeneSet.faa.m8 -o ./ \n";
print ANNO "source deactivate py2.7\n";
print ANNO "/lustre/sdb/zengl/tools/localperl/bin/perl /lustre/sdb/zengl/bin/module/Meta/4-9.ko_anno.pl kobas.anno.list ./\n";
print ANNO "/lustre/sdb/zengl/tools/localperl/bin/perl /lustre/sdb/zengl/bin/module/Meta/4-10.ko_profile.pl kegg.category.xls $outdir/07.Profile/$runid.geneabundance.txt kegg.profile\n";
print ANNO "/lustre/sdb/zengl/tools/localperl/bin/perl /lustre/sdb/zengl/bin/module/Anno/getKEGGfromBlast2.pl -i pathway.txt -format kobas -o pathways\n";
print ANNO "rm $runid.geneSet.faa.old.m8\n";
#COG
print ANNO "cd $outdir/05.Anno/COG\n";
print ANNO "source activate py2.7\n";
print ANNO "python /lustre/sdb/zengl/bin/module/software/eggnog-mapper/emapper.py -i $outdir/04.PredictGeneSet/$runid.geneSet.faa --output uniqGeneSet.egg -m diamond --cpu $thread\n";
print ANNO "source deactivate py2.7\n";
print ANNO "perl /lustre/sdb/zengl/bin/module/Meta/4-2.eggNOG_anno.pl -i uniqGeneSet.egg.emapper.annotations -o uniqGeneSet.egg.xls\n";
print ANNO "perl /lustre/sdb/zengl/bin/module/Meta/4-3.cog_profile.pl -i uniqGeneSet.egg.xls -p $outdir/07.Profile/$runid.geneabundance.txt -o eggNOG.profile\n";
print ANNO "perl /lustre/sdb/zengl/bin/module/Meta/4-4.cog_bar_eachSam.pl eggNOG.profile.function.xls bar\n";
#CAZy
print ANNO "cd $outdir/05.Anno/CAZy/\n";
print ANNO "/lustre/sdb/zengl/tools/localperl/bin/perl /lustre/sdb/zengl/bin/module/Anno/4-10.cazy_anno.pl $outdir/04.PredictGeneSet/$runid.geneSet.faa cazy /lustre/sdb/zengl/bin/module/Anno/CAZyme/V5/\n";
print ANNO "/lustre/sdb/zengl/tools/localperl/bin/perl /lustre/sdb/zengl/bin/module/Anno/st_CAZy_anno.pl  cazy.anno.xls >cazy.anno.st.xls\n";
print ANNO "perl /lustre/sdb/zengl/bin/module/Meta/4-15.cazy_profile.pl -i cazy.anno.xls -p $outdir/07.Profile/$runid.geneabundance.txt -o cazy.profile\n";

#Merge abuandance profile
print SH5 "cd $outdir/07.Profile\n";
print SH5 "perl $Bin/Profile-Merge-BAM.pl bam.list $runid.geneabundance.txt\n";

close (OA,OB,CD,IN,SH5,ANNO);

#merge sample cmd
open SH, ">$outdir/shell/S1.STAT.$runid.sh" or die "can not open file: $outdir/shell/S1.STAT.$runid.sh\n";
print SH "cd $outdir/stat/\n";
print SH "$java -jar $fqstat -i $outdir/$fqlist > $outdir/stat/raw.$runid.stat.xls\n";
print SH "$java -jar $fqstat -i $outdir/01.QC/clean.fq.list > $outdir/stat/clean.$runid.stat.xls\n";
print SH "$java -jar $fqstat -i $outdir/01.QC/clip.fq.list > $outdir/stat/clip.$runid.stat.xls\n";
close SH;

#MAGs quant calculate
open SH, ">$outdir/shell/S10.MAGsQua.$runid.sh" or die "can not open file: $outdir/shell/S10.MAGsQua.$runid.sh\n";
print SH "cd $outdir/08.MAGsPorfile	\n";
print SH "source activate metawrap\n";
print SH "metawrap quant_bins -b BINs/ -o ALL/ -t $thread $cleanfastq\n";
close SH;

