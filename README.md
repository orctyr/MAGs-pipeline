MAGs-pipeline

This website contains the MAG construction pipeline for RMGMC Metagenome Project.
Including:
(1) Reads quality control and low-quality reads timming
(2) metagenome assembly using two tools
(3) merge multiple sources of assembles
(4) SNP, Indels correction by re-sequencing methods
(5) MAGs construction by metabat2, maxbin and concoct
(6) Integrating MAGs by DAStools 
(7) MAGs quality control and dereplication
(8) Taxonomy identification
(9) maximum-likelihood tree construction based on MAGs

1. Package Required
	(1) Trimmomatic(v.0.33), http://www.usadellab.org/cms/index.php?page=trimmomatic
	(2) BWA(v.0.7.17), https://github.com/lh3/bwa
	(3) MEGAHIT(v.1.1.1), https://github.com/voutcn/megahit
	(4) IDBA-UD(v.1.1.3), https://www.psc.edu/user-resources/software/idba-ud
	(5) AMOS(v.3.1.0), http://amos.sourceforge.net/wiki/index.php/AMOS
	(6) SAMtools(v.1.9), http://samtools.sourceforge.net/
	(7) MaxBin(v.2.2.4, under metaWRAP)
	(8) MetaBAT2(v.2.11.1 under metaWRAP)
	(9) CONCOCT(v.0.4.0 under metaWRAP)
	(10) DAS Tool(v.1.1.1), https://github.com/cmks/DAS_Tool
	(11) CheckM(v.1.0.7 under metaWRAP)
	(12) metaWRAP(v.1.3), https://github.com/bxlab/metaWRAP
	(13) GTDB-Tk(v.0.1.6), https://ecogenomics.github.io/GTDBTk
	(14) PhyloPhlAn(v.1.0), http://segatalab.cibio.unitn.it/tools/phylophlan/
	(15) Prodigal(v.2.6.3), https://github.com/hyattpd/Prodigal

2. Step1: shell script generation
	Use Yak metagenome data as an example.
	perl Meta-Pipeline.pl -fqlist Yak-fqlist.txt -host a -set Yak-Group.txt
	
3. Step2: run pipeline (fold shell/)
	S1.QC*sh: Reads quality control and remove host genomes
	S2.idba*sh: IDBA assembly
	S2.megahit*sh: MEGAHIT assembly
	S3.merge*sh: merge two assemble results from above two steps
	S4.correct1(2).sh: contigs correction
	S5.BIN*sh: run metabat2 as the first binning methods
	S6.predict*sh & S7.GeneSet*sh: used in RMGMC geneset 
	S8.Abund.*sh: gene abundance calculation
	
4. Step3: run DAStools to improve the MAG quality and quantity
	perl Meta-DAStool.pl -input Yak-metabat2.list
	scripts will be generated in shell/ and run all scripts simultaneously
	