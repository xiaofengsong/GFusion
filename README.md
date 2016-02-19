# GFusion

GFusion is a software package to discovery fusion genes using RNA-Seq data.
The first public release of Gfusion is now available for download. 
Usage:
GFusion.pl #GFusion accepts RNA-Seq data input file of fastq format.
Software Prerequisites:
To use GFusion, you need to install Tophat , Bowtie, Samtools, Bioperl module firstly. 
TopHat: a fast splice junction mapper for RNA-Seq reads;
Samtools: utilities for the Sequence Alignment/Map (SAM) format;
Bowtie1: an tool for aligning sequencing reads to long reference sequences. Download and extract Bowtie 1 releases;
Bioperl:  a community effort to produce Perl code which is useful in biology.
To run GFusion: 	
perl Gfusion.pl -o out_file -r 0 -p 6 -index /path/to/bowtie1_index -gtf /path/to/genes.gtf -1 PE_reads_1.fastq -2 PE_reads_2.fastq 
Options:
-o: output directory;
-r: <int> the expected (mean) inner distance between mate pairs(default = 0);
-n: the minimum number of split reads(default = 1);
-p: <int>threads (default = 4);
-index: genome index (use bowtie index);
-gtf: gene annotation file(*.gtf file);
-1: pair-end reads file, PE_reads_1.fq
-2: pair-end reads file, PE_reads_2.fq;
-fq: single-end reads file, SE_reads.fq;
-help: usage help;
-intron_length: the minimum length of adjacent genes (default = 50000); 
-pos_length: the anchor length (more than 20bp, default = 0.4*read_length);
The following files (genome index, gene annotation) could be download from TopHat. Website: http://ccb.jhu.edu/software/tophat/igenomes.shtml 
Example:
perl Gfusion.pl -o out_file -r 0 -p 6 -index hg19.fa -gtf genes.gtf -1 test_1.fastq -2 test_2.fastq
Example files download:
	The alignment index file (-index): BowtieIndex files (Homo_sapiens, hg19 version) can be downloaded from http://ccb.jhu.edu/software/tophat/igenomes.shtml;
	The genome annotation file (-gtf): gtf file;
	Input file:test_1.fastq, test_2.fastq: test.fastq;
	Output file:the reported fusion genes in the out_file/result file.
5’ gene	5’ chromosome	5’ position	3’ gene	3’ chromosome	3’ position	strand	Split 	Spanning
BCAS4	chr20	49411710	BCAS3	chr17	59445688	ff	43	101

Output file format:
GFusion will report the fusion genes in the out_file/result file. The information of fusion genes is listed as follow format:
1. 5’ gene
2. 5’ chromosome
3. 5’ position
4. 3’ gene
5. 3’ chromosome
6. 3’ position
7. strand
8. the number of split reads
9. the number of spanning reads

