# GFusion
GFusion is a software package to detect fusion genes using RNA-Seq data

Download the latest version of GFusion. Version: 1.0.0. 

Authors: Qi Chen, Xiaofeng Song (xfsong@nuaa.edu.cn)

Usage:

perl GFusion.pl -o out_file -r 0 -p 6 -i /path/to/bowtie1_index -g /path/to/genes.gtf -1 PE_reads_1.fastq -2 PE_reads_2.fastq 

Option

-o: output file directory;

-r: the expected (mean) inner distance between mate pairs(default = 0);

-n: the minimum number of split reads(default = 1);

-p: threads (default = 4);

-i: genome index (use bowtie index);

-g: gene annotation file(*.gtf file);

-1: pair-end reads file, PE_reads_1.fq;

-2: pair-end reads file, PE_reads_2.fq;

-f: single-end reads file, SE_reads.fq;

-help: usage help;

-d: the minimum length of adjacent genes (default = 50000); 

-L: the anchor length (more than 20bp, default = 0.4*read_length);

Note:

The following files (genome index file, gene annotation gtf file) could be download from TopHat. Website: http://ccb.jhu.edu/software/tophat/igenomes.shtml 

Example:
perl GFusion.pl -o output -r 0 -p 12 -i /path/hg19 -g /path/genes.gtf -1 /path/test_1.fastq -2 /path/test_2.fastq

Example files download:

1. The alignment index file (-index): BowtieIndex files (Homo_sapiens, hg19 version) can be downloaded from http://ccb.jhu.edu/software/tophat/igenomes.shtml;

2. The genome annotation file (-gtf): GTF file can be downloaded from http://ccb.jhu.edu/software/tophat/igenomes.shtml;

3. Input file :The fastq format file: test_1.fastq, test_2.fastq can be downloaded from http://bioinfo.nuaa.edu.cn/GFusion/test_fastq.zip 

4. Output file: fusion_gene_list.txt can be downloaded from http://bioinfo.nuaa.edu.cn/GFusion/fusion_gene_list.txt

Software Prerequisites:
The following three software should be installed in your cluster or computer before running the GFusion.pl

TopHat: a fast splice junction mapper for RNA-Seq reads;

Bowtie1: an tool for aligning sequencing reads to long reference sequences. Download and extract Bowtie 1 releases;

Samtools: utilities for the Sequence Alignment/Map (SAM) format (Version: 0.1.19 recommended); Note: new version samtools cannot be used here.

Bioperl: a community effort to produce Perl code which is useful in biology.

Output file format:

GFusion will report the fusion genes in the output/result file. The information of fusion genes is listed as follow format:

1. 5’ gene name 

2. 5’ chromosome 

3. 5’ position

4. 3’ gene name 

5. 3’ chromosome 

6. 3’ position 

7. strand 

8. the number of split reads 

9. the number of spanning reads 

Reference:
Qi Chen, Ping Han, Xiaofeng Song. GFusion: a Novel Algorithm to Identify Fusion Genes from Cancer RNA-Seq Data. (under review in peer-reviewed journal)

Copyright (C) 2016 Xiaofeng Song.

