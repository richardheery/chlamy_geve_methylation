#!/bin/bash

dorado basecaller \
	sup,5mCG_5hmCG \
	pod5_files/15A \
	--reference ../genome_files/CC2937_T2T.fa | \
	samtools sort -@ 8 --write-index -o bam_files/chlamy_15A_5mCG_5hmCG_sorted.bam

dorado basecaller \
	sup,5mCG_5hmCG \
	pod5_files/18A \
	--reference ../genome_files/CC2937_T2T.fa | \
	samtools sort -@ 8 --write-index -o bam_files/chlamy_18A_5mCG_5hmCG_sorted.bam

dorado basecaller \
	sup,5mCG_5hmCG \
	pod5_files/7A \
	--reference ../genome_files/CC2937_T2T.fa | \
	samtools sort -@ 8 --write-index -o bam_files/chlamy_7A_5mCG_5hmCG_sorted.bam

dorado basecaller \
	sup,5mCG_5hmCG \
	pod5_files/8A \
	--reference ../genome_files/CC2937_T2T.fa | \
	samtools sort -@ 8 --write-index -o bam_files/chlamy_8A_5mCG_5hmCG_sorted.bam

dorado basecaller \
	sup,6mA \
	pod5_files/15A \
	--reference ../genome_files/CC2937_T2T.fa | \
	samtools sort -@ 8 --write-index -o bam_files/chlamy_15A_6mA_sorted.bam

dorado basecaller \
	sup,6mA \
	pod5_files/18A \
	--reference ../genome_files/CC2937_T2T.fa | \
	samtools sort -@ 8 --write-index -o bam_files/chlamy_18A_6mA_sorted.bam

dorado basecaller \
	sup,6mA \
	pod5_files/7A \
	--reference ../genome_files/CC2937_T2T.fa | \
	samtools sort -@ 8 --write-index -o bam_files/chlamy_7A_6mA_sorted.bam

dorado basecaller \
	sup,6mA \
	pod5_files/8A \
	--reference ../genome_files/CC2937_T2T.fa | \
	samtools sort -@ 8 --write-index -o bam_files/chlamy_8A_6mA_sorted.bam
