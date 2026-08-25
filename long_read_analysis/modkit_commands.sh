#!/bin/bash

mkdir bedMethyl_files
mkdir modkit_calls

# Create bedMethyl file for CG for chlamy_8A
modkit pileup -t 10 --preset traditional --with-header \
	-r ../genome_files/CC2937_T2T.fa \
	bam_files/chlamy_8A_5mCG_5hmCG_sorted.bam \
	- | gzip > bedMethyl_files/chlamy_8A_5mCG_5hmCG.bedMethyl.gz

# Create bedMethyl file for CG for chlamy_7A
modkit pileup -t 10 --preset traditional --with-header \
	-r ../genome_files/CC2937_T2T.fa \
	bam_files/chlamy_7A_5mCG_5hmCG_sorted.bam \
	- | gzip > bedMethyl_files/chlamy_7A_5mCG_5hmCG.bedMethyl.gz

# Create bedMethyl file for CG for chlamy_15A
modkit pileup -t 10 --preset traditional --with-header \
	-r ../genome_files/CC2937_T2T.fa \
	bam_files/chlamy_15A_5mCG_5hmCG_sorted.bam \
	- | gzip > bedMethyl_files/chlamy_15A_5mCG_5hmCG.bedMethyl.gz

# Create bedMethyl file for CG for chlamy_18A
modkit pileup -t 10 --preset traditional --with-header \
	-r ../genome_files/CC2937_T2T.fa \
	bam_files/chlamy_18A_5mCG_5hmCG_sorted.bam \
	- | gzip > bedMethyl_files/chlamy_18A_5mCG_5hmCG.bedMethyl.gz
	
# Create bedMethyl file for AT for chlamy_8A 
modkit pileup -t 10 --motif AT 0 --mod-thresholds a:0.995 --combine-strands --with-header \
	-r ../genome_files/CC2937_T2T.fa \
	bam_files/chlamy_8A_6mA_sorted.bam \
	- | gzip > bedMethyl_files/chlamy_8A_6mA.bedMethyl.gz

# Create bedMethyl file for AT for chlamy_7A 
modkit pileup -t 10 --motif AT 0 --mod-thresholds a:0.995 --combine-strands --with-header \
	-r ../genome_files/CC2937_T2T.fa \
	bam_files/chlamy_7A_6mA_sorted.bam \
	- | gzip > bedMethyl_files/chlamy_7A_6mA.bedMethyl.gz

# Create bedMethyl file for AT for chlamy_15A 
modkit pileup -t 10 --motif AT 0 --mod-thresholds a:0.995 --combine-strands --with-header \
	-r ../genome_files/CC2937_T2T.fa \
	bam_files/chlamy_15A_6mA_sorted.bam \
	- | gzip > bedMethyl_files/chlamy_15A_6mA.bedMethyl.gz

# Create bedMethyl file for AT for chlamy_18A 
modkit pileup -t 10 --motif AT 0 --mod-thresholds a:0.995 --combine-strands --with-header \
	-r ../genome_files/CC2937_T2T.fa \
	bam_files/chlamy_18A_6mA_sorted.bam \
	- | gzip > bedMethyl_files/chlamy_18A_6mA.bedMethyl.gz
	
# Extract per-read CG calls for chromosome 15 from chlamy_15A
modkit extract calls \
	--reference ../genome_files/CC2937_T2T.fa \
	--region chromosome_15 \
	-t 10 --cpg --pass-only --mapped-only --bgzf \
	bam_files/chlamy_15A_5mCG_5hmCG_sorted.bam \
	modkit_calls/chlamy_15A_5mCG_chr15_modkit_calls.tsv.gz

# Extract per-read AT calls for chromosome 15 from chlamy_15A
modkit extract calls \
	--reference ../genome_files/CC2937_T2T.fa \
	--region chromosome_15 \
	-t 10 --motif AT 0 --pass-only --mapped-only --bgzf --filter-threshold A:0.995 \
	bam_files/chlamy_15A_6mA_sorted.bam \
	modkit_calls/chlamy_15A_6mAT_chr15_modkit_calls.tsv.gz
	
# Extract per-read CG calls for chromosome 15 from Run3
modkit extract calls \
	--reference ../genome_files/CC2937_T2T.fa \
	--region chromosome_15 \
	-t 10 --cpg --pass-only --mapped-only --bgzf \
	bam_files/ThirdRun_GEVEs_5mCG_calls.sorted_remap_CC2937_T2T.bam \
	modkit_calls/ThirdRun_GEVEs_5mCG_modkit_calls.tsv.gz

# Extract per-read AT calls for chromosome 15 from Run3
modkit extract calls \
	--reference ../genome_files/CC2937_T2T.fa \
	--region chromosome_15 \
	-t 10 --motif AT 0 --pass-only --mapped-only --bgzf --filter-threshold A:0.995 \
	bam_files/ThirdRun_GEVEs_6mA_calls.sorted_remap_CC2937_T2T.bam \
	modkit_calls/ThirdRun_GEVEs_6mAT_modkit_calls.tsv.gz
