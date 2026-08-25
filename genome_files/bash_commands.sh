#!/bin/bash

# Create a BED file with the location of the GEVE and relicts in the CC2937 T2T assembly
seqkit subseq --bed Erazo_2025_GEVE_and_relicts.bed Erazo_2025.fasta > GEVE_and_relicts.fa
minimap2 -x asm5 -secondary=no -a CC2937_T2T.fa GEVE_and_relicts.fa | samtools view -o GEVE_and_relicts_T2T.bam
bedtools bamtobed -i GEVE_and_relicts_T2T.bam | cut -f 1-3 > GEVE_and_relicts_T2T.bed
grep chromosome_15 GEVE_and_relicts_T2T.bed > GEVE_T2T.bed

# Liftover annotation to CC2937 T2T assembly
conda create -n liftoff -c conda-forge -c bioconda --strict-channel-priority python=3.10 liftoff=1.6.3 libsqlite=3.48.0
conda run -n liftoff liftoff CC2937_T2T.fa Erazo_2025.fasta -g updated_chlamy_annotation.gff3 -p 10 -polish -o updated_chlamy_annotation_CC2937_T2T_liftoff.gff3
gffread updated_chlamy_annotation_CC2937_T2T_liftoff.gff3_polished -T --keep-genes -o updated_chlamy_annotation_CC2937_T2T_liftoff_polished.gtf

# Calculate depth across the GEVE
samtools depth -a -H -@ 15 \
	-b GEVE_T2T.bed \
	../long_read_analysis/bam_files/chlamy_7A_5mCG_5hmCG_sorted_remap_CC2937_T2T.bam \
	../long_read_analysis/bam_files/chlamy_8A_5mCG_5hmCG_sorted_remap_CC2937_T2T.bam \
	../long_read_analysis/bam_files/chlamy_15A_5mCG_5hmCG_sorted_remap_CC2937_T2T.bam \
	../long_read_analysis/bam_files/chlamy_18A_5mCG_5hmCG_sorted_remap_CC2937_T2T.bam > geve_depth.txt
	
# Create a chain to liftover from the contig assembly for the original runs to the T2T assembly
minimap2 -cx asm5 Erazo_2025.fasta CC2937_T2T.fa > Erazo_2025_to_T2T_aln.paf
paf2chain -i Erazo_2025_to_T2T_aln.paf > Erazo_2025_to_T2T.chain
sed -i '/chain/ s/\t/ /g' Erazo_2025_to_T2T.chain

# Create indices for CC2937 T2T assembly
STAR \
  --runThreadN 16 \
  --runMode genomeGenerate \
  --genomeFastaFiles CC2937_T2T.fa \
  --sjdbGTFfile updated_chlamy_annotation_CC2937_T2T_liftoff_polished.gtf \
  --sjdbOverhang 149 \
  --genomeDir CC2937_T2T_star_index 
  
bowtie2-build CC2937_T2T.fa CC2937_T2T_bowtie2_index/index

cat CC2937_T2T.fa chrP_and_chrL.fa > CC2937_T2T_chrP_and_chrL.fa
bs_seeker2-build.py -f CC2937_T2T_chrP_and_chrL.fa --aligner=bowtie2 -d CC2937_T2T_chrP_and_chrL_bsseeker_index

TElocal_indexer.py --afile CC2937_T2T.TEs_v3_8.bed.gtf --itype TE
