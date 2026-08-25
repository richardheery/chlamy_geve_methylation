#!/bin/bash

mkdir -p bam_files coverage_bigwig_files normalized_bigwigs

### Map single-end FASTQ files in fastq_files/PRJNA255778 with Bowtie2 to CC2937_T2T.fa
for fq in fastq_files/PRJNA255778/*.fastq.gz; do

	sample=$(basename "$fq" .fastq.gz)
	output_bam=bam_files/${sample}_sorted.bam
	echo "Mapping $sample"

  bowtie2 \
    --end-to-end \
    --very-sensitive \
    -k 1 \
    -p 10 \
    -x CC2937_T2T_bowtie2_index/index \
    -U "$fq" | \
    samtools sort -@ 10 --write-index -o "$output_bam"
  
done

### Map paired-end FASTQ files in fastq_files/PRJNA1029152 with Bowtie2 to CC2937_T2T.fa
while read -r fq1 fq2; do

	sample=$(basename "$fq1" _1.fastq.gz)
	output_bam=bam_files/${sample}_sorted.bam
	echo "Mapping $sample"

  bowtie2 \
    --end-to-end \
    --very-sensitive \
    -k 1 \
    -p 10 \
    -x CC2937_T2T_bowtie2_index/index \
    -1 "$fq1" -2 "$fq2" | \
    samtools sort -@ 10 --write-index -o "$output_bam" 

done < fastq_files/PRJNA1029152/fastq_pairs.txt

### Convert BAM files to CPM-normalized bigWig files
for bam in bam_files/*.bam; do

	sample=$(basename "$bam" .bam)
	bw=coverage_bigwig_files/${sample}.bw

	bamCoverage -b $bam --normalizeUsing CPM --numberOfProcessors 10 -o $bw
  
done

### Generate log2 ChIP/input-normalized bigWig files
while read -r chip input; do

    base=$(basename "$chip" .bw)
    out="normalized_bigwigs/${base}_log2_vs_input.bw"

    echo "Processing $chip"

    bigwigCompare \
        -b1 "$chip" \
        -b2 "$input" \
        --operation log2 \
        --pseudocount 1 \
        --skipZeroOverZero \
        --binSize 50 \
        --numberOfProcessors 10 \
        -o "$out"

done < coverage_bigwig_files/bigwig_pairs.txt
