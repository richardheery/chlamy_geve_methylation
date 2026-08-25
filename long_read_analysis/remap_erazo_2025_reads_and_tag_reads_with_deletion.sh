#!/bin/bash

### Remap reads to CC2937-T2T assembly

cd bam_files

for i in *5mCG_calls*.bam; do
	sample=$(basename "$i" .bam)
	samtools fastq -T MM,ML "${i}" | \
	minimap2 -y -ax map-ont -t 4 "../genome_files/CC2937_T2T.fa" - | \
	samtools view -u - | samtools sort -@ 4 --write-index -o "${sample}_remap_CC2937_T2T.bam"
done

### Filter BAM files for reads overlapping 100 bp at the 3' end of the deletion and tag reads based on whether they support the deletion or not

# Create BED files for the 50-bp regions immediately upstream and downstream of the 3' deletion breakpoint
echo -e "chromosome_15\t1499435\t1499485" > upstream_region_3_prime_deletion_breakpoint.bed
echo -e "chromosome_15\t1499485\t1499535" > downstream_region_3_prime_deletion_breakpoint.bed

### First run

# Extract reads overlapping the region downstream of the deletion which do not overlap the deletion and tag them
bedtools intersect -a FirstRun_GEVEs_5mCG_calls.sorted_remap_CC2937_T2T.bam -b downstream_region_3_prime_deletion_breakpoint.bed | \
	bedtools intersect -a - -b upstream_region_3_prime_deletion_breakpoint.bed -v | \
	samtools view -h | \
	awk 'BEGIN {OFS="\t"} /^@/ {print; next} {print $0, "HP:i:1"}' | samtools sort -o FirstRun_GEVEs_5mCG_calls_downstream_deletion1.sorted_remap_CC2937_T2T.bam --write-index

# Find reads overlapping the deletion and tag them
samtools view -h -L upstream_region_3_prime_deletion_breakpoint.bed FirstRun_GEVEs_5mCG_calls.sorted_remap_CC2937_T2T.bam  | \
	awk 'BEGIN {OFS="\t"} /^@/ {print; next} {print $0, "HP:i:2"}' | samtools sort -o FirstRun_GEVEs_5mCG_calls_downstream_deletion2.sorted_remap_CC2937_T2T.bam --write-index

# Merge both sets of reads 
samtools merge FirstRun_GEVEs_5mCG_calls_downstream_deletion_tagged.bam \
	FirstRun_GEVEs_5mCG_calls_downstream_deletion1.sorted_remap_CC2937_T2T.bam \
	FirstRun_GEVEs_5mCG_calls_downstream_deletion2.sorted_remap_CC2937_T2T.bam \
	--write-index -f

### Second run

# Extract reads overlapping the region downstream of the deletion which do not overlap the deletion and tag them
bedtools intersect -a SecondRun_GEVEs_5mCG_calls.sorted_remap_CC2937_T2T.bam -b downstream_region_3_prime_deletion_breakpoint.bed | \
	bedtools intersect -a - -b upstream_region_3_prime_deletion_breakpoint.bed -v | \
	samtools view -h | \
	awk 'BEGIN {OFS="\t"} /^@/ {print; next} {print $0, "HP:i:1"}' | samtools sort -o SecondRun_GEVEs_5mCG_calls_downstream_deletion1.sorted_remap_CC2937_T2T.bam --write-index

# Find reads overlapping the deletion and tag them
samtools view -h -L upstream_region_3_prime_deletion_breakpoint.bed SecondRun_GEVEs_5mCG_calls.sorted_remap_CC2937_T2T.bam  | \
	awk 'BEGIN {OFS="\t"} /^@/ {print; next} {print $0, "HP:i:2"}' | samtools sort -o SecondRun_GEVEs_5mCG_calls_downstream_deletion2.sorted_remap_CC2937_T2T.bam --write-index

# Merge both sets of reads 
samtools merge SecondRun_GEVEs_5mCG_calls_downstream_deletion_tagged.bam \
	SecondRun_GEVEs_5mCG_calls_downstream_deletion1.sorted_remap_CC2937_T2T.bam \
	SecondRun_GEVEs_5mCG_calls_downstream_deletion2.sorted_remap_CC2937_T2T.bam \
	--write-index -f
	
### Third run

# Extract reads overlapping the region downstream of the deletion which do not overlap the deletion and tag them
bedtools intersect -a ThirdRun_GEVEs_5mCG_calls.sorted_remap_CC2937_T2T.bam -b downstream_region_3_prime_deletion_breakpoint.bed | \
	bedtools intersect -a - -b upstream_region_3_prime_deletion_breakpoint.bed -v | \
	samtools view -h | \
	awk 'BEGIN {OFS="\t"} /^@/ {print; next} {print $0, "HP:i:1"}' | samtools sort -o ThirdRun_GEVEs_5mCG_calls_downstream_deletion1.sorted_remap_CC2937_T2T.bam --write-index

# Find reads overlapping the deletion and tag them
samtools view -h -L upstream_region_3_prime_deletion_breakpoint.bed ThirdRun_GEVEs_5mCG_calls.sorted_remap_CC2937_T2T.bam  | \
	awk 'BEGIN {OFS="\t"} /^@/ {print; next} {print $0, "HP:i:2"}' | samtools sort -o ThirdRun_GEVEs_5mCG_calls_downstream_deletion2.sorted_remap_CC2937_T2T.bam --write-index

# Merge both sets of reads, subsampling 25% of the reads
samtools merge - \
	ThirdRun_GEVEs_5mCG_calls_downstream_deletion1.sorted_remap_CC2937_T2T.bam \
	ThirdRun_GEVEs_5mCG_calls_downstream_deletion2.sorted_remap_CC2937_T2T.bam | \
	samtools view -s 123.25 -o ThirdRun_GEVEs_5mCG_calls_downstream_deletion_tagged.bam  --write-index

# Remove temporary files
rm FirstRun_GEVEs_5mCG_calls_downstream_deletion1.sorted_remap_CC2937_T2T.bam* FirstRun_GEVEs_5mCG_calls_downstream_deletion2.sorted_remap_CC2937_T2T.bam*
rm SecondRun_GEVEs_5mCG_calls_downstream_deletion1.sorted_remap_CC2937_T2T.bam* SecondRun_GEVEs_5mCG_calls_downstream_deletion2.sorted_remap_CC2937_T2T.bam*
rm ThirdRun_GEVEs_5mCG_calls_downstream_deletion1.sorted_remap_CC2937_T2T.bam* ThirdRun_GEVEs_5mCG_calls_downstream_deletion2.sorted_remap_CC2937_T2T.bam*
