#!/bin/bash
set -e

# Set parameters
reads1_f="$1"
reads2_f="$2"
reads1_r="$3"
reads2_r="$4"
name="$5"
genome="$6"
index="$7"
nthreads="$8"

# Check that all packages are available
command -v bs_seeker2-align.py 
command -v fastp 
command -v bbmerge.sh 
command -v sambamba 
command -v samtools 
command -v bam 
command -v cgmaptools

# Check that genome contains chrL and chrP
if ! zgrep -q ">chrL" $genome; then
	echo "Error: chrL not found in $genome" >&2
	exit 1
fi

if ! zgrep -q ">chrP" $genome; then
	echo "Error: chrP not found in $genome" >&2
	exit 1
fi

# Combine FASTQ files for forward and reverse reads
echo Combining reads...
cat $reads1_f $reads2_f > "${name}_combined_f.fastq.gz"
cat $reads1_r $reads2_r > "${name}_combined_r.fastq.gz"

# Trim reads
echo Trimming reads...
fastp -i "${name}_combined_f.fastq.gz" -o "${name}_f_trimmed.fastq.gz" \
	-I "${name}_combined_r.fastq.gz" -O "${name}_r_trimmed.fastq.gz" \
	--html "${name}_fastp_report.html" --json /dev/null \
	-w $nthreads 
rm "${name}_combined_f.fastq.gz" "${name}_combined_r.fastq.gz"

# Merge overlapping read pairs 
echo Merging overlapping reads...
bbmerge.sh in1="${name}_f_trimmed.fastq.gz" \
	in2="${name}_r_trimmed.fastq.gz" \
	qtrim=r -Xmx10g \
	out="${name}_merged.fastq.gz" outu1="${name}_f_unmerged.fastq.gz" outu2="${name}_r_unmerged.fastq.gz" 
rm "${name}_f_trimmed.fastq.gz" "${name}_r_trimmed.fastq.gz"

# Paired-end alignment of unmerged reads
echo Aligning unmerged reads...
bs_seeker2-align.py \
	--aligner=bowtie2 --bt2--end-to-end --bt2-p $nthreads -e 300 -X 3000 -m 4 \
	-1 "${name}_f_unmerged.fastq.gz" -2 "${name}_r_unmerged.fastq.gz" -o "${name}_unmerged.bam" \
	-d $index \
	-g $genome
rm "${name}_f_unmerged.fastq.gz" "${name}_r_unmerged.fastq.gz"

# Single-end alignment of merged reads
echo Aligning merged reads...
bs_seeker2-align.py \
	--aligner=bowtie2 --bt2--end-to-end --bt2-p $nthreads -e 400 -m 6 \
	-i "${name}_merged.fastq.gz" -o "${name}_merged.bam" \
	-d $index \
	-g $genome
rm "${name}_merged.fastq.gz"

# Sort the output bam files 
echo Sorting BAMs...
sambamba sort -t $nthreads "${name}_unmerged.bam" -o "${name}_unmerged.sorted.bam"
rm "${name}_unmerged.bam" 
sambamba sort -t $nthreads "${name}_merged.bam" -o "${name}_merged.sorted.bam"
rm "${name}_merged.bam"

# Softclip any remaining overlaps 
echo Clipping overlapping reads in unmerged BAM...
bam clipOverlap --in "${name}_unmerged.sorted.bam"  --out "${name}_unmerged.sorted.clipped.bam" 
rm "${name}_unmerged.sorted.bam" "${name}_unmerged.sorted.bam.bai"

# Merge both bam files
echo Merging BAM files...
sambamba merge -t $nthreads "${name}_combined.bam" "${name}_merged.sorted.bam" "${name}_unmerged.sorted.clipped.bam" 
rm "${name}_merged.sorted.bam" "${name}_merged.sorted.bam.bai" "${name}_unmerged.sorted.clipped.bam"

# Mark duplicates 
echo Marking duplicates...
sambamba markdup -r -t $nthreads "${name}_combined.bam" "${name}_combined.dedup.bam"
rm "${name}_combined.bam" "${name}_combined.bam.bai"

# Create md5 sum for output bam file
echo Creating MD5...
md5sum "${name}_combined.dedup.bam" > "${name}_combined.dedup.bam.md5" 

# Create a CGmap file from deduplicated BAM
echo Creating CGmap file from BAM...
cgmaptools convert bam2cgmap --bam "${name}_combined.dedup.bam" --genome $genome -o ${name}
	
echo Done

# Perform checks of methylation levels
echo Performing checks of methylation levels...

summary=$(zcat "${name}.CGmap.gz" | awk '$4 == "CG" && $1 != "chrL" && $1 != "chrP" && $1 != "chrM"' | awk '{cov+=$8; mC+=$7;} END{print "mCG", cov, mC, mC/cov*100;}')
echo $name $summary > "${name}_mCG_global_stats.txt"

summary=$(zcat "${name}.CGmap.gz" | awk '$1 == "chrL"' | awk '{cov+=$8; mC+=$7;} END{print "mCG", cov, mC, mC/cov*100;}')
echo $name $summary > "${name}_lambda_global_stats.txt" 

summary=$(zcat "${name}.CGmap.gz" | awk '$1 == "chrP" && $4 == "CG"' | awk '{cov+=$8; mC+=$7;} END{print "mCG", cov, mC, mC/cov*100;}')
echo $name $summary > "${name}_pUC19_mCG_global_stats.txt"
