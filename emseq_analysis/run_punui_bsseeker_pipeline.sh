#!/bin/bash

mkdir cgmap_files
cd cgmap_files

bash bsseeker_pipeline.sh \
	../Punui2_reads/Punui2_EKDL250032980-1A_23C5GMLT4_L2_1.fq.gz /dev/null \
	../Punui2_reads/Punui2_EKDL250032980-1A_23C5GMLT4_L2_2.fq.gz /dev/null \
	Punui2 \
	../../chlamy_genome/CC2937_T2T_chrP_and_chrL.fa \
	../../chlamy_genome/CC2937_T2T_chrP_and_chrL_bsseeker_index \
	8
