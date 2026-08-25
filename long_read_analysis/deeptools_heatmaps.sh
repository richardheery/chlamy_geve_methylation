# Create heatmaps for 5mC and 6mA across GEVE and host genes

computeMatrix scale-regions \
	-R /home/richardheery/chlamydomonas_project/chlamy_genome/t2t/GEVE_DVP_genes.bed /home/richardheery/chlamydomonas_project/chlamy_genome/t2t/GEVE_BVP_genes.bed \
	-b 3000 -a 3000 -m 5000 --binSize 100 \
	-S $(ls bigWigs/*5mC*) \
	--numberOfProcessors 10 -o plots/geve_gene_5mC_matrix.gz
	
plotHeatmap -m plots/geve_gene_5mC_matrix.gz \
	--regionsLabel DVP BVP \
	--colorMap RdBu_r \
	-o plots/geve_gene_5mC_heatmap.png
	
computeMatrix scale-regions \
	-R /home/richardheery/chlamydomonas_project/chlamy_genome/t2t/GEVE_DVP_genes.bed /home/richardheery/chlamydomonas_project/chlamy_genome/t2t/GEVE_BVP_genes.bed \
	-b 3000 -a 3000 -m 5000 --binSize 100 \
	-S $(ls bigWigs/*6mA*) \
	--numberOfProcessors 10 -o plots/geve_gene_6mA_matrix.gz
	
plotHeatmap -m plots/geve_gene_6mA_matrix.gz \
	--regionsLabel DVP BVP \
	--colorMap RdBu_r \
	-o plots/geve_gene_6mA_heatmap.png
	
computeMatrix scale-regions \
	-R /home/richardheery/chlamydomonas_project/chlamy_genome/t2t/chlamy_host_genes.bed /home/richardheery/chlamydomonas_project/chlamy_genome/t2t/GEVE_all_genes.bed \
	-b 3000 -a 3000 -m 5000 --binSize 100 \
	-S $(ls bigWigs/*5mC*) \
	--numberOfProcessors 10 -o plots/host_vs_geve_gene_5mC_matrix.gz
	
plotHeatmap -m plots/host_vs_geve_gene_5mC_matrix.gz \
	--regionsLabel Host GEVE \
	--colorMap RdBu_r \
	-o plots/host_vs_geve_gene_5mC_heatmap.png
	
computeMatrix scale-regions \
	-R /home/richardheery/chlamydomonas_project/chlamy_genome/t2t/chlamy_host_genes.bed /home/richardheery/chlamydomonas_project/chlamy_genome/t2t/GEVE_all_genes.bed \
	-b 3000 -a 3000 -m 5000 --binSize 100 \
	-S $(ls bigWigs/*6mA*) \
	--numberOfProcessors 10 -o plots/host_vs_geve_gene_6mA_matrix.gz
	
plotHeatmap -m plots/host_vs_geve_gene_6mA_matrix.gz \
	--regionsLabel Host GEVE \
	--colorMap RdBu_r \
	-o plots/host_vs_geve_gene_6mA_heatmap.png
	
computeMatrix scale-regions \
	-R all_mCG_genes.bed \
	-b 3000 -a 3000 -m 5000 --binSize 100 \
	-S $(ls bigWigs/*5mC*) \
	--numberOfProcessors 10 -o plots/all_mCG_gene_5mC_matrix.gz
	
plotHeatmap -m plots/all_mCG_gene_5mC_matrix.gz \
	--kmeans 3 \
	--colorMap RdBu_r \
	-o plots/all_mCG_gene_5mC_heatmap.png	
	
###

computeMatrix scale-regions \
	-R /home/richardheery/chlamydomonas_project/chlamy_genome/t2t/chlamy_host_genes.bed /home/richardheery/chlamydomonas_project/chlamy_genome/t2t/GEVE_all_genes.bed \
	-b 3000 -a 3000 -m 5000 --binSize 100 \
	-S bigWigs/chlamy_15A_5mC.bw bigWigs/chlamy_15A_6mA.bw \
	--numberOfProcessors 10 -o plots/chlamy_15A_geve_vs_host_5mC_6mA_matrix.gz
	
plotHeatmap -m plots/chlamy_15A_geve_vs_host_5mC_6mA_matrix.gz \
	--regionsLabel Host GEVE \
	--colorMap RdBu_r \
	-o final_figures/figure1/chlamy_15A_geve_vs_host_5mC_6mA_heatmap.png
	
computeMatrix scale-regions \
	-R /home/richardheery/chlamydomonas_project/chlamy_genome/t2t/GEVE_DVP_genes.bed /home/richardheery/chlamydomonas_project/chlamy_genome/t2t/GEVE_BVP_genes.bed \
	-b 3000 -a 3000 -m 5000 --binSize 100 \
	-S bigWigs/chlamy_15A_5mC.bw bigWigs/chlamy_15A_6mA.bw \
	--numberOfProcessors 10 -o plots/chlamy_15A_bvp_vs_dvp_5mC_6mA_matrix.gz
	
plotHeatmap -m plots/chlamy_15A_bvp_vs_dvp_5mC_6mA_matrix.gz \
	--regionsLabel DVP BVP \
	--colorMap RdBu_r \
	-o final_figures/figure1/chlamy_15A_bvp_vs_dvp_5mC_6mA_heatmap.png
