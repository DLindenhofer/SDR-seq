


#Downloaded bam file containing mixed NIH3T3 cells from 

export PATH=/Users/dominik.lindenhofer/samtools_installed:$PATH



#Change smaller file

samtools view -h 500_hgmm_3p_LT_Chromium_X_possorted_genome_bam.bam | sed 's/mm10___//' | samtools view -bS - > 500_hgmm_3p_LT_Chromium_X_possorted_genome_bam_edit.bam
samtools index 500_hgmm_3p_LT_Chromium_X_possorted_genome_bam_edit.bam
samtools view -h 500_hgmm_3p_LT_Chromium_X_possorted_genome_bam_edit.bam




