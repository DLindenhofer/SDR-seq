
#Install packages

#Required sinto and samtools
pip list
pip install sinto
conda install -c bioconda samtools

#Subsample bam file

#Use package sinto - run on jupyterhup

#gDNA

Outdir_SDRranger=<insert_path>
Outdir_Subset=<insert_path>

#60
sinto filterbarcodes -b $Outdir_SDRranger/gDNA/60/gDNA_with_bc.sorted.bam \
                     -c $Outdir_Subset/SDRranger_gDNA_60/SDR002_60_meta_gDNA_CB.txt \
                     --outdir $Outdir_Subset/SDRranger_gDNA_60

#60fail

sinto filterbarcodes -b $Outdir_SDRranger/gDNA/60fail/gDNA_with_bc.sorted.bam \
                     -c $Outdir_Subset/SDRranger_gDNA_60fail/SDR002_60fail_meta_gDNA_CB.txt \
                     --outdir $Outdir_Subset/SDRranger_gDNA_60fail


#120

sinto filterbarcodes -b $Outdir_SDRranger/gDNA/120/gDNA_with_bc.sorted.bam \
                     -c $Outdir_Subset/SDRranger_gDNA_120/SDR002_120_meta_gDNA_CB.txt \
                     --outdir $Outdir_Subset/SDRranger_gDNA_120


#240

sinto filterbarcodes -b $Outdir_SDRranger/gDNA/240/gDNA_with_bc.sorted.bam \
                     -c $Outdir_Subset/SDRranger_gDNA_240/SDR002_240_meta_gDNA_CB.txt \
                     --outdir $Outdir_Subset/SDRranger_gDNA_240



#Check number of reads


samtools view -c $Outdir_Subset/SDRranger_gDNA_60/SDR002_60.bam 
samtools view -c $Outdir_Subset/SDRranger_gDNA_60fail/SDR002_60fail.bam
samtools view -c $Outdir_Subset/SDRranger_gDNA_120/SDR002_120.bam
samtools view -c $Outdir_Subset/SDRranger_gDNA_240/SDR002_240.bam


#cDNA -

#Need to write both CB and SB into a specific tag for filtering

#Add both to the read name 

#Write back into a tag

#60

samtools view -h $Outdir_SDRranger/cDNA/60/RNA_with_bc_umi.sorted.bam \
  | sinto tagtoname --tag SB -b - \
  | sinto tagtoname --tag CB -b - \
  | sinto nametotag --tag CS --barcode_regex "[A-Za-z]+\.[A-Za-z]+:[A-Za-z]+" -b - \
  | samtools view -b - \
  > $Outdir_Subset/SDRranger_cDNA_60/RNA_with_bc_umi_edit.bam


#60fail

samtools view -h $Outdir_SDRranger/cDNA/60fail/RNA_with_bc_umi.sorted.bam \
  | sinto tagtoname --tag SB -b - \
  | sinto tagtoname --tag CB -b - \
  | sinto nametotag --tag CS --barcode_regex "[A-Za-z]+\.[A-Za-z]+:[A-Za-z]+" -b - \
  | samtools view -b - \
  > $Outdir_Subset/SDRranger_cDNA_60fail/RNA_with_bc_umi_edit.bam





#120

samtools view -h $Outdir_SDRranger/cDNA/120/RNA_with_bc_umi.sorted.bam \
  | sinto tagtoname --tag SB -b - \
  | sinto tagtoname --tag CB -b - \
  | sinto nametotag --tag CS --barcode_regex "[A-Za-z]+\.[A-Za-z]+:[A-Za-z]+" -b - \
  | samtools view -b - \
  > $Outdir_Subset/SDRranger_cDNA_120/RNA_with_bc_umi_edit.bam




#240

samtools view -h $Outdir_SDRranger/cDNA/240/RNA_with_bc_umi.sorted.bam \
  | sinto tagtoname --tag SB -b - \
  | sinto tagtoname --tag CB -b - \
  | sinto nametotag --tag CS --barcode_regex "[A-Za-z]+\.[A-Za-z]+:[A-Za-z]+" -b - \
  | samtools view -b - \
  > $Outdir_Subset/SDRranger_cDNA_240/RNA_with_bc_umi_edit.bam





#Filter out barcodes based on new tag

#Load samtools - index bam files




samtools index $Outdir_Subset/SDRranger_cDNA_60/RNA_with_bc_umi_edit.bam
samtools index $Outdir_Subset/SDRranger_cDNA_60fail/RNA_with_bc_umi_edit.bam
samtools index $Outdir_Subset/SDRranger_cDNA_120/RNA_with_bc_umi_edit.bam
samtools index $Outdir_Subset/SDRranger_cDNA_240/RNA_with_bc_umi_edit.bam



#60

sinto filterbarcodes -b $Outdir_Subset/SDRranger_cDNA_60/RNA_with_bc_umi_edit.bam \
					      -c $Outdir_Subset/SDRranger_cDNA_60/SDR002_60_meta_cDNA_CB.txt \
					      --barcodetag "CS" \
					      --outdir $Outdir_Subset/SDRranger_cDNA_60



#60 fail 

sinto filterbarcodes -b $Outdir_Subset/SDRranger_cDNA_60fail/RNA_with_bc_umi_edit.bam \
                     -c /$Outdir_Subset/SDRranger_cDNA_60fail/SDR002_60fail_meta_cDNA_CB.txt \
                     --barcodetag "CS" \
                     --outdir $Outdir_Subset/SDRranger_cDNA_60fail



#120

sinto filterbarcodes -b $Outdir_Subset/SDRranger_cDNA_120/RNA_with_bc_umi_edit.bam \
                     -c $Outdir_Subset/SDRranger_cDNA_120/SDR002_120_meta_cDNA_CB.txt \
                     --barcodetag "CS" \
                     --outdir $Outdir_Subset/SDRranger_cDNA_120



#240

sinto filterbarcodes -b $Outdir_Subset/SDRranger_cDNA_240/RNA_with_bc_umi_edit.bam \
                     -c $Outdir_Subset/SDRranger_cDNA_240/SDR002_240_meta_cDNA_CB.txt \
                     --barcodetag "CS" \
                     --outdir $Outdir_Subset/SDRranger_cDNA_240



#Check number of reads

samtools view -c $Outdir_Subset/SDRranger_cDNA_60/SDR002_60.bam 
samtools view -c $Outdir_Subset/SDRranger_cDNA_60fail/SDR002_60fail.bam 
samtools view -c $Outdir_Subset/SDRranger_cDNA_120/SDR002_120.bam 
samtools view -c $Outdir_Subset/SDRranger_cDNA_240/SDR002_240.bam 




#Subsample reads from bam file

#gDNA

#60
samtools view -s 0.39309785 -b $Outdir_Subset/SDRranger_gDNA_60/SDR002_60.bam \
                             > $Outdir_Subset/SDRranger_gDNA_60/SDR002_60_sub.bam

#60fail


#120 - unchanged

#240
samtools view -s 0.293196998 -b $Outdir_Subset/SDRranger_gDNA_240/SDR002_240.bam \
                             > $Outdir_Subset/SDRranger_gDNA_240/SDR002_240_sub.bam





#cDNA

#60

samtools view -s 0.591620461 -b $Outdir_Subset/SDRranger_cDNA_60/SDR002_60.bam  \
                             > $Outdir_Subset/SDRranger_cDNA_60/SDR002_60_sub.bam 


#60fail

samtools view -s 0.484370983 -b $Outdir_Subset/SDRranger_cDNA_60fail/SDR002_60fail.bam  \
                             > $Outdir_Subset/SDRranger_cDNA_60fail/SDR002_60fail_sub.bam 


#120 - unchanged



#240

samtools view -s 0.812812501 -b $Outdir_Subset/SDRranger_cDNA_240/SDR002_240.bam  \
                             > $Outdir_Subset/SDRranger_cDNA_240/SDR002_240_sub.bam 



#Check number of reads



#Check number of reads

#gDNA

samtools view -c $Outdir_Subset/SDRranger_gDNA_60/SDR002_60_sub.bam
samtools view -c $Outdir_Subset/SDRranger_gDNA_60fail/SDR002_60fail_sub.bam
samtools view -c $Outdir_Subset/SDRranger_gDNA_240/SDR002_240_sub.bam


#cDNA


samtools view -c $Outdir_Subset/SDRranger_cDNA_60fail/SDR002_60fail_sub.bam
samtools view -c $Outdir_Subset/SDRranger_cDNA_120/SDR002_120_sub.bam
samtools view -c $Outdir_Subset/SDRranger_cDNA_240/SDR002_240_sub.bam


#############
#_____________________________________________________

#Write tags needed for count_matrix function from SDRranger - GX and GN
#############
#gDNA


#60
samtools view -h $Outdir_Subset/SDRranger_gDNA_60/SDR002_60_sub.bam | \
awk 'BEGIN{FS="\t";OFS="\t"} {if($0 ~ /^@/) {print $0} else {print $0 "\tGX:Z:" $3 "\tGN:Z:" $3}}' | \
samtools view -b \
> $Outdir_Subset/SDRranger_gDNA_60/SDR002_60_sub_cmatin.bam 

samtools index $Outdir_Subset/SDRranger_gDNA_60/SDR002_60_sub_cmatin.bam


#60fail

samtools view -h $Outdir_Subset/SDRranger_gDNA_60fail/SDR002_60fail_sub.bam | \
awk 'BEGIN{FS="\t";OFS="\t"} {if($0 ~ /^@/) {print $0} else {print $0 "\tGX:Z:" $3 "\tGN:Z:" $3}}' | \
samtools view -b \
> $Outdir_Subset/SDRranger_gDNA_60fail/SDR002_60fail_sub_cmatin.bam 

samtools index $Outdir_Subset/SDRranger_gDNA_60fail/SDR002_60fail_sub_cmatin.bam


#120
samtools view -h $Outdir_Subset/SDRranger_gDNA_120/SDR002_120.bam | \
awk 'BEGIN{FS="\t";OFS="\t"} {if($0 ~ /^@/) {print $0} else {print $0 "\tGX:Z:" $3 "\tGN:Z:" $3}}' | \
samtools view -b \
> $Outdir_Subset/SDRranger_gDNA_120/SDR002_120_cmatin.bam 

samtools index $Outdir_Subset/SDRranger_gDNA_120/SDR002_120_cmatin.bam

#240
samtools view -h $Outdir_Subset/SDRranger_gDNA_240/SDR002_240_sub.bam | \
awk 'BEGIN{FS="\t";OFS="\t"} {if($0 ~ /^@/) {print $0} else {print $0 "\tGX:Z:" $3 "\tGN:Z:" $3}}' | \
samtools view -b \
> $Outdir_Subset/SDRranger_gDNA_240/SDR002_240_sub_cmatin.bam 

samtools index $Outdir_Subset/SDRranger_gDNA_240/SDR002_240_sub_cmatin.bam


#############
#cDNA

#60
samtools view -h $Outdir_Subset/SDRranger_cDNA_60/SDR002_60_sub.bam | \
awk 'BEGIN{FS="\t";OFS="\t"} {if($0 ~ /^@/) {print $0} else {print $0 "\tGX:Z:" $3 "\tGN:Z:" $3}}' | \
samtools view -b \
> $Outdir_Subset/SDRranger_cDNA_60/SDR002_60_sub_cmatin.bam 

samtools index $Outdir_Subset/SDRranger_cDNA_60/SDR002_60_sub_cmatin.bam



#60fail

samtools view -h $Outdir_Subset/SDRranger_cDNA_60fail/SDR002_60fail_sub.bam | \
awk 'BEGIN{FS="\t";OFS="\t"} {if($0 ~ /^@/) {print $0} else {print $0 "\tGX:Z:" $3 "\tGN:Z:" $3}}' | \
samtools view -b \
> $Outdir_Subset/SDRranger_cDNA_60fail/SDR002_60fail_sub_cmatin.bam

samtools index $Outdir_Subset/SDRranger_cDNA_60fail/SDR002_60fail_sub_cmatin.bam

#120
samtools view -h $Outdir_Subset/SDRranger_cDNA_120/SDR002_120.bam | \
awk 'BEGIN{FS="\t";OFS="\t"} {if($0 ~ /^@/) {print $0} else {print $0 "\tGX:Z:" $3 "\tGN:Z:" $3}}' | \
samtools view -b \
> $Outdir_Subset/SDRranger_cDNA_120/SDR002_120_cmatin.bam 

samtools index $Outdir_Subset/SDRranger_cDNA_120/SDR002_120_cmatin.bam


#240
samtools view -h $Outdir_Subset/SDRranger_cDNA_240/SDR002_240_sub.bam | \
awk 'BEGIN{FS="\t";OFS="\t"} {if($0 ~ /^@/) {print $0} else {print $0 "\tGX:Z:" $3 "\tGN:Z:" $3}}' | \
samtools view -b \
> $Outdir_Subset/SDRranger_cDNA_240/SDR002_240_sub_cmatin.bam 

samtools index $Outdir_Subset/SDRranger_cDNA_240/SDR002_240_sub_cmatin.bam





#Run count_matrix script from SDRranger












