


#Use this to make FASTQ and GTF files for building STAR ref and annotations



# Build gDNA


python 07_Species_Mix_buildReference.py \
-c 07_Species_Mix_gDNA_targets_comb.csv \
										-o 07_Species_Mix_gDNA_REF


nohup sbatch 07_Species_Mix_gDNA_REF_Build_ref.sh &



