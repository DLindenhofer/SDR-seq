


#Use this to make FASTQ and GTF files for building STAR ref and annotations



# Build RNA


python 07_Species_Mix_buildReference.py \
-c 07_Species_Mix_cDNA_TAP_REF_seqs_comb.csv \
										-o 07_Species_Mix_cDNA_TAP_REF-V2


sbatch 07_Species_Mix_cDNA_TAP_REF_Build_ref.sh




