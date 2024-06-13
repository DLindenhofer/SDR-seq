


#Use this to make FASTQ and GTF files for building STAR ref and annotations



# Build CROP-seq references 


python CRISPRi_CROP_buildReference.py -i CRISPRi_CROP_sgRNAs.csv \
							-c CRISPRI_CROP_vector_seq_mod.csv \
							-o CRISPRi_CROP_REF

