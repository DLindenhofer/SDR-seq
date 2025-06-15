
#General preparations

#Create reference dict to be compatible with GATK - NEEDED!

#For samtools installation
source activate samtools

#Index fasta reference file
#Same reference used as in CRISPRi experiment

samtools faidx Panel_size_gDNA_REF.fa

module load GATK
#GATK v4.2.3.0
gatk CreateSequenceDictionary -R Panel_size_gDNA_REF.fa

