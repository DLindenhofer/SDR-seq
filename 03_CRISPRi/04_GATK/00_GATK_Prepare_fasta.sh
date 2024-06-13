
#General preparations

#Create reference dict to be compatible with GATK - NEEDED!

#For samtools installation
source activate samtools

#Index fasta reference file

samtools faidx CRISPRi_gDNA_REF.fa

module load GATK
#GATK v4.2.3.0
gatk CreateSequenceDictionary -R CRISPRi_gDNA_REF.fa

