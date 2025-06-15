#!/bin/bash
#SBATCH -A lsteinme               # group to which you belong
#SBATCH -N 1                       # number of nodes
#SBATCH -n 8                        # number of cores
#SBATCH --mem 64G                  # memory pool for all cores
#SBATCH -t 1-00:00:00                   # runtime limit (D-HH:MM:SS)
#SBATCH -o /07_Species_Mix_cDNA_TAP_REF-V2/temp.out          # STDOUT
#SBATCH -e /07_Species_Mix_cDNA_TAP_REF-V2/temp.err          # STDERR
#SBATCH --mail-type=END,FAIL        # notifications for job done & fail
#SBATCH --mail-user=dominik.lindenhofer@embl.de # send-to address

module add STAR
module load STAR

STAR --runThreadN 1 \
--runMode genomeGenerate \
--genomeDir 07_Species_Mix_cDNA_TAP_REF_index-V2 \
--genomeFastaFiles 07_Species_Mix_cDNA_TAP_REF-V2.fa \
--sjdbGTFfile 07_Species_Mix_cDNA_TAP_REF-V2.gtf \
--sjdbOverhang 99 \
--genomeSAindexNbases 4