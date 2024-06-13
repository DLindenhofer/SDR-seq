#!/bin/bash
#SBATCH -A lsteinme               # group to which you belong
#SBATCH -N 1                       # number of nodes
#SBATCH -n 8                        # number of cores
#SBATCH --mem 64G                  # memory pool for all cores
#SBATCH -t 1-00:00:00                   # runtime limit (D-HH:MM:SS)
#SBATCH --mail-type=END,FAIL        # notifications for job done & fail
#SBATCH --mail-user=dominik.lindenhofer@embl.de # send-to address

module add STAR
module load STAR

STAR --runThreadN 1 \
--runMode genomeGenerate \
--genomeDir 05_BE_REF_index \
--genomeFastaFiles 05_BE_gDNA_REF.fa \
--sjdbGTFfile 05_BE_gDNA_REF.gtf \
--sjdbOverhang 99 \
--genomeSAindexNbases 4