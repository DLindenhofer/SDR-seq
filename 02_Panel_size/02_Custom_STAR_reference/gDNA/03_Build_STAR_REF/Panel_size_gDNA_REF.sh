#!/bin/bash
#SBATCH -A lsteinme               # group to which you belong
#SBATCH -N 1                       # number of nodes
#SBATCH -n 8                        # number of cores
#SBATCH --mem 128G                  # memory pool for all cores
#SBATCH -t 1-00:00:00                   # runtime limit (D-HH:MM:SS)
#SBATCH -o temp.out          # STDOUT
#SBATCH -e temp.err          # STDERR
#SBATCH --mail-type=END,FAIL        # notifications for job done & fail
#SBATCH --mail-user=dominik.lindenhofer@embl.de # send-to address

module add STAR
module load STAR

STAR --runThreadN 1 \
--runMode genomeGenerate \
--genomeDir Panel_size_gDNA_REF_index \
--genomeFastaFiles Panel_size_gDNA_REF.fa \
--sjdbGTFfile Panel_size_gDNA_REF.gtf \
--sjdbOverhang 99 \
--genomeSAindexNbases 4