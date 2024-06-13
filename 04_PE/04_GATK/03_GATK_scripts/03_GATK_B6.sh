#!/bin/bash
#SBATCH -A lsteinme              # group to which you belong
#SBATCH -N 1                        # number of nodes
#SBATCH -n 1                        # number of cores
#SBATCH --mem 4G                    # memory pool for all cores
#SBATCH -t 0-24:00:00                   # runtime limit (D-HH:MM:SS)
#SBATCH --mail-type=BEGIN,END,FAIL        # notifications for job done & fail
#SBATCH --mail-user=dominik.lindenhofer@embl.de # send-to address
#SBATCH --array=1-1000%100

module load GATK


# Specify the path to the config file
config="/02_Input_GATK/Batch_6.txt"

prefix=$(awk -F '\t' -v ArrayTaskID=$SLURM_ARRAY_TASK_ID 'NR==ArrayTaskID {print $1}' $config)
filename=$(awk -F '\t' -v ArrayTaskID=$SLURM_ARRAY_TASK_ID 'NR==ArrayTaskID {print $2}' $config)

ref="CRISPRi_gDNA_REF.fa"
out="04_GATK_Out"


gatk HaplotypeCaller --java-options "-Xmx4g" \
                     --max-reads-per-alignment-start 0 \
                         -R "$ref" \
                         -I "$filename" \
                         -O $out/${prefix}.vcf \
                         -ploidy 2












