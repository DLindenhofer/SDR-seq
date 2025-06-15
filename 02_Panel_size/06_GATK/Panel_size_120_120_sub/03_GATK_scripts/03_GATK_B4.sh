#!/bin/bash
#SBATCH -A lsteinme              # group to which you belong
#SBATCH -N 1                        # number of nodes
#SBATCH -n 1                        # number of cores
#SBATCH --mem 4G                    # memory pool for all cores
#SBATCH -t 0-24:00:00                   # runtime limit (D-HH:MM:SS)
#SBATCH -o /g/steinmetz/projects/SDR_Shared_John_Dominik/SDR002/20240918_GATK_Subset/SDR002_120_120_sub/04_GATK_Slurm_Out/slurm.%N.%j.out          # STDOUT
#SBATCH -e /g/steinmetz/projects/SDR_Shared_John_Dominik/SDR002/20240918_GATK_Subset/SDR002_120_120_sub/04_GATK_Slurm_Out/slurm.%N.%j.err          # STDERR
#SBATCH --mail-type=BEGIN,END,FAIL        # notifications for job done & fail
#SBATCH --mail-user=dominik.lindenhofer@embl.de # send-to address
#SBATCH --array=1-1000%100

module load GATK


# Specify the path to the config file
config="/g/steinmetz/projects/SDR_Shared_John_Dominik/SDR002/20240918_GATK_Subset/SDR002_120_120_sub/02_Input_GATK/Batch_4.txt"

prefix=$(awk -F '\t' -v ArrayTaskID=$SLURM_ARRAY_TASK_ID 'NR==ArrayTaskID {print $1}' $config)
filename=$(awk -F '\t' -v ArrayTaskID=$SLURM_ARRAY_TASK_ID 'NR==ArrayTaskID {print $2}' $config)

ref="/g/steinmetz/projects/SDR_Shared_John_Dominik/SDR_Genome_REFs/SDR002_gDNA_REF-V2/SDR002_gDNA_REF-V2.fa"
out="/g/steinmetz/projects/SDR_Shared_John_Dominik/SDR002/20240918_GATK_Subset/SDR002_120_120_sub/04_GATK_Out"


gatk HaplotypeCaller --java-options "-Xmx4g" \
                     --max-reads-per-alignment-start 0 \
                         -R "$ref" \
                         -I "$filename" \
                         -O $out/${prefix}.vcf \
                         -ploidy 2












