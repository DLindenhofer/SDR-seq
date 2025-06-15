#!/bin/bash
#SBATCH -A lsteinme              # group to which you belong
#SBATCH -N 1                        # number of nodes
#SBATCH -n 8                        # number of cores
#SBATCH --mem 64G                    # memory pool for all cores
#SBATCH -t 0-24:00:00                   # runtime limit (D-HH:MM:SS)
#SBATCH -o slurm.%N.%j.out          # STDOUT
#SBATCH -e slurm.%N.%j.err          # STDERR
#SBATCH --mail-type=BEGIN,END,FAIL        # notifications for job done & fail
#SBATCH --mail-user=dominik.lindenhofer@embl.de # send-to address
#SBATCH --array=1-17


# Specify the path to the config file
config=/g/steinmetz/projects/SDR_Shared_John_Dominik/SDR002/20240918_GATK_Subset/SDR002_60_60fail_sub/01_config_file_sinto.txt

prefix=$(awk -F '\t' -v ArrayTaskID=$SLURM_ARRAY_TASK_ID 'NR==ArrayTaskID {print $1}' $config)
filename=$(awk -F '\t' -v ArrayTaskID=$SLURM_ARRAY_TASK_ID 'NR==ArrayTaskID {print $2}' $config)

Input_bam="/g/steinmetz/projects/SDR_Shared_John_Dominik/SDR002/20240329_Subset/SDRranger_gDNA_60fail/SDR002_60fail_sub_cmatin.bam"
out_dir="/g/steinmetz/projects/SDR_Shared_John_Dominik/SDR002/20240918_GATK_Subset/SDR002_60_60fail_sub/02_Input_cells_bam"



sinto filterbarcodes -b $Input_bam \
					 -c "$filename" \
					 --outdir $out_dir



