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
#SBATCH --array=1-5


# Specify the path to the config file
#Change array size according to number of batches to submit

config=<INSERT_PATH>

prefix=$(awk -F '\t' -v ArrayTaskID=$SLURM_ARRAY_TASK_ID 'NR==ArrayTaskID {print $1}' $config)
filename=$(awk -F '\t' -v ArrayTaskID=$SLURM_ARRAY_TASK_ID 'NR==ArrayTaskID {print $2}' $config)

#Use gDNA bam file from SDRranger output as input 
Input_bam="gDNA_with_bc.sorted.bam"
out_dir="02_Input_cells_bam"



sinto filterbarcodes -b $Input_bam \
					 -c "$filename" \
					 --outdir $out_dir



