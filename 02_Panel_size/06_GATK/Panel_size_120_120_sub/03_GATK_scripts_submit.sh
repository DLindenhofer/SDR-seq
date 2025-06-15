#!/bin/bash
#SBATCH -A lsteinme              # group to which you belong
#SBATCH -N 1                        # number of nodes
#SBATCH -n 1                        # number of cores
#SBATCH --mem 1G                    # memory pool for all cores
#SBATCH -t 0-24:00:00                   # runtime limit (D-HH:MM:SS)
#SBATCH -o slurm.%N.%j.out          # STDOUT
#SBATCH -e slurm.%N.%j.err          # STDERR
#SBATCH --mail-type=BEGIN,END,FAIL        # notifications for job done & fail
#SBATCH --mail-user=dominik.lindenhofer@embl.de # send-to address

# Path to the folder containing the job scripts
job_script_folder="/g/steinmetz/projects/SDR_Shared_John_Dominik/SDR002/20240918_GATK_Subset/SDR002_120_120_sub/03_GATK_scripts"

# Maximum number of jobs to submit at the same time
max_jobs=7

# Loop through all job scripts in the folder
for job_script in "$job_script_folder"/*.sh; do
    # Check how many jobs are currently in the queue
    while [ "$(squeue -u lindenho | grep -E ' R | PD ' | wc -l)" -ge "$max_jobs" ]; do
        # Wait for 60 seconds before checking again
        sleep 60
    done
    
    # Submit the job script
    sbatch "$job_script"
    
    # Optionally, add a small delay between submissions to avoid flooding the scheduler
    sleep 60
done
