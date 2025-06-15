#!/bin/bash
#SBATCH -A lsteinme              # group to which you belong
#SBATCH -N 1                        # number of nodes
#SBATCH -n 1                        # number of cores
#SBATCH --mem 12G                    # memory pool for all cores
#SBATCH -t 0-24:00:00                   # runtime limit (D-HH:MM:SS)
#SBATCH -o slurm.%N.%j.out          # STDOUT
#SBATCH -e slurm.%N.%j.err          # STDERR
#SBATCH --mail-type=BEGIN,END,FAIL        # notifications for job done & fail
#SBATCH --mail-user=dominik.lindenhofer@embl.de # send-to address

#For samtools
module load SAMtools/1.17-GCC-12.2.0
Loop_dir="/g/steinmetz/projects/SDR_Shared_John_Dominik/SDR002/20240918_GATK_Subset/SDR002_240_240_sub/02_Input_cells_bam"

for file in $Loop_dir/*.bam; do
    # Extract the basename for use in SM tag
    base=$(basename "$file" .bam)
    echo "Processing $base..."

    # Intermediate filenames
    file_withRG="${file%.bam}_withRG.bam"
    file_withRG_MAPQ="${file%.bam}_withRG_MAPQ.bam"

    # Step 1: Add Read Groups
    samtools addreplacerg -r "@RG\tID:1\tLB:lib1\tPL:illumina\tPU:unit1\tSM:$base" \
        -o "$file_withRG" "$file"

    # Step 2: Modify MAPQ scores, then sort and write back to original filename
    samtools view -h "$file_withRG" | \
    awk 'BEGIN{FS="\t";OFS="\t"} {if($5=="255") $5="60"; print}' | \
    samtools view -b | \
    samtools sort -o "$file"

    # Step 3: Index the now overwritten original BAM file
    samtools index "$file"

    # Optionally, remove intermediate files (if any)
    rm -f "$file_withRG" "$file_withRG_MAPQ"

done



