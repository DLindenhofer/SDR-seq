
#samtools required
source activate /samtools



Loop_dir="02_Input_cells_bam"

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



