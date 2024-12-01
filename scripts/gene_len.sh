#!/bin/bash
set -e  # Exit immediately if a command exits with a non-zero status

# Variables
ANNOTATION_GTF="/mnt/d/Jiahe/Microsoft/miRBase_Data_Scraper/human_reference/annotations/gencode.v47.chr_patch_hapl_scaff.annotation.gtf"
OUTPUT_CSV="/mnt/d/Jiahe/Microsoft/miRBase_Data_Scraper/gh38_genes.csv"

# Function to check if a command exists
command_exists () {
    command -v "$1" >/dev/null 2>&1
}

# Check for required commands
for cmd in awk; do
    if ! command_exists "$cmd"; then
        echo "Error: Required command '$cmd' is not installed."
        exit 1
    fi
done

# Check if the annotation GTF file exists
if [ ! -f "$ANNOTATION_GTF" ]; then
    echo "Error: Annotation GTF file not found at $ANNOTATION_GTF."
    exit 1
fi

echo "Using annotation file: $ANNOTATION_GTF"

# Write header to CSV
echo "gene_id,start,end,length" > "$OUTPUT_CSV"

# Parse the GTF file and extract required fields
echo "Processing $ANNOTATION_GTF to generate $OUTPUT_CSV..."

awk '
BEGIN { FS = "\t"; OFS = "," }
$3 == "gene" {
    # Extract gene_id from the attributes column
    match($9, /gene_id "([^"]+)"/, arr)
    gene_id = arr[1]
    start = $4
    end = $5
    gene_length = end - start + 1
    print gene_id, start, end, gene_length
}
' "$ANNOTATION_GTF" >> "$OUTPUT_CSV"

echo "CSV file generated at $OUTPUT_CSV."
