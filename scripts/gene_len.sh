#!/bin/bash
set -e  # Exit immediately if a command exits with a non-zero status

# Variables
# Uncomment and set the appropriate paths for human and E. coli annotations
# For Human
# ANNOTATION_GTF="/mnt/d/Jiahe/Microsoft/miRBase_Data_Scraper/human_reference/annotations/gencode.v47.chr_patch_hapl_scaff.annotation.gtf"
# OUTPUT_CSV="/mnt/d/Jiahe/Microsoft/miRBase_Data_Scraper/gh38_genes.csv"

# For E. coli
ANNOTATION_GFF="/mnt/d/Jiahe/Microsoft/miRBase_Data_Scraper/ecoli_reference/annotations/ecoli_annotations.gff"
OUTPUT_CSV="/mnt/d/Jiahe/Microsoft/miRBase_Data_Scraper/ecoli_genes.csv"

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

# Check if the annotation file exists
if [ ! -f "$ANNOTATION_GFF" ]; then
    echo "Error: Annotation file not found at $ANNOTATION_GFF."
    exit 1
fi

echo "Using annotation file: $ANNOTATION_GFF"

# Write header to CSV
echo "gene_id,gene_name,start,end,length" > "$OUTPUT_CSV"

# Parse the annotation file and extract required fields
echo "Processing $ANNOTATION_GFF to generate $OUTPUT_CSV..."

awk '
BEGIN { FS = "\t"; OFS = "," }
$3 == "gene" {
    gene_id = ""
    gene_name = ""

    # Attempt to extract gene_id from gene_id "ENSG..."
    if (match($9, /gene_id "([^"]+)"/, arr)) {
        gene_id = arr[1]
    }
    # Else, attempt to extract gene_id from ID=gene-b0001
    else if (match($9, /ID=gene-([^;]+)/, arr)) {
        gene_id = arr[1]
    }

    # Attempt to extract gene_name from gene_name "BRCA1"
    if (match($9, /gene_name "([^"]+)"/, arr)) {
        gene_name = arr[1]
    }
    # Else, attempt to extract gene_name from Name=thrL
    else if (match($9, /Name=([^;]+)/, arr)) {
        gene_name = arr[1]
    }

    # Calculate gene length
    start = $4
    end = $5
    gene_length = end - start + 1  # Inclusive of both start and end

    # Print to CSV
    print gene_id, gene_name, start, end, gene_length
}
' "$ANNOTATION_GFF" >> "$OUTPUT_CSV"

echo "CSV file generated at $OUTPUT_CSV."
