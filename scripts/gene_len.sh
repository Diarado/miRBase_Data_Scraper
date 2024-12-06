#!/bin/bash
set -e  # Exit immediately if a command exits with a non-zero status
set -o pipefail  # Fail a pipeline if any command errors

# Add debugging output
exec 1> >(tee -a "debug_output.log")
exec 2>&1

# Directory to store downloaded and processed data
DATA_DIR="/mnt/d/Jiahe/Microsoft/miRBase_Data_Scraper/reference"
mkdir -p "$DATA_DIR"
cd "$DATA_DIR"

# Output CSV file
OUTPUT_CSV="$DATA_DIR/human_genes.csv"

# NCBI Accession for Homo sapiens
ACCESSION="GCF_000001405.40"

echo "Script started"
echo "Working directory: $(pwd)"
echo "DATA_DIR: $DATA_DIR"

# ----------------------------
# Function to check if a command exists
# ----------------------------
command_exists () {
    command -v "$1" >/dev/null 2>&1
}

# ----------------------------
# Check for Required Commands
# ----------------------------

REQUIRED_COMMANDS=("datasets" "awk" "bedtools" "unzip" "gunzip" "wget" "samtools")

for cmd in "${REQUIRED_COMMANDS[@]}"; do
    if ! command_exists "$cmd"; then
        echo "Error: Required command '$cmd' is not installed or not in PATH."
        echo "Please install '$cmd' before running this script."
        exit 1
    fi
done

# ----------------------------
# Step 1: Download and Extract
# ----------------------------

echo "Downloading human genome and annotations..."
if [ ! -f human_genome.zip ]; then
    datasets download genome accession "$ACCESSION" \
        --include gff3,rna,cds,protein,genome,seq-report \
        --filename human_genome.zip
fi

echo "Extracting files..."
if [ -f human_genome.zip ]; then
    unzip -o human_genome.zip -d human_genome
else
    echo "Error: human_genome.zip not found"
    exit 1
fi

# ----------------------------
# Step 2: Locate Files
# ----------------------------

# Locate the GFF3 file
GFF3_FILE=$(find human_genome -type f \( -name "*.gff3" -o -name "*.gff" \) | head -n 1)

# Verify GFF3 file
if [ -z "$GFF3_FILE" ]; then
    echo "Error: GFF3 file not found"
    echo "Contents of human_genome directory:"
    ls -R human_genome
    exit 1
fi

echo "Found GFF3 file: $GFF3_FILE"
echo "GFF3 file size: $(ls -lh "$GFF3_FILE")"

# ----------------------------
# Step 3: Process FASTA Files
# ----------------------------

echo "Combining all FASTA files into a single genome FASTA..."

# Find all genomic FASTA files
ALL_FASTA_FILES=$(find human_genome -type f \( -name "*genomic.fna.gz" -o -name "*genomic.fna" \))

# Temporary combined FASTA file
COMBINED_FASTA="$DATA_DIR/combined_genomic.fna"

# Initialize the combined FASTA file
> "$COMBINED_FASTA"

# Iterate over all FASTA files and append their contents
for fasta in $ALL_FASTA_FILES; do
    if [[ "$fasta" == *.gz ]]; then
        echo "Unzipping and adding $fasta to combined FASTA..."
        gunzip -c "$fasta" >> "$COMBINED_FASTA"
    else
        echo "Adding $fasta to combined FASTA..."
        cat "$fasta" >> "$COMBINED_FASTA"
    fi
done

echo "All FASTA files have been combined."
echo "Combined FASTA size: $(ls -lh "$COMBINED_FASTA")"

# ----------------------------
# Step 4: Index FASTA
# ----------------------------

echo "Indexing the combined genome FASTA file..."
samtools faidx "$COMBINED_FASTA"

# ----------------------------
# Step 5: Process GFF3
# ----------------------------

echo "Processing GFF3 file..."
BED_FILE="$DATA_DIR/human_genes.bed"

# Write header to CSV
echo "gene_id,gene_name,start,end,length,sequence,summary" > "$OUTPUT_CSV"

# Create BED file with improved attribute parsing
awk -v BED_FILE="$BED_FILE" '
BEGIN { 
    FS="\t"; 
    OFS="\t";
    print "Starting GFF3 processing..." > "/dev/stderr"
}
$3 == "gene" {
    print "Processing gene entry at line " NR > "/dev/stderr"
    
    # Initialize variables
    gene_id = "NA"
    gene_name = "NA"
    description = "NA"
    summary = "NA"
    
    # Parse attributes
    n = split($9, attrs, ";")
    for (i = 1; i <= n; i++) {
        gsub(/^[ \t]+/, "", attrs[i])
        
        # Handle Dbxref attributes specially
        if (attrs[i] ~ /^Dbxref=/) {
            split(attrs[i], dbxrefs, ",")
            for (j in dbxrefs) {
                if (dbxrefs[j] ~ /^GeneID:/) {
                    split(dbxrefs[j], temp, ":")
                    gene_id = temp[2]
                }
            }
            continue
        }
        
        # Parse key-value pairs
        split(attrs[i], pair, "=")
        key = pair[1]
        value = pair[2]
        
        if (key == "ID") {
            if (gene_id == "NA") gene_id = value
        } else if (key == "Name") {
            gene_name = value
        } else if (key == "description") {
            description = value
        } else if (key == "Note" && value ~ /Summary:/) {
            # Extract the summary from the Note field
            summary = value
            # Remove "Summary:" prefix if present
            sub(/^Summary:[ ]*/, "", summary)
        }
    }
    
    # If no specific summary found, use description
    if (summary == "NA") {
        summary = description
    }
    
    # Print BED format with summary included
    if (gene_id != "NA") {
        print $1, $4-1, $5, gene_id, ($5 - $4 + 1), $7, gene_name, summary
    }
}
END {
    print "GFF3 processing complete" > "/dev/stderr"
}' "$GFF3_FILE" > "$BED_FILE"

# Verify BED file creation
if [ ! -s "$BED_FILE" ]; then
    echo "Error: BED file is empty or was not created"
    exit 1
fi

echo "BED file created successfully"
echo "First few lines of BED file:"
head -n 5 "$BED_FILE"

# ----------------------------
# Step 6: Extract Sequences
# ----------------------------

echo "Extracting gene sequences..."
GENE_FASTA="$DATA_DIR/human_genes_sequences.fasta"
bedtools getfasta -fi "$COMBINED_FASTA" -bed "$BED_FILE" -s -name -fo "$GENE_FASTA"

if [ ! -s "$GENE_FASTA" ]; then
    echo "Error: FASTA file is empty or was not created"
    exit 1
fi

echo "Sequence extraction complete"

# ----------------------------
# Step 7: Create Final CSV
# ----------------------------

echo "Creating final CSV..."
TEMP_SEQ_FILE="$DATA_DIR/temp_sequences.txt"

# Extract sequences to temporary file
awk '
/^>/ {
    gsub(/^>/, "", $0)
    split($0, header, ":")
    gene_id = header[1]
    getline
    print gene_id "\t" $0
}' "$GENE_FASTA" > "$TEMP_SEQ_FILE"

# Combine information and create final CSV
awk '
BEGIN {
    FS="\t"
    OFS=","
    while ((getline < "'$TEMP_SEQ_FILE'") > 0) {
        split($0, fields, "\t")
        sequences[fields[1]] = fields[2]
    }
}
{
    gene_id = $4
    gene_name = $7
    start = $2 + 1
    end = $3
    len = $5
    sequence = sequences[gene_id]
    if (sequence == "") sequence = "NA"
    
    # Handle the summary field
    summary = $8
    gsub(/"/, "\"\"", summary)
    summary = "\"" summary "\""
    
    print gene_id, gene_name, start, end, len, sequence, summary
}' "$BED_FILE" > "$OUTPUT_CSV"

# Clean up
rm -f "$TEMP_SEQ_FILE"

echo "Processing complete. Output saved to $OUTPUT_CSV"
echo "Final CSV file size: $(ls -lh "$OUTPUT_CSV")"
echo "First few lines of CSV:"
head -n 5 "$OUTPUT_CSV"