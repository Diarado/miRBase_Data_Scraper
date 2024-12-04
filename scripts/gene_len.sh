#!/bin/bash
set -e  # Exit immediately if a command exits with a non-zero status

# ----------------------------
# Variables
# ----------------------------

# Directory to store downloaded and processed data
DATA_DIR="$HOME/ncbi_human_genome"
mkdir -p "$DATA_DIR"
cd "$DATA_DIR"

# Output CSV file
OUTPUT_CSV="$DATA_DIR/human_genes.csv"

# NCBI Taxonomy ID for Homo sapiens
TAXON_ID="9606"

# ----------------------------
# Function to check if a command exists
# ----------------------------
command_exists () {
    command -v "$1" >/dev/null 2>&1
}

# ----------------------------
# Check for Required Commands
# ----------------------------

REQUIRED_COMMANDS=("datasets" "awk" "bedtools" "unzip" "gunzip")

for cmd in "${REQUIRED_COMMANDS[@]}"; do
    if ! command_exists "$cmd"; then
        echo "Error: Required command '$cmd' is not installed."
        exit 1
    fi
done

# ----------------------------
# Step 1: Download Human Genome and Annotations from NCBI
# ----------------------------

echo "Downloading human genome and annotations from NCBI Datasets CLI..."

# Download the human genome dataset with annotations
datasets download genome taxon "$TAXON_ID" --filename human_genome.zip

# Unzip the downloaded dataset
echo "Unzipping human_genome.zip..."
unzip -o human_genome.zip -d human_genome

# Paths to important files
GFF3_FILE=$(find human_genome -type f -name "*.gff3" | head -n 1)
FASTA_FILE=$(find human_genome -type f -name "*.fna.gz" | head -n 1)

# Check if GFF3 and FASTA files are found
if [ -z "$GFF3_FILE" ]; then
    echo "Error: GFF3 annotation file not found in the downloaded dataset."
    exit 1
fi

if [ -z "$FASTA_FILE" ]; then
    echo "Error: Genome FASTA file not found in the downloaded dataset."
    exit 1
fi

echo "Found GFF3 file: $GFF3_FILE"
echo "Found FASTA file: $FASTA_FILE"

# Unzip the FASTA file if it's gzipped
if [[ "$FASTA_FILE" == *.gz ]]; then
    echo "Unzipping FASTA file..."
    gunzip -k "$FASTA_FILE"  # -k to keep the original gzipped file
    FASTA_FILE="${FASTA_FILE%.gz}"
    echo "Fasta file unzipped to: $FASTA_FILE"
fi

# ----------------------------
# Step 2: Index the Genome FASTA with bedtools
# ----------------------------

echo "Indexing the genome FASTA file with bedtools..."

bedtools faidx "$FASTA_FILE"

# ----------------------------
# Step 3: Process the GFF3 File to Extract Gene Information
# ----------------------------

echo "Processing GFF3 file to extract gene information..."

# Write header to CSV
echo "gene_id,gene_name,start,end,length,sequence,summary" > "$OUTPUT_CSV"

# Temporary BED file for gene coordinates
BED_FILE="$DATA_DIR/human_genes.bed"

# Extract gene entries and convert to BED format
awk '
BEGIN { FS = "\t" }
$3 == "gene" {
    # Initialize variables
    gene_id = ""
    gene_name = ""
    description = ""
    chrom = $1
    start = $4
    end = $5

    # Parse attributes
    split($9, attrs, ";")
    for(i in attrs) {
        split(attrs[i], pair, "=")
        key = pair[1]
        value = pair[2]
        if(key == "gene_id") {
            gene_id = value
        }
        if(key == "gene_name") {
            gene_name = value
        }
        if(key == "description") {
            description = value
        }
    }

    # Calculate gene length
    gene_length = end - start + 1

    # Print BED format: chrom, start-1 (0-based), end, gene_id, gene_length, strand
    # bedtools uses 0-based start, so subtract 1
    print chrom "\t" (start - 1) "\t" end "\t" gene_id "\t" gene_length "\t" $7 "\t" gene_name "\t" description
}' "$GFF3_FILE" > "$BED_FILE"

# ----------------------------
# Step 4: Extract Gene Sequences Using bedtools
# ----------------------------

echo "Extracting gene sequences using bedtools..."

# Extract sequences and save to FASTA
GENE_FASTA="$DATA_DIR/human_genes_sequences.fasta"
bedtools getfasta -fi "$FASTA_FILE" -bed "$BED_FILE" -name -fo "$GENE_FASTA"

# ----------------------------
# Step 5: Compile Gene Information into CSV
# ----------------------------

echo "Compiling gene information into CSV..."

# Create associative array to store sequences
declare -A gene_sequences

# Read the FASTA file and store sequences in the array
while read -r line; do
    if [[ "$line" == ">"* ]]; then
        gene_id_seq=$(echo "$line" | sed 's/>//')
        current_gene="$gene_id_seq"
        gene_sequences["$current_gene"]=""
    else
        gene_sequences["$current_gene"]+="$line"
    fi
done < "$GENE_FASTA"

# Process BED file and append information to CSV
while IFS=$'\t' read -r chrom start end gene_id gene_length strand gene_name description; do
    # Retrieve sequence from the associative array
    sequence="${gene_sequences[$gene_id]}"
    
    # Handle cases where description might contain quotes or commas
    # Enclose the description in double quotes and escape any existing double quotes
    summary=$(echo "$description" | sed 's/"/""/g')
    summary="\"$summary\""

    # Append to CSV
    echo "$gene_id,$gene_name,$start,$end,$gene_length,$sequence,$summary" >> "$OUTPUT_CSV"
done < "$BED_FILE"

echo "CSV file generated at $OUTPUT_CSV."
