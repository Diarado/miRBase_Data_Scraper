#!/usr/bin/env bash
set -euo pipefail

DATA_DIR="./reference_data"
ACCESSION="GCF_000001405.40"

mkdir -p "$DATA_DIR"
cd "$DATA_DIR"

OUTPUT_CSV="human_genes_with_summaries.csv"

echo "Starting test pipeline for first 5 genes..."
# ----------------------------
# Check for required tools
# ----------------------------
REQUIRED_COMMANDS=("datasets" "awk" "bedtools" "unzip" "gunzip" "wget" "samtools" "xmllint" "efetch")

for cmd in "${REQUIRED_COMMANDS[@]}"; do
    if ! command -v "$cmd" > /dev/null; then
        echo "Error: Required command '$cmd' not found. Please install it."
        exit 1
    fi
done

# ----------------------------
# Step 1: Download human genome and annotations if not present
# ----------------------------
if [ ! -f "human_genome.zip" ]; then
    echo "Downloading human genome and annotations..."
    datasets download genome accession "$ACCESSION" \
        --include gff3,rna,cds,protein,genome,seq-report \
        --filename human_genome.zip
else
    echo "human_genome.zip already exists. Skipping download."
fi

# ----------------------------
# Step 2: Extract downloaded data if needed
# ----------------------------
GFF3_FILE=""
if [ ! -d "human_genome" ] || [ -z "$(find human_genome -type f -name '*.gff*' 2>/dev/null)" ]; then
    echo "Extracting downloaded data..."
    unzip -o human_genome.zip -d human_genome
    echo "Extraction complete."
else
    echo "human_genome directory with GFF3 files already present. Skipping extraction."
fi

GFF3_FILE=$(find human_genome -type f -name "*.gff*" | head -n 1)
if [ -z "$GFF3_FILE" ]; then
    echo "Error: No GFF3 file found after extraction."
    exit 1
fi
echo "Found GFF3 file: $GFF3_FILE"

# ----------------------------
# Step 3: Combine FASTA and index only if not done
# ----------------------------
COMBINED_FASTA="combined_genomic.fna"

if [ ! -f "$COMBINED_FASTA" ]; then
    echo "Combining FASTA files..."
    > "$COMBINED_FASTA"
    ALL_FASTA_FILES=$(find human_genome -type f -name "*genomic.fna*" )
    for f in $ALL_FASTA_FILES; do
        if [[ "$f" == *.gz ]]; then
            gunzip -c "$f" >> "$COMBINED_FASTA"
        else
            cat "$f" >> "$COMBINED_FASTA"
        fi
    done
    echo "FASTA combination complete."
else
    echo "$COMBINED_FASTA already exists. Skipping FASTA combination."
fi

if [ ! -f "$COMBINED_FASTA.fai" ]; then
    echo "Indexing FASTA..."
    samtools faidx "$COMBINED_FASTA"
    echo "FASTA indexing complete."
else
    echo "$COMBINED_FASTA.fai already exists. Skipping indexing."
fi

# ----------------------------
# Step 4: Parse GFF3 to extract gene annotations if not done
# ----------------------------
BED_FILE="human_genes.bed"
if [ ! -s "$BED_FILE" ]; then
    echo "Extracting gene annotations from GFF3..."
    awk -F"\t" -v OFS="\t" '
    $3 == "gene" {
        gene_id="NA"; gene_name="NA"; desc="NA";
        n = split($9, attrs, ";")
        for(i=1; i<=n; i++){
            split(attrs[i], pair, "=")
            key=pair[1]; value=pair[2]
            gsub(/^ +| +$/, "", key)
            gsub(/^ +| +$/, "", value)
            
            if(key=="Name"){
                gene_name=value
            }
            
            # Find GeneID in Dbxref
            if(key=="Dbxref"){
                split(value, dbx, ",")
                for(j in dbx){
                    if(dbx[j] ~ /^GeneID:/){
                        split(dbx[j], t, ":")
                        gene_id=t[2]
                    }
                }
            }
            
            if(key=="description"){
                desc=value
            }
        }
        if(gene_id != "NA") {
            gene_len=($5-$4+1)
            print $1, $4-1, $5, gene_id, gene_len, $7, gene_name, desc
        }
    }' "$GFF3_FILE" > "$BED_FILE"

    if [ ! -s "$BED_FILE" ]; then
        echo "No gene entries found in GFF3."
        exit 1
    fi
    echo "Gene annotations extracted."
else
    echo "$BED_FILE already exists. Skipping GFF3 parsing."
fi

# ----------------------------
# Step 5: Extract gene sequences if not done
# ----------------------------
GENE_FASTA="human_genes_sequences.fasta"
if [ ! -s "$GENE_FASTA" ]; then
    echo "Extracting gene sequences..."
    bedtools getfasta -fi "$COMBINED_FASTA" -bed "$BED_FILE" -s -name -fo "$GENE_FASTA"
    echo "Gene sequences extracted."
else
    echo "$GENE_FASTA already exists. Skipping sequence extraction."
fi

# ----------------------------
# Step 6: Create initial CSV if not done
# ----------------------------
if [ ! -f "initial_genes.csv" ]; then
    echo "Creating initial CSV..."
    awk '
    /^>/ {
        header=substr($0,2)
        split(header,h,":")
        gene_id=h[1]
        getline seq
        seqs[gene_id]=seq
        next
    }
    END {
        # Nothing printed here, just storing sequences
    }' "$GENE_FASTA" > /dev/null

    awk '
    BEGIN{
        FS="\t"; OFS=","
        print "gene_id,gene_name,start,end,gene_len,sequence,description"
    }
    FNR==NR {
        # reading FASTA sequences
        if($0 ~ /^>/) {
           gene_id=substr($0,2)
           split(gene_id,a,":")
           gene_id=a[1]
           getline seq
           seqs[gene_id]=seq
        }
        next
    }
    NR>FNR {
        # reading BED_FILE
        gene_id=$4
        gene_name=$7
        start=$2+1
        end=$3
        gene_len=$5
        desc=$8
        seq="NA"
        if(gene_id in seqs) {
            seq=seqs[gene_id]
        }
        print gene_id,gene_name,start,end,gene_len,seq,desc
    }
    ' "$GENE_FASTA" "$BED_FILE" > initial_genes.csv
    echo "Initial CSV created."
else
    echo "initial_genes.csv already exists. Skipping initial CSV creation."
fi

# ----------------------------
# Limit to first 5 genes for testing
# ----------------------------
# echo "Filtering to first 5 genes..."
# head -n 3 initial_genes.csv > test_genes.csv
# echo "Test file created: test_genes.csv"

# ----------------------------
# Step 7: Fetch Official Full Name and Summary from NCBI for ALL genes
# ----------------------------
echo "Fetching official full name and summary for all genes from NCBI..."
TEMP_IDS=$(mktemp)
cut -d',' -f1 initial_genes.csv | tail -n +2 | sort | uniq > "$TEMP_IDS"

TEMP_XML=$(mktemp)
TEMP_PARSED=$(mktemp)

echo "gene_id,official_full_name,official_summary" > "$TEMP_PARSED"

while read -r GENE_ID; do
    [ -z "$GENE_ID" ] && continue
    echo "Processing GENE_ID: $GENE_ID"

    # Use timeout to avoid hanging forever. If it takes more than 30s, kill it.
    if ! timeout 30 efetch -db gene -id "$GENE_ID" -format docsum > "$TEMP_XML" 2>/dev/null; then
        echo "Error: efetch timed out or failed for $GENE_ID"
        OFFICIAL_NAME="NA"
        OFFICIAL_SUMMARY=""
    else
        OFFICIAL_NAME=$(xmllint --xpath 'string(//DocumentSummary/Item[@Name="OfficialFullName"])' "$TEMP_XML" 2>/dev/null || true)
        [ -z "$OFFICIAL_NAME" ] && OFFICIAL_NAME=$(xmllint --xpath 'string(//DocumentSummary/Item[@Name="Description"])' "$TEMP_XML" 2>/dev/null || true)
        OFFICIAL_SUMMARY=$(xmllint --xpath 'string(//DocumentSummary/Item[@Name="Summary"])' "$TEMP_XML" 2>/dev/null || true)

        # If no summary found via docsum, attempt HTML scraping with timeout
        if [ -z "$OFFICIAL_SUMMARY" ]; then
            echo "No summary in docsum, fetching HTML page..."

            # Use timeout for wget as well, --timeout sets a network timeout, but we can wrap in 'timeout' too
            if ! timeout 30 wget -qO gene_page.html "https://www.ncbi.nlm.nih.gov/gene/${GENE_ID}"; then
                echo "Error: wget timed out or failed for $GENE_ID"
                OFFICIAL_SUMMARY=""
            else
                OFFICIAL_SUMMARY=$(xmllint --html --xpath 'string(//dt[normalize-space(text())="Summary"]/following-sibling::dd[1])' gene_page.html 2>/dev/null || true)
                if [ -z "$OFFICIAL_SUMMARY" ]; then
                    OFFICIAL_SUMMARY=""
                    echo "No HTML summary found for GENE_ID: $GENE_ID"
                else
                    echo "HTML summary found for GENE_ID: $GENE_ID"
                fi
            fi
        else
            echo "Summary found in docsum XML for GENE_ID: $GENE_ID"
        fi
    fi

    [ -z "$OFFICIAL_NAME" ] && OFFICIAL_NAME="NA"

    OFFICIAL_NAME=$(echo "$OFFICIAL_NAME" | sed 's/"/""/g')
    OFFICIAL_SUMMARY=$(echo "$OFFICIAL_SUMMARY" | sed 's/"/""/g')

    echo "$GENE_ID,\"$OFFICIAL_NAME\",\"$OFFICIAL_SUMMARY\"" >> "$TEMP_PARSED"
    echo "Finished processing GENE_ID: $GENE_ID"
    sleep 0.3
done < "$TEMP_IDS"

echo "Merging summary info into final CSV..."
awk -F',' 'NR==FNR {
    if(NR>1){
        gene_id=$1; fullname=$2; summary=$3
        names[gene_id]=fullname
        sums[gene_id]=summary
    }
    next
}
NR>FNR {
    gene_id=$1
    if(gene_id in names){
        print $0","names[gene_id]","sums[gene_id]
    } else {
        print $0",NA,NA"
    }
}' "$TEMP_PARSED" initial_genes.csv > "$OUTPUT_CSV"

echo "Processing complete. Output for all genes: $OUTPUT_CSV"

rm -f "$TEMP_IDS" "$TEMP_XML" "$TEMP_PARSED" gene_page.html
echo "Script finished."
