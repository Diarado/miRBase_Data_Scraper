#!/usr/bin/env python3
# Enhanced version with logging, retries, email parameter, and partial summary logging.
import requests
import gzip
import os
from pyfaidx import Fasta
import csv
import time
import logging

# Set up logging
LOG_DIR = "./human_reference"
os.makedirs(LOG_DIR, exist_ok=True)
LOG_FILE = os.path.join(LOG_DIR, "human_reference.log")

logging.basicConfig(
    filename=LOG_FILE,
    filemode='a',
    format='%(asctime)s [%(levelname)s] %(message)s',
    level=logging.INFO
)

console = logging.StreamHandler()
console.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s [%(levelname)s] %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
console.setFormatter(formatter)
logging.getLogger('').addHandler(console)

# Constants
DATA_DIR = "./human_reference"
ACCESSION = "GCF_000001405.40"  # Human GRCh38.p14
BASE_URL = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/"
GFF3_FILE_NAME = ACCESSION + "_GRCh38.p14_genomic.gff.gz"
FASTA_FILE_NAME = ACCESSION + "_GRCh38.p14_genomic.fna.gz"

os.makedirs(DATA_DIR, exist_ok=True)

def download_file(url, dest):
    if not os.path.exists(dest):
        logging.info(f"Downloading {url} ...")
        r = requests.get(url, stream=True)
        r.raise_for_status()
        with open(dest, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
        logging.info(f"Downloaded {dest}")
    else:
        logging.info(f"{dest} already exists. Skipping download.")

gff3_url = BASE_URL + GFF3_FILE_NAME
fasta_url = BASE_URL + FASTA_FILE_NAME

gff3_path = os.path.join(DATA_DIR, GFF3_FILE_NAME)
fasta_path = os.path.join(DATA_DIR, FASTA_FILE_NAME)

# ----------------------------
# Step 1: Download GFF3 and FASTA
# ----------------------------
download_file(gff3_url, gff3_path)
download_file(fasta_url, fasta_path)

# ----------------------------
# Step 2: Parse GFF3 to extract gene annotations
# We'll extract only "gene" features and find GeneID and Name
# ----------------------------
genes = []
logging.info("Parsing GFF3...")
with gzip.open(gff3_path, 'rt') as f:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.strip().split('\t')
        if len(parts) < 9:
            continue
        seqid, source, feature_type, start, end, score, strand, phase, attributes = parts
        if feature_type == "gene":
            attr_dict = {}
            for attr in attributes.split(';'):
                if '=' in attr:
                    key, val = attr.split('=', 1)
                    attr_dict[key] = val
            
            # Extract gene_id from Dbxref if present
            gene_id = "NA"
            if 'Dbxref' in attr_dict:
                dbxrefs = attr_dict['Dbxref'].split(',')
                for dbx in dbxrefs:
                    if dbx.startswith("GeneID:"):
                        gene_id = dbx.split(':')[1]
                        break
            
            gene_name = attr_dict.get('Name', 'NA')
            start_pos = int(start)
            end_pos = int(end)
            gene_len = end_pos - start_pos + 1
            strand_symbol = strand
            
            if gene_id != "NA":
                genes.append({
                    "gene_id": gene_id,
                    "gene_name": gene_name,
                    "seqid": seqid,
                    "start": start_pos,
                    "end": end_pos,
                    "gene_len": gene_len,
                    "strand": strand_symbol
                })

logging.info(f"Found {len(genes)} genes.")

# debug usage
# genes = genes[:50]
# logging.info("Processing only the first 50 genes for testing purposes.")

# ----------------------------
# Step 3: Extract sequences from FASTA
# ----------------------------
logging.info("Indexing FASTA with pyfaidx...")
decompressed_fasta = os.path.join(DATA_DIR, "genome.fna")
if not os.path.exists(decompressed_fasta):
    logging.info("Decompressing FASTA...")
    with gzip.open(fasta_path, 'rb') as fin, open(decompressed_fasta, 'wb') as fout:
        fout.write(fin.read())
    logging.info("Decompression complete.")

logging.info("Loading FASTA...")
genome = Fasta(decompressed_fasta)

logging.info("Extracting sequences for each gene...")
for g in genes:
    seq = genome[g["seqid"]][g["start"]-1:g["end"]].seq
    # Reverse complement if strand is '-'
    if g["strand"] == '-':
        seq = seq[::-1].translate(str.maketrans('ACGTacgtNn','TGCAtgcaNn'))
    g["sequence"] = seq

# ----------------------------
# Step 4: Fetch summary for each gene using E-utilities with retry
# Include email, exponential backoff, log URL, and partial summary logging.
# ----------------------------
def fetch_gene_summary(gene_id, max_retries=5):
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    params = {
        "db": "gene",
        "id": gene_id,
        "retmode": "json",
        "email": "jt828@cornell.edu"  # using requested email
    }

    for attempt in range(max_retries):
        resp = requests.get(base_url, params=params)
        req_url = resp.url  # get the final requested URL for logging
        if resp.status_code == 200:
            data = resp.json()
            result = data.get("result", {})
            gene_data = result.get(gene_id, {})
            summary = gene_data.get("summary", "")
            summary = summary.replace("\n", " ")
            # Log first 10 words of summary if available
            words = summary.split()
            first_10_words = ' '.join(words[:10]) if words else "(no summary)"
            logging.info(f"Successfully fetched summary for {gene_id} from {req_url}. First 10 words: {first_10_words}")
            return summary
        elif resp.status_code == 429:
            # Too many requests: wait and retry
            wait_time = 2 ** attempt
            logging.warning(f"429 Too Many Requests for {gene_id} at {req_url}. Waiting {wait_time} seconds before retry {attempt+1}...")
            time.sleep(wait_time)
        else:
            # Some other error
            logging.warning(f"Could not fetch summary for {gene_id} at {req_url}: {resp.status_code} {resp.reason}")
            return ""

    # If all retries fail
    logging.warning(f"Exhausted retries for {gene_id} after {max_retries} attempts. No summary available.")
    return ""

logging.info("Fetching summaries from NCBI...")
for i, g in enumerate(genes):
    g["summary"] = fetch_gene_summary(g["gene_id"])

    # Add a small delay every 10 requests to reduce rate-limiting
    if i % 10 == 0 and i > 0:
        time.sleep(1)

# ----------------------------
# Step 5: Write output as CSV
# ----------------------------
OUTPUT_CSV = os.path.join(DATA_DIR, "human_genes_with_summaries.csv")
logging.info(f"Writing output to {OUTPUT_CSV}...")
with open(OUTPUT_CSV, 'w', newline='', encoding='utf-8') as outfile:
    writer = csv.writer(outfile)  # Default delimiter is comma
    writer.writerow(["gene_id", "gene_name", "start", "end", "gene_len", "sequence", "summary"])
    for g in genes:
        writer.writerow([g["gene_id"], g["gene_name"], g["start"], g["end"], g["gene_len"], g["sequence"], g.get("summary", "")])

logging.info("Processing complete.")
