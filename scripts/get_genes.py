#!/usr/bin/env python3
import requests
import gzip
import os
from pyfaidx import Fasta
import csv
import time
import logging
from typing import Dict, List

# Toggle between human and E. coli by commenting/uncommenting:
# ORGANISM = "human"
ORGANISM = "ecoli"

# Configuration dictionary for different organisms
ORGANISM_CONFIG = {
    "human": {
        "log_dir": "./human_reference",
        "accession": "GCF_000001405.40",
        "base_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/",
        "assembly": "GRCh38.p14",
        "gff3_suffix": "_genomic.gff.gz",
        "fasta_suffix": "_genomic.fna.gz",
        "output_prefix": "human_genes",
        "email": "jt828@cornell.edu"
    },
    "ecoli": {
        "log_dir": "./ecoli_reference",
        "accession": "GCF_000005845.2",
        "base_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/",
        "assembly": "ASM584v2",
        "gff3_suffix": "_genomic.gff.gz",
        "fasta_suffix": "_genomic.fna.gz",
        "output_prefix": "ecoli_genes",
        "email": "jt828@cornell.edu"
    }
}

def setup_logging(config: Dict) -> None:
    """Set up logging configuration"""
    log_dir = config["log_dir"]
    os.makedirs(log_dir, exist_ok=True)
    log_file = os.path.join(log_dir, f"{config['output_prefix']}.log")

    logging.basicConfig(
        filename=log_file,
        filemode='a',
        format='%(asctime)s [%(levelname)s] %(message)s',
        level=logging.INFO
    )

    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s [%(levelname)s] %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    console.setFormatter(formatter)
    logging.getLogger('').addHandler(console)

def get_file_paths(config: Dict) -> tuple:
    """Generate file paths based on configuration"""
    accession = config["accession"]
    assembly = config["assembly"]
    
    gff3_file_name = f"{accession}_{assembly}{config['gff3_suffix']}"
    fasta_file_name = f"{accession}_{assembly}{config['fasta_suffix']}"
    
    gff3_url = config["base_url"] + gff3_file_name
    fasta_url = config["base_url"] + fasta_file_name
    
    gff3_path = os.path.join(config["log_dir"], gff3_file_name)
    fasta_path = os.path.join(config["log_dir"], fasta_file_name)
    
    return gff3_url, fasta_url, gff3_path, fasta_path

def download_file(url: str, dest: str) -> None:
    """Download file with progress logging"""
    if not os.path.exists(dest):
        logging.info(f"Downloading {url} ...")
        try:
            r = requests.get(url, stream=True)
            r.raise_for_status()
            total_size = int(r.headers.get('content-length', 0))
            block_size = 8192
            downloaded = 0
            
            with open(dest, 'wb') as f:
                for chunk in r.iter_content(chunk_size=block_size):
                    f.write(chunk)
                    downloaded += len(chunk)
                    if total_size > 0:
                        percent = (downloaded / total_size) * 100
                        logging.info(f"Download progress: {percent:.1f}%")
                        
            logging.info(f"Downloaded {dest}")
        except Exception as e:
            logging.error(f"Error downloading {url}: {str(e)}")
            raise
    else:
        logging.info(f"{dest} already exists. Skipping download.")

def parse_gff3(gff3_path: str) -> List[Dict]:
    """Parse GFF3 file and extract gene information"""
    genes = []
    logging.info("Parsing GFF3...")
    try:
        with gzip.open(gff3_path, 'rt') as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue
                    
                seqid, source, feature_type, start, end, score, strand, phase, attributes = parts
                
                if feature_type == "gene":
                    attr_dict = dict(attr.split('=', 1) for attr in attributes.split(';') if '=' in attr)
                    
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
                    
                    if gene_id != "NA":
                        genes.append({
                            "gene_id": gene_id,
                            "gene_name": gene_name,
                            "seqid": seqid,
                            "start": start_pos,
                            "end": end_pos,
                            "gene_len": gene_len,
                            "strand": strand
                        })
                        
        logging.info(f"Found {len(genes)} genes.")
        return genes
    except Exception as e:
        logging.error(f"Error parsing GFF3: {str(e)}")
        raise

def extract_sequences(genes: List[Dict], fasta_path: str, data_dir: str) -> List[Dict]:
    """Extract sequences for genes from FASTA file"""
    logging.info("Processing FASTA file...")
    decompressed_fasta = os.path.join(data_dir, "genome.fna")
    
    if not os.path.exists(decompressed_fasta):
        logging.info("Decompressing FASTA...")
        with gzip.open(fasta_path, 'rb') as fin, open(decompressed_fasta, 'wb') as fout:
            fout.write(fin.read())
        logging.info("Decompression complete.")

    try:
        genome = Fasta(decompressed_fasta)
        logging.info("Extracting sequences for each gene...")
        
        for g in genes:
            seq = genome[g["seqid"]][g["start"]-1:g["end"]].seq
            if g["strand"] == '-':
                seq = seq[::-1].translate(str.maketrans('ACGTacgtNn','TGCAtgcaNn'))
            g["sequence"] = seq
            
        return genes
    except Exception as e:
        logging.error(f"Error processing FASTA: {str(e)}")
        raise

def fetch_gene_summary(gene_id: str, email: str, max_retries: int = 5) -> str:
    """Fetch gene summary from NCBI E-utilities with retry logic"""
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    params = {
        "db": "gene",
        "id": gene_id,
        "retmode": "json",
        "email": email
    }

    for attempt in range(max_retries):
        try:
            resp = requests.get(base_url, params=params)
            
            if resp.status_code == 200:
                data = resp.json()
                result = data.get("result", {})
                gene_data = result.get(gene_id, {})
                summary = gene_data.get("summary", "")
                summary = summary.replace("\n", " ")
                
                words = summary.split()
                first_10_words = ' '.join(words[:10]) if words else "(no summary)"
                logging.info(f"Fetched summary for {gene_id}. First 10 words: {first_10_words}")
                return summary
                
            elif resp.status_code == 429:
                wait_time = 2 ** attempt
                logging.warning(f"Rate limit hit for {gene_id}. Waiting {wait_time}s...")
                time.sleep(wait_time)
            else:
                logging.warning(f"Error fetching {gene_id}: {resp.status_code}")
                time.sleep(1)
                
        except Exception as e:
            logging.error(f"Error in fetch_gene_summary for {gene_id}: {str(e)}")
            time.sleep(1)

    logging.warning(f"Failed to fetch summary for {gene_id} after {max_retries} attempts")
    return ""

def write_output(genes: List[Dict], config: Dict) -> None:
    """Write results to CSV file"""
    output_file = os.path.join(config["log_dir"], f"{config['output_prefix']}_with_summaries.csv")
    logging.info(f"Writing output to {output_file}...")
    
    try:
        with open(output_file, 'w', newline='', encoding='utf-8') as outfile:
            writer = csv.writer(outfile)
            writer.writerow(["gene_id", "gene_name", "start", "end", "gene_len", "sequence", "summary"])
            for g in genes:
                writer.writerow([
                    g["gene_id"],
                    g["gene_name"],
                    g["start"],
                    g["end"],
                    g["gene_len"],
                    g["sequence"],
                    g.get("summary", "")
                ])
        logging.info("Output file written successfully.")
    except Exception as e:
        logging.error(f"Error writing output: {str(e)}")
        raise

def main():
    """Main execution function"""
    # Get configuration for selected organism
    if ORGANISM not in ORGANISM_CONFIG:
        raise ValueError(f"Unsupported organism: {ORGANISM}")
    
    config = ORGANISM_CONFIG[ORGANISM]
    
    # Setup logging
    setup_logging(config)
    logging.info(f"Starting processing for {ORGANISM}")
    
    try:
        # Create output directory
        os.makedirs(config["log_dir"], exist_ok=True)
        
        # Get file paths
        gff3_url, fasta_url, gff3_path, fasta_path = get_file_paths(config)
        
        # Download files
        download_file(gff3_url, gff3_path)
        download_file(fasta_url, fasta_path)
        
        # Parse GFF3
        genes = parse_gff3(gff3_path)
        
        # Extract sequences
        genes = extract_sequences(genes, fasta_path, config["log_dir"])
        
        # test
        # genes = genes[:5]
        
        # Fetch summaries
        logging.info("Fetching gene summaries...")
        for i, g in enumerate(genes):
            g["summary"] = fetch_gene_summary(g["gene_id"], config["email"])
            if i % 10 == 0 and i > 0:
                time.sleep(1)
        
        # Write output
        write_output(genes, config)
        
        logging.info("Processing complete.")
        
    except Exception as e:
        logging.error(f"Error in main execution: {str(e)}")
        raise

if __name__ == "__main__":
    main()