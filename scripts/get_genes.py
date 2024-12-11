#!/usr/bin/env python3
import requests
import gzip
import os
from pyfaidx import Fasta
import csv
import time
import logging
from typing import Dict, List
import math
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry

# Toggle between human and E. coli by commenting/uncommenting:
ORGANISM = "human"
# ORGANISM = "ecoli"

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
        "email": "jt828@cornell.edu",
        "batch_size": 5000
    },
    "ecoli": {
        "log_dir": "./ecoli_reference",
        "accession": "GCF_000005845.2",
        "base_url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/",
        "assembly": "ASM584v2",
        "gff3_suffix": "_genomic.gff.gz",
        "fasta_suffix": "_genomic.fna.gz",
        "output_prefix": "ecoli_genes",
        "email": "jt828@cornell.edu",
        "batch_size": 5000
    }
}

def create_session_with_retries() -> requests.Session:
    """Create a session with retry logic"""
    session = requests.Session()
    retries = Retry(
        total=5,  # number of retries
        backoff_factor=1,  # wait 1, 2, 4, 8, 16 seconds between retries
        status_forcelist=[500, 502, 503, 504]  # retry on these HTTP status codes
    )
    session.mount('https://', HTTPAdapter(max_retries=retries))
    return session

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
            session = create_session_with_retries()
            r = session.get(url, stream=True)
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
                        if downloaded % (5 * block_size) == 0:  # Log every 5 blocks
                            logging.info(f"Download progress: {percent:.1f}%")
                        
            logging.info(f"Downloaded {dest}")
        except Exception as e:
            logging.error(f"Error downloading {url}: {str(e)}")
            if os.path.exists(dest):
                os.remove(dest)  # Clean up partial download
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
    
    session = create_session_with_retries()

    for attempt in range(max_retries):
        try:
            resp = session.get(base_url, params=params)
            
            if resp.status_code == 200:
                data = resp.json()
                result = data.get("result", {})
                gene_data = result.get(gene_id, {})
                summary = gene_data.get("summary", "")
                summary = summary.replace("\n", " ")
                
                words = summary.split()
                first_10_words = ' '.join(words[:10]) if words else "(no summary)"
                logging.info(f"Successfully fetched summary for {gene_id} from {base_url}. First 10 words: {first_10_words}")
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

def write_batch_output(genes: List[Dict], batch_num: int, config: Dict) -> None:
    """Write results for a single batch to CSV file"""
    output_file = os.path.join(config["log_dir"], f"{config['output_prefix']}_batch_{batch_num}.csv")
    logging.info(f"Writing batch {batch_num} to {output_file}...")
    
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
        logging.info(f"Batch {batch_num} written successfully.")
    except Exception as e:
        logging.error(f"Error writing batch {batch_num}: {str(e)}")
        raise

def process_gene_batch(genes: List[Dict], batch_num: int, config: Dict) -> None:
    """Process a batch of genes including fetching summaries and writing output"""
    logging.info(f"Processing batch {batch_num} with {len(genes)} genes...")
    
    # Fetch summaries for the batch
    for i, g in enumerate(genes):
        g["summary"] = fetch_gene_summary(g["gene_id"], config["email"])
        if i % 10 == 0 and i > 0:
            time.sleep(1)  # Rate limiting
    
    # Write batch output
    write_batch_output(genes, batch_num, config)

def combine_batch_files(config: Dict, total_batches: int) -> None:
    """Combine all batch files into a single output file"""
    final_output = os.path.join(config["log_dir"], f"{config['output_prefix']}_combined.csv")
    logging.info(f"Combining {total_batches} batch files into {final_output}...")
    
    try:
        # Write header to combined file
        with open(final_output, 'w', newline='', encoding='utf-8') as outfile:
            writer = csv.writer(outfile)
            writer.writerow(["gene_id", "gene_name", "start", "end", "gene_len", "sequence", "summary"])
            
            # Append each batch file
            for batch_num in range(total_batches):
                batch_file = os.path.join(config["log_dir"], f"{config['output_prefix']}_batch_{batch_num}.csv")
                if os.path.exists(batch_file):
                    with open(batch_file, 'r', encoding='utf-8') as batch:
                        next(batch)  # Skip header
                        for line in batch:
                            outfile.write(line)
                    
                    # Optionally remove batch file after combining
                    os.remove(batch_file)
                else:
                    logging.warning(f"Batch file {batch_file} not found")
                
        logging.info("All batch files combined successfully.")
    except Exception as e:
        logging.error(f"Error combining batch files: {str(e)}")
        raise

def main():
    """Main execution function with batch processing"""
    if ORGANISM not in ORGANISM_CONFIG:
        raise ValueError(f"Unsupported organism: {ORGANISM}")
    
    config = ORGANISM_CONFIG[ORGANISM]
    setup_logging(config)
    logging.info(f"Starting batch processing for {ORGANISM}")
    
    try:
        # Create output directory
        os.makedirs(config["log_dir"], exist_ok=True)
        
        # Get and download files
        gff3_url, fasta_url, gff3_path, fasta_path = get_file_paths(config)
        download_file(gff3_url, gff3_path)
        download_file(fasta_url, fasta_path)
        
        # Parse GFF3 and extract sequences
        genes = parse_gff3(gff3_path)
        genes = extract_sequences(genes, fasta_path, config["log_dir"])
        
        # Calculate number of batches
        batch_size = config["batch_size"]
        total_batches = math.ceil(len(genes) / batch_size)
        logging.info(f"Processing {len(genes)} genes in {total_batches} batches of {batch_size}")
        
        # Process each batch
        for batch_num in range(total_batches):
            start_idx = batch_num * batch_size
            end_idx = min((batch_num + 1) * batch_size, len(genes))
            batch_genes = genes[start_idx:end_idx]
            
            # Skip if batch file already exists
            batch_file = os.path.join(config["log_dir"], f"{config['output_prefix']}_batch_{batch_num}.csv")
            if os.path.exists(batch_file):
                logging.info(f"Batch {batch_num} already exists, skipping...")
                continue
            
            process_gene_batch(batch_genes, batch_num, config)
            
            # Add delay between batches
            if batch_num < total_batches - 1:
                time.sleep(2)
        
        # Combine all batch files
        combine_batch_files(config, total_batches)
        
        logging.info("All processing complete.")
        
    except Exception as e:
        logging.error(f"Error in main execution: {str(e)}")
        raise

if __name__ == "__main__":
    main()