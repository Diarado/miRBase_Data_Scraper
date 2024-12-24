# get amino acid sequence and description, just like get_genes
# We only collect all data after 01/01/2024
from Bio import Entrez
from Bio import SeqIO
import time
from datetime import datetime
import csv
import pandas as pd
import logging
import sys
from pathlib import Path

def setup_logging(timestamp):
    """
    Set up logging configuration
    """
    Path('logs').mkdir(exist_ok=True)
    
    log_filename = f'logs/protein_fetch_{timestamp}.log'
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_filename),
            logging.StreamHandler(sys.stdout)
        ]
    )
    
    logging.info("=== Script Started ===")
    logging.info(f"Log file created at: {log_filename}")
    
    return log_filename

def search_recent_proteins(start=0):
    """
    Search for proteins added or updated after the specified date with pagination
    """
    logging.info(f"Searching for proteins modified after {START_DATE} - batch starting at {start}")
    Entrez.email = EMAIL
    
    try:
        search_handle = Entrez.esearch(
            db="protein",
            term=f"{START_DATE}[MDAT]:3000[MDAT]",
            retmax=BATCH_SIZE,
            retstart=start,
            sort="modification date"
        )
        
        search_results = Entrez.read(search_handle)
        search_handle.close()
        
        total_count = int(search_results["Count"])
        protein_ids = search_results["IdList"]
        
        logging.info(f"Found {len(protein_ids)} proteins in current batch (Total: {total_count})")
        return protein_ids, total_count
    
    except Exception as e:
        logging.error(f"Error in search: {str(e)}")
        return [], 0

def process_protein_batch(protein_ids, existing_data=None):
    """
    Process a batch of protein IDs
    """
    protein_data_list = []
    successful_fetches = 0
    failed_fetches = 0
    
    for protein_id in protein_ids:
        time.sleep(0.5)  # NCBI rate limit compliance
        
        # Skip if already processed
        if existing_data is not None and any(d['protein_id'] == protein_id for d in existing_data):
            logging.info(f"Skipping already processed protein ID: {protein_id}")
            continue
            
        protein_data = get_protein_info(protein_id)
        if protein_data:
            protein_data_list.append(protein_data)
            successful_fetches += 1
        else:
            failed_fetches += 1
    
    return protein_data_list, successful_fetches, failed_fetches

def get_protein_info(protein_id):
    """
    Get detailed information for a specific protein
    """
    logging.info(f"Fetching information for protein ID: {protein_id}")
    
    try:
        handle = Entrez.efetch(
            db="protein",
            id=protein_id,
            rettype="gb",
            retmode="text"
        )
        
        record = SeqIO.read(handle, "genbank")
        handle.close()
        
        features_list = []
        for feature in record.features:
            if feature.type in ["Region", "Domain", "Site"]:
                if 'note' in feature.qualifiers:
                    features_list.append(feature.qualifiers['note'][0])
        
        function = record.annotations.get('function', '')
        if not function:
            for feature in record.features:
                if 'function' in feature.qualifiers:
                    function = feature.qualifiers['function'][0]
                    break
        
        protein_data = {
            'protein_id': protein_id,
            'sequence': str(record.seq),
            'description': record.description,
            'date_modified': record.annotations.get('date', 'Not available'),
            'organism': record.annotations.get('organism', 'Not available'),
            'length': len(record.seq),
            'function': function,
            'features': '; '.join(features_list),
            'gene_name': record.annotations.get('gene_name', ''),
            'taxonomy': '; '.join(record.annotations.get('taxonomy', [])),
            'accessions': '; '.join(record.annotations.get('accessions', [])),
        }
        
        return protein_data
        
    except Exception as e:
        logging.error(f"Error fetching protein {protein_id}: {str(e)}")
        return None

def save_intermediate_results(protein_data_list, timestamp):
    """
    Save intermediate results to prevent data loss
    """
    intermediate_file = f'output/protein_data_intermediate_{timestamp}.csv'
    try:
        if protein_data_list:
            df = pd.DataFrame(protein_data_list)
            df.to_csv(intermediate_file, index=False)
            logging.info(f"Saved intermediate results to {intermediate_file}")
    except Exception as e:
        logging.error(f"Error saving intermediate results: {str(e)}")

def load_intermediate_results(timestamp):
    """
    Load intermediate results if they exist
    """
    intermediate_file = f'output/protein_data_intermediate_{timestamp}.csv'
    try:
        if Path(intermediate_file).exists():
            df = pd.read_csv(intermediate_file)
            return df.to_dict('records')
    except Exception as e:
        logging.error(f"Error loading intermediate results: {str(e)}")
    return []

def save_final_results(protein_data_list, timestamp):
    """
    Save final results to CSV file
    """
    csv_filename = f'output/protein_data_{timestamp}.csv'
    try:
        df = pd.DataFrame(protein_data_list)
        columns_order = [
            'protein_id', 'description', 'gene_name', 'organism',
            'length', 'function', 'features', 'taxonomy',
            'date_modified', 'accessions', 'sequence'
        ]
        df = df[columns_order]
        df.to_csv(csv_filename, index=False)
        logging.info(f"Saved final CSV to {csv_filename}")
        return csv_filename
    except Exception as e:
        logging.error(f"Error saving final CSV: {str(e)}")
        return None

# Configuration
BATCH_SIZE = 5000  # Number of proteins to process in each batch
EMAIL = "your_email@example.com"  # Your email for NCBI
START_DATE = "2024/01/01"  # Start date for protein search

def main():
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_filename = setup_logging(timestamp)
    Path('output').mkdir(exist_ok=True)
    
    try:
        start_position = 0
        total_successful_fetches = 0
        total_failed_fetches = 0
        all_protein_data = load_intermediate_results(timestamp)
        
        # Get initial batch and total count
        protein_ids, total_count = search_recent_proteins(start=start_position)
        
        while protein_ids:
            logging.info(f"\n=== Processing batch starting at position {start_position} ===")
            
            # Process current batch
            batch_data, successful, failed = process_protein_batch(protein_ids, all_protein_data)
            
            # Update totals
            total_successful_fetches += successful
            total_failed_fetches += failed
            all_protein_data.extend(batch_data)
            
            # Save intermediate results
            save_intermediate_results(all_protein_data, timestamp)
            
            # Progress report
            progress = (start_position + BATCH_SIZE) / total_count * 100
            logging.info(f"Progress: {progress:.1f}% ({start_position + BATCH_SIZE}/{total_count})")
            
            # Get next batch
            start_position += BATCH_SIZE
            protein_ids, _ = search_recent_proteins(start=start_position)
            
            if not protein_ids or start_position >= total_count:
                break
        
        # Save final results
        save_final_results(all_protein_data, timestamp)
        
        # Log summary statistics
        logging.info("\n=== Final Summary ===")
        logging.info(f"Total proteins processed: {total_successful_fetches + total_failed_fetches}")
        logging.info(f"Successful fetches: {total_successful_fetches}")
        logging.info(f"Failed fetches: {total_failed_fetches}")
        logging.info(f"Success rate: {(total_successful_fetches/(total_successful_fetches + total_failed_fetches))*100:.1f}%")
        
    except Exception as e:
        logging.error(f"Critical error in main execution: {str(e)}")
    finally:
        logging.info("=== Script Finished ===")

if __name__ == "__main__":
    main()