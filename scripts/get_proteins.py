# get amino acid sequence and description, just like get_genes
# We only collect data from 1.1, 2024 to today
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
    # Create logs directory if it doesn't exist
    Path('logs').mkdir(exist_ok=True)
    
    # Configure logging
    log_filename = f'logs/protein_fetch_{timestamp}.log'
    
    # Set up logging to both file and console
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

def search_recent_proteins(start_date="2024/01/01"):
    """
    Search for proteins added or updated after the specified date
    """
    logging.info(f"Searching for proteins modified after {start_date}")
    Entrez.email = "your_email@example.com"  
    
    try:
        search_handle = Entrez.esearch(
            db="protein",
            term=f"{start_date}[MDAT]:3000[MDAT]",
            retmax=10, # TODO
            sort="modification date"
        )
        
        search_results = Entrez.read(search_handle)
        search_handle.close()
        
        protein_count = len(search_results["IdList"])
        logging.info(f"Found {protein_count} proteins")
        return search_results["IdList"]
    
    except Exception as e:
        logging.error(f"Error in search: {str(e)}")
        return []

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
        
        # Extract features and annotations
        features_list = []
        for feature in record.features:
            if feature.type in ["Region", "Domain", "Site"]:
                if 'note' in feature.qualifiers:
                    features_list.append(feature.qualifiers['note'][0])
        
        # Get protein function from annotations
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
        
        logging.info(f"Successfully retrieved data for protein ID: {protein_id}")
        return protein_data
        
    except Exception as e:
        logging.error(f"Error fetching protein {protein_id}: {str(e)}")
        return None

def save_sequences_file(protein_data_list, timestamp):
    """
    Save just the sequences in a separate file
    """
    fasta_filename = f'protein_sequences.fasta'
    logging.info(f"Saving sequences to FASTA file: {fasta_filename}")
    
    try:
        with open(fasta_filename, 'w') as f:
            for data in protein_data_list:
                if data:
                    f.write(f">{data['protein_id']} {data['description']}\n")
                    sequence = data['sequence']
                    for i in range(0, len(sequence), 60):
                        f.write(sequence[i:i+60] + '\n')
        
        logging.info(f"Successfully saved {len(protein_data_list)} sequences to FASTA file")
        return fasta_filename
    except Exception as e:
        logging.error(f"Error saving FASTA file: {str(e)}")
        return None

def save_to_csv(protein_data_list, timestamp):
    """
    Save protein data to CSV file
    """
    csv_filename = f'protein_data.csv'
    logging.info(f"Saving data to CSV file: {csv_filename}")
    
    try:
        df = pd.DataFrame(protein_data_list)
        
        # Reorder columns for better readability
        columns_order = [
            'protein_id', 'description', 'gene_name', 'organism',
            'length', 'function', 'features', 'taxonomy',
            'date_modified', 'accessions', 'sequence'
        ]
        
        df = df[columns_order]
        df.to_csv(csv_filename, index=False)
        
        logging.info(f"Successfully saved {len(protein_data_list)} entries to CSV file")
        return csv_filename
    except Exception as e:
        logging.error(f"Error saving CSV file: {str(e)}")
        return None

def main():
    # Generate timestamp for file names
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Set up logging
    log_filename = setup_logging(timestamp)
    
    # Create output directory if it doesn't exist
    Path('output').mkdir(exist_ok=True)
    
    try:
        # Search for proteins
        protein_ids = search_recent_proteins()
        
        if not protein_ids:
            logging.error("No proteins found or search failed")
            return
        
        # List to store all protein data
        protein_data_list = []
        successful_fetches = 0
        failed_fetches = 0
        
        # Process each protein
        total_proteins = len(protein_ids)
        for i, protein_id in enumerate(protein_ids, 1):
            logging.info(f"Processing protein {i} of {total_proteins} ({(i/total_proteins)*100:.1f}%)")
            
            time.sleep(0.5)  # Respect NCBI rate limits
            
            protein_data = get_protein_info(protein_id)
            if protein_data:
                protein_data_list.append(protein_data)
                successful_fetches += 1
            else:
                failed_fetches += 1
        
        # Save results
        csv_file = save_to_csv(protein_data_list, timestamp)
        # fasta_file = save_sequences_file(protein_data_list, timestamp)
        
        # Log summary statistics
        logging.info("\n=== Summary Statistics ===")
        logging.info(f"Total proteins processed: {total_proteins}")
        logging.info(f"Successful fetches: {successful_fetches}")
        logging.info(f"Failed fetches: {failed_fetches}")
        logging.info(f"Success rate: {(successful_fetches/total_proteins)*100:.1f}%")
        
        logging.info("\n=== Output Files ===")
        logging.info(f"CSV file: {csv_file}")
        # logging.info(f"FASTA file: {fasta_file}")
        logging.info(f"Log file: {log_filename}")
        
    except Exception as e:
        logging.error(f"Critical error in main execution: {str(e)}")
    finally:
        logging.info("=== Script Finished ===")

if __name__ == "__main__":
    main()