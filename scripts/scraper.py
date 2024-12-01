import requests
from bs4 import BeautifulSoup
import csv
import logging
import time
from urllib.parse import urljoin

# Configure logging
logging.basicConfig(
    filename='scraper_with_sequence.log',
    filemode='a',
    format='%(asctime)s - %(levelname)s - %(message)s',
    level=logging.INFO
)

def get_summary_and_sequence(detail_url, session, headers):
    """
    Fetches the summary and sequence from the miRNA detail page.
    Empty summary is considered valid.
    
    Parameters:
        detail_url (str): The full URL to the miRNA detail page
        session (requests.Session): The session object for persistent connections
        headers (dict): The headers to use for the request
    
    Returns:
        tuple: (summary_text, sequence_text)
    """
    try:
        logging.info(f"Fetching detail page: {detail_url}")
        response = session.get(detail_url, headers=headers, timeout=10)
        response.raise_for_status()
        soup = BeautifulSoup(response.text, 'html.parser')
        
        # Initialize variables
        summary_text = ""
        sequence_text = ""
        
        # Try to find the description - it's okay if we don't find it
        description_div = soup.find('div', style=lambda value: value and 'margin-left:15px' in value)
        if description_div:
            pre_tag = description_div.find('pre')
            if pre_tag:
                summary_text = pre_tag.get_text(separator=' ', strip=True)
                logging.info(f"Summary found: {summary_text[:50]}...")
        
        # Try to find the sequence
        sequence_div = soup.find('div', {'id': 'hairpinSequence'})
        if sequence_div:
            sequence_span = sequence_div.find('span', {'class': 'text-monospace'})
            if sequence_span:
                sequence_text = sequence_span.get_text(strip=True)
                logging.info(f"Sequence found: {sequence_text[:50]}...")
        
        return summary_text, sequence_text
            
    except requests.exceptions.Timeout:
        logging.error(f"Timeout while fetching {detail_url}")
        return "", ""
        
    except requests.exceptions.HTTPError as http_err:
        logging.error(f"HTTP error while fetching {detail_url}: {http_err}")
        if http_err.response.status_code == 429:  # Too Many Requests
            logging.info("Rate limited, waiting 60 seconds...")
            time.sleep(60)
            return get_summary_and_sequence(detail_url, session, headers)
        return "", ""
        
    except Exception as err:
        logging.error(f"Error while fetching {detail_url}: {err}")
        return "", ""

def scrape_mirbase_with_sequence(url, output_csv='miRNA_human_with_sequence.csv', max_rows=3000):
    """
    Scrapes Start, End, Name, Sequence, and Summary information from miRBase and saves to a CSV file.
    """
    headers = {
        'User-Agent': (
            'Mozilla/5.0 (Windows NT 10.0; Win64; x64) '
            'AppleWebKit/537.36 (KHTML, like Gecko) '
            'Chrome/58.0.3029.110 Safari/537.3'
        )
    }
    
    base_url = "https://mirbase.org"
    session = requests.Session()
    session.headers.update(headers)

    try:
        logging.info(f"Sending GET request to {url}")
        response = session.get(url)
        response.raise_for_status()
        
        soup = BeautifulSoup(response.text, 'html.parser')
        table = soup.find('table', {'id': 'results-table'})
        
        if not table:
            logging.error("Could not find the results table on the page.")
            return
        
        thead = table.find('thead')
        if not thead:
            logging.error("The table does not contain a <thead> section.")
            return
        
        headers_row = thead.find_all('th')
        header_titles = [header.get_text(strip=True) for header in headers_row]
        
        try:
            start_idx = header_titles.index('Start')
            end_idx = header_titles.index('End')
            name_idx = header_titles.index('Name')
            chr_idx = header_titles.index('Chromosome')
            strand_idx = header_titles.index('Strand')
        except ValueError as e:
            logging.error(f"Could not find required columns: {e}")
            return
        
        tbody = table.find('tbody')
        if not tbody:
            logging.error("The table does not contain a <tbody> section.")
            return
        
        rows = tbody.find_all('tr')
        logging.info(f"Found {len(rows)} rows in the table.")
        
        data_list = []
        
        for idx, row in enumerate(rows, start=1):
            if idx > max_rows:
                logging.info(f"Reached the maximum limit of {max_rows} rows.")
                break
            
            cols = row.find_all('td')
            if len(cols) < max(start_idx, end_idx, name_idx) + 1:
                logging.warning(f"Row {idx} does not have enough columns. Skipping.")
                continue
            
            start = cols[start_idx].get_text(strip=True)
            end = cols[end_idx].get_text(strip=True)
            chromosome = cols[chr_idx].get_text(strip=True)
            strand = cols[strand_idx].get_text(strip=True)
            
            
            name_cell = cols[name_idx]
            name_link = name_cell.find('a')
            if name_link and 'href' in name_link.attrs:
                name = name_link.get_text(strip=True)
                relative_link = name_link['href']
                detail_url = urljoin(base_url, relative_link)
            else:
                name = name_cell.get_text(strip=True)
                detail_url = ""
                logging.warning(f"Row {idx} does not have a valid link for Name: {name}")
            
            summary = ""
            sequence = ""
            if detail_url:
                summary, sequence = get_summary_and_sequence(detail_url, session, headers)
                time.sleep(1)  # Polite delay between requests
            
            seq_len = len(sequence)
            sum_len = len(summary)
            
            data_list.append({
                'Name': name,
                'Start': start,
                'End': end,
                'Sequence': sequence,
                'Seq_len': seq_len,
                'Summary': summary,
                'Summary_len': sum_len,
                'Chr': chromosome,
                'Strand': strand
            })
            
            print(f"Processed {idx} / {max_rows} rows.")
        
        if not data_list:
            logging.warning("No data extracted from the table.")
            return
        
        with open(output_csv, mode='w', newline='', encoding='utf-8') as csvfile:
            fieldnames = ['Name', 'Start', 'End', 'Sequence', 'Seq_len', 'Summary', 'Summary_len', 'Chr', 'Strand']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            
            writer.writeheader()
            for entry in data_list:
                writer.writerow(entry)
        
        logging.info(f"Data successfully saved to '{output_csv}'.")
        print(f"Data successfully saved to '{output_csv}'.")
    
    except requests.exceptions.HTTPError as http_err:
        logging.error(f"HTTP error occurred: {http_err}")
    except Exception as err:
        logging.error(f"An error occurred: {err}")

if __name__ == "__main__":
    mirbase_url = "https://mirbase.org/browse/results/?organism=hsa"
    scrape_mirbase_with_sequence(mirbase_url)