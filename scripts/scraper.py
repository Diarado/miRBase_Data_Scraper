import requests
from bs4 import BeautifulSoup
import csv
import logging
import time
from urllib.parse import urljoin

# Configure logging
logging.basicConfig(
    filename='scraper_with_limited_summary.log',
    filemode='a',
    format='%(asctime)s - %(levelname)s - %(message)s',
    level=logging.INFO
)

def get_summary(detail_url, session, headers, word_limit=50):
    """
    Fetches the summary from the miRNA detail page and limits it to the first `word_limit` words.
    
    Parameters:
        detail_url (str): The full URL to the miRNA detail page.
        session (requests.Session): The session object for persistent connections.
        headers (dict): The headers to use for the request.
        word_limit (int): Number of words to retain in the summary.
    
    Returns:
        str: The extracted and truncated summary text, or an empty string if not found.
    """
    try:
        logging.info(f"Fetching detail page: {detail_url}")
        response = session.get(detail_url, headers=headers)
        response.raise_for_status()
        soup = BeautifulSoup(response.text, 'html.parser')
        
        # Find the <pre> tag containing the description
        # Find by the containing div's style
        description_div = soup.find('div', style=lambda value: value and 'margin-left:15px' in value)
        
        if description_div:
            pre_tag = description_div.find('pre')
            if pre_tag:
                summary_text = pre_tag.get_text(separator=' ', strip=True)
                # Split the summary into words and take the first `word_limit` words
                summary_words = summary_text.split()[:word_limit]
                summary_short = ' '.join(summary_words)
                logging.info(f"Summary found: {summary_short[:60]}...")  # Log first 60 chars
                return summary_short
            else:
                logging.warning(f"<pre> tag not found within the description div for URL: {detail_url}")
                return ""
        else:
            logging.warning(f"Description div not found for URL: {detail_url}")
            return ""
    
    except requests.exceptions.HTTPError as http_err:
        logging.error(f"HTTP error while fetching {detail_url}: {http_err}")
    except Exception as err:
        logging.error(f"Error while fetching {detail_url}: {err}")
    
    return ""

def scrape_mirbase_start_end_summary(url, output_csv='miRNA_human.csv', max_rows=1917):
    """
    Scrapes Start, End, Name, and Summary (first five words) information from miRBase and saves to a CSV file.
    
    Parameters:
        url (str): The miRBase results page URL.
        output_csv (str): The name of the output CSV file.
        max_rows (int): The maximum number of rows to scrape and write to the CSV.
    """
    # Headers to mimic a browser visit
    headers = {
        'User-Agent': (
            'Mozilla/5.0 (Windows NT 10.0; Win64; x64) '
            'AppleWebKit/537.36 (KHTML, like Gecko) '
            'Chrome/58.0.3029.110 Safari/537.3'
        )
    }
    
    # Base URL to construct full URLs for detail pages
    base_url = "https://mirbase.org"

    # Use a session for persistent connections
    session = requests.Session()
    session.headers.update(headers)

    try:
        logging.info(f"Sending GET request to {url}")
        response = session.get(url)
        response.raise_for_status()
        
        # Parse the HTML content
        soup = BeautifulSoup(response.text, 'html.parser')
        
        # Find the table containing the results by ID
        table = soup.find('table', {'id': 'results-table'})
        
        if not table:
            logging.error("Could not find the results table on the page.")
            print("Could not find the results table on the page.")
            return
        
        # Extract table headers to identify the column indices for Start, End, and Name
        thead = table.find('thead')
        if not thead:
            logging.error("The table does not contain a <thead> section.")
            print("The table does not contain a <thead> section.")
            return
        
        headers_row = thead.find_all('th')
        header_titles = [header.get_text(strip=True) for header in headers_row]
        
        try:
            start_idx = header_titles.index('Start')
            end_idx = header_titles.index('End')
            name_idx = header_titles.index('Name')
        except ValueError as e:
            logging.error(f"Could not find required columns in the table headers: {e}")
            print("Could not find 'Start', 'End', or 'Name' columns in the table headers.")
            return
        
        # Iterate over each row in the table body
        tbody = table.find('tbody')
        if not tbody:
            logging.error("The table does not contain a <tbody> section.")
            print("The table does not contain a <tbody> section.")
            return
        
        rows = tbody.find_all('tr')
        logging.info(f"Found {len(rows)} rows in the table.")
        
        # List to store Name, Start, End, and Summary information
        data_list = []
        
        for idx, row in enumerate(rows, start=1):
            # **LIMITING TO FIRST `max_rows` ROWS**
            if idx > max_rows:
                logging.info(f"Reached the maximum limit of {max_rows} rows. Stopping the scraper.")
                break  # Stop processing after max_rows
            
            cols = row.find_all('td')
            # Ensure there are enough columns
            if len(cols) < max(start_idx, end_idx, name_idx) + 1:
                logging.warning(f"Row {idx} does not have enough columns. Skipping.")
                print(f"Row {idx}: Not enough columns. Skipping.")
                continue  # Skip rows that do not have enough columns
            
            # Extract Start and End
            start = cols[start_idx].get_text(strip=True)
            end = cols[end_idx].get_text(strip=True)
            
            # Extract Name and link
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
            
            # Fetch Summary from detail page if link exists
            summary = ""
            if detail_url:
                summary = get_summary(detail_url, session, headers)
                # Optional: Delay between requests to be polite to the server
                time.sleep(1)  # Sleep for 1 second
            else:
                logging.info(f"No detail URL for Name: {name}. Summary will be empty.")
            
            # Append the data to the list
            data_list.append({
                'Name': name,
                'Start': start,
                'End': end,
                'Summary': summary
            })
            
            # Print progress
            print(f"Processed {idx} / {max_rows} rows.")
        
        if not data_list:
            logging.warning("No data extracted from the table.")
            print("No data extracted from the table.")
            return
        
        # Write the data to a CSV file
        with open(output_csv, mode='w', newline='', encoding='utf-8') as csvfile:
            fieldnames = ['Name', 'Start', 'End', 'Summary']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            
            writer.writeheader()
            for entry in data_list:
                writer.writerow(entry)
        
        logging.info(f"Data successfully saved to '{output_csv}'.")
        print(f"Data successfully saved to '{output_csv}'.")
    
    except requests.exceptions.HTTPError as http_err:
        logging.error(f"HTTP error occurred: {http_err}")
        print(f"HTTP error occurred: {http_err}")  # HTTP error
    except Exception as err:
        logging.error(f"An error occurred: {err}")
        print(f"An error occurred: {err}")  # Other errors

if __name__ == "__main__":
    mirbase_url = "https://mirbase.org/browse/results/?organism=hsa"
    scrape_mirbase_start_end_summary(mirbase_url)
