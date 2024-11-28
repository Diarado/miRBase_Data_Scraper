import requests
from bs4 import BeautifulSoup
import csv
import logging

# Configure logging
logging.basicConfig(
    filename='scraper.log',
    filemode='a',
    format='%(asctime)s - %(levelname)s - %(message)s',
    level=logging.INFO
)

def scrape_mirbase_start_end(url, output_csv='start_end_data.csv'):
    # Headers to mimic a browser visit
    headers = {
        'User-Agent': (
            'Mozilla/5.0 (Windows NT 10.0; Win64; x64) '
            'AppleWebKit/537.36 (KHTML, like Gecko) '
            'Chrome/58.0.3029.110 Safari/537.3'
        )
    }

    try:
        logging.info(f"Sending GET request to {url}")
        response = requests.get(url, headers=headers)
        response.raise_for_status()  # Raise exception for HTTP errors

        # Parse the HTML content
        soup = BeautifulSoup(response.text, 'html.parser')

        # Find the table containing the results by ID
        table = soup.find('table', {'id': 'results-table'})

        if not table:
            logging.error("Could not find the results table on the page.")
            print("Could not find the results table on the page.")
            return

        # Extract table headers to identify the column indices for Start and End
        thead = table.find('thead')
        if not thead:
            logging.error("The table does not contain a <thead> section.")
            print("The table does not contain a <thead> section.")
            return

        headers = thead.find_all('th')
        header_titles = [header.get_text(strip=True) for header in headers]

        try:
            start_idx = header_titles.index('Start')
            end_idx = header_titles.index('End')
        except ValueError:
            logging.error("Could not find 'Start' or 'End' columns in the table headers.")
            print("Could not find 'Start' or 'End' columns in the table headers.")
            return

        # Iterate over each row in the table body
        tbody = table.find('tbody')
        if not tbody:
            logging.error("The table does not contain a <tbody> section.")
            print("The table does not contain a <tbody> section.")
            return

        rows = tbody.find_all('tr')

        # List to store Start and End information
        start_end_list = []

        for row in rows:
            cols = row.find_all('td')
            # Ensure there are enough columns
            if len(cols) < max(start_idx, end_idx) + 1:
                continue  # Skip rows that do not have enough columns

            # Extract text from the Start and End columns
            start = cols[start_idx].get_text(strip=True)
            end = cols[end_idx].get_text(strip=True)
            start_end_list.append({'Start': start, 'End': end})

        if not start_end_list:
            logging.warning("No Start and End data found.")
            print("No Start and End data found.")
            return

        # Write the data to a CSV file
        with open(output_csv, mode='w', newline='', encoding='utf-8') as csvfile:
            fieldnames = ['Start', 'End']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

            writer.writeheader()
            for entry in start_end_list:
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
    scrape_mirbase_start_end(mirbase_url)
