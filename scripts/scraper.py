import requests
from bs4 import BeautifulSoup
import csv
import logging
import time
from urllib.parse import urljoin
import subprocess
import pandas as pd
import numpy as np
from io import StringIO

# Configure logging
logging.basicConfig(
    filename='scraper_with_sequence.log',
    filemode='a',
    format='%(asctime)s - %(levelname)s - %(message)s',
    level=logging.INFO  
)

def get_gene_summary(gene_ids):
    if isinstance(gene_ids, list):
        gene_ids = ','.join([str(gid) for gid in gene_ids])
    
    cmd = (
        f"esearch -db gene -query {gene_ids} | "
        f"efetch -format docsum | "
        f"xtract -pattern DocumentSummary -element Id Name ScientificName TaxID CommonName Summary"
    )
    result = subprocess.run(cmd, stdout=subprocess.PIPE, shell=True, text=True)
    
    data = StringIO(result.stdout)
    names = ["gene_id", "gene_name", "org_name", "org_id", "common_org_name", "summary"]
    df = pd.read_csv(data, sep="\t", header=None, names=names)
    
    df['tag'] = np.where(
        (df.summary.isnull()) | (df.summary.astype(str).str.startswith('DISCONTINUED')),
        'no_summary',
        'has_summary'
    )
    df = df.loc[:, ['tag', 'gene_id', 'gene_name', 'org_name', 'common_org_name', 'summary']]
    
    return df

def get_dict():
    ref = {
        'Alveolata': ['Symbiodinium microadriaticum'],
        'Chromalveolata': {
            'Heterokontophyta': {}
        },
        'Metazoa': {
            'Bilateria': {
                'Deuterostoma': {
                    'Xenoturbella bocki': ['Xenoturbella bocki'],
                    'Chordata': {
                        'Cephalochordata': ['Branchiostoma belcheri', 'Branchiostoma floridae'],
                        'Urochordata': ['Ciona intestinalis', 'Ciona savignyi', 'Oikopleura dioica'],
                        'Vertebrata': {
                            'Agnathostomata': ['Petromyzon marinus'],
                            'Amphibia': ['Xenopus laevis', 'Xenopus tropicalis'],
                            'Aves': ['Anas platyrhynchos', 'Columba livia', 'Gallus gallus', 'Taeniopygia guttata'],
                            'Mammalia': {
                                'Carnivora': ['Canis familiaris'],
                                'Cingulata': ['Dasypus novemcinctus'],
                                'Lagomorpha': ['Oryctolagus cuniculus'],
                                'Laurasiatheria': ['Artibeus jamaicensis', 'Eptesicus fuscus', 'Equus caballus', 'Pteropus alecto'],
                                'Metatheria': ['Macropus eugenii', 'Monodelphis domestica', 'Sarcophilus harrisii'],
                                'Primates': {
                                    'Atelidae': ['Ateles geoffroyi', 'Lagothrix lagotricha'],
                                    'Cebidae': ['Callithrix jacchus', 'Saguinus labiatus', 'Saimiri boliviensis'],
                                    'Cercopithecidae': ['Macaca mulatta', 'Macaca nemestrina', 'Papio hamadryas', 'Pygathrix bieti'],
                                    'Cheirogaleidae': ['Microcebus murinus'],
                                    'Daubentoniidae': ['Daubentonia madagascariensis'],
                                    'Galagidae': ['Otolemur garnettii'],
                                    'Hominidae': ['Gorilla gorilla', 'Homo sapiens', 'Pan paniscus', 'Pan troglodytes', 'Pongo pygmaeus', 'Symphalangus syndactylus'],
                                    'Hylobatidae': ['Nomascus leucogenys'],
                                    'Lemuridae': ['Lemur catta']
                                },
                                'Prototheria': ['Ornithorhynchus anatinus'],
                                'Rodentia': ['Cavia porcellus', 'Cricetulus griseus', 'Mus musculus', 'Rattus norvegicus'],
                                'Ruminantia': ['Bos taurus', 'Capra hircus', 'Ovis aries'],
                                'Scandentia': ['Tupaia chinensis'],
                                'Suina': ['Sus scrofa']
                            },
                            'Sauria': ['Alligator mississippiensis', 'Anolis carolinensis', 'Chrysemys picta', 'Ophiophagus hannah', 'Python bivittatus'],
                            'Teleostei': [
                                'Astatotilapia burtoni', 'Cyprinus carpio', 'Danio rerio', 'Electrophorus electricus',
                                'Fugu rubripes', 'Gadus morhua', 'Hippoglossus hippoglossus', 'Ictalurus punctatus',
                                'Metriaclima zebra', 'Neolamprologus brichardi', 'Oreochromis niloticus',
                                'Oryzias latipes', 'Paralichthys olivaceus', 'Pundamilia nyererei',
                                'Salmo salar', 'Tetraodon nigroviridis'
                            ]
                        }
                    },
                    'Echinodermata': ['Lytechinus variegatus', 'Patiria miniata', 'Strongylocentrotus purpuratus'],
                    'Hemichordata': ['Saccoglossus kowalevskii']
                },
                'Ecdysozoa': {
                    'Arthropoda': {
                        'Chelicerata': ['Ixodes scapularis', 'Parasteatoda tepidariorum', 'Rhipicephalus microplus', 'Tetranychus urticae'],
                        'Crustacea': ['Daphnia pulex', 'Marsupenaeus japonicus', 'Triops cancriformis'],
                        'Hexapoda': [
                            'Acyrthosiphon pisum', 'Aedes aegypti', 'Anopheles gambiae', 'Apis mellifera',
                            'Bactrocera dorsalis', 'Biston betularia', 'Bombyx mori', 'Culex quinquefasciatus',
                            'Dinoponera quadriceps', 'Drosophila ananassae', 'Drosophila erecta',
                            'Drosophila grimshawi', 'Drosophila melanogaster', 'Drosophila mojavensis',
                            'Drosophila persimilis', 'Drosophila pseudoobscura', 'Drosophila sechellia',
                            'Drosophila simulans', 'Drosophila virilis', 'Drosophila willistoni',
                            'Drosophila yakuba', 'Heliconius melpomene', 'Locusta migratoria', 'Manduca sexta',
                            'Nasonia giraulti', 'Nasonia longicornis', 'Nasonia vitripennis', 'Plutella xylostella',
                            'Polistes canadensis', 'Spodoptera frugiperda', 'Tribolium castaneum'
                        ],
                        'Mandibulata': ['Strigamia maritima']
                    },
                    'Nematoda': [
                        'Ascaris suum', 'Brugia malayi', 'Caenorhabditis brenneri', 'Caenorhabditis briggsae',
                        'Caenorhabditis elegans', 'Caenorhabditis remanei', 'Haemonchus contortus',
                        'Heligmosomoides polygyrus', 'Panagrellus redivivus', 'Pristionchus pacificus',
                        'Strongyloides ratti'
                    ]
                },
                'Lophotrochozoa': {
                    'Annelida': ['Capitella teleta'],
                    'Brachiopoda': ['Glottidia pyramidata', 'Terebratulina retusa'],
                    'Mollusca': ['Haliotis rufescens', 'Lottia gigantea', 'Melibe leonina'],
                    'Nemertea': ['Cerebratulus lacteus'],
                    'Platyhelminthes': [
                        'Echinococcus granulosus', 'Echinococcus multilocularis', 'Fasciola hepatica',
                        'Gyrodactylus salaris', 'Mesocestoides corti', 'Schistosoma japonicum',
                        'Schistosoma mansoni', 'Schmidtea mediterranea'
                    ]
                }
            }
        },
        'Cnidaria': ['Hydra magnipapillata', 'Nematostella vectensis'],
        'Porifera': ['Amphimedon queenslandica', 'Leucosolenia complicata', 'Sycon ciliatum'],
        'Mycetozoa': ['Dictyostelium discoideum'],
        'Viridiplantae': {
            'Chlorophyta': ['Chlamydomonas reinhardtii'],
            'Coniferophyta': ['Cunninghamia lanceolata', 'Picea abies', 'Pinus densata', 'Pinus taeda'],
            'Embryophyta': {
                'species': ['Physcomitrella patens', 'Selaginella moellendorffii'],
                'Magnoliophyta': {
                    'species': ['Amborella trichopoda'],
                    'eudicotyledons': {
                        'Amaranthaceae': ['Salicornia europaea'],
                        'Araliaceae': ['Panax ginseng'],
                        'Asteraceae': [
                            'Cynara cardunculus', 'Helianthus annuus', 'Helianthus argophyllus',
                            'Helianthus ciliaris', 'Helianthus exilis', 'Helianthus paradoxus',
                            'Helianthus petiolaris', 'Helianthus tuberosus'
                        ],
                        'Brassicaceae': [
                            'Arabidopsis lyrata', 'Arabidopsis thaliana', 'Brassica napus',
                            'Brassica oleracea', 'Brassica rapa', 'Camelina sativa'
                        ],
                        'Caricaceae': ['Carica papaya'],
                        'Cucurbitaceae': ['Cucumis melo', 'Cucumis sativus'],
                        'Euphorbiaceae': ['Hevea brasiliensis', 'Manihot esculenta', 'Ricinus communis'],
                        'Fabaceae': [
                            'Acacia auriculiformis', 'Acacia mangium', 'Arachis hypogaea', 'Glycine max',
                            'Glycine soja', 'Lotus japonicus', 'Medicago truncatula', 'Phaseolus vulgaris',
                            'Vigna unguiculata'
                        ],
                        'Lamiales': [
                            'Avicennia marina', 'Digitalis purpurea', 'Rehmannia glutinosa',
                            'Salvia miltiorrhiza', 'Salvia sclarea'
                        ],
                        'Linaceae': ['Linum usitatissimum'],
                        'Malvaceae': [
                            'Gossypium arboreum', 'Gossypium herbaceum', 'Gossypium hirsutum',
                            'Gossypium raimondii', 'Theobroma cacao'
                        ],
                        'Myrtaceae': ['Eugenia uniflora'],
                        'Paeoniaceae': ['Paeonia lactiflora'],
                        'Ranunculaceae': ['Aquilegia caerulea'],
                        'Rhizophoraceae': ['Bruguiera cylindrica', 'Bruguiera gymnorhiza'],
                        'Rosaceae': ['Fragaria vesca', 'Malus domestica', 'Prunus persica'],
                        'Rutaceae': [
                            'Citrus clementina', 'Citrus reticulata', 'Citrus sinensis',
                            'Citrus trifoliata'
                        ],
                        'Salicaceae': ['Populus euphratica', 'Populus trichocarpa'],
                        'Solanaceae': ['Nicotiana tabacum', 'Solanum lycopersicum', 'Solanum tuberosum'],
                        'Vitaceae': ['Vitis vinifera']
                    },
                    'monocotyledons': [
                        'Aegilops tauschii', 'Asparagus officinalis', 'Brachypodium distachyon',
                        'Elaeis guineensis', 'Festuca arundinacea', 'Hordeum vulgare', 'Oryza sativa',
                        'Saccharum officinarum', 'Saccharum sp.', 'Sorghum bicolor',
                        'Triticum aestivum', 'Triticum turgidum', 'Vriesea carinata', 'Zea mays'
                    ]
                }
            }
        },
        'Viruses': [
            'Bandicoot papillomatosis carcinomatosis virus type 1',
            'Bandicoot papillomatosis carcinomatosis virus type 2',
            'BK polyomavirus', 'Bovine foamy virus', 'Bovine herpesvirus 1',
            'Bovine herpesvirus 5', 'Bovine leukemia virus', 'Duck enteritis virus',
            'Epstein Barr virus', 'Gorilla gorilla gorilla polyomavirus 1',
            'Herpes B virus', 'Herpes Simplex Virus 1', 'Herpes Simplex Virus 2',
            'Herpesvirus of turkeys', 'Herpesvirus saimiri strain A11',
            'Human cytomegalovirus', 'Human herpesvirus 6B', 'Human immunodeficiency virus 1',
            'Infectious laryngotracheitis virus', 'JC polyomavirus',
            'Kaposi sarcoma-associated herpesvirus', "Marek's disease virus type 1",
            "Marek's disease virus type 2", 'Merkel cell polyomavirus',
            'Mouse cytomegalovirus', 'Mouse gammaherpesvirus 68',
            'Pan troglodytes verus polyomavirus 2a', 'Pseudorabies virus',
            'Raccoon polyomavirus', 'Rhesus lymphocryptovirus', 'Rhesus monkey rhadinovirus',
            'Simian foamy virus', 'Simian virus 40', 'Torque teno virus'
        ]
    }
    
    return ref

def get_summary_and_sequence(detail_url, session, headers):
    """
    Fetches the summary, sequence, and structure from the miRNA detail page.
    Empty summary is considered valid.
    
    Parameters:
        detail_url (str): The full URL to the miRNA detail page
        session (requests.Session): The session object for persistent connections
        headers (dict): The headers to use for the request
    
    Returns:
        tuple: (summary_text, sequence_text, structure_text)
    """
    try:
        logging.info(f"Fetching detail page: {detail_url}")
        response = session.get(detail_url, headers=headers, timeout=30)
        response.raise_for_status()
        soup = BeautifulSoup(response.text, 'html.parser')
        
        # Initialize variables
        summary_text = ""
        sequence_text = ""
        structure_text = ""
        
        # Try to find the description
        description_div = soup.find('div', style=lambda value: value and 'margin-left:15px' in value)
        if description_div:
            pre_tag = description_div.find('pre')
            if pre_tag:
                summary_text = pre_tag.get_text(separator=' ', strip=True)
                logging.info(f"Summary found: {summary_text[:50]}...")
        
        # Try to find the sequence and structure
        sequence_div = soup.find('div', {'id': 'hairpinSequence'})
        if sequence_div:
            # Find all text-monospace spans
            monospace_spans = sequence_div.find_all('span', {'class': 'text-monospace'})
            
            # First span contains the sequence
            if len(monospace_spans) > 0:
                sequence_text = monospace_spans[0].get_text(strip=True)
                logging.info(f"Sequence found: {sequence_text[:50]}...")
            
            # Second span contains the structure
            if len(monospace_spans) > 1:
                structure_text = monospace_spans[1].get_text(strip=True)
                logging.info(f"Structure found: {structure_text[:50]}...")
        
        return summary_text, sequence_text, structure_text
            
    except requests.exceptions.Timeout:
        logging.error(f"Timeout while fetching {detail_url}")
        return "", "", ""
        
    except requests.exceptions.HTTPError as http_err:
        logging.error(f"HTTP error while fetching {detail_url}: {http_err}")
        if http_err.response.status_code == 429:  # Too Many Requests
            logging.info("Rate limited, waiting 60 seconds...")
            time.sleep(60)
            return get_summary_and_sequence(detail_url, session, headers)
        return "", "", ""
        
    except Exception as err:
        logging.error(f"Error while fetching {detail_url}: {err}")
        return "", "", ""

# output_csv='miRNA_human_with_sequence.csv'
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
        # this for Arab
        # table = soup.find('table', {'id': 'ath-results-table'})
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
                summary, sequence, structure = get_summary_and_sequence(detail_url, session, headers)
                time.sleep(1)  # Polite delay between requests
            
            seq_len = len(sequence)
            sum_len = len(summary)
            
            data_list.append({
                'Name': name,
                'Start': start,
                'End': end,
                'Sequence': sequence,
                'Seq_structure': structure,
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
            fieldnames = ['Name', 'Start', 'End', 'Sequence', 'Seq_structure', 'Seq_len', 'Summary', 'Summary_len', 'Chr', 'Strand']
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
    # for human
    mirbase_url = "https://mirbase.org/browse/results/?organism=hsa"
    
    # for mouse
    # mirbase_url = "https://mirbase.org/browse/results/?organism=mmu"
    
    # for rat
    # mirbase_url = "https://mirbase.org/browse/results/?organism=rno"
    
    # for fly
    # mirbase_url = "https://mirbase.org/browse/results/?organism=dme"
    
    # for worm
    # mirbase_url = "https://mirbase.org/browse/results/?organism=cel"
    
    # for Arabidopsis
    # mirbase_url = "https://mirbase.org/browse/results/?organism=ath"
    
    scrape_mirbase_with_sequence(mirbase_url)       