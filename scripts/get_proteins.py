# get amino acid sequence and description, just like get_genes

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

