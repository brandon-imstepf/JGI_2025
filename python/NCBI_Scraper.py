import os
import requests
from bs4 import BeautifulSoup

base_url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/"
output_dir = "NCBI_Datasets"
os.makedirs(output_dir, exist_ok=True)

DEBUG = True  # Set to True for debug output

def safe_get(url, timeout=120):
    try:
        return requests.get(url, timeout=timeout)
    except Exception as e:
        print(f"Error fetching {url}: {e}")
        return None

def get_taxid_from_assembly_report(report_url):
    response = safe_get(report_url)
    if not response or response.status_code != 200:
        if DEBUG:
            print(f"Failed to fetch assembly report: {report_url}")
        return None
    for line in response.text.splitlines():
        if line.startswith('# Taxid:') or line.startswith('# TaxID:'):
            parts = line.split(':', 1)
            if len(parts) > 1:
                taxid = parts[1].strip()
                if DEBUG:
                    print(f"Found TaxID: {taxid} in {report_url}")
                return taxid
    if DEBUG:
        print(f"No TaxID found in {report_url}")
    return None

def is_prokaryote(taxid):
    tax_url = f"https://taxonomy.jgi.doe.gov/sc/id/{taxid}"
    response = safe_get(tax_url)
    if not response or response.status_code != 200:
        if DEBUG:
            print(f"Failed to fetch taxonomy for TaxID {taxid}")
        return False
    tax_text = response.text.lower()
    if DEBUG:
        print(f"Taxonomy for {taxid}: {tax_text}")
    if any(x in tax_text for x in ["d:bacteria", "d:archaea", "prokaryote", "prokaryota"]):
        return True
    else:
        print(f"[DEBUG] Taxonomy for TaxID {taxid} did not match prokaryote: {tax_text}")
        return False

def download_files_from_dir(dir_url, rel_path, depth, max_depth=4, max_dirs=5):
    response = safe_get(dir_url)
    if not response:
        return
    soup = BeautifulSoup(response.text, 'html.parser')
    subdirs = []
    files = []
    for link in soup.find_all('a'):
        href = link.get('href')
        if href and href.endswith('/') and href != '../':
            subdirs.append(href)
        elif href:
            files.append(href)
    if depth < max_depth and subdirs:
        for i, subdir in enumerate(subdirs[:max_dirs] if depth == 1 else subdirs):
            if DEBUG:
                print(f"{'  '*depth}[{i+1}/{len(subdirs)}] Descending into {dir_url+subdir}")
            download_files_from_dir(dir_url + subdir, os.path.join(rel_path, subdir.strip('/')), depth+1, max_depth, max_dirs)
    else:
        # Find required files
        fna_file = next((f for f in files if f.endswith("genomic.fna.gz") and not f.endswith("from_genomic.fna.gz")), None)
        gff_file = next((f for f in files if f.endswith("genomic.gff.gz")), None)
        report_file = next((f for f in files if f.endswith("_assembly_report.txt")), None)
        if fna_file and gff_file and report_file:
            report_url = dir_url + report_file
            taxid = get_taxid_from_assembly_report(report_url)
            if taxid and is_prokaryote(taxid):
                for file in [fna_file, gff_file]:
                    file_url = dir_url + file
                    out_path = os.path.join(output_dir, file)
                    if os.path.exists(out_path):
                        if DEBUG:
                            print(f"File {out_path} already exists, skipping.")
                        continue
                    print(f"{'  '*depth}Downloading {file_url} -> {out_path}")
                    file_response = safe_get(file_url)
                    if file_response:
                        with open(out_path, 'wb') as f:
                            f.write(file_response.content)
            else:
                if DEBUG:
                    print(f"Skipping {dir_url} due to non-prokaryote taxonomy.")
        else:
            if DEBUG:
                print(f"Skipping {dir_url}: Required files not found.")

# Start recursive traversal from base_url, limit top-level dirs for testing
MAX_TOP_DIRS = 999999
response = safe_get(base_url)
if not response:
    print("Failed to fetch base URL.")
    exit(1)
soup = BeautifulSoup(response.text, 'html.parser')
gcf_dirs = [link.get('href') for link in soup.find_all('a') if link.get('href') and link.get('href').endswith('/')]
for i, dir in enumerate(gcf_dirs[:MAX_TOP_DIRS]):
    print(f"[{i+1}/{min(MAX_TOP_DIRS, len(gcf_dirs))}] Processing {base_url+dir}")
    download_files_from_dir(base_url + dir, dir.strip('/'), depth=1, max_depth=4, max_dirs=5)
