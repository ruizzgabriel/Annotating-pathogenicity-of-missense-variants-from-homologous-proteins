### This script collects the PFAM codes corresponding to the non-pathogenic genes

import pandas as pd
import requests
import time

non_patho = pd.read_csv('./Domains/non_patho_unique_genes.txt', sep='\t')
print(non_patho.head())


def get_pfam_from_uniprot(gene):
    try:
        # UniProt API endpoint
        url = f"https://rest.uniprot.org/uniprotkb/search?query=gene:{gene}+organism_id:9606&fields=accession,gene_names,xref_pfam"
        response = requests.get(url)
        
        if response.status_code == 200:
            data = response.json()
            if data['results']:
                # Get the first reviewed (Swiss-Prot) result if available
                for result in data['results']:
                    if result['entryType'] == "UniProtKB reviewed (Swiss-Prot)":
                        # Extract Pfam IDs from cross-references
                        pfam_ids = []
                        if 'uniProtKBCrossReferences' in result:
                            for xref in result['uniProtKBCrossReferences']:
                                if xref['database'] == 'Pfam':
                                    pfam_id = xref['id']
                                    pfam_ids.append(pfam_id)
                        
                        # Return semicolon-separated Pfam IDs or None if no Pfam domains
                        return ';'.join(pfam_ids) if pfam_ids else None
                # If no reviewed entry found, try the first result
                result = data['results'][0]
                pfam_ids = []
                if 'uniProtKBCrossReferences' in result:
                    for xref in result['uniProtKBCrossReferences']:
                        if xref['database'] == 'Pfam':
                            pfam_id = xref['id']
                            pfam_ids.append(pfam_id)
                return ';'.join(pfam_ids) if pfam_ids else None
        return None
    
    except Exception as e:
        print(f"Error retrieving data for {gene}: {str(e)}")
        return None

# Add Pfam column by applying the function to each gene
non_patho['PFAM'] = non_patho['Gene(s)'].apply(lambda x: get_pfam_from_uniprot(x))

# Add small delay between requests to avoid overwhelming the API
time.sleep(1)

# Display the resulting DataFrame
print(non_patho)

# Optionally save to file
non_patho.to_csv('./Domains/genes_with_pfam.txt', sep='\t', index=False)