import pandas as pd
import subprocess
import os

### Here we're filtering a file called AlphaMissense_aa_substitutions.tsv that has been retrieved from:
# https://alphamissense.hegelab.org/
# https://zenodo.org/records/8208688
# https://zenodo.org/records/8208688/files/AlphaMissense_aa_substitutions.tsv.gz?download=1




# Load the ENST codes from the file
accession = pd.read_csv('./Alphamissense/Unique_accession_from_pairs.txt')

# Define the output directory and create it if it doesn't exist
output_dir = "./Alphamissense/Alphamissense_by_accessions/"
os.makedirs(output_dir, exist_ok=True)

# Loop through each ENST code in the 'Entry' column
for code in accession['Entry']:
    # Skip any header-like entry if necessary
    if code != 'Entry':
        # Define the output filename based on the ENST code
        output_filename = f"{output_dir}{code}.txt"
        
        # Run the grep command to extract lines containing the ENST code
        with open(output_filename, 'w') as output_file:
            subprocess.run(
                ["grep", code, "./Alphamissense/AlphaMissense_aa_substitutions.tsv"],
                stdout=output_file
            )

        print(f"Extracted data for {code} into {output_filename}")
