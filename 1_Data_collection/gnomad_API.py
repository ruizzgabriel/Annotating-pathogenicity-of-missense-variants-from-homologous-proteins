# Retrieve data from Gnomad database about all human reviewed proteins in UniProt

import requests
import os
import pathlib
import pandas as pd
import time
import random

# get the directory path of the current file
dir_path = pathlib.Path(__file__).parent.absolute()

# set the working directory to the directory path of the current file
os.chdir(dir_path)

# Set constants for the API URL, dataset, and reference genome.
GNOMAD_API_URL = "https://gnomad.broadinstitute.org/api"
DATASET_LIST = ["gnomad_r4"]    # "gnomad_r3", "gnomad_r2_1"]
REFERENCE_GENOME = ["GRCh38"]  # , "GRCh37"]


def fetch(jsondata, url=GNOMAD_API_URL):
    """
    Make a POST request to the GNOMAD API and return the response as JSON.

    Args:
        jsondata (dict): The JSON payload to send in the request.
        url (str, optional): The API endpoint URL. Defaults to GNOMAD_API_URL.

    Raises:
        Exception: If there are errors in the response JSON.

    Returns:
        dict: The JSON response from the API.
    """
    headers = {"Content-Type": "application/json"}
    response = requests.post(url, json=jsondata, headers=headers)
    response.raise_for_status()  # raise an exception if the request fails
    json = response.json()
    if "errors" in json:
        raise Exception(str(json["errors"]))
    return json


def get_variant_list(gene_id, genomad_dataset, ref_genome):
    """
    Query the GNOMAD API for a list of variants for the given gene.

    Args:
        gene_id (str): The Ensembl gene ID to query.

    Returns:
        list: A list of variant objects, each containing consequence, ref, alt, pos,
              transcript_id, hgvsc, hgvsp, rsid, and variant_id fields.
    """
    # Define the GraphQL query string using f-strings to interpolate the gene ID, dataset, and reference genome.
    query = f"""
        {{
            gene(gene_id: "{gene_id}", reference_genome: {ref_genome}) {{
                variants(dataset: {genomad_dataset}) {{
                    chrom
                    pos
                    rsid
                    ref
                    alt
                    transcript_id
                    hgvsc
                    hgvsp
                    consequence
                    gene_symbol
                    variant_id: variantId
                    genome {{
                        af
                        homozygote_count
                        hemizygote_count
                    }}
                    exome {{
                        af
                        homozygote_count
                        hemizygote_count
                    }}  
                    }}
                clinvar_variants {{
                    clinical_significance
                    clinvar_variation_id
                    }}       
            }}
        }}
    """
    # Create the JSON payload to send in the request.
    req_variantlist = {"query": query, "variables": {"withFriends": False}}

    # ----
    # Open a file for writing --> Only to check the output format, it can be removed
    # response = fetch(req_variantlist)
    # with open("output.txt", "w") as f:
    # Redirect print output to the file
    #    print(response, file=f)
    # print(type(response))
    # ----

    # Send the request and return the list of variants from the response.
    return fetch(req_variantlist)


# Load the pandas DataFrame with the Ensembl codes and gene names
uniprot = pd.read_csv("data/uniprot_db_human_2023_04_06_with_NM.tsv", sep=",")
#print(uniprot)
# Create a list to store the ENSG codes that give errors
error_list = []

# Create a counter to control how many searches have been done
count = 0
# Loop through each row in the DataFrame
for i, row in uniprot.iterrows():
    # Wait 3 seconds after each search, to avoid 'Too many requests' errors
    time.sleep(3)
    # Increment the value of the count variable by 1
    count += 1
    # Retrieve the ENSG code from "Bgee" column
    ensg = row["Bgee"]

    # Print counter and ENSG code
    #print(count)
    #print(ensg)

    # Check if the count is multiple of 10
    if count % 10 == 0:
        # If yes, generate a random integer between 13 and 53 and pause the script for that number of seconds, to avoid 'Too many requests' errors
        time.sleep(random.randint(13, 53))
    
    if count % 100 == 0:
        print(f"{count} processed out of {len(uniprot)} UniProt entries.")

    # Retrieve the gene name from "Gene Names" column
    gene_name = row["Gene Names"]
    
    # Check that ensg variable is not a NaN
    if pd.isna(ensg):
        # If yes, skip this row
        continue
    try:
        # For every gnomAD version in DATASET_LIST do:
        for dataset in DATASET_LIST:
            # For every genome version in REFERENCE_GENOME do:
            for genome in REFERENCE_GENOME:
                # Check not to search version gnomad_r3 for genome version CRCh37 because is not available
                if dataset == "gnomad_r3" and genome == "GRCh37":
                    continue
                # Set the names of the output files
                variants_file = f"data/{dataset}_{genome}_variants.txt"
                error_file = f"data/{dataset}_{genome}_error.txt"

                # Call the get_variant_list function with the ENSG code
                variant_list = get_variant_list(ensg, dataset, genome)

                # get the list of variants
                variants = variant_list["data"]["gene"]["variants"]

                # get the list of clinical significances
                clinvar_variants = variant_list["data"]["gene"]["clinvar_variants"]

                # create a dataframe for variants
                variants_df = pd.DataFrame(variants)

                # create a dataframe for clinical significances
                clinvar_variants_df = pd.DataFrame(clinvar_variants)

                # concatenate the two dataframes horizontally (variants and clinvar)
                result_df = pd.concat([variants_df, clinvar_variants_df], axis=1)

                # Open the variants file for appending
                with open(variants_file, "a") as f:
                    # Write the result dataframe to the file
                    result_df.to_csv(f, header=f.tell() == 0, index=False)

    except Exception as e:
        # if an error occurs, add the gene name to the error list
        error_list.append((ensg, gene_name))
        print(f"Error: {e}")

# Write the error list to the error file with both ENSG code and gene name
with open(error_file, "w") as f:
    for ensg, gene_name in error_list:
        f.write(f"{ensg}\t{gene_name}\n")
