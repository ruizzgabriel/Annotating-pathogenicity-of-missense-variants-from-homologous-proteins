# Some API searches returned an error. Here, we retry them.

import pandas as pd
import glob
import time
import random
import requests
import os

# Set constant for the API URL
GNOMAD_API_URL = "https://gnomad.broadinstitute.org/api"


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
    return fetch(req_variantlist)


def retry_errors(errors_file, gnomad_v, ref_genome):
    error_list = []
    # Set a variable used to introduce pauses to avoid 'Too many requests' errors
    time_count = 0
    for i, error in errors_file.iterrows():
        # Increment the value of the count variable by 1
        time_count += 1

        # Set the ENSG code
        ensg = error[0]

        print(f"{ensg}_________count__{time_count}")

        # Set the gene name
        gene_name = error[1].split()[0]

        # Wait 4 seconds after each search, to avoid 'Too many requests' errors
        time.sleep(4)

        # Check if the count is multiple of 10
        if time_count % 10 == 0:
            # If yes, generate a random integer between 13 and 53 and pause the script for that number of seconds, to avoid 'Too many requests' errors
            time.sleep(random.randint(30, 40))

        # Check that ensg variable is not a NaN
        if pd.isna(ensg):
            # If yes, skip this row
            continue
        try:
            # Check not to search version gnomad_r3 for genome version CRCh37 because is not available
            if gnomad_v == "gnomad_r3" and ref_genome == "GRCh37":
                continue

            # Set the name of the output file
            variants_file = f"data/{gnomad_v}_{ref_genome}_variants_rescued.txt"

            # Call the get_variant_list function with the ENSG code
            variant_list = get_variant_list(ensg, gnomad_v, ref_genome)

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
            error_list.append((ensg, gene_name, e))
            print(f"Error: {e}")
    return error_list


# errors = pd.read_csv("data/gnomad_r2_1_GRCh37_error.txt", sep="\t", header=None)

# Define the folder path
folder_path = "./data"

# Define the file pattern
file_pattern = "gnomad_r*_error.txt"

# Save the file list matching the folder path and file pattern
file_list_pattern = f"{folder_path}/{file_pattern}"
file_list_pattern = os.path.normpath(file_list_pattern)
file_list = glob.glob(file_list_pattern)

# For every file in the file list do:
for file in file_list:
    print(f"------------------WE START WITH THE FILE: {file} ------------------")
    # Open the file of the corresponding gnomAD and genome versions
    errors_file = pd.read_csv(file, sep="\t", header=None)
    # Check which is the gnomAD and genome version of the current file
    if "gnomad_r3" in file:
        gnomad_v = "gnomad_r3"
    elif "gnomad_r2_1" in file:
        gnomad_v = "gnomad_r2_1"
    elif "gnomad_r4" in file:
        gnomad_v = "gnomad_r4"
    if "GRCh37" in file:
        ref_genome = "GRCh37"
    elif "GRCh38" in file:
        ref_genome = "GRCh38"

    # Retry the search for the entries that we could not find results in the first execution
    error_list = retry_errors(errors_file, gnomad_v, ref_genome)

    # Set the name of the new error file
    error_file = f"data/{gnomad_v}_{ref_genome}_error_notFoundGenes.txt"
    error_file = os.path.normpath(error_file)
    # Write the error list to the error file with both ENSG code and gene name from those that are not found in gnomAD
    with open(error_file, "w") as f:
        for ensg, gene_name, error in error_list:
            f.write(f"{ensg}\t{gene_name}\t{error}\n")
