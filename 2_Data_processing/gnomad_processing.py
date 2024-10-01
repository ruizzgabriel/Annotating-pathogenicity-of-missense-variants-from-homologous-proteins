# Process Gnomad file, obtained in the previous step via API. 

import glob
import re
import csv
import pandas as pd
import numpy as np
import os
import time
from decimal import Decimal


def write_df(output_df, new_row):
    """
    Write an entry into a DataFrame

    Args:
        output_df(DataFrame): Where the information have to be introduced
        new_row(list): A new entry in the DataFrame

    Returns:
        output_df(DataFrame): The same as before but with the new entry
    """
    output_df = output_df.append(
        pd.Series(new_row, index=output_df.columns), ignore_index=True
    )
    return output_df


def biomart_dic():
    """
    Extract ENST and NM information from biomart file and create a dictionary relating them.

    Returns:
        biomart_d(Dictionary): relating each ENST(key) with it/its NM(value/s)
    """
    # Set the path to the biomart file
    biomart_file = os.path.normpath("../1_Data_collection/data/biomart_human_genes_2023_05_08.txt")
    # Read data extracted from biomart database using pandas package
    biomart = pd.read_csv(biomart_file, delimiter="\t", header="infer")
    # Filter the DataFrame for the columns of interest and drop duplicates
    filtered_biomart = biomart[
        ["Transcript stable ID", "RefSeq mRNA ID"]
    ].drop_duplicates(keep="first")

    # Initialize the dictionary
    biomart_d = {}
    # Iterate over the biomart DataFrame
    for index, value in filtered_biomart.iterrows():
        # Save ENST code
        ENST = value["Transcript stable ID"]
        # Save NM code
        NM = value["RefSeq mRNA ID"]
        # Check if both ENST and NM codes are not NAs
        if type(ENST) != float and type(NM) != float:
            # Check if the key element NM is already in the dictionary
            if ENST in biomart_d:
                # If yes, add the ENST code to a list including the different ENST codes for the NM
                biomart_d[ENST].append(NM)
            else:
                # If not, create the entry in the dictionary
                biomart_d[ENST] = [NM]

    return biomart_d


def entry_processing(entry, biomart_d):
    """
    Process the entry and divide the elements that we are interested in into different columns.

    Args:
        entry(list): containing the whole entry in the dataset

    Returns:
        new_row(list): new row of the output file
    """
    # Obtain the different informations from raw data
    rsid = entry[2]
    enst = entry[5]
    # Obtain NM code from ENST code using the biomart dictionary
    NM = biomart_d[enst] if enst in biomart_d else np.nan
    prot_change = entry[7]
    consequence = entry[8]
    gene = entry[9]
    # Obtain clinical significance if is present
    clinical_significance = entry[13] if len(entry) > 13 else np.nan
    # Extract initial aminoacid from the protein change
    initial_aa = prot_change[2:5]
    # Extract the position from the protein change
    position = re.findall(r"\d+", prot_change[3:10])
    position = position[0] if position else None
    # Extract final aminoacid from the protein change
    final_aa = match[1] if (match := re.search(r"\d+(.*)", prot_change)) else np.nan
    # Extract the genome allele frequence
    genome_af = (
        entry[11].split(",")[0].split(":")[1] if len(entry) > 11 and entry[11] else None
    )
    # If the genome allele frequence is not None, we transform it to a Decimal object
    if genome_af is not None:
        genome_af = Decimal(genome_af)

    # Extract the exome allele frequence
    exome_af = (
        entry[12].split(",")[0].split(":")[1] if len(entry) > 11 and entry[12] else None
    )

    # If the exome allele frequence is not None, we transform it to a Decimal object
    if exome_af is not None:
        exome_af = Decimal(exome_af)

    return [
        gene,
        enst,
        NM,
        rsid,
        prot_change,
        position,
        initial_aa,
        final_aa,
        consequence,
        clinical_significance,
        genome_af,
        exome_af,
    ]


def data_extraction(entries_gnomad, gnomad_v, genome):
    """
    Data extraction and processment from the raw data files of gnomAD

    Args:
        entries_gnomad(csv.reader): Contains the information of a csv file
        gnomad_v(string): Indicates the gnomAD version
        genome(string): Indicates the genome version

    Returns:
        error_list(list): Containing the non-processed entries and the error obtained
    """
    output_df = pd.DataFrame(
        columns=[
            "gene",
            "enst",
            "NM",
            "rsid",
            "prot_change",
            "position",
            "initial_aa",
            "final_aa",
            "consequence",
            "clinical_significance",
            "genome_af",
            "exome_af",
        ]
    )
    # Create dictionary relating ENST codes with NM codes
    biomart_d = biomart_dic()

    header = True
    error_list = []
    count = 0
    s_time = time.time()
    output_file_path = os.path.normpath("./data")
    for entry in entries_gnomad:
        if entry:
            count += 1
            try:
                new_row = entry_processing(entry, biomart_d)
                output_df = write_df(output_df, new_row)
                # print(output_df)
                if count % 1000 == 0:
                    with open(
                        f"{output_file_path}/{gnomad_v}_{genome}_processed.csv",
                        "a",
                        newline="",
                    ) as gnomad_file:
                        # Condition: if is the first time we write in the file, we add the header
                        if header == True:
                            output_df.to_csv(gnomad_file, header=True, index=False)
                            header = False
                        # The rest of times we do not include another header each time we write
                        else:
                            output_df.to_csv(gnomad_file, header=False, index=False)

                    # Reinitialize the dataframe with the columns
                    output_df = pd.DataFrame(
                        columns=[
                            "gene",
                            "enst",
                            "NM",
                            "rsid",
                            "prot_change",
                            "position",
                            "initial_aa",
                            "final_aa",
                            "consequence",
                            "clinical_significance",
                            "genome_af",
                            "exome_af",
                        ]
                    )

                if count % 1000000 == 0:
                    f_time = time.time()
                    minutes = (f_time - s_time) / 60
                    print(
                        f"The first {count} entries of {gnomad_v}_{genome} has been processed. The last milion entries processing taked {minutes:.2f} minutes."
                    )
                    s_time = time.time()

            except Exception as e:
                # if an error occurs, add a tuple containing the entry to the error list and the error type
                error_list.append((entry, e))
                print(f"Error: {e}")

    return error_list


# Set the path to data directory
folder_path = "../1_Data_collection/data"
folder_path = os.path.normpath(folder_path)

# Set the files pattern of gnomAD
#file_pattern = "gnomad_r*_variants.txt"
file_pattern = "gnomad_r*_variants_rescued.txt"

# List the files matching the pattern
file_list = glob.glob(f"{folder_path}/{file_pattern}")
print(file_list)

for file in file_list:
    with open(file, "r") as csvfile:
        entries_gnomad = csv.reader(csvfile, delimiter=",")
        next(entries_gnomad)

        file_name = os.path.basename(file)
        if "r2" in file_name:
            gnomad_v = "gnomad_r2_1"
            genome = file_name.split("_")[3].split(".")[0]
        elif "r3" in file_name:
            gnomad_v = "gnomad_r3"
            genome = file_name.split("_")[2].split(".")[0]
        
        elif "r4" in file_name:
            gnomad_v = "gnomad_r4"
            genome = file_name.split("_")[2].split(".")[0]


        print(f"-----We start to process {gnomad_v}_{genome}-----")
        start_time = time.time()
        error_list = data_extraction(entries_gnomad, gnomad_v, genome)
        error_file_path = os.path.normpath("./data")
        error_file = (
            f"{error_file_path}/{gnomad_v}_{genome}_errors_in_processing_step.txt"
        )
        with open(error_file, "w") as f:
            for entry, error_type in error_list:
                f.write(f"/{entry}\t{error_type}\n")
        end_time = time.time()
        # calculate the elapsed time
        elapsed_time = (end_time - start_time) / 3600
        print(f"Elapsed time {gnomad_v}_{genome}: {elapsed_time:.2f} hours")