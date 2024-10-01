# Get equivalent positions for Gnomad variants using Pfam alignments
# This code was adapted from Sergi Soldevila in order to consider all Pfam codes, not only the first.

import pandas as pd
import os
import glob
import numpy as np
import contextlib
import time
import argparse
from colorama import Fore, Style
from Bio.Data.IUPACData import protein_letters_3to1


""" # Create the parser
parser = argparse.ArgumentParser()
# Add an argument
parser.add_argument('--version_pattern', type=str, help='Version pattern')
#parser.add_argument('--up_low_case', type=str, help='Keep pfam vs upperCase')
# Parse the arguments
args = parser.parse_args()
"""

full_vs_seed = input("Write the desired alignment type. (Mix condition prioritizes seed alignment and uses full ones when seed is not found)\nYou can choose between:\n\tfull\n\tseed\n\tseed+full(mix)\nPlease, write the version as indicated: ")


up_low_case = 'keep_pfam'

def pos_align(gnomad_df, full_vs_seed):
    """
    Add different information to the gnomAD DataFrame:
        - Alignment Pfam position from the protein change position
        - In which Pfam alignment was found this position
        - Check if the initial aminoacid in the protein change matches with the aa in the alignment position
    
    Args:
        gnomad_df (DataFrame): Containing the filtered gnomAD file
    
    Returns:
        gnomad_df (DataFrame): The input DataFrame with the added information.
    """
    # Set the path where Pfam alignments are located
    if full_vs_seed == 'seed':
        pfam_path = "../1_Data_collection/data/Pfam_SEED_HUMAN_2024_03_20"
    elif full_vs_seed == 'full':
        pfam_path = "../1_Data_collection/data/Pfam_HUMAN_InterPro_2023-07-03"

    print(f'We are using Pfams with the path... {pfam_path}')
    pfam_path = os.path.normpath(pfam_path)
    # Set a boolean to indicate if the position is already found in the row to avoid checking more Pfam alignments
    switch = False
    # Set a counter to check the cases for which the alignment was not found
    cases = 0
    processed_rows = 0
    start_time = time.time()
    # Iterate over the rows of gnomAD DataFrame
    for index, row in gnomad_df.iterrows():
        processed_rows += 1
        if processed_rows % 100000 == 0:
            elapsed_time = time.time() - start_time
            print(f"Processed {processed_rows} rows in {round(elapsed_time / 60, 2)} mins.")

        # Save the list of Pfam codes in the current row
        pfam_list = row["Pfam"]
        # Iterate over the list of Pfam codes
        for pfam_code in pfam_list:
            # If the position has not been retrieved for this row yet
            if switch == False:
                # Handle the possibility that the corresponding Pfam alignment file is not found
                with contextlib.suppress(FileNotFoundError):
                    # Set the name of the file from the Pfam code
                    pfam_file = pfam_path + "/" + str(pfam_code) + "_HUMAN.txt"
                    pfam_file = os.path.normpath(pfam_file)
                    # print(pfam_file)

                    # Open the corresponding Pfam alignment file
                    with open(pfam_file) as f:
                        # Save the content of the file
                        file_content = f.readlines()
                    # Iterate over the contents of the file
                    for text in file_content:
                        # Remove the newline characters
                        text = text.replace("\n", "")
                        # Save the UniProt entry code and transform it to string
                        uniprot = row["Uniprot_entry_name"]
                        uniprot = str(uniprot)
                        # Check if the text starts by the UniProt code
                        if text.startswith(uniprot):
                            # If yes, save the rang, initial and final position variables from the text
                            rang = text.split("/")[1].split("-", maxsplit=1)
                            pos_i = float(rang[0])
                            pos_f = float(rang[1].split()[0])
                            # Save the number indicating the position of the aminoacidic change and transform it to float
                            aa_change_pos = row["position"]
                            #aa_change_pos = aa_change_pos.replace('[', '').replace(']', '').replace("'", "").strip()
                            aa_change_pos = float(aa_change_pos)
                            # Check if the aminoacidic change position fits inside the range of the alignment position
                            if aa_change_pos >= pos_i and aa_change_pos <= pos_f:
                                # If yes, initialize two counters to 0.
                                compt_aa = 0
                                compt = 0
                                # Save the part of the text corresponding to the alignment
                                alignment = text.split()[1]

                                # Iterate over each position of the alignment
                                for align_pos in alignment:
                                    # Add one to the counter each time
                                    compt = compt + 1
                                    # If the position is an aminoacid, add one to the other counter
                                    if align_pos.isalpha() == True:
                                        compt_aa = compt_aa + 1
                                        # Check if the counter of aminoacids is equal to this expression: position of the aminoacid change - initial position in the pfam alignment + 1
                                        if compt_aa == (aa_change_pos - pos_i + 1):
                                            # If yes, save the Pfam code where the alignment was found
                                            gnomad_df.at[index, "Pfam_found"] = pfam_code
                                            # Add the counter of positions in the alignment to the Pos_align column in the corresponding row
                                            gnomad_df.at[index, "Pos_align"] = compt

                                            # Extract the initial aa and transform to one-letter code using protein_letters_3to1 dictionary from Bio
                                            initial_aa = row["initial_aa"]
                                            if initial_aa in protein_letters_3to1:
                                                initial_aa = protein_letters_3to1[initial_aa]
                                                # Check if the aminoacid in the alignment is the same as the initial one, compute it depending on the chosen case (considering only uppercase or all residues)
                                                if (up_low_case == 'upperCase' and align_pos.upper() == initial_aa
                                                    or up_low_case == 'keep_pfam' and align_pos == initial_aa):
                                                    # If yes, introduce True to Initial_aa_coincidence column
                                                    gnomad_df.at[index, "Initial_aa_coincidence"] = True
                                                else:
                                                    # If not, introduce False
                                                    gnomad_df.at[index, "Initial_aa_coincidence"] = False
                                            else:
                                                gnomad_df.at[index, "Initial_aa_coincidence"] = np.nan


                                            # change the switch boolean to True, meaning that it was found the alignment and the position is added
                                            switch = True
                                            # Break the for loop
                                            break

                                # Reinitialize both counters
                                compt = 0
                                compt_aa = 0

        # Check if the switch boolean is still false, meaning that after checking all pfam codes, the alignment was not found and any position was added
        if switch == False:
            # If yes, add 0 as the position of alignment
            gnomad_df.at[index, "Pos_align"] = 0
            # Add one to the cases counter, that indicates the amount of cases for which the alignment was not found
            cases += 1
        # Set the switch boolean to False, before start checking the next row
        switch = False

    print(f"\n -The {full_vs_seed} alignment, was not found in {cases:,} cases.\n\n")
    return gnomad_df


# Ask the user which version of gnomAD files is desired
version_files = input("Write the desired gnomAD files version. You can choose between:\n\tfreq-2\n\tfreq-3\n\tfreq-4\n\tfreq-5\n\tfreq-6\n\tfreq-7\n\tfreq-8\n\tno_freq\nPlease, write the version as indicated: ")
version_pattern = '' if version_files == 'no_freq' else '_' + version_files

# Ask the user if wants to compute only upper case residues from Pfam or turn all of them to upper case
#case_files = input("Pfam alignments have lower case residues indicating less conservation. Write the residues that you want to compute. You can change all residues to upper case (up) or keep Pfam structure (pfam):\n\tupper_case(up)\n\tkeep_pfam(pfam)\nPlease, write the version as indicated between parenthesis: ")
#up_low_case = 'upperCase' if case_files == 'up' else 'keep_pfam'

# Read biomart file
folder_path = "../1_Data_collection/data"
biomart_file = f"{folder_path}/biomart_human_genes_2023_05_23.txt"
biomart = pd.read_csv(biomart_file, delimiter="\t", header="infer")

# Initialize an empty dictionary of biomart file
pfam_d = {}

# Group the DataFrame by the Transcript stable ID column
biomart_transcripts = biomart.groupby("Transcript stable ID")

# Iterate over the transcripts
for transcript, group in biomart_transcripts:
    # Save the different Pfam codes related to every transcript
    pfam_codes = group["Pfam ID"].unique().tolist()
    # Save the relation between the transcript and its Pfam codes in the dictionary
    pfam_d[transcript] = pfam_codes

# Save the UniProt code related to every transcript
biomart_uniprot_d = dict(zip(biomart['Transcript stable ID'], biomart['UniProtKB/Swiss-Prot ID']))

# Filter those transcripts in the dictionary that have associated an NA
pfam_d = {key: values for key, values in pfam_d.items() if values != [np.nan]}
biomart_uniprot_d = {key: value for key, value in biomart_uniprot_d.items() if not pd.isna(value)}
#print(biomart_uniprot_d)

# Read the UniProt file
uniprot_file = f"{folder_path}/uniprot_db_human_2023_04_06_with_NM.tsv"
uniprot = pd.read_csv(uniprot_file, delimiter=",", header="infer")
#print(uniprot.columns)
uniprot_d = dict(zip(uniprot['Entry'], uniprot['Entry Name']))


# Read the filtered gnomAD files
folder_path = "../3_Data_filtering/data"
folder_path = os.path.normpath(folder_path)

# Set the path to the gnomAD file
#file_pattern = f"gnomad_r*_AD_missense_genes{args.version_pattern}_joined.txt"
file_pattern = f"gnomad_r*_AD_missense_genes{version_pattern}_joined.txt"
#print(file_pattern)

# Search the files matching the pattern
gnomad_file_list = glob.glob(f"{folder_path}/{file_pattern}")
print(gnomad_file_list)

# Add selenocysteine to protein_letters_3to1 dictionary from Bio because it was found in gnomAD
protein_letters_3to1["Sec"] = "U"

# Iterate over the different gnomAD files
for gnomad_file in gnomad_file_list:
    print(f"We start processing... {os.path.basename(gnomad_file)}")
    # Read the gnomAD file as a DataFrame
    gnomad = pd.read_csv(gnomad_file, delimiter=",", header="infer")
    #print(gnomad)
    gnomad["Pfam"] = gnomad["enst"].map(pfam_d)
    gnomad["Uniprot_entry"] = gnomad["enst"].map(biomart_uniprot_d)
    gnomad['Uniprot_entry_name'] = gnomad["Uniprot_entry"].map(uniprot_d)
    
    gnomad.to_csv(f'./data/gnomad{version_pattern}_with_pfams_no_align.txt', sep="\t", index=False)

    #print(gnomad[["enst", 'Uniprot_entry', 'Uniprot_entry_name', "Pfam"]])
    gnomad = gnomad.dropna(subset=['Pfam', 'Uniprot_entry_name'])
    #print(gnomad)
    print(Fore.MAGENTA + '--Start aligning positions--')
    print(Style.RESET_ALL)

    if full_vs_seed == 'mix':
        initial_length = len(gnomad)
        # First call with full_vs_seed = 'seed'
        gnomad = pos_align(gnomad, full_vs_seed='seed')

        # Filter out rows where Pos_align is already filled
        seed_rows = gnomad[gnomad['Pos_align'] != 0]
        #print(f'FILLED ROWS SEED: {len(seed_rows)}')

        # Filter out rows where Pos_align is not filled
        unfilled_rows = gnomad[gnomad['Pos_align'] == 0]
        #print(f'UNFILLED ROWS: {len(unfilled_rows)}')

        # Second call with full_vs_seed = 'full' for unfilled rows
        #print(f'ROWS trying FULL: {len(unfilled_rows)}')
        full_rows = pos_align(unfilled_rows, full_vs_seed='full')
        full_rows = full_rows[full_rows['Pos_align'] != 0]
        #print(f'ROWS FILLED FULL: {len(full_rows)}')

        # Concatenate the filled and unfilled rows back into a single DataFrame
        gnomad = pd.concat([seed_rows, full_rows])
        #print(f'LENGTH AFTER MIX: {len(gnomad)}')
        # This won't remove anything but just to assure that all is correct:
        gnomad = gnomad[gnomad['Pos_align'] != 0]

        # Sort the DataFrame by index to match the original order
        gnomad = gnomad.sort_index()
        print(f'Rows with alignment: {len(gnomad)}')
        print(f'We could not find seed nor full alignments in: {initial_length - len(gnomad)} rows.')

    else:
        gnomad = pos_align(gnomad, full_vs_seed)

    # Remove 0s
    filtered_gnomad = gnomad[gnomad["Pos_align"] != 0]

    if 'r2' in gnomad_file:
        gnomad_v = 'r2_1'
    elif 'r3' in gnomad_file:
        gnomad_v = 'r3'
    elif 'r4' in gnomad_file:
        gnomad_v = 'r4'
    
    # Check if data directory exists
    output_path = "./data"
    if not os.path.exists(output_path):
        # If does not exist, create it
        os.makedirs(output_path)

    
    if full_vs_seed == 'full':
        #output_file = f"{output_path}/gnomad_{gnomad_v}_pos_align_FULL_PFAM_interpro{args.version_pattern}.txt"
        output_file = f"{output_path}/gnomad_{gnomad_v}_pos_align_FULL_PFAM_interpro{version_pattern}.txt"
    elif full_vs_seed == 'seed':
        #output_file = f"{output_path}/gnomad_{gnomad_v}_pos_align_SEED_PFAM_interpro{args.version_pattern}.txt"
        output_file = f"{output_path}/gnomad_{gnomad_v}_pos_align_SEED_PFAM_interpro{version_pattern}.txt"
    elif full_vs_seed == 'mix':
        output_file = f"{output_path}/gnomad_{gnomad_v}_pos_align_SEED+FULL_PFAM_interpro{version_pattern}.txt"

    output_file = os.path.normpath(output_file)
    filtered_gnomad.to_csv(output_file, sep="\t", index=False)
