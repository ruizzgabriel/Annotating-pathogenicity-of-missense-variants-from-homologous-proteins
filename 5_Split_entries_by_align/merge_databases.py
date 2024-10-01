# Merge variants from ClinVar and Gnomad, but knowing which ones come from one database and which ones from the other

import os
import pandas as pd
import glob
import argparse
from colorama import Fore, Style

""" # Create the parser
parser = argparse.ArgumentParser()
# Add an argument
parser.add_argument('--version_pattern', type=str, help='Version pattern')
# Parse the arguments
args = parser.parse_args() """


# Ask the user which version of gnomAD files is desired
version_files = input("Write the desired gnomAD files version. You can choose between:\n\tfreq-2\n\tfreq-3\n\tfreq-4\n\tfreq-5\n\tfreq-6\n\tfreq-7\n\tfreq-8\n\tno_freq\nPlease, write the version as indicated: ")
version_pattern = '' if version_files == 'no_freq' else '_' + version_files

full_vs_seed = input("Write the desired alignment type. (Mix condition prioritizes seed alignment and uses full ones when seed is not found)\nYou can choose between:\n\tfull\n\tseed\n\tseed+full(mix)\nPlease, write the version as indicated: ")


# Set data path
data_path = "../4_Align_pos/data"

##### ClinVar #####
# Set ClinVar file path
""" if up_low_case == "upperCase":
    clinvar_path = f"{data_path}/clinvar_pos_align_ALL_PFAM_interpro_{up_low_case}.txt"
else:
    clinvar_path = f"{data_path}/clinvar_pos_align_ALL_PFAM_interpro.txt"
 """

# Read ClinVar file
if full_vs_seed == 'full':
    clinvar_path = f"{data_path}/clinvar_pos_align_FULL_PFAM_interpro.txt" 
elif full_vs_seed == 'seed':
    clinvar_path = f"{data_path}/clinvar_pos_align_SEED_PFAM_interpro.txt"
elif full_vs_seed == 'mix':
    clinvar_path = f"{data_path}/clinvar_pos_align_SEED+FULL_PFAM_interpro.txt"

clinvar_path = os.path.normpath(clinvar_path)
clinvar_df = pd.read_csv(clinvar_path, sep="\t", header=0)

# Simplify ClinVar file to avoid non-used columns
columns_to_remove = [
    "Name",
    "Protein change",
    "GRCh37Chromosome",
    "GRCh37Location",
    "GRCh38Chromosome",
    "GRCh38Location",
    "NM",
    "Pfam",
]
# Use the drop() method to remove the specified columns
clinvar_df = clinvar_df.drop(columns_to_remove, axis=1)

# Add a column indicating that all these rows are from ClinVar
clinvar_df["database"] = "ClinVar"
# print(clinvar_df.columns)

# Change column names in order to be equal to gnomAD
# Set a dictionary with old names and new ones.
column_name_changes = {
    "Gene(s)": "gene",
    "Clinical significance (Last reviewed)": "clinical_significance",
    "NM_no_version": "NM",
    "Prot_change": "prot_change",
    "Initial aa": "initial_aa",
    "Position": "position",
    "Final aa": "final_aa",
    "dbSNP ID": "rsid",
}
# Rename the columns using the dictionary
clinvar_df = clinvar_df.rename(columns=column_name_changes)
# print(clinvar_df.columns)


##### gnomAD #####

# Set the path to the gnomAD file
if full_vs_seed == 'full':
    #gnomad_pattern = f"{data_path}/gnomad_r*_pos_align_ALL_PFAM_interpro{args.version_pattern}.txt"
    gnomad_pattern = f"{data_path}/gnomad_r*_pos_align_FULL_PFAM_interpro{version_pattern}.txt"
elif full_vs_seed == 'seed':
    gnomad_pattern = f"{data_path}/gnomad_r*_pos_align_SEED_PFAM_interpro{version_pattern}.txt"
elif full_vs_seed == 'mix':
    gnomad_pattern = f"{data_path}/gnomad_r*_pos_align_SEED+FULL_PFAM_interpro{version_pattern}.txt"

gnomad_pattern = os.path.normpath(gnomad_pattern)

print(f'We are using:\n\t{clinvar_path}\n\t{gnomad_pattern}')
# Search the files matching the pattern
gnomad_file_list = glob.glob(gnomad_pattern)
# print(gnomad_file_list)

# Initialize a new list to save all gnomAD DataFrames
gnomad_df_list = []
# Iterate over the gnomAD files
for gnomad_file in gnomad_file_list:
    # Read gnomad_file
    gnomad_df = pd.read_csv(gnomad_file, delimiter="\t", header="infer")

    # Remove consequence column
    remove_cols_gnomad = ["consequence", "Pfam"]
    gnomad_df = gnomad_df.drop(remove_cols_gnomad, axis=1)

    # Save the gnomAD version
    #gnomad_v = "gnomAD_r2_1" if "r2" in gnomad_file else "gnomAD_r3"

    # Add a column indicating that all these rows are from gnomAD and which version
    gnomad_df["database"] = "gnomAD"     # gnomad_v
    # Add the df to the list
    gnomad_df_list.append(gnomad_df)

# Concatenate the gnomAD DataFrames in the list
all_gnomad = pd.concat(gnomad_df_list)


##### Concatenate databases #####

# Concatenate the ClinVar DataFrame with the concatenated gnomAD DataFrame
all_db_entries = pd.concat([clinvar_df, all_gnomad], ignore_index=True, axis=0)

print(Fore.GREEN + '\n\nSTATISTICS:\n')
print(Style.RESET_ALL)
print(f' -Only merge: {len(all_db_entries.index):,}')


##### Remove variants in which initial aa is not the same as in the Pfam alignment #####

# Filter merged DataFrame to avoid cases in which the initial aa was not equal in Pfam and the protein change, meaning that Pfam numeration was not equal to NM numeration
all_db_entries = all_db_entries[all_db_entries['Initial_aa_coincidence'] == True]
print(f' -Filtering equal initial amino acid: {len(all_db_entries.index):,}')


##### Remove entries with unclear consequences --> those appearing in both databases (patho and non-patho) #####

# Create an empty mask with all False (we will set True later for the rows we want to keep)
mask = pd.Series(False, index=all_db_entries.index)

# Group the df for detecting variants with unclear consequences
grouped_df = all_db_entries.groupby(['gene' ,'initial_aa', 'position', 'final_aa'])

for group_name, group in grouped_df:
    if len(group) > 1 and group['database'].nunique() > 1:
        #if len(group) > 2:
        #    print(group[['gene' ,'initial_aa', 'position', 'final_aa', 'database']])
        # Don't change the mask (keep it False) for these rows
        continue
    # Otherwise, set the mask to True for these rows
    mask[group.index] = True

# Use the mask to index all_db_entries and get a DataFrame with only the rows we want to keep
all_db_entries = all_db_entries[mask]
#print(all_db_entries)
print(f' -After the removal of non-clear pathogenicity annotation: {len(all_db_entries.index):,}')


##### Quantify number of variants of gnomAD and ClinVar #####
print(Fore.YELLOW + '\n\nQUANTIFY VARIANTS:\n')
print(Style.RESET_ALL)

# Assuming 'df' is your DataFrame
count_clin = all_db_entries['database'].value_counts().get('ClinVar', 0)
print(f" -ClinVar: There are {count_clin:,} pathogenic variants from ClinVar.")

count_gnom = all_db_entries['database'].value_counts().get('gnomAD', 0)
print(f" -gnomAD: There are {count_gnom:,} non-pathogenic variants from gnomAD.\n\n")


# Check if data directory exists
output_path = "./data"
if not os.path.exists(output_path):
    # If does not exist, create it
    os.makedirs(output_path)

# Set output file name and normalize its path
if full_vs_seed == 'full':
    #output_file = f"{output_path}/merged_db_interpro{args.version_pattern}.txt"
    output_file = f"{output_path}/FULL_align/merged_db_interpro{version_pattern}.txt"

elif full_vs_seed == 'seed':
    #output_file = f"{output_path}/merged_db_interpro_SEED_{args.version_pattern}.txt"
    output_file = f"{output_path}/SEED_align/merged_db_interpro_SEED{version_pattern}.txt"
elif full_vs_seed == 'mix':
    output_file = f"{output_path}/SEED+FULL_align/merged_db_interpro_SEED+FULL{version_pattern}.txt"

output_file = os.path.normpath(output_file)
# Save the concatenated DataFrame containing the entries from the three files
all_db_entries.to_csv(output_file, sep="\t", index=False)
