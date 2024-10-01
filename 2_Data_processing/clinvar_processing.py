# Process ClinVar file that was retrieved directly from its webserver.

import pandas as pd, numpy as np
import re
import os

def rescue_clinvar_empty_rows(row):
    """
    Adds Uniprot_entry_name and Pfam codes from UniProt to a row in ClinVar DataFrame. This function is used as a second option to rescue the empty entries for these information.
    Since only one NM code was introduced in UniProt file, we here try to relate those NM that are described in UniProt database but we have not in our UniProt file.

    Args:
        row(pandas.Series): A row from clinvar DataFrame

    Returns:
        row(pandas.Series): The input row with added Uniprot_entry_name and Pfam codes

    """
    try:
        # Check if the row is empty for the rows Uniprot_entry_name and Pfam
        if pd.isnull(row["Pfam"]) and pd.isnull(row["Uniprot_entry_name"]):
            # Get the NM value for the current row
            NM = row["NM_no_version"]

            # Find the corresponding UniProt entry in BioMart DataFrame
            entry_uniprot = biomart.loc[
                biomart["RefSeq mRNA ID"] == NM, "UniProtKB/Swiss-Prot ID"
            ].values[0]

            """ # Check if the entry_uniprot is not empty and s
            if len(entry_uniprot) > 0:
                entry_uniprot = entry_uniprot[0] """

            # Check if the entry_uniprot is not empty
            if type(entry_uniprot) != float:
                # Find the corresponding UniProt entry name and Pfam codes in UniProt DataFrame
                entry_name_uniprot = uniprot.loc[
                    uniprot["Entry"] == entry_uniprot, "Entry Name"
                ]

                # Check if the entry_uniprot is not empty
                if type(entry_name_uniprot) != float:
                    # Add the entry_uniprot to the ClinVar DataFrame
                    row["Uniprot_entry_name"] = entry_name_uniprot.iloc[0]

                # Get the Pfam code/s for the current row
                pfam = uniprot.loc[uniprot["Entry"] == entry_uniprot, "Pfam"]
                if type(pfam) != float:
                    # Add the Pfam code/s to the ClinVar DataFrame
                    row["Pfam"] = pfam.iloc[0]

    except Exception as e:
        print(f"{NM} __________Error: {e}")

    return row


# Set the path to raw data directory
folder_path = "../1_Data_collection/data"

# Set the path to the clinvar file
clinvar_file = f"{folder_path}/clinvar_homosapiens_2023_11_14.txt"
clinvar_file = os.path.normpath(clinvar_file)

# Read data extracted from clinvar database using pandas package
clinvar = pd.read_csv(clinvar_file, delimiter="\t", header="infer")

# Drop the empty column
#clinvar = clinvar.drop("Unnamed: 15", axis=1)


# Create new columns extracting information from column 'Name' in 'clinvar' DataFrame
try:
    # Create a new column 'NM'. The lambda function splits the values in the 'Name' column by whitespace, takes the first element of the resulting list, and splits it by '(' to remove any additional information after the gene name.
    clinvar["NM"] = clinvar["Name"].apply(lambda x: x.split()[0].split("(")[0])

    # Create a new column 'Prot_change' with the 3-letter code aminoacid change. Almost all cases were found to be separated by whitespace, but a single rare case was found without it. To handle it, it is also considered that can be separated with ']'.
    clinvar["Prot_change"] = clinvar["Name"].apply(
        lambda x: re.split(r"[ ]|\]", x)[-1].replace("(", "").replace(")", "")
        if "(p." in re.split(r"[ ]|\]", x)[-1]
        else None
    )

    # Create a new column 'Initial aa'. The lambda function splits the values in the 'Name' column by whitespace and takes the last element of the resulting list, which is the protein change. Extracts the first 3-letter amino acid code between '(p.' and ')', if present. If the pattern '(p.' is not found, the lambda function returns None, and the new column value is set to None.
    clinvar["Initial aa"] = clinvar["Name"].apply(
        lambda x: x.split()[-1][3:6] if "(p." in x.split()[-1] else None
    )
    # Create a new column 'Position'. Extracts the numeric position from the protein change and save it as an integer. If the protein change is None, it introduces an NAType value.
    clinvar["Position"] = pd.to_numeric(
        clinvar["Prot_change"].str.extract(r"(\d+)")[0], errors="coerce"
    ).astype("Int64")

    # Create a new column 'Final aa'. The lambda function obtains the protein change, as before. It searches the final 3-letter amino acid code in the protein change. If 'fs' is found, it is introduced as the value. If the length of the protein change is longer than 15, it is supposed to be an strange case including other confusing information and lambda function returns None.
    clinvar["Final aa"] = clinvar["Name"].apply(
        lambda x: x.split()[-1][-4:-1]
        if "(p." in x.split()[-1]
        and "fs" not in x.split()[-1]
        and len(x.split()[-1]) < 16
        else (
            "fs" if "fs" in x.split()[-1][-4:-1] and len(x.split()[-1]) < 16 else None
        )
    )
# If any error is found, add 1 to the error counter
except Exception as e:
    print(e)


### ADD UNIPROT CODES ###
# Read uniprot file
folder_path = "../1_Data_collection/data"
uniprot_file = f"{folder_path}/uniprot_db_human_2023_04_06_with_NM.tsv"
uniprot_file = os.path.normpath(uniprot_file)
uniprot = pd.read_csv(uniprot_file, delimiter=",", header="infer")
# Add a column to UniProt with a short NM code, without considering the number version
uniprot["NM_no_version"] = uniprot["NM_code"].str.split(".").str[0]

# Create a dictionary relating NM codes with UniProt entries
nm_uniprot_dict = uniprot.set_index("NM_no_version")["Entry Name"].to_dict()

# Add a column to ClinVar with a short NM code, without considering the number version
clinvar["NM_no_version"] = clinvar["NM"].str.split(".").str[0]

# Add UniProt entries to ClinVar based on NM codes
clinvar["Uniprot_entry_name"] = clinvar["NM_no_version"].map(nm_uniprot_dict)


### ADD PFAM CODES ###

print(f"Before merge: {len(clinvar)}")
# Add Pfam entries from UniProt to ClinVar. 
clinvar = pd.merge(
    clinvar, uniprot[["NM_no_version", "Pfam"]], on="NM_no_version", how="left"
)

# We drop duplicated lines
print(f"After merge: {len(clinvar)}")
clinvar = clinvar.drop_duplicates()
print(f'After drop duplicates:  {len(clinvar)}')

# There are still duplicated lines, but with different Pfam codes. We join those duplicated rows joining the Pfam codes to conserve all of them
# Convert 'Pfam' column to string
clinvar['Pfam'] = clinvar['Pfam'].astype(str)
clinvar['Pfam'] = clinvar.groupby(['Name'])['Pfam'].transform(''.join)
print(f'After join: {len(clinvar)}')
# After this join, there are duplicated lines with the joined Pfam codes, so we can drop duplicates again
clinvar = clinvar.drop_duplicates()
print(f'After drop duplicates:  {len(clinvar)}')

# Convert 'nan' string values with null values
clinvar['Pfam'] = clinvar['Pfam'].where(clinvar['Pfam'].str.startswith('PF'), np.nan)

print(f"We have {clinvar['Uniprot_entry_name'].isnull().sum()} empty rows for the Uniprot entry name and {clinvar['Pfam'].isnull().sum()} empty rows for Pfam column.")


### Some entries have not been related with an UniProt code yet. We will retry it using another database: BioMart
# Read biomart file
biomart_file = f"{folder_path}/biomart_human_genes_2023_05_23.txt"
biomart_file = os.path.normpath(biomart_file)

biomart = pd.read_csv(biomart_file, delimiter="\t", header="infer")

# Apply the add_info function to each row in df_clinvar
clinvar = clinvar.apply(rescue_clinvar_empty_rows, axis=1)
print(f"After trying to rescue more entries, we have {clinvar['Uniprot_entry_name'].isnull().sum()} empty rows for the Uniprot entry name and {clinvar['Pfam'].isnull().sum()} empty rows for Pfam column.")

# Set output file name
output_file = "./data/clinvar_homosapiens_processed.txt"
output_file = os.path.normpath(output_file)

# Write the modified DataFrame to a txt file
clinvar.to_csv(
    output_file,
    sep="\t",
    index=False,
)