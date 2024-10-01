import pandas as pd, os
from Bio.Data.IUPACData import protein_letters_3to1
import blosum as bl

######################################################
#################### Normal pairs ####################
######################################################

def group_df_normal_pairs(alignment, version_pattern, full_vs_seed):
    # sourcery skip: extract-method
    """
    Iterate over all files in which we divided the entries, extract paired entries for columns initial_aa, final_aa and Pos_align and save them in a single DataFrame.
    Args:
        pfam_list (list): containing the name of the files where the entries are saved
    
    """
    # Initialize an empty DataFrame
    pairs_df = pd.DataFrame()

    # Open the file as a DataFrame
    pfam = pd.read_csv(alignment, sep="\t", header=0)
    if len(pfam) > 1:
        # Save the pfam code in order to save an output file with its name
        pfam_code = pfam['Pfam_found'][0]
        # Group entries by the desired columns and filter by groups larger than 1 entry
        pfam = pfam.groupby(['initial_aa', 'final_aa', 'Pos_align']).filter(lambda x: len(x) > 1)
    
        # Append pairs to the whole DataFrame
        #pairs_df = pairs_df.append(pfam, ignore_index=True)
        pairs_df = pfam.reset_index()

        # Remove duplicates from the whole DataFrame
        pairs_df = pairs_df.drop_duplicates()

        if len(pairs_df) > 0:
            # Check if data directory exists
            if full_vs_seed == 'full':
                output_path = f"./data/normal_pairs/FULL_align/pfam_pairs_AD_missense{version_pattern}/pairs_by_pfam_code"
            elif full_vs_seed == 'seed':
                output_path = f"./data/normal_pairs/SEED_align/pfam_pairs_AD_missense{version_pattern}/pairs_by_pfam_code"
            else:
                output_path = f"./data/normal_pairs/SEED+FULL_align/pfam_pairs_AD_missense{version_pattern}/pairs_by_pfam_code"

            if not os.path.exists(output_path):
                # If does not exist, create it
                os.makedirs(output_path)

            # Set the output file name
            output_file = f"{output_path}/{pfam_code}_pairs{version_pattern}.txt"
            output_file = os.path.normpath(output_file)

            # Save the DataFrame in a txt file
            pairs_df.to_csv(output_file, sep="\t", index=False)
    
            return pairs_df.groupby(['initial_aa', 'Pos_align', 'final_aa'])
        else:
            return None





#######################################################################################################################################
def get_positive_scores_blosum62(final_aa_list):
    """
    Create a dictionary with the 'final_aa' present in a DataFrame group. It relates the aa with its scores in the BLOSUM 62 matrix.
    It is created for each group, to only retrieve the information of the aminoacids (aa) present in the group.
    Only positive scores are saved, because they correspond to those aminoacidic changes that we consider equivalent.

    Args:
        final_aa_list (list): Containing the final aminoacids for the treated group
    
    Returns:
        eq_aa_change_d (dict): Containing the scores between the aminoacids using the BLOSUM 62 matrix
    """
    # Save BLOSUM 62 matrix from blosum package
    blosum62 = bl.BLOSUM(62)
    # Initialize an empty dictionary to save the scores
    eq_aa_change_d = {}
    # Iterate over the aas in the list
    for aa_1 in final_aa_list:              
        # Check that the aa is not in the dict
        if aa_1 not in eq_aa_change_d:
            # If not, save it in the dict
            eq_aa_change_d[aa_1] = {}
            # Iterate over the aas in the list
            for aa_2 in final_aa_list:          
                # Compute the score using the BLOSUM 62 matrix and both aas
                score = blosum62[aa_1][aa_2]
                # Save only the positive scores
                if score > 0:
                    eq_aa_change_d[aa_1][aa_2] = score

    return eq_aa_change_d


def get_equivalent_initial_aa(group):
    """
    Extract the different 'initial_aa' for a group and compute the score using function get_positive_scores_blosum62.
    Obtain the possible aa substitutions based on the BLOSUM 62 scores.

    Args:
        group (tuple): Containing different rows of the merged DataFrame (from gnomAD and ClinVar) that match in 'final_aa' and 'Pos_align'.
    
    Returns:
        group (tuple): Same as in the input but with a new column (Equivalent_initial_aa) indicating the possible substitutions based on the BLOSUM62 positive scores 
        considering only the aas present inside the group
    """
    
    # Obtain the aas present in the protein change's 'initial_aa'
    aa_list = group['initial_aa'].tolist()

    # Transform the 'initial_aa' 3-letter codes to 1-letter code using protein_letters_3to1 dictionary from Bioconductor
    aa_list = [protein_letters_3to1[aa] for aa in aa_list if aa in protein_letters_3to1]
    #print(aa_list)
    # Get positive scores dictionary between possible equivalent aminoacidic changes in initial_aa
    eq_aa_change_d = get_positive_scores_blosum62(aa_list)
    #print(eq_aa_change_d)

    # Create a list to hold the indices to be dropped
    indices_to_drop = []

    group = group.reset_index(drop = True)
    # Iterate over the rows of the group
    for index, entry in group.iterrows():
        # Transform the 'initial_aa' 3-letter code to 1-letter code to be able to compare it
        if entry['initial_aa'] in protein_letters_3to1:
            initial_aa = protein_letters_3to1[entry['initial_aa']]
            # Save the 1-letter code in the 'initial_aa' column
            group.loc[index, 'initial_aa'] = initial_aa

            # Get the dictionary keys of the possible equivalent aas, which are 1-letter codes, in a list format
            equivalent_aa = list(eq_aa_change_d[initial_aa])
            # Save the possible equivalent aas in string format without separator
            group.at[index, 'Equivalent_initial_aa'] = ''.join(equivalent_aa)
        else:
            #group.at[index, 'Equivalent_initial_aa'] = ''
            # Add the index to the list of indices to be dropped
            indices_to_drop.append(index)
        # Drop the rows from the group
    group = group.drop(indices_to_drop)

    return group

def get_equivalent_final_aa(group):
    """
    Extract the different 'final_aa' for a group and compute the score using function get_positive_scores_blosum62.
    Obtain the possible aa substitutions based on the BLOSUM 62 scores.

    Args:
        group (tuple): Containing different rows of the merged DataFrame (from gnomAD and ClinVar) that match in 'initial_aa' and 'Pos_align'.
    
    Returns:
        group (tuple): Same as in the input but with a new column (Equivalent_final_aa) indicating the possible substitutions based on the BLOSUM62 positive scores 
        considering only the aas present inside the group
    """
    
    # Obtain the aas present in the protein change's 'final_aa'
    aa_list = group['final_aa'].tolist()

    # Transform the 'final_aa' 3-letter codes to 1-letter code using protein_letters_3to1 dictionary from Bioconductor
    aa_list = [protein_letters_3to1[aa] for aa in aa_list if aa in protein_letters_3to1]
    #print(aa_list)
    # Get positive scores dictionary between possible equivalent aminoacidic changes in final_aa
    eq_aa_change_d = get_positive_scores_blosum62(aa_list)
    #print(eq_aa_change_d)

    group = group.reset_index(drop = True)
    # Iterate over the rows of the group
    for index, entry in group.iterrows():
        # Transform the 'final_aa' 3-letter code to 1-letter code to be able to compare it
        final_aa = protein_letters_3to1[entry['final_aa']]
        # Save the 1-letter code in the 'final_aa' column
        group.loc[index, 'final_aa'] = final_aa

        # Get the dictionary keys of the possible equivalent aas, which are 1-letter codes, in a list format
        equivalent_aa = list(eq_aa_change_d[final_aa])
        # Save the possible equivalent aas in string format without separator
        group.at[index, 'Equivalent_final_aa'] = ''.join(equivalent_aa)

    return group

def get_pairs(group, aa_2_compute):  # sourcery skip: use-contextlib-suppress
    """
    Retrieve the pairs for a group of merged DataFrame that contains ClinVar and gnomAD information. Paired rows are considered those with\n
    - the same 'initial_aa', 'Pos_align' and equivalent 'final_aa'\n
    - the same 'final_aa', 'Pos_align' and equivalent 'initial_aa'\n
    It will be computed differently depending on the user's choice.\n
    
    Possible equivalent aas for every row are indicated in Equivalent_final_aa or Equivalent_initial_aa column.

    Args:
        group (tuple): Containing different rows of the merged DataFrame (from gnomAD and ClinVar) that match in 'initial_aa' and 'Pos_align'.
    
    Returns: 
        group (tuple): Containing the paired rows considering equivalent 'final_aa'
    
    """
    # Set new column 'pertain_to_pair' as False in all group 
    group['pertain_to_pair'] = False
    # Iterate over the entries of the group
    group = group.reset_index(drop = True)
    for index, entry in group.iterrows():
        #final_aa = entry['final_aa']
        #print(final_aa, '//////////', equivalent_aas)
        equivalent_aas = ''
        try:
            # Obtain the possible equivalent final_aa or initial_aa from its column (in string format)
            equivalent_aas = entry[f'Equivalent_{aa_2_compute}']
        except KeyError:
            continue     #group[group['pertain_to_pair'] == True]
        
        if equivalent_aas:
            # Iterate again over the entries of the group
            for index2, entry2 in group.iterrows():
                #print(index, '::::::::::::::::.', index2)

                # Iterate over the equivalent aas in the row of the first iteration
                for aa in equivalent_aas:
                    #print(aa, '........', entry2['final_aa'])

                    # Check if the aa matches with the final_aa of the entry in the second iteration, and if both indices does not match, meaning that is not comparing the same entry
                    if aa == entry2[aa_2_compute] and index != index2:
                        #print(aa, '-----', entry2['final_aa'])
                        # If conditions are accomplished, set pertain_to_pair column for both rows as True
                        group.at[index, 'pertain_to_pair'] = True
                        group.at[index2, 'pertain_to_pair'] = True
                        # Break the loop
                        break
    # Return only the rows of the group that pertain to a pair
    return group[group['pertain_to_pair'] == True]


def get_pairs_both_aa(group):  # sourcery skip: use-contextlib-suppress
    """
    Retrieve the pairs for a group of merged DataFrame that contains ClinVar and gnomAD information. Paired rows are considered those with\n
    - the same 'initial_aa', 'Pos_align' and equivalent 'final_aa'\n
    - the same 'final_aa', 'Pos_align' and equivalent 'initial_aa'\n
    It will be computed differently depending on the user's choice.\n
    
    Possible equivalent aas for every row are indicated in Equivalent_final_aa or Equivalent_initial_aa column.

    Args:
        group (tuple): Containing different rows of the merged DataFrame (from gnomAD and ClinVar) that match in 'initial_aa' and 'Pos_align'.
    
    Returns: 
        group (tuple): Containing the paired rows considering equivalent 'final_aa'
    
    """
    # Set new column 'pertain_to_pair' as False in all group 
    group['pertain_to_pair'] = False
    # Iterate over the entries of the group
    group = group.reset_index(drop = True)
    for index, entry in group.iterrows():
        #final_aa = entry['final_aa']
        #print(final_aa, '//////////', equivalent_aas)
        equivalent_initial_aas = ''
        equivalent_final_aas = ''
        try:
            # Obtain the possible equivalent final_aa or initial_aa from its column (in string format)
            equivalent_initial_aas = entry['Equivalent_initial_aa']
            equivalent_final_aas = entry['Equivalent_final_aa']
        except KeyError:
            continue     #group[group['pertain_to_pair'] == True]
        
        if equivalent_initial_aas:
            # Iterate again over the entries of the group
            for index2, entry2 in group.iterrows():
                #print(index, '::::::::::::::::.', index2)

                # Iterate over the equivalent aas in the row of the first iteration
                for aa in equivalent_initial_aas:
                    #print(aa, '........', entry2['final_aa'])

                    # Check if the aa matches with the final_aa of the entry in the second iteration, and if both indices does not match, meaning that is not comparing the same entry
                    if aa == entry2['initial_aa'] and index != index2:
                        #print(aa, '-----', entry2['final_aa'])
                        # If conditions are accomplished, set pertain_to_pair column for both rows as True
                        group.at[index, 'pertain_to_pair'] = True
                        group.at[index2, 'pertain_to_pair'] = True
                        # Break the loop
                        break

        if equivalent_final_aas:
            # Iterate again over the entries of the group
            for index2, entry2 in group.iterrows():
                #print(index, '::::::::::::::::.', index2)

                # Iterate over the equivalent aas in the row of the first iteration
                for aa in equivalent_final_aas:
                    #print(aa, '........', entry2['final_aa'])

                    # Check if the aa matches with the final_aa of the entry in the second iteration, and if both indices does not match, meaning that is not comparing the same entry
                    if aa == entry2['final_aa'] and index != index2:
                        #print(aa, '-----', entry2['final_aa'])
                        # If conditions are accomplished, set pertain_to_pair column for both rows as True
                        group.at[index, 'pertain_to_pair'] = True
                        group.at[index2, 'pertain_to_pair'] = True
                        # Break the loop
                        break
    # Return only the rows of the group that pertain to a pair
    return group[group['pertain_to_pair'] == True]


def group_df_chosen_aa(pfam_df, aa_2_compute):
    """
    Group the Pfam DataFrame depending on the aa that has to be computed by BLOSUM62.

    Args:
        pfam_df (DataFrame): DataFrame containing the information of the entries in the corresponding Pfam alignment
        aa_2_compute (string): Describing which aa has to be computed by BLOSUM62 to find similar aa

    Returns:
        grouped_df (DataFrame): Input DataFrame grouped by columns of interest and filtered to obtain groups bigger than 1 entry
    """
    # Check if BLOSUM62 has to be computed for initial_aa
    if aa_2_compute == 'initial_aa':
        # Group entries by the final_aa and Pos_align columns and filter by groups larger than 1 entry
        pfam_df = pfam_df.groupby(['final_aa', 'Pos_align']).filter(lambda x: len(x) > 1)
        # Group again by the same columns to be able to iterate over them
        return pfam_df.groupby(['final_aa', 'Pos_align'])

    # Check if BLOSUM62 has to be computed for final_aa
    elif aa_2_compute == 'final_aa':
        # Group entries by the initial_aa and Pos_align columns and filter by groups larger than 1 entry
        pfam_df = pfam_df.groupby(['initial_aa', 'Pos_align']).filter(lambda x: len(x) > 1)
        # Group again by the same columns to be able to iterate over them
        return pfam_df.groupby(['initial_aa', 'Pos_align'])

    # Check if BLOSUM62 has to be computed for both initial_aa and final_aa
    else:
        # Group entries by the Pos_align column and filter by groups larger than 1 entry
        pfam_df = pfam_df.groupby(['Pos_align']).filter(lambda x: len(x) > 1)
        # Group again by the same columns to be able to iterate over them
        return pfam_df.groupby(['Pos_align'])


def extract_blosum_pairs(alignment, aa_2_compute, version_pattern, full_vs_seed):
    """
    Iterate over all Pfam files in which we divided the entries and extract pairs with similar initial_aa and/or final_aa and Pos_align. Finally, save the filtered pairs in separated txt files for each Pfam alignment.
    
    Args:
        pfam_list (list): containing the name of the files where the entries are saved
    
    """
    #print(f'------------------- Start with {alignment}')
    # Open the file as a DataFrame
    pfam = pd.read_csv(alignment, sep="\t", header=0)
    if len(pfam) > 1:
        # Save the pfam code in order to save an output file with its name
        pfam_code = pfam['Pfam_found'][0]

        # Conserve the rows that match the initial aa in the protein change with the aa present in its pfam alignment
        pfam = pfam[pfam['Initial_aa_coincidence'] == True]

        # Group the Pfam df depending on the aa that has to be computed by BLOSUM62
        pfam_grouped = group_df_chosen_aa(pfam, aa_2_compute)      
        
        # Initialize groups list
        groups_list = []
        # Iterate over the groups in grouped Pfam DataFrame
        #############################################################################################################
        for group_name, group in pfam_grouped:
            if aa_2_compute == 'initial_aa':
                # Obtain the possible equivalent initial_aa changes to obtain pairs within the group
                pfam_eq_aas = get_equivalent_initial_aa(group)
                # Obtain pairs within the group
                pfam_pairs = get_pairs(pfam_eq_aas, aa_2_compute)
                # Check that the group is not empty
                if not pfam_pairs.empty:
                    # Append the pairs to groups list
                    groups_list.append(pfam_pairs)
            
            elif aa_2_compute == 'final_aa':
                # Obtain the possible equivalent final_aa changes to obtain pairs within the group
                pfam_eq_aas = get_equivalent_final_aa(group)
                # Obtain pairs within the group
                pfam_pairs = get_pairs(pfam_eq_aas, aa_2_compute)
                # Check that the group is not empty
                if not pfam_pairs.empty:
                    # Append the pairs to groups list
                    groups_list.append(pfam_pairs)
            
            else:
                # This condition considers the option 'both_aa'
                # Obtain the possible equivalent initial_aa changes to obtain pairs within the group
                pfam_eq_aas = get_equivalent_initial_aa(group)
                # Obtain the possible equivalent final_aa changes to obtain pairs within the group
                pfam_eq_aas = get_equivalent_final_aa(pfam_eq_aas)

                # Get pairs
                pfam_pairs = get_pairs_both_aa(pfam_eq_aas)
                #pfam_pairs = get_pairs(pfam_eq_aas, 'initial_aa')
                #pfam_pairs = get_pairs(pfam_pairs, 'final_aa')
                # Check that the group is not empty
                if not pfam_pairs.empty:
                    # Append the pairs to groups list
                    groups_list.append(pfam_pairs)

        if groups_list:
            # Concatenate groups to obtain a single DataFrame for all pairs in a Pfam alignment file
            pfam_pairs = pd.concat(groups_list)

            # Remove duplicates from the whole DataFrame
            pfam_pairs.drop_duplicates(inplace=True)

            # Set the output file name
            if full_vs_seed == 'full':
                output_path = f"./data/homologous_pairs/FULL_align/blosum62_{aa_2_compute}/pfam_pairs_interpro_+score{version_pattern}/pairs_by_pfam_code"
            elif full_vs_seed == 'seed':
                output_path = f"./data/homologous_pairs/SEED_align/blosum62_{aa_2_compute}/pfam_pairs_interpro_+score{version_pattern}/pairs_by_pfam_code"
            else:
                output_path = f"./data/homologous_pairs/SEED+FULL_align/blosum62_{aa_2_compute}/pfam_pairs_interpro_+score{version_pattern}/pairs_by_pfam_code"

            output_path = os.path.normpath(output_path)
            
            # Check if the output data directory exists
            if not os.path.exists(output_path):
                # If does not exist, create it
                os.makedirs(output_path)
            
            output_file = f"{output_path}/{pfam_code}_pairs_interpro_+score{version_pattern}.txt"
            output_file = os.path.normpath(output_file)
            #print(f'{pfam_code} done')
            # Save the DataFrame in a txt file
            pfam_pairs.to_csv(output_file, sep="\t", index=False)

            return pfam_pairs
        
        else:
            return None

