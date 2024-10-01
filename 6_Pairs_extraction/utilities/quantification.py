
########################################################
#################### QUANTIFICATION ####################
########################################################

def extract_groups(pairs_df, pair_type, aa_2_compute=None):
    """
    Read a Pfam file and group it by the columns:
        - Normal pairs: 'initial_aa', 'Pos_align' and 'final_aa'
        - Homologous pairs: 'initial_aa', 'Pos_align' and/or 'final_aa'
    
    Args:
        pfam_file (String): name of a Pfam alignment file
        pair_type(str): indicating if we want to extract normal pairs or homologous pairs computed by BLOSUM62
        aa_2_compute(str): 

    Returns:
        group (tuple): Containing the grouped DataFrame by the desired columns
    """
    # Read the pfam file with pairs
    #pairs_df = pd.read_csv(pfam_file, sep="\t", header=0)

    if pair_type == 'blosum' and aa_2_compute is not None:
        if aa_2_compute == 'final_aa':
            # Return the grouped DataFrame
            return pairs_df.groupby(['initial_aa', 'Pos_align'])
        elif aa_2_compute == 'initial_aa':
            # Return the grouped DataFrame
            return pairs_df.groupby(['final_aa', 'Pos_align'])
        else:
            # Return the grouped DataFrame
            return pairs_df.groupby(['Pos_align'])
        


def extract_pairs(group, pair_type, N_N, P_P, N_P, aa_2_compute=None):
    """
    Extract pairs from a DataFrame group. Calls function save_pairs() to save the pairs in the correct variable.

    Args:
        pair_type(str): indicating if we want to extract normal pairs or homologous pairs computed by BLOSUM62
        if pair_type == normal:
            group (tuple): Containing the grouped Pfam DataFrame by the 'initial_aa', 'Pos_align' and 'final_aa' columns
        else:
            group (tuple): Containing the grouped Pfam DataFrame by the 'initial_aa' and 'Pos_align' columns
    """
    # Iterate over the entries of the group
    group = group.reset_index(drop = True)
    for index, entry in group.iterrows():
        # Save the variables of interest for this row
        pfam_code = entry['Pfam_found']
        initial_aa = entry['initial_aa']
        pos_align = entry['Pos_align']
        final_aa = entry['final_aa']
        gene_name = entry['gene']
        pos_prot_change = entry['position']
        #variation_id = entry['VariationID']
        db_1 = entry['database']
        uniprot1 = entry['Uniprot_entry_name']
        # Iterate over the entries of the group
        for index2, entry2 in group.loc[index+1:].iterrows():
            # Check that we are not comparing one row to itself
            if index != index2:
                # Save the variables of interest for this other row. Some of the previous saved variables are supposed to be common in both, that is why are not saved twice
                gene_name2 = entry2['gene']
                pos_prot_change2 = entry2['position']
                #variation_id2 = entry2['VariationID']
                db_2 = entry2['database']
                uniprot2 = entry2['Uniprot_entry_name']
                if pair_type == 'normal':
                    # Prepare the pair row format with the desired information in list format
                    if pos_prot_change > pos_prot_change2:
                        pair_row = [pfam_code, initial_aa, pos_align, final_aa, gene_name, gene_name2, pos_prot_change, pos_prot_change2, db_1, db_2, uniprot1, uniprot2]
                    elif pos_prot_change2 > pos_prot_change:
                        pair_row = [pfam_code, initial_aa, pos_align, final_aa, gene_name2, gene_name, pos_prot_change2, pos_prot_change, db_2, db_1, uniprot2, uniprot1]
                    else:
                        if gene_name >= gene_name2:
                            pair_row = [pfam_code, initial_aa, pos_align, final_aa, gene_name, gene_name2, pos_prot_change, pos_prot_change2, db_1, db_2, uniprot1, uniprot2]
                        else:
                            pair_row = [pfam_code, initial_aa, pos_align, final_aa, gene_name2, gene_name, pos_prot_change2, pos_prot_change, db_2, db_1, uniprot2, uniprot1]

                else:
                    final_aa2 = entry2['final_aa']
                    initial_aa2 = entry2['initial_aa']
                    # Prepare the pair row format with the desired information in list format
                    if aa_2_compute == 'final_aa':
                        if final_aa2 in entry['Equivalent_final_aa']:
                            pair_row = [pfam_code, initial_aa, pos_align, final_aa, final_aa2, gene_name, gene_name2, pos_prot_change, pos_prot_change2, db_1, db_2, uniprot1, uniprot2]
                        else:
                            pair_row = ''
                    elif aa_2_compute == 'initial_aa':
                        if initial_aa2 in entry['Equivalent_initial_aa']:
                            pair_row = [pfam_code, initial_aa, pos_align, final_aa, final_aa2, gene_name, gene_name2, pos_prot_change, pos_prot_change2, db_1, db_2, uniprot1, uniprot2]
                        else:
                            pair_row = ''
                    else:
                        if (initial_aa2 in entry['Equivalent_initial_aa']) and (final_aa2 in entry['Equivalent_final_aa']):
                            pair_row = [pfam_code, initial_aa, pos_align, final_aa, final_aa2, gene_name, gene_name2, pos_prot_change, pos_prot_change2, db_1, db_2, uniprot1, uniprot2]
                        else:
                            pair_row = ''
                #print('yes')
                # Call the save_pairs() function to save the pair row in the correct variable, which determines the pair type
                save_pairs(pair_row, db_1, db_2, N_N, P_P, N_P)

def save_pairs(pair_row, db_1, db_2, N_N, P_P, N_P):
    """
    Save a pair in the correct pair type.

    Args:
        - pair_row(list): Containing the necessary information to describe the pair
        - db_1 (string): Indicates if the first entry of the pair pertains to gnomAD or ClinVar
        - db_2 (string): Indicates if the second entry of the pair pertains to gnomAD or ClinVar
    """
    if pair_row:
        #print(pair_row)
        # Check if databases in both entries are equal
        if db_1 == db_2:
            # If yes, check if the database is ClinVar
            if db_1 == 'ClinVar':
                # If yes, add the pair entry to P_P (Pathogenic/Pathogenic) list
                P_P.append(pair_row)
            elif db_1 == 'gnomAD':
                # If is gnomAD, add the pair entry to N_N (Non-pathogenic/Non-pathogenic) list
                N_N.append(pair_row)
        else:
            # If are not equal, add the pair entry to N_P (Non-pathogenic/Pathogenic) list
            N_P.append(pair_row)
