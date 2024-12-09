from Bio.Data.IUPACData import protein_letters_3to1, protein_letters_1to3
import os
import pandas as pd
import json
from subprocess import Popen, PIPE
from io import StringIO
from sqlalchemy import text
import checking_tools


def search_var_annotation(variants_df, gene_name, ref_aa, pos, mut_aa):
    """
    Search if the variant is previously annotated in ClinVar or gnomAD. Use the file with our restrictions and filters.
    """

    patho_l = ['Likely pathogenic', 'Likely risk allele', 'Pathogenic', 'Pathogenic/Likely pathogenic', 'Pathogenic/Likely risk allele', 'Pathogenic; risk factor']
    neutral_l = ['Neutral', 'Benign', 'Benign/Likely benign', 'Likely benign']
    
    existing_vars = ((variants_df['gene'] == gene_name) & (variants_df['initial_aa'] == ref_aa) & (variants_df['final_aa'] == mut_aa) & (variants_df['position'] == pos))
    filtered_df = variants_df[existing_vars]

    if not filtered_df.empty:
        annot_list = filtered_df['clinical_significance'].tolist()
        
        # Ensure filtered_df is explicitly copied before modification
        filtered_df = filtered_df.copy()

        filtered_df.drop(columns=['initial_aa', 'position', 'final_aa'], inplace=True)
        filtered_df.rename(columns={'gene':'Gene', 'clinical_significance': 'Clinical significance', 'prot_change': 'Protein change', 'NM': 'RefSeq NM code', 'database': 'Database'}, inplace=True)
        pd.set_option('future.no_silent_downcasting', True)
        filtered_df = filtered_df.fillna('-')
        # Convert lists in the 'RefSeq NM code' column to strings
        filtered_df['RefSeq NM code'] = filtered_df['RefSeq NM code'].apply(lambda x: ', '.join(x) if isinstance(x, list) else str(x))

        
        if (any(item in neutral_l for item in annot_list)) and (any(item in patho_l for item in annot_list)):
            warn_user = 'This variant has <b>uncertain pathogenicity</b> because is annotated as disease-causing in ClinVar and as non-disease causing in gnomAD.'
            return warn_user, filtered_df
        elif any(item in neutral_l for item in annot_list):
            warn_user = 'This variant is annotated as <b>non-pathogenic</b> in gnomAD.'
            return warn_user, filtered_df
        elif any(item in patho_l for item in annot_list):
            warn_user = 'This variant is annotated as <b>pathogenic</b> in ClinVar.'
            return warn_user, filtered_df
    else:
        return False, None
    
    
def search_annot_clinvar_gnomad(gene_name, uniprot_entry_name, ref_aa, pos, mut_aa):
    """
    Search annotation in ClinVar and gnomAD global files.
    """
    patho_l = ['Likely pathogenic', 'Likely risk allele', 'Pathogenic', 'Pathogenic/Likely pathogenic', 'Pathogenic/Likely risk allele', 'Pathogenic; risk factor']
    neutral_l = ['Neutral', 'Benign', 'Benign/Likely benign', 'Likely benign']

    clinvar_df = pd.DataFrame()
    gnomad_df = pd.DataFrame()
    # Check if exists in ClinVar
    original_clinvar = pd.read_csv('/var/www/homologous/data/clinvar_homosapiens_processed.txt', sep='\t', usecols=['Gene(s)', 'Uniprot_entry_name', 'Prot_change', 'Clinical significance (Last reviewed)', 'Initial aa', 'Position', 'Final aa', 'VariationID'])
    existing_clinvar = ((original_clinvar['Uniprot_entry_name'] == uniprot_entry_name) & (original_clinvar['Initial aa'] == ref_aa) & (original_clinvar['Final aa'] == mut_aa) & (original_clinvar['Position'] == pos))
    clinvar_df = original_clinvar[existing_clinvar]
    
    # Check if exists in gnomAD
    gnomad_file = '/var/www/homologous/data/gnomad_r4_GRCh38_missense_vars.txt'
    # Construct the awk command
    """ awk_command = (
        f"awk -F',' '$1 == \"{gene_name}\" && $6 == \"{pos}\" && "
        f"$7 == \"{ref_aa}\" && $8 == \"{mut_aa}\"' {gnomad_file}"
    ) """

    grep_command = f"grep '{gene_name}' {gnomad_file}"
    awk_command = (
        f"awk -F',' -v pos='{pos}' -v ref_aa='{ref_aa}' -v mut_aa='{mut_aa}' "
        "'$6 == pos && $7 == ref_aa && $8 == mut_aa'"
    )

    command = f"{grep_command} | {awk_command}"
    
    # Run the command with Popen
    with Popen(command, shell=True, stdout=PIPE, stderr=PIPE, text=True) as process:
        output, error = process.communicate()

        if process.returncode == 0 and output:
            header = ['gene','enst','NM','rsid','prot_change','position','initial_aa','final_aa','consequence','clinical_significance','genome_af','exome_af']
            gnomad_df = pd.read_csv(StringIO(output), sep=",", names=header)
            #print(gnomad_df)
            if output.strip(): # check the output is not empty
                for line in output.splitlines():
                    if line.strip():  # check if the line is not empty
                        line = str(line)
                        for patho_v in patho_l:
                            if patho_v in line:
                                annot_list.append('Pathogenic')
                                break
                        for non_patho_v in neutral_l:
                            if non_patho_v in line:
                                annot_list.append('Neutral')
                                break
        else:
            print("Error:\n", error)

    # Query to select specific columns and read directly into a DataFrame
    """ query = f"SELECT * FROM homolvar_gnomad_no_freq WHERE Gene = '{gene_name}' AND position = {pos} AND initial_aa = '{ref_aa}' AND final_aa = '{mut_aa}';"
    gnomad_df = pd.read_sql_query(query, engine) """


    annot_list = []
    if not clinvar_df.empty:
        annot_list_raw = clinvar_df['Clinical significance (Last reviewed)'].tolist()
        for annot in annot_list_raw:
            annot = annot.split('(')[0]
            annot_list.append(annot)
        ######   !!!!!!
        clinvar_df.rename(columns={'Gene(s)':'Gene', 'Clinical significance (Last reviewed)': 'Clinical significance', 'Prot_change':'Protein change'}, inplace=True)
        clinvar_df.drop(columns=['Uniprot_entry_name', 'Initial aa', 'Position', 'Final aa'], inplace=True)
        clinvar_df['Database'] = 'ClinVar'
        #print(clinvar_df)
    if not gnomad_df.empty:
        ######   !!!!!!
        #print(gnomad_df.columns)
        # SQL APPROACH
        #gnomad_df.rename(columns={'RefSeq_NM_code':'RefSeq NM code', 'Clinical_significance':'Clinical significance', 'Protein_change':'Protein change'}, inplace=True)
        # BASH QUERY APPROACH
        gnomad_df.rename(columns={'gene':'Gene', 'clinical_significance': 'Clinical significance', 'NM':'RefSeq NM code', 'Clinical_significance':'Clinical significance', 'prot_change':'Protein change'}, inplace=True)

        annot_list_raw = gnomad_df['Clinical significance'].tolist()
        for annot in annot_list_raw:
            annot_list.append(annot)

        # SQL APPROACH
        #gnomad_df.drop(columns=['id','initial_aa', 'position', 'final_aa'], inplace=True)
        # BASH QUERY APPROACH
        gnomad_df.drop(columns=['initial_aa', 'position', 'final_aa', 'enst', 'genome_af', 'exome_af'], inplace=True)
        gnomad_df['Database'] = 'gnomAD'
        #print(gnomad_df)
    
    if (not clinvar_df.empty) and (not gnomad_df.empty):
        output_df = pd.concat([clinvar_df, gnomad_df], ignore_index=True, sort=False)
    elif not clinvar_df.empty:
        output_df = clinvar_df
    elif not gnomad_df.empty:
        output_df = gnomad_df
    else:
        return None, None
    
    warn_user = ''
    # Check the list of annotations filled by ClinVar and gnomAD
    if (any(item in neutral_l for item in annot_list)) and (any(item in patho_l for item in annot_list)):
        annot = 'Uncertain'
        warn_user = 'This variant has <b>uncertain pathogenicity</b> because is annotated as disease-causing in ClinVar and as non-disease causing in gnomAD.'
    elif any(item in neutral_l for item in annot_list):
        #print(annot_list)
        warn_user = 'This variant is annotated as <b>non-pathogenic</b> in gnomAD.'
        annot = 'Neutral'
    elif any(item in patho_l for item in annot_list):
        annot = 'Patho'
        warn_user = 'This variant is annotated as <b>pathogenic</b> in ClinVar.'
    else:
        db_list = output_df['Database'].tolist()
        #print(db_list)
        if 'ClinVar' in db_list and 'gnomAD' in db_list:
            warn_user = 'This variant is annotated in ClinVar and gnomAD.'
        elif 'ClinVar' in db_list:
            warn_user = 'This variant is annotated in ClinVar.'
        elif 'gnomAD' in db_list:
            warn_user = 'This variant is annotated in gnomAD.'
            
    return warn_user, output_df
    

def get_different_codes(users_input, prot_entries_df):
    """
    Map the input that gives the user, and retrieve the different information about this gene/protein. If is not found, it returns None values.
    """
    users_input = users_input.upper().strip()

    if users_input in prot_entries_df['Gene_name(HGNC)'].values:
        gene_name = users_input
        entry = prot_entries_df[prot_entries_df['Gene_name(HGNC)'] == users_input]['Uniprot_entry'].values[0]
        entry_name = prot_entries_df[prot_entries_df['Gene_name(HGNC)'] == users_input]['Uniprot_entry_name'].values[0]
        np_code = prot_entries_df[prot_entries_df['Gene_name(HGNC)'] == users_input]['RefSeq_NP_code'].values[0]
        return gene_name, entry, entry_name, np_code

    elif users_input in prot_entries_df['Uniprot_entry'].values:
        entry = users_input
        gene_name = prot_entries_df[prot_entries_df['Uniprot_entry'] == users_input]['Gene_name(HGNC)'].values[0]
        entry_name = prot_entries_df[prot_entries_df['Uniprot_entry'] == users_input]['Uniprot_entry_name'].values[0]
        np_code = prot_entries_df[prot_entries_df['Uniprot_entry'] == users_input]['RefSeq_NP_code'].values[0]
        return gene_name, entry, entry_name, np_code
    
    elif users_input in prot_entries_df['Uniprot_entry_name'].values:
        entry_name = users_input
        gene_name = prot_entries_df[prot_entries_df['Uniprot_entry_name'] == users_input]['Gene_name(HGNC)'].values[0]
        entry = prot_entries_df[prot_entries_df['Uniprot_entry_name'] == users_input]['Uniprot_entry'].values[0]
        np_code = prot_entries_df[prot_entries_df['Uniprot_entry_name'] == users_input]['RefSeq_NP_code'].values[0]
        return gene_name, entry, entry_name, np_code

    elif users_input in prot_entries_df['RefSeq_NP_code'].values:
        np_code = users_input
        gene_name = prot_entries_df[prot_entries_df['RefSeq_NP_code'] == users_input]['Gene_name(HGNC)'].values[0]
        entry = prot_entries_df[prot_entries_df['RefSeq_NP_code'] == users_input]['Uniprot_entry'].values[0]
        entry_name = prot_entries_df[prot_entries_df['RefSeq_NP_code'] == users_input]['Uniprot_entry_name'].values[0]
        return gene_name, entry, entry_name, np_code

    return None, None, None, None
    

def search_for_pairs(uniprot_code, prot_align_df, variants_df, ref_aa, pos, mut_aa, strict_vs_homologous, curated_vs_pfam):
    """
    Seek for homologous variants to the one indicated by the user.
    """
    # Initialize this variable as None
    strict_pairs_df = None
    
    warn_user, eq_pos, align_id, align_path = get_eq_pos_and_check_pfam(uniprot_code, prot_align_df, pos, curated_vs_pfam)     ########### !!!!
    #print(eq_pos, pfam)
    # If it is still all OK, it continues:
    if warn_user == False:
        warn_user, prots_w_same_align = get_proteins_from_align(align_path, uniprot_code)
        #print(prots_w_same_align)
        
        # If it is still all OK, it continues:
        if warn_user == False:
            eq_proteins_d = get_similar_vars_prot_aligns(uniprot_code, prots_w_same_align, eq_pos, pos, align_id, curated_vs_pfam)      ########### !!!!

            #print(eq_proteins_d)
            if eq_proteins_d:
                # Check if we want to search strict pairs
                if strict_vs_homologous == 'strict':
                    #print('strict')

                    strict_pairs_df = get_strict_pairs(variants_df, eq_proteins_d, mut_aa, ref_aa)
                    #print('STRICT')
                    #print(strict_pairs_df.columns)
                    if (strict_pairs_df is not None) and (len(strict_pairs_df) > 0):
                        strict_pairs_df.drop(columns=['Pos_align', 'initial_aa', 'position', 'final_aa', 'Uniprot_entry_name', 'database'], inplace=True)
                        #print(strict_pairs_df.columns)
                        # 'gene', 'clinical_significance', 'prot_change', 'initial_aa', 'position', 'final_aa', 'Uniprot_entry_name', 'Uniprot_link', 'NM_link', 'Pos_align', 'database', 'Database_link'                        
                        strict_pairs_df.rename(columns={'gene':'Gene', 'clinical_significance': 'Clinical significance', 'prot_change': 'Protein change', 'NM_link': 'RefSeq NM code', 'Database_link': 'Database', 'Uniprot_link':'UniProt entry'}, inplace=True)
                        strict_pairs_df['Pair type'] = 'Strict'
                        strict_pairs_df['Alignment type'] = curated_vs_pfam
                        """ if curated_vs_pfam == 'CURATED':
                            strict_pairs_df['Alignment type'] = 'CURATED'
                        else:
                            strict_pairs_df['Alignment type'] = 'PFAM' """

                        return strict_pairs_df, False
                    warn_user = 'There are not strict pairs for this variant.'
                    return None, warn_user

                # Else, we want to look for homologous pairs:
                else:
                    # Get similar aas to the mut_aa (included) based on a positive score in BLOSUM62 matrix
                    similar_mut_aa = get_similar_aas(mut_aa)
                    #print(similar_mut_aa)
                    homolog_pairs_df = get_homologous_pairs(variants_df, eq_proteins_d, similar_mut_aa, ref_aa)
                    #print('HOMOLOGOUS')
                    #print(homolog_pairs_df['prot_change'])
                    if homolog_pairs_df is not None:
                        if not homolog_pairs_df.empty:
                            homolog_pairs_df.drop(columns=['Pos_align', 'initial_aa', 'position', 'final_aa', 'Uniprot_entry_name', 'database'], inplace=True)
                            #print(homolog_pairs_df.columns)
                            homolog_pairs_df.rename(columns={'gene':'Gene', 'clinical_significance': 'Clinical significance', 'prot_change': 'Protein change', 'NM_link': 'RefSeq NM code', 'Database_link': 'Database', 'Uniprot_link':'UniProt entry'}, inplace=True)
                            homolog_pairs_df['Pair type'] = 'Similar'
                            homolog_pairs_df['Alignment type'] = curated_vs_pfam
                            """ if curated_vs_pfam == 'CURATED':
                                homolog_pairs_df['Alignment type'] = 'CURATED'
                            else:
                                homolog_pairs_df['Alignment type'] = 'PFAM' """
                            
                            return homolog_pairs_df, False
                    warn_user = 'There are not strict nor homologous pairs for this variant.'
                    return None, warn_user

    return None, warn_user


def get_similar_aas(mut_aa):
    """
    Get similar amino acids to the one given by the user. Those are the amino acids with a positive score in the BLOSUM62 scoring substitution matrix.
    """
    mut_aa = protein_letters_3to1[mut_aa]
    positive_matches = []
    # Load the BLOSUM62 dictionary from a JSON file
    with open('/var/www/homologous/data/blosum62.json', 'r') as json_file:
        blosum62 = json.load(json_file)
    
    # Iterate through the matrix and find residues with positive scores
    for other_residue, score in blosum62[mut_aa].items():
        if other_residue != mut_aa:
            if score > 0:
                if other_residue not in ['B', 'Z', 'X', 'J', '*']:
                    other_residue = protein_letters_1to3[other_residue]
                    positive_matches.append(other_residue)
    
    return positive_matches


def get_eq_pos_and_check_pfam(uniprot_code, prot_align_df, pos, curated_vs_pfam):
    """
    Get the equivalent position and the alignment path.
    """
    # Retrieve equivalent position in the alignment and the corresponding Pfam alignment
    #print(curated_vs_pfam)
    if 'Pfam' in curated_vs_pfam:
        align_id = prot_align_df[prot_align_df['prot_pos'] == pos]['pfam_align'].values[0]
    else:
        align_id = prot_align_df[prot_align_df['prot_pos'] == pos]['curated_align'].values[0]

    eq_pos = prot_align_df[prot_align_df['prot_pos'] == pos]['eq_pos_align'].values[0]

    # Try to find the SEED alignment of the protein indicated by the user
    if pd.notnull(align_id):
        if 'PF' in align_id and 'SEED' in curated_vs_pfam:
            align_path = f'/var/www/homologous/data/Pfam_SEED_HUMAN_2024_03_20/{align_id}_HUMAN.txt'
        elif 'PF' in align_id and 'FULL' in curated_vs_pfam:
            align_path = f'/var/www/homologous/data/Pfam_FULL_HUMAN_2023_07_03/{align_id}_HUMAN.txt'

        else:
            with open('/var/www/homologous/data/CURATED_aligns_ids.json', 'r') as json_file:
                curated_aligns_ids = json.load(json_file)

            file_pattern = curated_aligns_ids[uniprot_code]
            align_path = f'/var/www/homologous/data/aligns_manually_curated/{file_pattern}'
        
        warn_user = False
        return warn_user, eq_pos, align_id, align_path
    else:
        warn_user = 'We could not look for homologous variants because Pfam has not a seed alignment for this region.'
        return warn_user, None, None, None


def get_proteins_from_align(align_path, uniprot_code):
    """
    Get the proteins that appear with the protein of interest in the alignment. These are the proteins that share at least one domain with it, so we can compare these regions.
    """
    prots_with_same_align = []
    try:
        with open(align_path) as f:
            # Save the content of the file
            file_content = f.readlines()
            # Iterate over the contents of the file
            for line in file_content:
                # Remove the newline characters
                line = line.replace("\n", "")
                prot_uniprot = line.split('/')[0]
                #if prot_uniprot != uniprot_code:
                prots_with_same_align.append(prot_uniprot)
    except FileNotFoundError:
        return 'There is not an alignment for this Pfam.', None
    
    if prots_with_same_align:
        #print(prots_with_same_pfam)
        return False, prots_with_same_align
    else:
        return 'There are not homologous proteins from which extrapolate the annotation of this variant.', None


def get_similar_vars_prot_aligns(uniprot_code, prots_w_same_align, eq_pos, pos, align_id, curated_vs_pfam):
    """
    Construct a dictionary in which we add the possible strict pairs that can be matched with the user's variant.
    The dictionary contains the original positions and ref_aas, for the protein in the same align as the original, based on the equivalent position.
    """
    eq_proteins_d = {}
    #print(pfam, type(align_id))

    # Iterate over the prots that have the same align as the region in which is present the user's variant
    for prot in prots_w_same_align:
        if 'SEED' in curated_vs_pfam:
            prot_align_path = f'/var/www/homologous/data/seed_align_by_prots/{prot}_align_seed_pfam.txt'
        elif 'FULL' in curated_vs_pfam:
            prot_align_path = f'/var/www/homologous/data/full_align_by_prots/{prot}_align_full_pfam.txt'
        else:
            prot_align_path = f'/var/www/homologous/data/curated_aligns_by_prots/{prot}_align_CURATED.txt'

        #print(prot_align_path)
        if os.path.exists(prot_align_path):
            prot_align_df = pd.read_csv(prot_align_path, sep=',')
        else:
            continue
        
        # Filter the df to get only the variants that match the same align
        if 'Pfam' in curated_vs_pfam:
            sub_df_same_align = prot_align_df[prot_align_df['pfam_align'] == align_id]
        else:
            sub_df_same_align = prot_align_df[prot_align_df['curated_align'] == align_id]

        if not sub_df_same_align.empty:
            # Get the original positions and ref_aas for this protein (w same align as the original) based on the equivalent position
            pos_prot_w_same_align = sub_df_same_align[sub_df_same_align['eq_pos_align'] == float(eq_pos)]['prot_pos'].tolist()
            res_prot_w_same_align = sub_df_same_align[sub_df_same_align['eq_pos_align'] == float(eq_pos)]['residue'].tolist()
            
            if pos_prot_w_same_align and res_prot_w_same_align:
                # Save in the dictionary of eq proteins, the pos and residue in the equivalent position of the user's choice
                eq_proteins_d[prot] = {'pos':pos_prot_w_same_align, 'res': [protein_letters_1to3[res.upper()] for res in res_prot_w_same_align]}

    return eq_proteins_d


def get_homologous_pairs(variants_df, eq_proteins_d, similar_mut_aa, ref_aa):
    """
    Search homologous variants based on the previously obtained information:
    - variants_df
    - eq_proteins_d: with the original positions and ref_aas for the protein in the same align as the original
    - similar_mut_aa: similar amino acids to the one given by the user
    - ref_aa
    """
    hom_pairs_general_df = pd.DataFrame()
    variants_df = variants_df[variants_df['Uniprot_entry_name'].isin(list(eq_proteins_d.keys()))]
    for prot in eq_proteins_d:
        positions = eq_proteins_d[prot]['pos']
        residues = eq_proteins_d[prot]['res']
        #print(residues)
        for sim_aa in similar_mut_aa:
            for i in range(len(positions)):
                pos = positions[i]
                res = residues[i]
                
                # Check that the reference amino acid is equal to the one annotated in the alignment
                # In addition, check that the variant is not non-synonymous, just a final check that we only display missense variants
                if res == ref_aa and res != sim_aa:
                    hom_vars_df = variants_df[(variants_df['Uniprot_entry_name'] == prot) &
                        (variants_df['position'] == float(pos)) &
                        (variants_df['initial_aa'] == res) &
                        (variants_df['final_aa'] == sim_aa)
                        ]
                    if hom_vars_df is not None:  # Ensure that strict_pairs is not None or empty
                        hom_pairs_general_df = pd.concat([hom_pairs_general_df, hom_vars_df], ignore_index=True)

    if len(hom_pairs_general_df) > 0:
        return hom_pairs_general_df
    else:
        return None


def get_strict_pairs(variants_df, eq_proteins_d, mut_aa, ref_aa):
    """
    Get strict pairs, which are those variants with the same protein change in homologous proteins.
    """
    strict_pairs_general_df = pd.DataFrame()
    variants_df = variants_df[variants_df['Uniprot_entry_name'].isin(list(eq_proteins_d.keys()))]
    for prot in eq_proteins_d:
        positions = eq_proteins_d[prot]['pos']
        residues = eq_proteins_d[prot]['res']
        for i in range(len(positions)):
            pos = positions[i]
            res = residues[i]

            # Check that the reference amino acid is equal to the one annotated in the alignment
            # In addition, check that the variant is not non-synonymous, just a final check that we only display missense variants
            if res == ref_aa and res != mut_aa:
        
                prot_vars_df = variants_df[(variants_df['Uniprot_entry_name'] == prot) &
                    (variants_df['position'] == float(pos)) &
                    (variants_df['initial_aa'] == res) &
                    (variants_df['final_aa'] == mut_aa)
                    ]

                if len(prot_vars_df) > 0:  # Ensure that strict_pairs is not None or empty
                    strict_pairs_general_df = pd.concat([strict_pairs_general_df, prot_vars_df], ignore_index=True)
    
    #print(strict_pairs_general_df)
    if len(strict_pairs_general_df) > 0:
        return strict_pairs_general_df
    else:
        return None


def format_exist_annot(existing_annot_df, app):
    exist_annot_df_curated = None
    if existing_annot_df is not None:
        prot_entries = os.path.join(app.root_path, 'data', 'prot_entry_types_mapped.txt')
        prot_entries_df = pd.read_csv(prot_entries, sep='\t')
        exist_annot_df_curated = checking_tools.format_vars_df(existing_annot_df, prot_entries_df)
        #print(exist_annot_df_curated.columns)
        # Order columns of final df
        custom_order = ['Gene', 'Clinical significance', 'Protein change', 'UniProt entry', 'RefSeq NM code', 'Database']
        exist_annot_df_curated = exist_annot_df_curated[custom_order]
        return exist_annot_df_curated
    return None
    