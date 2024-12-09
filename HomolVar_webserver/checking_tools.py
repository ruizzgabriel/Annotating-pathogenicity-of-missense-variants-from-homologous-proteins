from Bio.Data.IUPACData import protein_letters_3to1, protein_letters_1to3
from markupsafe import Markup
import os
import pandas as pd

def check_pos(pos):
    """
    Check the position indicated by the user.
    """
    # Check if pos is a Markup object
    if isinstance(pos, Markup):
        pos = pos.striptags()  # Strip HTML tags if it's a Markup object

    # Ensure pos is not empty and convert it to an integer
    if pos and pos != 'None':
        try:
            pos = int(pos)  # Convert to integer
        except ValueError:
            pos = None  # Handle cases where conversion fails
    else:
        pos = None  # Handle None or empty cases
    return pos


def get_aa(aa):
    """
    Check that the aa given is correct. If it is given in 1-letter format, we save it like 3-letter format code. If it does not exist, returns None value.
    """
    aa_dict = {
        'Alanine': ('Ala', 'A'),
        'Arginine': ('Arg', 'R'),
        'Asparagine': ('Asn', 'N'),
        'Aspartic acid': ('Asp', 'D'),
        'Cysteine': ('Cys', 'C'),
        'Glutamine': ('Gln', 'Q'),
        'Glutamic acid': ('Glu', 'E'),
        'Glycine': ('Gly', 'G'),
        'Histidine': ('His', 'H'),
        'Isoleucine': ('Ile', 'I'),
        'Leucine': ('Leu', 'L'),
        'Lysine': ('Lys', 'K'),
        'Methionine': ('Met', 'M'),
        'Phenylalanine': ('Phe', 'F'),
        'Proline': ('Pro', 'P'),
        'Serine': ('Ser', 'S'),
        'Threonine': ('Thr', 'T'),
        'Tryptophan': ('Trp', 'W'),
        'Tyrosine': ('Tyr', 'Y'),
        'Valine': ('Val', 'V')
    }
    if aa.capitalize() in protein_letters_3to1:
        return aa.capitalize()
    elif aa.upper() in protein_letters_1to3:
        return protein_letters_1to3[aa.upper()]
    elif aa.capitalize() in aa_dict:
        return aa_dict[aa.capitalize()][0]
    elif ('*' in aa) or ('Ter' in aa):
        return 'Ter'
    elif ('fs' in aa):
        return 'fs'
    else:
        return None


def check_protein_type(uniprot_code, entry, gene, prot_var, prot_entries_df):
    """
    Here, we'll check if the protein given by the user is Autosomal Dominant, Autosomal Recessive or is a non-patho protein (non-vital protein or involved in Complex Disease)
    """
    
    # Convert 'omim' to numeric, forcing errors to NaN
    prot_entries_df['omim'] = pd.to_numeric(prot_entries_df['omim'], errors='coerce')
    # Filter to keep only rows with numeric values in 'omim'
    prot_entries_df = prot_entries_df[prot_entries_df['omim'].notna()]    

    # Get the unique protein names from those ones that have OMIM id
    omim_prots = prot_entries_df['Uniprot_entry_name'].dropna()
    omim_prots = list(set(omim_prots))

    if uniprot_code in omim_prots:
        #if mim_num != False:
        mim_num = prot_entries_df[prot_entries_df['Uniprot_entry_name'] == uniprot_code]['omim'].iloc[0]
        #print(mim_num)
        mim_num = str(mim_num).split('.')[0]

        form_data = {
            'output_user': f'This protein is involved in a mendelian disorder:  OMIM <a href="https://www.omim.org/entry/{mim_num}" target="_blank">{mim_num}</a>',   # https://www.omim.org/search?index=entry&start=1&search=AARS1
            'omim': mim_num,
            'patho_prot': True,
            'protein':uniprot_code,
            'uniprot_entry': entry,
            'prot_var': prot_var,
            'gene_name' : gene
        }
        
        #print(form_data)
    
    else:     # Complex disease or non-vital proteins  -> non-patho proteins
        form_data = {
            'output_user' : 'This protein is not considered as autosomal dominant or recessive in OMIM, so we consider this protein as not disease-causing. Maybe its function is not vital for the organism or maybe this protein is involved in complex disease.',
            'omim': False,
            'patho_prot': False,
            'protein': uniprot_code,
            'uniprot_entry': entry,
            'prot_var': prot_var,
            'gene_name' : gene
        }
    return form_data


def check_users_input(uniprot_code, ref_aa, pos):
    """
    This function checks if we have an alignment for the protein chosen and if the reference amino acid given by the user is the same as in the UniProt's sequence.
    We check if we have a curated alignment and a Pfam alignment for the protein given
    If something is not correct, it will warn the user.
    """
    curated_align_df = pd.DataFrame()
    pfam_align_df = pd.DataFrame()
    
    # Search the original protein align file
    pfam_align_path = f'/var/www/homologous/data/seed_align_by_prots/{uniprot_code}_align_seed_pfam.txt'
    pfam_FULL_align_path = f'/var/www/homologous/data/full_align_by_prots/{uniprot_code}_align_full_pfam.txt'
    if os.path.exists(pfam_align_path):
        pfam_align_df = pd.read_csv(pfam_align_path, sep=',')
        align_type = 'Pfam SEED'
    elif os.path.exists(pfam_FULL_align_path):
        align_type = 'Pfam FULL'
        pfam_align_df = pd.read_csv(pfam_FULL_align_path, sep=',')

    curated_align_path = f'/var/www/homologous/data/curated_aligns_by_prots/{uniprot_code}_align_CURATED.txt'
    if os.path.exists(curated_align_path):
        align_type = 'CURATED'
        curated_align_df = pd.read_csv(curated_align_path, sep=',')

    # If we have alignment curated and from Pfam
    if (len(curated_align_df) > 0) and (len(pfam_align_df) > 0):
        warn_user_curated, curated_align_df = check_users_pos_in_align(curated_align_df, pos, ref_aa)
        warn_user_pfam, pfam_align_df = check_users_pos_in_align(pfam_align_df, pos, ref_aa)

        if warn_user_pfam and warn_user_curated:
            # It may be the same error in both cases so we return one of them. We return the curated_align_df to plot it, because it will be much better than the alignment from Pfam
            return warn_user_curated, curated_align_df, None, align_type
        elif warn_user_curated:
            return False, None, pfam_align_df, align_type
        elif warn_user_pfam:
            return False, curated_align_df, None, align_type
        # No warn_users found so we return both dfs to seek for pairs and plot the curated align
        return False, curated_align_df, pfam_align_df, align_type
    
    # If we only have curated alignment
    elif len(curated_align_df) > 0:
        warn_user, curated_align_df = check_users_pos_in_align(curated_align_df, pos, ref_aa)
        if warn_user:
            # If warn_user means that the position or the ref_aa is not correct in the alignment, we will notify the user but we return the align df anyway to perform the alignment plot
            return warn_user, curated_align_df, None, align_type
        else:
            # All OK so we will seek for pairs in curated align
            return False, curated_align_df, None, align_type
    
    # If we only have Pfam alignment
    elif len(pfam_align_df) > 0:
        warn_user, pfam_align_df = check_users_pos_in_align(pfam_align_df, pos, ref_aa)
        if warn_user:
            # If warn_user means that the position or the ref_aa is not correct in the alignment, we will notify the user but we return the align df anyway to perform the alignment plot
            return warn_user, None, pfam_align_df, align_type
        else: 
            # All OK so we will seek for pairs in Pfam align
            return False, None, pfam_align_df, align_type

    # If we don't have any alignment, return the user that the chosen protein is not in our database
    else:
        warn_user = 'no_alignment'
        # We will not be able to do the plot so we return dfs=None
        return warn_user, None, None, None


def check_users_pos_in_align(align_df, pos, ref_aa):
    """
    Check if the reference amino acid given by the user matches with the one in the alignment in the same position.
    """
    # Check that the position given is in the protein length
    if pos in align_df['prot_pos']:

        # Check that the reference aa given by the user corresponds to the correct aa in the given position in the UniProt sequence.
        res = align_df[align_df['prot_pos'] == pos]['residue'].values[0]
        
        # If it is different, set an error to show
        if protein_letters_1to3[res] != ref_aa:
            warn_user = f'The provided residue ({ref_aa}) does not match the residue present in the UniProt sequence at the given position ({protein_letters_1to3[res]}).\nPlease verify that your variant is correctly indicated.'
            # We return also the df for performing the plot
            return warn_user, align_df
        
        # If it is correct, retrieve the equivalent position and the Pfam alignment in which is found
        else:
            warn_user = False
            return warn_user, align_df
    else:
        warn_user = f'The position indicated is out of the protein length.'
        return warn_user, align_df
    

def check_consequence_of_pairs(pairs_df):
    """
    Check if the pairs obtained are damaging, non-damaging or differ between them leading to uncertain pathogenesis.
    """
    #print(pairs_df.columns)
    for i, row in pairs_df.iterrows():
        if type(row['Clinical significance']) == float:
            if 'gnomAD' in row['Database']:
                pairs_df.loc[i, 'Clinical significance'] = 'Neutral'
            if 'ClinVar' in row['Database']:
                pairs_df.loc[i, 'Clinical significance'] = 'Pathogenic'

    # Check if there are results for CURATED alignment
    curated_annot = None
    all_annot = None
    curated_pairs_df = pairs_df[pairs_df['Alignment type'] == 'CURATED']

    patho_l = ['Likely pathogenic', 'Likely risk allele', 'Pathogenic', 'Pathogenic/Likely pathogenic', 'Pathogenic/Likely risk allele', 'Pathogenic; risk factor']
    neutral_l = ['Neutral', 'Benign', 'Benign/Likely benign', 'Likely benign']

    # If we have curated align:
    if len(curated_pairs_df) > 0:
        annot_list = curated_pairs_df['Clinical significance'].tolist()
        if (any(item in neutral_l for item in annot_list)) and (any(item in patho_l for item in annot_list)):
            curated_annot = 'Uncertain consequence'
        elif any(item in neutral_l for item in annot_list):
            curated_annot = 'NON-DAMAGING'
        elif any(item in patho_l for item in annot_list):
            curated_annot = 'DAMAGING'
        
    annot_list = pairs_df['Clinical significance'].tolist()
    #print(pairs_df.columns)
    #print(annot_list)
    if (any(item in neutral_l for item in annot_list)) and (any(item in patho_l for item in annot_list)):
        all_annot = 'Uncertain consequence'
    elif any(item in neutral_l for item in annot_list):
        all_annot = 'NON-DAMAGING'
    elif any(item in patho_l for item in annot_list):
        all_annot = 'DAMAGING'
    
    # If there is a curated alignment
    if curated_annot:
        # If pathogenicity annotation is the same as the obtained in the Pfam aligns:
        if curated_annot == all_annot:
            # Return all pairs
            return curated_annot, pairs_df
        # If not, we trust the curated annotation and only show the pairs obtained from the curated align
        else:
            return curated_annot, curated_pairs_df
    # If no curated alignment is available, we return all pairs
    else:
        return all_annot, pairs_df


def format_vars_df(vars_df, prot_entries_df):
    """
    Format variants found previously annotated in ClinVar or gnomAD. Add links to some values, etc.
    """
    #print(vars_df.columns)
    for i, row in vars_df.iterrows():

        clin_sig = row['Clinical significance']
        if pd.isna(clin_sig) and row['Database'] == 'gnomAD':
            vars_df.loc[i, 'Clinical significance'] = 'Neutral'
        elif pd.isna(clin_sig) and row['Database'] == 'ClinVar':
            vars_df.loc[i, 'Clinical significance'] = 'Pathogenic'

        nm_code = row['RefSeq NM code']
        nm_code = nm_code.replace("]", '').replace("[", '').replace("'", '').split(',')
        nm_l = []
        for nm in nm_code:
            if '.' in nm:
                nm = nm.split('.')[0]
            nm = f"<a href='https://www.ncbi.nlm.nih.gov/nuccore/{nm}' target='_blank'>{nm}</a>"
            nm_l.append(nm)
        #print(nm_l)
        # Join the formatted codes and update the DataFrame
        vars_df.loc[i, 'RefSeq NM code'] = ', '.join(nm_l)

        entry = prot_entries_df[prot_entries_df['Gene_name(HGNC)'] == row['Gene']]['Uniprot_entry'].values[0]
        vars_df.loc[i, 'UniProt entry'] = f"<a href='https://www.uniprot.org/uniprotkb/{entry}/entry/' target='_blank'>{entry}</a>"

        # Link gnomAD
        #print(vars_df.columns)
        rsid = row['rsid']
        db = row['Database']

        if not pd.isna(rsid) and db == 'gnomAD':
            vars_df.loc[i, 'Database'] = f"<a href='https://gnomad.broadinstitute.org/variant/{rsid}' target='_blank'>gnomAD</a>"
        elif db == 'gnomAD':
            vars_df.loc[i, 'Database'] = f"<a href='https://gnomad.broadinstitute.org/' target='_blank'>gnomAD</a>"

        elif 'VariationID' in row:
            var_id = row['VariationID']
            if var_id is not None and db == 'ClinVar':
                if '.0' in var_id:
                    var_id = var_id.replace('.0', '')
                var_id = int(var_id)
                vars_df.loc[i, 'Database'] = f"<a href='https://www.ncbi.nlm.nih.gov/clinvar/variation/{var_id}/' target='_blank'>ClinVar</a>"
        elif db == 'ClinVar':
            vars_df.loc[i, 'Database'] = f"<a href='https://www.ncbi.nlm.nih.gov/clinvar/' target='_blank'>ClinVar</a>"

        if db != 'ClinVar' and db != 'gnomAD':
            if db == 'CFERV':
                vars_df.loc[i, 'Database'] = f"<a href='https://webapp2.pharm.emory.edu/cferv/' target='_blank'>CFERV</a>"
            elif db == 'CFERV':
                vars_df.loc[i, 'Database'] = f"<a href='https://alf06.uab.es/grindb/grindbdata' target='_blank'>BCN-GRIN</a>"
            elif db == 'GRIN-Leipzig':
                vars_df.loc[i, 'Database'] = f"<a href='https://grin-database.de/gen_table' target='_blank'>GRIN-Leipzig</a>"
            elif db == 'LOVD':
                vars_df.loc[i, 'Database'] = f"<a href='https://www.lovd.nl/' target='_blank'>LOVD</a>"
            elif db == 'UniProt':
                vars_df.loc[i, 'Database'] = f"<a href='https://www.uniprot.org/uniprotkb/{entry}/entry' target='_blank'> UniProt</a>"

    return vars_df


def format_pairs_df(pairs_df, prot_var, entry):

    # Show first results for strict pairs
    pairs_df = pairs_df.sort_values(by='Pair type', ascending=False)

    # Check that we did not found as homologous, the same variant the user gave us
    pairs_df = pairs_df[pairs_df['Protein change'] != prot_var]

    for i, row in pairs_df.iterrows():
        if row['Alignment type'] == 'Pfam FULL':
            pairs_df.loc[i, 'Alignment type'] = f"<a href='https://www.ebi.ac.uk/interpro/protein/UniProt/{entry}/entry/pfam/#table' target='_blank'> Pfam FULL</a>"
        elif row['Alignment type'] == 'Pfam SEED':
            pairs_df.loc[i, 'Alignment type'] = f"<a href='https://www.ebi.ac.uk/interpro/protein/UniProt/{entry}/entry/pfam/#table' target='_blank'> Pfam SEED</a>"
    
    return pairs_df