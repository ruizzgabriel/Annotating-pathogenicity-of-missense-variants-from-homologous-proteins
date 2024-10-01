# Quantificate pairs filtered by sequence identity %

import pandas as pd
import numpy as np
from utilities.input_files import obtain_pair_files


def filter_by_identity(pairs_df, identity_per):
    """
    Filter file of pairs by the desired identity percentage. Returns the same pairs DataFrame with pairs that have a higher identity than the threshold set.

    Args:
        - pairs_df(DataFrame): containing pairs of variants N-N, P-P or N-P
        - identity_per(int): percentage of identity desired to filter
    """
    if type(identity_per) != float:
        identity_per = float(identity_per)
    
    if len(pairs_df) > 0:
        pairs_df = pairs_df[pairs_df['seq_identity(%)'] > identity_per]

    return len(pairs_df)


def read_pair_file(pairs_file):

    pairs_file = pd.read_csv(pairs_file, sep="\t", header=0)
    if len(pairs_file) > 0:
        #pairs_file['seq_identity(%)'] = pairs_file['seq_identity(%)'].str.rstrip('%').astype(float)
        pairs_file['seq_identity(%)'] = pairs_file['seq_identity(%)'].astype(str).str.rstrip('%').astype(float)
    
    return pairs_file

def extract_num_pfams(NN, PP, NP, identity_per=None):
    """
    Extract the total number of Pfam codes that formed pairs in this filter.
    """
    pairs_files = [NN, PP, NP]
    pfams_list = []
    for pairs_df in pairs_files:
        if len(pairs_df) > 0:
            if identity_per != None:
                pairs_df = pairs_df[pairs_df['seq_identity(%)'] > identity_per]
            pfams_list.append(pairs_df['Pfam_code'].unique().tolist())

    #print(pfams_list)
    all_pfam_list = [pfam for sublist in pfams_list for pfam in sublist]
    num_pfams = len(set(all_pfam_list))
    print(num_pfams)
    return num_pfams


def search_variants(variants_l, variant_a, variant_b=None):
    v1 = False
    if variant_b:
        v2 = False
    else:
        v2 = None

    for c_variant in variants_l:
        if variant_a == c_variant[0]:
            c_variant[1] += 1
            v1 = True
        if variant_b is not None and variant_b == c_variant[0]:
            c_variant[1] += 1
            v2 = True
        
        if v1 == v2 and v1 == True:
            break
        if v2 is None and v1 == True:
            break

    if v1 == False:
        variants_l.append([variant_a, 1])
    if v2 == False:
        variants_l.append([variant_b, 1])
    return variants_l


def extract_patho_variants(PP, NP, pair_type, identity_per=None):
    """
    Pfam_code	initial_aa	pos_align	final_aa	gene_name	gene_name2	pos_prot_change	pos_prot_change2	db1	db2	uniprot1	
    uniprot2	score_substmat	seq_identity(%)

    Pfam_code	initial_aa	pos_align	final_aa	final_aa2	gene_name	gene_name2	pos_prot_change	pos_prot_change2	
    db1	db2	uniprot1	uniprot2	score_substmat	seq_identity(%)
    """
    if identity_per:
        PP = PP[PP['seq_identity(%)'] > identity_per]
        NP = NP[NP['seq_identity(%)'] > identity_per]

    patho_variants = []
    for pairs_file in [PP, NP]:
        for index, row in pairs_file.iterrows():
            # save info from variant 1 and for variants 2 without loosing the reference of which one is
            pfam_code = row['Pfam_code']
            i_aa = row['initial_aa']
            pos_align = row['pos_align']
            f_aa = row['final_aa']
            gene_name = row['gene_name']
            gene_name2 = row['gene_name2']
            pos_prot_change = row['pos_prot_change']
            pos_prot_change2 = row['pos_prot_change2']
            db1 = row['db1']
            db2 = row['db2']
            uniprot1 = row['uniprot1']
            uniprot2 = row['uniprot2']
            #score_substmat = row['score_substmat']
            #seq_id = row['seq_identity(%)']

            if pair_type == 'blosum':
                f_aa2 = row['final_aa2']
                
                variant1 = [pfam_code, i_aa, pos_align, f_aa, gene_name, pos_prot_change, db1, uniprot1]
                variant2 = [pfam_code, i_aa, pos_align, f_aa2, gene_name2, pos_prot_change2, db2, uniprot2]
            else:
                variant1 = [pfam_code, i_aa, pos_align, f_aa, gene_name, pos_prot_change, db1, uniprot1]
                variant2 = [pfam_code, i_aa, pos_align, f_aa, gene_name2, pos_prot_change2, db2, uniprot2]
            
            if db1 == db2 and db1 == 'ClinVar':
                patho_variants = search_variants(variants_l=patho_variants, variant_a=variant1, variant_b=variant2)
            if db1 == 'ClinVar':
                patho_variants = search_variants(variants_l=patho_variants, variant_a=variant1)
            elif db2 == 'ClinVar':
                patho_variants = search_variants(variants_l=patho_variants, variant_a=variant2)

    return len(patho_variants)


def extract_non_patho_variants(NN, NP, pair_type, identity_per=None):
    if identity_per:
        NN = NN[NN['seq_identity(%)'] > identity_per]
        NP = NP[NP['seq_identity(%)'] > identity_per]

    non_patho_variants = []
    for pairs_file in [NN, NP]:
        for index, row in pairs_file.iterrows():
            # save info from variant 1 and for variants 2 without loosing the reference of which one is
            pfam_code = row['Pfam_code']
            i_aa = row['initial_aa']
            pos_align = row['pos_align']
            f_aa = row['final_aa']
            gene_name = row['gene_name']
            gene_name2 = row['gene_name2']
            pos_prot_change = row['pos_prot_change']
            pos_prot_change2 = row['pos_prot_change2']
            db1 = row['db1']
            db2 = row['db2']
            uniprot1 = row['uniprot1']
            uniprot2 = row['uniprot2']
            #score_substmat = row['score_substmat']
            #seq_id = row['seq_identity(%)']

            if pair_type == 'blosum':
                f_aa2 = row['final_aa2']
                
                variant1 = [pfam_code, i_aa, pos_align, f_aa, gene_name, pos_prot_change, db1, uniprot1]
                variant2 = [pfam_code, i_aa, pos_align, f_aa2, gene_name2, pos_prot_change2, db2, uniprot2]
            else:
                variant1 = [pfam_code, i_aa, pos_align, f_aa, gene_name, pos_prot_change, db1, uniprot1]
                variant2 = [pfam_code, i_aa, pos_align, f_aa, gene_name2, pos_prot_change2, db2, uniprot2]
            
            if db1 == db2 and db1 == 'gnomAD':
                non_patho_variants = search_variants(variants_l=non_patho_variants, variant_a=variant1, variant_b=variant2)
            if db1 == 'gnomAD':
                non_patho_variants = search_variants(variants_l=non_patho_variants, variant_a=variant1)
            elif db2 == 'gnomAD':
                non_patho_variants = search_variants(variants_l=non_patho_variants, variant_a=variant2)

    return len(non_patho_variants)


normal_or_blosum = ['normal', 'blosum']
freqs = ['no_freq', 'freq-8', 'freq-7', 'freq-6', 'freq-5', 'freq-4', 'freq-3', 'freq-2']

results_table_25 = [['BLOSUM_AA', 'NORMAL/BLOSUM', 'AF_FILTER', 'N-N_pairs', 'P-P_pairs', 'N-P_pairs', 'TOTAL_PAIRS']]
results_table_30 = [['BLOSUM_AA', 'NORMAL/BLOSUM', 'AF_FILTER', 'Num_Pfams', 'Patho_variants', 'Non-patho_variants', 'N-N_pairs', 'P-P_pairs', 'N-P_pairs', 'TOTAL_PAIRS']]
results_table_40 = [['BLOSUM_AA', 'NORMAL/BLOSUM', 'AF_FILTER', 'N-N_pairs', 'P-P_pairs', 'N-P_pairs', 'TOTAL_PAIRS']]

full_vs_seed = input("Write the desired alignment type. (Mix condition prioritizes seed alignment and uses full ones when seed is not found)\nYou can choose between:\n\tfull\n\tseed\n\tseed+full(mix)\nPlease, write the version as indicated: ")

for pair_type in normal_or_blosum:
    if pair_type == 'normal':
        for AF_freq in freqs:
            # Get dfs one for each pair type file
            N_P_file, N_N_file, P_P_file = obtain_pair_files(pair_type, AF_freq, full_vs_seed)
            print(pair_type, AF_freq)
            
            N_N_file = read_pair_file(N_N_file)
            P_P_file = read_pair_file(P_P_file)
            N_P_file = read_pair_file(N_P_file)     


            # Filter files of pairs by the indicated identity threshold and save its lengths
            NN_25 = filter_by_identity(N_N_file, 25)
            NN_30 = filter_by_identity(N_N_file, 30)
            NN_40 = filter_by_identity(N_N_file, 40)

            PP_25 = filter_by_identity(P_P_file, 25)
            PP_30 = filter_by_identity(P_P_file, 30)
            PP_40 = filter_by_identity(P_P_file, 40)

            NP_25 = filter_by_identity(N_P_file, 25)
            NP_30 = filter_by_identity(N_P_file, 30)
            NP_40 = filter_by_identity(N_P_file, 40)

            num_pfams_30 = extract_num_pfams(N_N_file, P_P_file, N_P_file, 30)
            #num_patho_vars = '-'
            num_patho_vars_30 = 'not_calculated'
            num_non_patho_vars_30 = 'not_calculated'
            if AF_freq == 'freq-6':
                #num_patho_vars = extract_patho_variants(P_P_file, N_P_file, pair_type)
                num_patho_vars_30 = extract_patho_variants(P_P_file, N_P_file, pair_type, 30)
                num_non_patho_vars_30 = extract_non_patho_variants(N_N_file, N_P_file, pair_type, 30)


            # Prepare row of results to append to the results df
            results_row_25 = ['-', pair_type, AF_freq, NN_25, PP_25, NP_25, (NN_25 + PP_25 + NP_25)]
            results_table_25.append(results_row_25)

            results_row_30 = ['-', pair_type, AF_freq, num_pfams_30, num_patho_vars_30, num_non_patho_vars_30, NN_30, PP_30, NP_30, (NN_30 + PP_30 + NP_30)]
            results_table_30.append(results_row_30)

            results_row_40 = ['-', pair_type, AF_freq, NN_40, PP_40, NP_40, (NN_40 + PP_40 + NP_40)]
            results_table_40.append(results_row_40)
            
    else:
        aa_blosum = ['final_aa', 'both_aa', 'initial_aa']
        for aa_2_compute in aa_blosum:
            for AF_freq in freqs:
                print(pair_type, AF_freq, aa_2_compute)
                # Get dfs one for each pair type file
                N_P_file, N_N_file, P_P_file = obtain_pair_files(pair_type, AF_freq, full_vs_seed, aa_2_compute)

                N_N_file = read_pair_file(N_N_file)
                P_P_file = read_pair_file(P_P_file)
                N_P_file = read_pair_file(N_P_file)

                # Filter files of pairs by the indicated identity threshold
                NN_25 = filter_by_identity(N_N_file, 25)
                NN_30 = filter_by_identity(N_N_file, 30)
                NN_40 = filter_by_identity(N_N_file, 40)

                PP_25 = filter_by_identity(P_P_file, 25)
                PP_30 = filter_by_identity(P_P_file, 30)
                PP_40 = filter_by_identity(P_P_file, 40)

                NP_25 = filter_by_identity(N_P_file, 25)
                NP_30 = filter_by_identity(N_P_file, 30)
                NP_40 = filter_by_identity(N_P_file, 40)

                num_pfams_30 = extract_num_pfams(N_N_file, P_P_file, N_P_file, 30)
                num_patho_vars_30 = 'not_calculated'
                num_non_patho_vars_30 = 'not_calculated'
                if aa_2_compute == 'final_aa' and AF_freq == 'freq-6':
                    #extract_patho_variants(P_P_file, N_P_file, pair_type)
                    num_patho_vars_30 = extract_patho_variants(P_P_file, N_P_file, pair_type, 30)
                    num_non_patho_vars_30 = extract_non_patho_variants(N_N_file, N_P_file, pair_type, 30)

                # Prepare row of results to append to the results df
                results_row_25 = [aa_2_compute, pair_type, AF_freq, NN_25, PP_25, NP_25, (NN_25 + PP_25 + NP_25)]
                results_table_25.append(results_row_25)

                results_row_30 = [aa_2_compute, pair_type, AF_freq, num_pfams_30, num_patho_vars_30, num_non_patho_vars_30, NN_30, PP_30, NP_30, (NN_30 + PP_30 + NP_30)]
                results_table_30.append(results_row_30)

                results_row_40 = [aa_2_compute, pair_type, AF_freq, NN_40, PP_40, NP_40, (NN_40 + PP_40 + NP_40)]
                results_table_40.append(results_row_40)


# Save the results dfs in txt files
results_table_25 = '\n'.join('\t'.join(str(item) for item in inner_list) for inner_list in results_table_25)
results_table_30 = '\n'.join('\t'.join(str(item) for item in inner_list) for inner_list in results_table_30)
results_table_40 = '\n'.join('\t'.join(str(item) for item in inner_list) for inner_list in results_table_40)

output_paths = [[f'./data/quanti_pairs_25_identity_{full_vs_seed}.txt', results_table_25],
                [f'./data/quanti_pairs_30_identity_{full_vs_seed}.txt', results_table_30],
                [f'./data/quanti_pairs_40_identity_{full_vs_seed}.txt', results_table_40]]

for path, results_table in output_paths:
    with open(path, 'w') as outfile:
        outfile.write(results_table)
