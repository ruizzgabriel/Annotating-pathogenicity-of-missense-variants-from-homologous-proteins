# Search envolut distances between variants in a pair --> calculate matrix score and seq identity from PFAM alignments
# CHANGES ARE GOING TO BE DONE IN SITU -> IN THE SAME FILES PRESENT IN 6_Pairs_extraction TO AVOID DUPLICATED DATA IN TWO FOLDERS

import os, pandas as pd
from colorama import Fore, Style
from utilities.pairs_alignment_features import (
    score_seqs,
    seq_identity,
    extract_sequences,
)

from utilities.input_files import user_chooses_condition, obtain_pair_files, save_file

def float_to_int(value):
    if value == '':
        return ''
    try:
        return int(value)
    except ValueError:
        return ''


def print_to_user(aa_2_compute, freq, full_vs_seed):
    if aa_2_compute != None:
        print(Fore.YELLOW + f"\n\nWe start computing HOMOLOGOUS pairs.\n\t -> {aa_2_compute}")
    else:
        print(Fore.YELLOW + f"\n\nWe'll compute STRICT pairs.\n")
    print(Style.RESET_ALL)

    print(Fore.GREEN + f"\n\nWe start extracting pairs with:\n\tFREQUENCE: {freq.split('_')[-1]}")
    print(Fore.MAGENTA + f"\n\tALIGNMENT: {full_vs_seed}")
    print(Style.RESET_ALL)


def extract_alignment(group, full_vs_seed):
    pfam = group["Pfam_code"].iloc[0]

    if full_vs_seed == 'full':
        pfam_file = (f"../1_Data_collection/data/Pfam_HUMAN_InterPro_2023-07-03/{pfam}_HUMAN.txt")
    elif full_vs_seed == 'seed':
        pfam_file = (f"../1_Data_collection/data/Pfam_SEED_HUMAN_2024_03_20/{pfam}_HUMAN.txt")

    # Open the corresponding Pfam alignment file
    try:
        with open(pfam_file) as f:
            # Save the content of the file
            file_content = f.readlines()
        return file_content
    except FileNotFoundError:
        return ''
    


def add_identity_n_score(pairs_df, group, pfam_align):
    for index, pair in group.iterrows():
        seq1, seq2 = extract_sequences(pair, pfam_align)
        """
        print(pfam)
        print(seq1, '\n\n', seq2)
        print('############### ############### ############### ############### ###############')
        """
        score_substmat = score_seqs(seq1, seq2, 1, "BLOSUM62")
        seq_id = seq_identity(seq1, seq2)

        if score_substmat != '':
            pairs_df.at[index, "score_substmat"] = score_substmat

        #else:
            #print(score_substmat)
            #print(f'{seq1}\n{seq2}\n\n')

        if seq_id:
            pairs_df.at[index, "seq_identity(%)"] = seq_id
    
    # pairs["score_substmat"] = pairs["score_substmat"].astype(int) # Does not work with empty strings that are also present
    if 'score_substmat' in pairs_df.columns:
        pairs_df.loc[:, "score_substmat"] = pairs_df["score_substmat"].apply(float_to_int)

    return pairs_df


def add_mix_stats(N_P_file, N_N_file, P_P_file):

    pair_files = [N_P_file, N_N_file, P_P_file]
    alignments = ['seed', 'full']

    for file in pair_files:

        pairs = pd.read_csv(file, sep="\t", header=0)

        # If seq_identity column exists, remove it to recalculate it from 0
        if 'seq_identity' in pairs.columns:
            pairs.drop('seq_identity', axis=1, inplace=True)

        print(f'LENGTH PAIRS:  {len(pairs)} pairs')

        if len(pairs) > 0:
            # print(pairs)
            pairs['score_substmat'] = ''
            pairs['seq_identity(%)'] = ''

            pairs_grouped = pairs.groupby("Pfam_code")
            for group_name, group in pairs_grouped:
                pfam_align = extract_alignment(group, full_vs_seed='seed')
                if pfam_align:
                    pairs = add_identity_n_score(pairs, group, pfam_align)

        seed_rows = pairs[(pairs['score_substmat'] != '') | (pairs['seq_identity(%)'] != '')]
        print(f'LENGTH PAIRS SEED:  {len(seed_rows)}')

        #unfilled_pairs = pairs.loc[(pairs['score_substmat'].isnull() | (pairs['score_substmat'] == '')) & (pairs['seq_identity(%)'].isnull() | (pairs['seq_identity(%)'] == ''))]
        unfilled_pairs = pairs[(pairs['score_substmat'] == '') & (pairs['seq_identity(%)'] == '')]
        print(f'LENGTH PAIRS TO TRY WITH FULL:  {len(unfilled_pairs)}')

        if len(unfilled_pairs) > 0:
            pairs_grouped = unfilled_pairs.groupby("Pfam_code")
            for group_name, group in pairs_grouped:
                pfam_align = extract_alignment(group, full_vs_seed='full')
                full_rows = add_identity_n_score(unfilled_pairs, group, pfam_align)

            # Concatenate the filled and unfilled rows back into a single DataFrame
            if len(seed_rows) > 0:
                pairs = pd.concat([seed_rows, full_rows])
            else:
                pairs = full_rows.copy()

            print(f'LENGTH PAIRS after concat:  {len(pairs)}')
            print(f"LENGTH PAIRS without identity nor score:  {len(pairs[(pairs['score_substmat'] == '') & (pairs['seq_identity(%)'] == '')])}")
            print("FILE: {} processed.".format(file.split('\\')[-1]))
            print('\n')

            # Set data path
            save_file(pairs, file)
        else:
            continue


def add_stats(N_P_file, N_N_file, P_P_file, full_vs_seed):

    pair_files = [N_P_file, N_N_file, P_P_file]
    #pair_files = [N_N_file]

    for file in pair_files:

        pairs = pd.read_csv(file, sep="\t", header=0)

        if 'seq_identity' in pairs.columns:
            pairs.drop('seq_identity', axis=1, inplace=True)

        if len(pairs) > 0:
            # print(pairs)

            pairs_grouped = pairs.groupby("Pfam_code")
            for group_name, group in pairs_grouped:
                pfam_align = extract_alignment(group, full_vs_seed)
                pairs = add_identity_n_score(pairs, group, pfam_align)

            ########## REMOVE ##########
            if (len(pairs[pairs["uniprot1"].isnull() | (pairs["uniprot1"] == "")]) != 0) or (len(pairs[pairs["uniprot2"].isnull() | (pairs["uniprot2"] == "")]) != 0):
                print(pair_type, file)
                print(f' - UNIPROT CODE 1, empty rows: {len(pairs[pairs["uniprot1"].isnull() | (pairs["uniprot1"] == "")])}')
                print(f' - UNIPROT CODE 2, empty rows: {len(pairs[pairs["uniprot2"].isnull() | (pairs["uniprot2"] == "")])}')
            if (len(pairs[pairs["score_substmat"].isnull() | (pairs["score_substmat"] == "")]) !=0) or (len(pairs[pairs["seq_identity(%)"].isnull() | (pairs["seq_identity(%)"] == "")]) !=0):
                print(f'- SCORES MATRIX, empty rows:  {len(pairs[pairs["score_substmat"].isnull() | (pairs["score_substmat"] == "")])}')
                print(f'- SEQUENCE IDENTITY, empty rows:  {len(pairs[pairs["seq_identity(%)"].isnull() | (pairs["seq_identity(%)"] == "")])}')
            ############################
            
            print("FILE: {} processed.".format(file.split('\\')[-1]))
            # print(pairs)
            #output_name = file.split("\\")[-1]
            # Set data path
            save_file(pairs, file)
        else:
            continue



######################## SELECT CONDITION ########################

# Ask the user if want to treat normal pairs or wants to compute BLOSUM62; the allele frequence threshold; and the amino acid to compute with BLOSUM, if necessary.
#N_P_file, N_N_file, P_P_file, pair_type, version_pattern, aa_2_compute = user_chooses_condition()


manually_vs_auto = input("Do you want to introduce all parameters manually? This will allow you execute several conditions in different terminals. If not, you'll only introduce the AF filter:\n\t-auto\n\t-manual\n\nWrite your decision as it appears in the options above please: ") 



if manually_vs_auto == 'auto':
    # Ask the user which version of pairs files is desired
    freq = input("Write the desired gnomAD files version. You can choose between:\n\tfreq-2\n\tfreq-3\n\tfreq-4\n\tfreq-5\n\tfreq-6\n\tfreq-7\n\tfreq-8\n\tno_freq\nPlease, write the version as indicated: ")

    #full_vs_seed = input("Write the desired alignment type. You can choose between:\n\tfull\n\tseed\nPlease, write the version as indicated: ")
    align_types = ['mix', 'full', 'seed']

    normal_vs_blosum = ['normal', 'blosum']

    for full_vs_seed in align_types:
        for pair_type in normal_vs_blosum:

            if pair_type == 'blosum':
            # Ask the user if wants to compute BLOSUM62 to initial_aa, final_aa or both
            #aa_2_compute = input("We are going to compute BLOSUM62 to obtain pairs of similar aminoacid changes. Write the aminoacid/s which you want to compute:\n\tinitial_aa\n\tfinal_aa\n\tboth_aa\nPlease, write your choice as indicated: ")
                blosum_aas = ['initial_aa', 'final_aa', 'both_aa']
                for aa_2_compute in blosum_aas:
                    print_to_user(aa_2_compute, freq, full_vs_seed)
                    N_P_file, N_N_file, P_P_file = obtain_pair_files(pair_type, freq, full_vs_seed, aa_2_compute)
                    
                    if full_vs_seed == 'mix':
                        add_mix_stats(N_P_file, N_N_file, P_P_file)
                    else:
                        add_stats(N_P_file, N_N_file, P_P_file, full_vs_seed)
            else:
                aa_2_compute = None
                print_to_user(aa_2_compute, freq, full_vs_seed)
                N_P_file, N_N_file, P_P_file = obtain_pair_files(pair_type, freq, full_vs_seed, aa_2_compute)
                print(N_P_file)
                
                if full_vs_seed == 'mix':
                        add_mix_stats(N_P_file, N_N_file, P_P_file)
                else:
                    add_stats(N_P_file, N_N_file, P_P_file, full_vs_seed)

else:
    # Ask the user if want to extract normal pairs or wants to compute BLOSUM62
    pair_type = input("Write the desired pair type. You can choose between:\n\tNormal pairs(normal)\n\tCompute BLOSUM62(blosum)\nPlease, write the version as indicated between parenthesis: ")

    # Ask the user which version of pairs files is desired
    freq = input("Write the desired gnomAD files version. You can choose between:\n\tfreq-2\n\tfreq-3\n\tfreq-4\n\tfreq-5\n\tfreq-6\n\tfreq-7\n\tfreq-8\n\tno_freq\nPlease, write the version as indicated: ")

    # full_vs_seed = input("Write the desired alignment type. You can choose between:\n\tfull\n\tseed\nPlease, write the version as indicated: ")
    align_types = ['mix', 'seed', 'full']

    if pair_type == 'blosum':
        # Ask the user if wants to compute BLOSUM62 to initial_aa, final_aa or both
        aa_2_compute = input("We are going to compute BLOSUM62 to obtain pairs of similar aminoacid changes. Write the aminoacid/s which you want to compute:\n\tinitial_aa\n\tfinal_aa\n\tboth_aa\nPlease, write your choice as indicated: ")
    else:
        aa_2_compute = None

    for full_vs_seed in align_types:
        print_to_user(aa_2_compute, freq, full_vs_seed)
        N_P_file, N_N_file, P_P_file = obtain_pair_files(pair_type, freq, full_vs_seed, aa_2_compute)

        if full_vs_seed == 'mix':
            add_mix_stats(N_P_file, N_N_file, P_P_file)
        else:
            add_stats(N_P_file, N_N_file, P_P_file, full_vs_seed)
