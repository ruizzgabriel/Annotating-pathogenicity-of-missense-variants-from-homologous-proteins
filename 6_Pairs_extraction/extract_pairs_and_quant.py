# We here detect pairs (similar or homologous) of variants, extract to different files and quantify them

import os, glob, pandas as pd
from colorama import Fore, Style
#from utilities.extraction import extract_normal_pairs, extract_blosum_pairs
from utilities.extraction import group_df_normal_pairs, extract_blosum_pairs
from utilities.quantification import extract_groups, extract_pairs
from utilities.handle_files import get_output_path, get_input_path, list_to_txt

def find_pairs(only_quanti, pair_type, pfam_list, aa_2_compute):
    # Initialize empty lists
        # Non-pathogenic/Non-pathogenic
    N_N = []
        # Non-pathogenic/Pathogenic
    N_P = []
        # Pathogenic/Pathogenic
    P_P = []

    g = 0
    if pair_type == 'blosum':
        # Extract filtered Pfam files with exact pairs and equivalent final_aa changes
        for alignment in pfam_list:
            # If only quanti is false means that pairs files by Pfam codes have also to be computed
            if only_quanti == False:
                pfam_pairs = extract_blosum_pairs(alignment, aa_2_compute, version_pattern, full_vs_seed)
            else:
                print('NOOO')
                pfam_pairs = pd.read_csv(alignment, sep="\t", header=0)

            # Quantification
            if pfam_pairs is not None:
                pfam_grouped = extract_groups(pfam_pairs, pair_type='blosum', aa_2_compute=aa_2_compute)
                for name, group in pfam_grouped:

                    if aa_2_compute == 'both_aa':
                        if ('Equivalent_final_aa' not in group.columns) or ('Equivalent_initial_aa' not in group.columns):
                            #print(group.columns)
                            g += 1
                            continue
                    #print(group)
                    #print(group['Equivalent_initial_aa'])
                    extract_pairs(group, pair_type='blosum', N_N=N_N, P_P=P_P, N_P=N_P, aa_2_compute=aa_2_compute)

        print(f'Errors in both_aa: {g}')

    else:
        # Iterate over the list of files
        for alignment in pfam_list:
            if only_quanti == False:
                pairs_df = group_df_normal_pairs(alignment, version_pattern, full_vs_seed)
            else:
                pairs_df = pd.read_csv(alignment, sep="\t", header=0)
                pairs_df = pairs_df.groupby(['initial_aa', 'Pos_align', 'final_aa'])

            if pairs_df is not None:
                for name, group in pairs_df:
                    extract_pairs(group, pair_type='normal', N_N=N_N, P_P=P_P, N_P=N_P)

    # Print the obtained pairs
    print(f"From the {(len(N_N) + len(N_P) + len(P_P)):,} total pairs:", '\n', f"We have {len(N_N):,} pairs corresponding to Non-pathogenic/Non-pathogenic", '\n', f"We have {len(P_P):,} pairs corresponding to Pathogenic/Pathogenic", '\n', f"We have {len(N_P):,} pairs corresponding to Pathogenic/Non-pathogenic")

    return N_N, N_P, P_P


def start_program(pair_type, only_quanti, full_vs_seed, version_pattern, aa_2_compute):
    N_P_file, N_N_file, P_P_file, output_path = get_output_path(pair_type, full_vs_seed, version_pattern, aa_2_compute)
    pfam_list, only_quanti, aa_2_compute = get_input_path(only_quanti, pair_type, full_vs_seed, version_pattern, aa_2_compute)

    if aa_2_compute != None:
        print(Fore.YELLOW + f"\n\nWe start computing HOMOLOGOUS pairs.\n\t -> {aa_2_compute}")
    else:
        print(Fore.YELLOW + f"\n\nWe'll compute STRICT pairs.\n")
    print(Style.RESET_ALL)

    print(Fore.GREEN + f"\n\nWe start extracting pairs with:\n\tFREQUENCE: {version_pattern.split('_')[-1]}")
    print(Fore.MAGENTA + f"\n\tALIGNMENT: {full_vs_seed}")
    print(Style.RESET_ALL)

    N_N, N_P, P_P = find_pairs(only_quanti, pair_type, pfam_list, aa_2_compute)

    # Save the different classified DataFrames with their corresponding file names
    #print(N_N[0])
    list_to_txt(N_N, N_N_file, pair_type)
    list_to_txt(P_P, P_P_file, pair_type)
    list_to_txt(N_P, N_P_file, pair_type)

    # Set output path name
    output_file_path = f"{output_path}/pairs_by_pfam_code/PF*.txt"
    output_file_path = os.path.normpath(output_file_path)
    #print(output_file_path)
    # List the files matching the pattern
    pfam_pairs_list = glob.glob(output_file_path)
    print(f'We found pairs in {len(pfam_pairs_list)} Pfam families.')


######################################################################
# ATTENTION: be careful with the commented inputs. Some of them are silenced because in this moment we are not interested on them

manually_vs_auto = input("Do you want to introduce all parameters manually? This will allow you execute several conditions in different terminals. If not, you'll only introduce the AF filter:\n\t-auto\n\t-manual\n\nWrite your decision as it appears in the options above please: ") 

if manually_vs_auto == 'auto':
    # Ask the user which version of pairs files is desired
    version_files = input("Write the desired gnomAD files version. You can choose between:\n\tfreq-2\n\tfreq-3\n\tfreq-4\n\tfreq-5\n\tfreq-6\n\tfreq-7\n\tfreq-8\n\tno_freq\nPlease, write the version as indicated: ")
    version_pattern = '' if version_files == 'no_freq' else '_' + version_files

    # Ask the user if already have the pair files or wants to compute the whole process. Useful when you have to include other information that was not included in the final pairs quanti files (N_P, P_P and N_N)
    #only_quanti = input("Do you have the pair files? If you only want to quantify them write 'quanti'. This is useful when you have to include other information that was not included in the final pairs quanti files (N_P, P_P and N_N).\n If you come from step 5 and have not still extracted the pairs by Pfam code write 'all':\n\tOnly quantification(quanti)\n\tCompute whole process(all)\nPlease, write the version as indicated between parenthesis: ")
    only_quanti = 'all'

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
                    start_program(pair_type, only_quanti, full_vs_seed, version_pattern, aa_2_compute)
            else:
                aa_2_compute = None
                start_program(pair_type, only_quanti, full_vs_seed, version_pattern, aa_2_compute)

else:
    # Ask the user if wants to extract normal pairs or wants to compute BLOSUM62 to search homologous pairs
    pair_type = input("Write the desired pair type. You can choose between:\n\tNormal pairs(normal)\n\tCompute BLOSUM62(blosum)\nPlease, write the version as indicated between parenthesis: ")

    # Ask the user which version of pairs files is desired
    version_files = input("Write the desired gnomAD files version. You can choose between:\n\tfreq-2\n\tfreq-3\n\tfreq-4\n\tfreq-5\n\tfreq-6\n\tfreq-7\n\tfreq-8\n\tno_freq\nPlease, write the version as indicated: ")
    version_pattern = '' if version_files == 'no_freq' else '_' + version_files

    # Ask the user if already have the pair files or wants to compute the whole process. Useful when you have to include other information that was not included in the final pairs quanti files (N_P, P_P and N_N)
    only_quanti = input("Do you have the pair files? If you only want to quantify them write 'quanti'. This is useful when you have to include other information that was not included in the final pairs quanti files (N_P, P_P and N_N).\n If you come from step 5 and have not still extracted the pairs by Pfam code write 'all':\n\tOnly quantification(quanti)\n\tCompute whole process(all)\nPlease, write the version as indicated between parenthesis: ")
    
    # full_vs_seed = input("Write the desired alignment type. You can choose between:\n\tfull\n\tseed\nPlease, write the version as indicated: ")
    align_types = ['mix', 'full', 'seed']

    if pair_type == 'blosum':
        # Ask the user if wants to compute BLOSUM62 to initial_aa, final_aa or both
        aa_2_compute = input("We are going to compute BLOSUM62 to obtain pairs of similar aminoacid changes. Write the aminoacid/s which you want to compute:\n\tinitial_aa\n\tfinal_aa\n\tboth_aa\nPlease, write your choice as indicated: ")
    else:
        aa_2_compute = None

    for full_vs_seed in align_types:
            start_program(pair_type, only_quanti, full_vs_seed, version_pattern, aa_2_compute)
