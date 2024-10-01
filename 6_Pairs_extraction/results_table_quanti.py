# Get a quantification table for all pair types (PP, NN, NP) in all conditions

import sys, pandas as pd, os, glob
import math


sys.path.insert(1, '../7_Stats_pairs/utilities/')
from input_files import obtain_pair_files


def read_pair_file(pairs_file):

    pairs_file = pd.read_csv(pairs_file, sep="\t", header=0)
    return pairs_file


def count_pfams(path_folder):
    
    directory_path = os.path.dirname(path_folder)
    path_pfams = os.path.join(directory_path, 'pairs_by_pfam_code')

    # Set output path name
    output_file_path = f"{path_pfams}/PF*.txt"
    output_file_path = os.path.normpath(output_file_path)

    # List the files matching the pattern
    pfam_pairs_list = glob.glob(output_file_path)

    #num_pfams = f'{len(pfam_pairs_list):,}'
    num_pfams = f'{len(pfam_pairs_list)}'

    return num_pfams

full_vs_seed = input("Write the desired alignment type. (Mix condition prioritizes seed alignment and uses full ones when seed is not found)\nYou can choose between:\n\tfull\n\tseed\n\tseed+full(mix)\nPlease, write the version as indicated: ")

""" 
if full_vs_seed != 'full':
    normal_or_blosum = ['normal', 'blosum']
else:
    normal_or_blosum = ['normal']
 """
normal_or_blosum = ['normal', 'blosum']
freqs = ['no_freq', 'freq-8', 'freq-7', 'freq-6', 'freq-5', 'freq-4', 'freq-3', 'freq-2']

######################
#normal_or_blosum = ['normal']#, 'blosum']
#freqs = ['no_freq']  #, 'freq-8', 'freq-7', 'freq-6', 'freq-5', 'freq-4', 'freq-3', 'freq-2']
######################

results_table = [['BLOSUM_AA', 'NORMAL/BLOSUM', 'AF_FILTER', 'PFAM_num', 'N-N_pairs', 'P-P_pairs', 'N-P_pairs', 'TOTAL_PAIRS', f'% error', 'PP + NN', 'log(PP+NN)', 'log(NP)', 'log(Total pairs)', 'NP/PP', 'NP/NN']]
print('\n\n')

for pair_type in normal_or_blosum:
    if pair_type == 'normal':
        for AF_freq in freqs:
            # Get dfs one for each pair type file
            N_P_file, N_N_file, P_P_file = obtain_pair_files(pair_type, AF_freq, full_vs_seed)

            # Get the number of Pfam families in this condition. We pass the path of one of the files, it is not important which of them, the path is extracted in the function.
            num_pfams = count_pfams(N_P_file)

            NN = read_pair_file(N_N_file)
            PP = read_pair_file(P_P_file)
            NP = read_pair_file(N_P_file)
            
            total_pairs = len(NN) + len(PP) + len(NP)
            per_error = f'{round((100*len(NP)) / total_pairs, 2)}%'

            # Calculate different parameters
            pp_plus_nn = len(PP) + len(NN)
            log_pp_nn = math.log10(pp_plus_nn)
            log_total = math.log10(total_pairs)
            if len(NP) > 0:
                log_np = math.log10(len(NP))
                np_vs_pp = len(NP) / len(PP)
                np_vs_nn = len(NP) / len(NN)
                results_row = ['-', pair_type, AF_freq, num_pfams, f'{len(NN)}', f'{len(PP)}', f'{len(NP)}', total_pairs, per_error, f'{pp_plus_nn}', round(log_pp_nn, 2), round(log_np, 2), round(log_total, 2), round(np_vs_pp, 2), round(np_vs_nn, 2)]
            else:
                log_np = '-'
                np_vs_pp = '-'
                np_vs_nn = '-'
                results_row = ['-', pair_type, AF_freq, num_pfams, f'{len(NN)}', f'{len(PP)}', f'{len(NP)}', total_pairs, per_error, f'{pp_plus_nn}', log_pp_nn, log_np, round(log_total, 2), np_vs_pp, np_vs_nn]

            #results_row = ['-', pair_type, AF_freq, num_pfams, f'{len(NN):,}', f'{len(PP):,}', f'{len(NP):,}', '{:,.0f}'.format(total_pairs), per_error, f'{pp_plus_nn:,}', round(log_pp_nn, 2), round(log_np, 2), round(log_total, 2), round(np_vs_pp, 2), round(np_vs_nn, 2)]

            results_table.append(results_row)

            print(f"{pair_type}, {AF_freq}    has been quantified.")

    else:
        aa_blosum = ['final_aa', 'both_aa', 'initial_aa']
        for aa_2_compute in aa_blosum:
            results_table.append(['','','','','','',''])

            for AF_freq in freqs:
                # Get dfs one for each pair type file
                N_P_file, N_N_file, P_P_file = obtain_pair_files(pair_type, AF_freq, full_vs_seed, aa_2_compute)

                # Get the number of Pfam families in this condition. We pass the path of one of the files, it is not important which of them, the path is extracted in the function.
                num_pfams = count_pfams(N_P_file)

                NN = read_pair_file(N_N_file)
                PP = read_pair_file(P_P_file)
                NP = read_pair_file(N_P_file)

                total_pairs = len(NN) + len(PP) + len(NP)
                per_error = f'{round((100*len(NP)) / total_pairs, 2)}%'

                # Calculate different parameters
                # Calculate different parameters
                pp_plus_nn = len(PP) + len(NN)
                log_pp_nn = math.log10(pp_plus_nn)
                log_total = math.log10(total_pairs)
                if len(NP) > 0:
                    log_np = math.log10(len(NP))
                    np_vs_pp = len(NP) / len(PP)
                    if len(NN) > 0:
                        np_vs_nn = len(NP) / len(NN)
                        results_row = [aa_2_compute, pair_type, AF_freq, num_pfams, f'{len(NN)}', f'{len(PP)}', f'{len(NP)}', total_pairs, per_error, f'{pp_plus_nn}', round(log_pp_nn, 2), round(log_np, 2), round(log_total, 2), round(np_vs_pp, 2), round(np_vs_nn, 2)]
                    else:
                        np_vs_nn = '-'
                        results_row = [aa_2_compute, pair_type, AF_freq, num_pfams, f'{len(NN)}', f'{len(PP)}', f'{len(NP)}', total_pairs, per_error, f'{pp_plus_nn}', round(log_pp_nn, 2), round(log_np, 2), round(log_total, 2), round(np_vs_pp, 2), np_vs_nn]

                else:
                    log_np = '-'
                    np_vs_pp = '-'
                    np_vs_nn = '-'
                    results_row = [aa_2_compute, pair_type, AF_freq, num_pfams, f'{len(NN)}', f'{len(PP)}', f'{len(NP)}', total_pairs, per_error, f'{pp_plus_nn}', log_pp_nn, log_np, round(log_total, 2), np_vs_pp, np_vs_nn]

                results_table.append(results_row)
                print(f"{pair_type}, {aa_2_compute}, {AF_freq}    has been quantified.")
print('\n\n')

# Save the results dfs in txt files
results_table = '\n'.join('\t'.join(str(item) for item in inner_list) for inner_list in results_table)

output_paths = [[f'./data/QUANTI_TABLE_PAIRS_{full_vs_seed}.txt', results_table]]

for path, results_table in output_paths:
    with open(path, 'w') as outfile:
        outfile.write(results_table)
