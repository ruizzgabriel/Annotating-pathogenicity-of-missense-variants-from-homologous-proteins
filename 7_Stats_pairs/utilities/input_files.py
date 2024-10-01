import os

def user_chooses_condition():
    """
    DEPRACATED --> revise paths
    """

    # Ask the user if want to treat normal pairs or wants to compute BLOSUM62
    pair_type = input("Write the desired pair type. You can choose between:\n\tNormal pairs(normal)\n\tCompute BLOSUM62(blosum)\nPlease, write the version as indicated between parenthesis: ")
    # Ask the user which version of pairs files is desired
    version_files = input("Write the desired gnomAD files version. You can choose between:\n\tfreq-2\n\tfreq-3\n\tfreq-4\n\tfreq-5\n\tfreq-6\n\tfreq-7\n\tfreq-8\n\tno_freq\nPlease, write the version as indicated: ")
    version_pattern = "" if version_files == "no_freq" else f"_{version_files}"

    # Initialize this variable as None because in case is not defined, it couldn't be returned. We need it for saving the file in the main program, if necessary.
    aa_2_compute = None

    if pair_type == "blosum":
        # Ask the user if wants to treat pairs with computed BLOSUM62 to initial_aa, final_aa or both
        aa_2_compute = input("We are going to compute BLOSUM62 to obtain pairs of similar aminoacid changes. Write the aminoacid/s which you want to compute:\n\tinitial_aa\n\tfinal_aa\n\tboth_aa\nPlease, write your choice as indicated: ")
        # Set data path
        data_path = f"../6_Pairs_extraction/data/homologous_pairs/blosum62_{aa_2_compute}/pfam_pairs_interpro_+score{version_pattern}"
        data_path = os.path.normpath(data_path)

        # Set the names of the output files
        N_P_file = f"{data_path}/N_P_pairs_interpro_+score_blosum62_{aa_2_compute}{version_pattern}.txt"
        N_N_file = f"{data_path}/N_N_pairs_interpro_+score_blosum62_{aa_2_compute}{version_pattern}.txt"
        P_P_file = f"{data_path}/P_P_pairs_interpro_+score_blosum62_{aa_2_compute}{version_pattern}.txt"

    else:
        # Set data path
        data_path = f"../6_Pairs_extraction/data/normal_pairs/pfam_pairs_AD_missense{version_pattern}"

        # Set the output file names
        N_P_file = f"{data_path}/N_P_pairs{version_pattern}.txt"
        N_N_file = f"{data_path}/N_N_pairs{version_pattern}.txt"
        P_P_file = f"{data_path}/P_P_pairs{version_pattern}.txt"

    N_P_file = os.path.normpath(N_P_file)
    N_N_file = os.path.normpath(N_N_file)
    P_P_file = os.path.normpath(P_P_file)
    
    return N_P_file, N_N_file, P_P_file, pair_type, version_pattern, aa_2_compute


def obtain_pair_files(pair_type, freq, full_vs_seed, aa_2_compute=None):

    freq = "" if freq == "no_freq" else f"_{freq}"

    if full_vs_seed == 'full':
        path_align = 'FULL_align'
    elif full_vs_seed == 'seed':
        path_align = 'SEED_align'
    elif full_vs_seed == 'mix':
        path_align = 'SEED+FULL_align'

    if pair_type == "blosum":
        # Ask the user if wants to treat pairs with computed BLOSUM62 to initial_aa, final_aa or both
        #aa_2_compute = input("We are going to compute BLOSUM62 to obtain pairs of similar aminoacid changes. Write the aminoacid/s which you want to compute:\n\tinitial_aa\n\tfinal_aa\n\tboth_aa\nPlease, write your choice as indicated: ")

        # Set data path
        data_path = f"../6_Pairs_extraction/data/homologous_pairs/{path_align}/blosum62_{aa_2_compute}/pfam_pairs_interpro_+score{freq}"
        data_path = os.path.normpath(data_path)

        # Set the names of the output files
        N_P_file = f"{data_path}/N_P_pairs_interpro_+score_blosum62_{aa_2_compute}{freq}.txt"
        N_N_file = f"{data_path}/N_N_pairs_interpro_+score_blosum62_{aa_2_compute}{freq}.txt"
        P_P_file = f"{data_path}/P_P_pairs_interpro_+score_blosum62_{aa_2_compute}{freq}.txt"

    else:
        # Set data path
        data_path = f"../6_Pairs_extraction/data/normal_pairs/{path_align}/pfam_pairs_AD_missense{freq}"

        # Set the output file names
        N_P_file = f"{data_path}/N_P_pairs{freq}.txt"
        N_N_file = f"{data_path}/N_N_pairs{freq}.txt"
        P_P_file = f"{data_path}/P_P_pairs{freq}.txt"

    N_P_file = os.path.normpath(N_P_file)
    N_N_file = os.path.normpath(N_N_file)
    P_P_file = os.path.normpath(P_P_file)
    
    return N_P_file, N_N_file, P_P_file


def save_file(pairs_file, file_path):
    #output_path = f"{file}/{output_name}"
    # print(output_path)

    # Save the pairs file with the new data to a CSV file
    pairs_file.to_csv(file_path, header=True, index=False, sep="\t")
