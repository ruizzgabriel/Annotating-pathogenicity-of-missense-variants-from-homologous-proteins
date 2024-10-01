import pandas as pd, os, glob

def get_output_path(pair_type, full_vs_seed, version_pattern, aa_2_compute):
    ################## OUPTUT FILE NAMES ##################
    if pair_type == 'blosum':
        
        # Set data path and the output data folder        
        if full_vs_seed == 'full':
            data_path = f"./data/homologous_pairs/FULL_align/blosum62_{aa_2_compute}/pfam_pairs_interpro_+score{version_pattern}/pairs_by_pfam_code"
            output_path = f"./data/homologous_pairs/FULL_align/blosum62_{aa_2_compute}/pfam_pairs_interpro_+score{version_pattern}"
        elif full_vs_seed == 'seed':
            data_path = f"./data/homologous_pairs/SEED_align/blosum62_{aa_2_compute}/pfam_pairs_interpro_+score{version_pattern}/pairs_by_pfam_code"
            output_path = f"./data/homologous_pairs/SEED_align/blosum62_{aa_2_compute}/pfam_pairs_interpro_+score{version_pattern}"
        elif full_vs_seed == 'mix':
            data_path = f"./data/homologous_pairs/SEED+FULL_align/blosum62_{aa_2_compute}/pfam_pairs_interpro_+score{version_pattern}/pairs_by_pfam_code"
            output_path = f"./data/homologous_pairs/SEED+FULL_align/blosum62_{aa_2_compute}/pfam_pairs_interpro_+score{version_pattern}"

        print(data_path)

        # Set the files pattern of pairs
        #pfam_pattern = f"{data_path}/PF*_pairs_interpro_+score_{up_low_case}{version_pattern}.txt"
        pfam_pattern = f"{data_path}/PF*_pairs_interpro_+score{version_pattern}.txt"
        pfam_pattern = os.path.normpath(pfam_pattern)

        if not os.path.exists(output_path):
            # If does not exist, create it
            os.makedirs(output_path)

        # Set the names of the output files
        N_P_file = f"{output_path}/N_P_pairs_interpro_+score_blosum62_{aa_2_compute}{version_pattern}.txt"
        N_N_file = f"{output_path}/N_N_pairs_interpro_+score_blosum62_{aa_2_compute}{version_pattern}.txt"
        P_P_file = f"{output_path}/P_P_pairs_interpro_+score_blosum62_{aa_2_compute}{version_pattern}.txt"


    else:
        # Set data path and the output data folder
        if full_vs_seed == 'full':
            #data_path = f"./data/normal_pairs/pfam_pairs_AD_missense_{up_low_case}{version_pattern}/pairs_by_pfam_code"
            data_path = f"./data/normal_pairs/FULL_align/pfam_pairs_AD_missense{version_pattern}/pairs_by_pfam_code"
            output_path = f"./data/normal_pairs/FULL_align/pfam_pairs_AD_missense{version_pattern}"
        elif full_vs_seed == 'seed':
            data_path = f"./data/normal_pairs/SEED_align/pfam_pairs_AD_missense{version_pattern}/pairs_by_pfam_code"
            output_path = f"./data/normal_pairs/SEED_align/pfam_pairs_AD_missense{version_pattern}"
        elif full_vs_seed == 'mix':
            data_path = f"./data/normal_pairs/SEED+FULL_align/pfam_pairs_AD_missense{version_pattern}/pairs_by_pfam_code"
            output_path = f"./data/normal_pairs/SEED+FULL_align/pfam_pairs_AD_missense{version_pattern}"

        # Set the files pattern of pairs *
        #pfam_pattern = f"{data_path}/PF*_pairs_{up_low_case}{version_pattern}.txt"
        pfam_pattern = f"{data_path}/PF*_pairs{version_pattern}.txt"
        pfam_pattern = os.path.normpath(pfam_pattern)

        if not os.path.exists(output_path):
                # If does not exist, create it
                os.makedirs(output_path)
                
        # Set the output file names *
        #N_P_file = f"{output_path}/N_P_pairs_{up_low_case}{version_pattern}.txt"
        N_P_file = f"{output_path}/N_P_pairs{version_pattern}.txt"

        #N_N_file = f"{output_path}/N_N_pairs_{up_low_case}{version_pattern}.txt"
        N_N_file = f"{output_path}/N_N_pairs{version_pattern}.txt"

        #P_P_file = f"{output_path}/P_P_pairs_{up_low_case}{version_pattern}.txt"
        P_P_file = f"{output_path}/P_P_pairs{version_pattern}.txt"


    N_P_file = os.path.normpath(N_P_file)
    N_N_file = os.path.normpath(N_N_file)
    P_P_file = os.path.normpath(P_P_file)
    
    return N_P_file, N_N_file, P_P_file, output_path


def get_input_path(only_quanti, pair_type, full_vs_seed, version_pattern, aa_2_compute):
    ################## INPUT FILE NAMES ##################
    if only_quanti != 'quanti':
        only_quanti = False
        # Set data path
        if full_vs_seed == 'full':
            #data_path = f"../5_Split_entries_by_align/data/pfam_groups_interpro_{up_low_case}{version_pattern}"
            data_path = f"../5_Split_entries_by_align/data/FULL_align/pfam_groups_interpro{version_pattern}"
        elif full_vs_seed == 'seed':
            data_path = f"../5_Split_entries_by_align/data/SEED_align/pfam_groups_interpro{version_pattern}"
        else:
            data_path = f"../5_Split_entries_by_align/data/SEED+FULL_align/pfam_groups_interpro{version_pattern}"


        # Set the files pattern of Pfam
        pfam_pattern = f"{data_path}/PF*.txt"
        pfam_pattern = os.path.normpath(pfam_pattern)

        # List the files matching the pattern
        pfam_list = glob.glob(pfam_pattern)
        print(f"In step 5, there are {len(pfam_list)} alignment files for the indicated frequency.\n\n")


    ################## Input file names only for quanti ##################
    # THIS PIECE OF CODE IS USEFUL IF YOU WANT ONLY TO QUANTIFY PAIRS. THAT IS, WHEN YOU ALREADY HAVE THE FILES OF PAIRS AND YOU WANT TO OBTAIN THE N_P, P_P AND N_N FINAL FILES.
    # THAT IS USEFUL IF YOU NEED TO RETRIEVE AN EXTRA INFORMATION FOR EACH ROW. FOR INSTANCE, AT FIRST WE DID NOT INCLUDE UNIPROT NAMES AND WE HAD TO IMPLEMENT THIS TO INCLUDE THEM.
    else: 
        only_quanti = True
        # Set the output file name
        if full_vs_seed == 'full':
            if pair_type == "blosum":
                data_path = f"./data/homologous_pairs/FULL_align/blosum62_{aa_2_compute}/pfam_pairs_interpro_+score{version_pattern}/pairs_by_pfam_code"
            else:
                data_path = f"./data/normal_pairs/FULL_align/pfam_pairs_AD_missense{version_pattern}/pairs_by_pfam_code"
        
        elif full_vs_seed == 'seed':
            if pair_type == "blosum":
                data_path = f"./data/homologous_pairs/SEED_align/blosum62_{aa_2_compute}/pfam_pairs_interpro_+score{version_pattern}/pairs_by_pfam_code"
            else:
                data_path = f"./data/normal_pairs/SEED_align/pfam_pairs_AD_missense{version_pattern}/pairs_by_pfam_code"
        else:
            if pair_type == "blosum":
                data_path = f"./data/homologous_pairs/SEED+FULL_align/blosum62_{aa_2_compute}/pfam_pairs_interpro_+score{version_pattern}/pairs_by_pfam_code"
            else:
                data_path = f"./data/normal_pairs/SEED+FULL_align/pfam_pairs_AD_missense{version_pattern}/pairs_by_pfam_code"
        data_path = os.path.normpath(data_path)

        # Set the files pattern of Pfam
        pfam_pattern = f"{data_path}/PF*_pairs*.txt"
        pfam_pattern = os.path.normpath(pfam_pattern)

        # List the files matching the pattern
        pfam_list = glob.glob(pfam_pattern)
        print(f"In step 6, there are {len(pfam_list)} alignment files for the indicated frequency. Let's quantify the pairs!\n\n")

    return pfam_list, only_quanti, aa_2_compute


def list_to_txt(rows_list, output_name, pair_type):
    """
        Transform a list of rows into a DataFrame.
        Write the DataFrame to a tab-separated text file.

        Args:
            rows_list (list): The list with the information to be written.
            output_name (str): The name of the output text file.
            pair_type(str): Indicating normal pairs computation or BLOSUM62 computation. (normal or blosum)
    """
    if pair_type == 'normal':
        header = ['Pfam_code', 'initial_aa', 'pos_align', 'final_aa', 'gene_name', 'gene_name2', 'pos_prot_change', 'pos_prot_change2', 'db1', 'db2', 'uniprot1', 'uniprot2']  # 'variation_id', 'variation_id2',
    else:
        header = ['Pfam_code', 'initial_aa', 'pos_align', 'final_aa', 'final_aa2', 'gene_name', 'gene_name2', 'pos_prot_change', 'pos_prot_change2', 'db1', 'db2', 'uniprot1', 'uniprot2']
    
    df = pd.DataFrame(rows_list, columns=header)
    
    df.to_csv(output_name, sep="\t", index=False)