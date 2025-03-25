import pandas as pd
from Bio.Data.IUPACData import protein_letters_3to1

# Read file including variants from pairs
vars_pairs_strict = pd.read_excel('./data/Vars_from_strict_pairs.xlsx')
# Read file with dbNSFP predictions
strict_pred = pd.read_csv('./results/server_uvic/strict_pairs_w_predictions_dbNSFP.txt', sep=',')
strict_pred.drop_duplicates(inplace=True)
print(strict_pred.shape)

# List of specific columns to keep
specific_columns = ["#chr", "pos(1-based)", "aaref", "aaalt", "aapos", "genename", 'HGVSp_VEP', 'Uniprot_entry', 'Uniprot_acc']

# Define the columns to keep
columns_to_keep = [
    "#chr", "pos(1-based)", "ref", "alt", "aaref", "aaalt", "rs_dbSNP", "aapos", "genename", "Ensembl_geneid", 
    "Ensembl_transcriptid",
    "Ensembl_proteinid", "Uniprot_acc", "Uniprot_entry", "HGVSc_snpEff", "HGVSp_snpEff", "HGVSc_VEP",
    "HGVSp_VEP", "APPRIS", "GENCODE_basic", "TSL", "VEP_canonical", "MANE", "cds_strand", "refcodon",
    "codonpos", "codon_degeneracy", "SIFT4G_pred", "Polyphen2_HVAR_pred", "MutationTaster_pred",
    "MutationAssessor_pred", "PROVEAN_pred", "VEST4_score", "VEST4_rankscore", "MetaSVM_pred",
    "MetaLR_pred", "MetaRNN_pred", "M-CAP_pred", "REVEL_score", "REVEL_rankscore", "MutPred_score",
    "MutPred_rankscore", "MVP_score", "MVP_rankscore", "gMVP_score", "gMVP_rankscore", "MPC_score",
    "MPC_rankscore", "PrimateAI_pred", "DEOGEN2_pred", "BayesDel_addAF_pred", "BayesDel_noAF_pred",
    "ClinPred_pred", "LIST-S2_pred", "VARITY_R_score", "VARITY_R_rankscore", "VARITY_ER_score",
    "VARITY_ER_rankscore", "VARITY_R_LOO_score", "VARITY_R_LOO_rankscore", "VARITY_ER_LOO_score",
    "VARITY_ER_LOO_rankscore", "ESM1b_pred", "AlphaMissense_pred", "PHACTboost_score",
    "PHACTboost_rankscore", "MutFormer_score", "MutFormer_rankscore", "MutScore_score",
    "MutScore_rankscore", "Aloft_pred", "CADD_raw", "CADD_raw_rankscore", "CADD_phred", "DANN_score",
    "DANN_rankscore", "fathmm-XF_coding_pred", "Eigen-raw_coding", "Eigen-raw_coding_rankscore",
    "Eigen-phred_coding", "Eigen-PC-raw_coding", "Eigen-PC-raw_coding_rankscore", "Eigen-PC-phred_coding",
    "GERP++_NR", "GERP++_RS", "GERP++_RS_rankscore", "GERP_91_mammals", "GERP_91_mammals_rankscore",
    "phyloP100way_vertebrate", "phyloP100way_vertebrate_rankscore", "phyloP470way_mammalian",
    "phyloP470way_mammalian_rankscore", "phyloP17way_primate", "phyloP17way_primate_rankscore",
    "phastCons100way_vertebrate", "phastCons100way_vertebrate_rankscore", "phastCons470way_mammalian",
    "phastCons470way_mammalian_rankscore", "phastCons17way_primate", "phastCons17way_primate_rankscore",
    "bStatistic", "bStatistic_converted_rankscore"
]  
# Remove the columns from dbNSFP ranskscores
strict_pred = strict_pred[[col for col in columns_to_keep if not 'rankscore' in col]]
#print(strict_pred.shape)

# Get a file with reduced number of columns
strict_pred.to_csv('./results/strict_pairs_w_preds_REDUCED_COLS.txt', index=False)