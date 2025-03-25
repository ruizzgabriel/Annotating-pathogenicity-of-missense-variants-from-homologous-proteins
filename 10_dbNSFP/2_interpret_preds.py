import pandas as pd
from utilities.pred_scores import *

# Read file with predictions for our variants forming pairs
strict_pred = pd.read_csv('./results/strict_pairs_w_preds_REDUCED_COLS.txt')


  ## Interpret predictions and get an output for each of one: D(amaging), T(olerated) or A(mbiguous)
# Mutation Taster
strict_pred['MutationTaster_output'] = strict_pred['MutationTaster_pred'].apply(mutationTaster)
# Mutation Assessor
strict_pred['MutationAssessor_output'] = strict_pred['MutationAssessor_pred'].apply(mutationAssessor)
# RPOVEAN
strict_pred['PROVEAN_output'] = strict_pred['PROVEAN_pred'].apply(PROVEAN)
# VEST4
strict_pred['VEST4_output'] = strict_pred['VEST4_score'].apply(VEST4)
# MetaSVM_pred
strict_pred['MetaSVM_output'] = strict_pred['MetaSVM_pred']
# MetaLR_pred
strict_pred['MetaLR_output'] = strict_pred['MetaLR_pred']
# MetaRNN_pred
strict_pred['MetaRNN_output'] = strict_pred['MetaRNN_pred'].apply(MetaRNN)
# M-CAP_pred
strict_pred['M-CAP_output'] = strict_pred['M-CAP_pred'].apply(MCAP)
# REVEL_score
strict_pred['REVEL_output'] = strict_pred['REVEL_score'].apply(REVEL)
# MVP_score
strict_pred['MVP_output'] = strict_pred['MVP_score'].apply(MVP_gMVP)
# gMVP_score
strict_pred['gMVP_output'] = strict_pred['gMVP_score'].apply(MVP_gMVP)
# PrimateAI_pred
strict_pred['PrimateAI_output'] = strict_pred['PrimateAI_pred']
# DEOGEN2_pred
strict_pred['DEOGEN2_output'] = strict_pred['DEOGEN2_pred'].apply(DEOGEN2)
# BayesDel_addAF_pred
strict_pred['BayesDel_addAF_output'] = strict_pred['BayesDel_addAF_pred']
# BayesDel_noAF_pred
strict_pred['BayesDel_noAF_output'] = strict_pred['BayesDel_noAF_pred']
# ClinPred_pred
strict_pred['ClinPred_output'] = strict_pred['ClinPred_pred']
# LIST-S2_pred
strict_pred['LIST-S2_output'] = strict_pred['LIST-S2_pred'].apply(LIST_S2)
# VARITY_R_score
strict_pred['VARITY_R_output'] = strict_pred['VARITY_R_score'].apply(VARITY)
# VARITY_ER_score
strict_pred['VARITY_ER_output'] = strict_pred['VARITY_ER_score'].apply(VARITY)
# VARITY_R_LOO_score
strict_pred['VARITY_R_LOO_output'] = strict_pred['VARITY_R_LOO_score'].apply(VARITY)
# VARITY_ER_LOO_score
strict_pred['VARITY_ER_LOO_output'] = strict_pred['VARITY_ER_LOO_score'].apply(VARITY)
# ESM1b_pred
strict_pred['ESM1b_output'] = strict_pred['ESM1b_pred'].apply(ESM1b)
# PHACTboost_score -> no clue
#strict_pred['PHACTboost_output'] = strict_pred['PHACTboost_score'].apply(PHACT)
# MutScore_score
strict_pred['MutScore_output'] = strict_pred['MutScore_score'].apply(MutScore)
# CADD_phred
strict_pred['CADD_output'] = strict_pred['CADD_phred'].apply(CADD)
# DANN_score
strict_pred['DANN_output'] = strict_pred['DANN_score'].apply(DANN)
# fathmm-XF_coding_pred
strict_pred['fathmm-XF_coding_output'] = strict_pred['fathmm-XF_coding_pred'].apply(fathmm_XF_coding)

# phyloP100way_vertebrate
strict_pred['phyloP100way_vertebrate_output'] = strict_pred['phyloP100way_vertebrate'].apply(phyloP)
# phyloP470way_mammalian
strict_pred['phyloP470way_mammalian_output'] = strict_pred['phyloP470way_mammalian'].apply(phyloP)
# phyloP17way_primate
strict_pred['phyloP17way_primate_output'] = strict_pred['phyloP17way_primate'].apply(phyloP)


# phastCons100way_vertebrate
strict_pred['phastCons100way_vertebrate_output'] = strict_pred['phastCons100way_vertebrate'].apply(phastCons)
# phastCons470way_mammalian
strict_pred['phastCons470way_mammalian_output'] = strict_pred['phastCons470way_mammalian'].apply(phastCons)
# phastCons17way_primate
strict_pred['phastCons17way_primate_output'] = strict_pred['phastCons17way_primate'].apply(phastCons)


# Filter columns
specific_columns = ["#chr", "pos(1-based)", "aaref", "aaalt", "aapos", "genename", 'HGVSp_VEP', 'Uniprot_entry', 'Uniprot_acc']
strict_pred = strict_pred[[col for col in strict_pred.columns if 'output' in col or col in specific_columns]]
strict_pred.drop_duplicates(inplace=True)
#print(strict_pred)

# Get the output file
strict_pred.to_csv('./results/strict_pairs_PREDICTIONS.txt', index=False)