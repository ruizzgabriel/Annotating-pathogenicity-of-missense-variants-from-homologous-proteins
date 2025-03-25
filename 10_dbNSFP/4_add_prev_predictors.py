import pandas as pd

def get_sift_score(value):
    if isinstance(value, str):
        value = value.strip()
        if value == 'Deleterious':
            return 'D'
        elif value == 'Neutral':
            return 'T'
    return 'A'

def get_poly_score(value):
    if isinstance(value, str):
        value = value.strip()
        if 'damaging' in value:
            return 'D'
        elif value == 'benign':
            return 'T'
    return 'A'

def get_alph_score(value):
    if isinstance(value, str):
        value = value.strip()
        if value == 'pathogenic':
            return 'D'
        elif value == 'benign':
            return 'T'
    return 'A'

def translate_score(row):
    alphamissense = row['ALPHAMISSENSE_prediction']
    sift = row['SIFT_prediction']
    polyphen = row['POLY_prediction']

    sift_score = get_sift_score(sift)
    poly_score = get_poly_score(polyphen)
    alph_score = get_alph_score(alphamissense)
    return alph_score, sift_score, poly_score


# File of the first predictors compared: AlphaMissense, SIFT and PolyPhen-2
prev_pred = pd.read_csv('./data/Evaluated_vars_from_strict_pairs(alphamissense).txt')
#print(prev_pred)

# Interpret the predictors to get the standard outputs: D(amaging), T(olerated) or A(mbiguous)
prev_pred[['AlphaMissense_score', 'SIFT_score', 'PolyPhen-2']] = prev_pred.apply(translate_score, axis=1, result_type='expand')

#print(prev_pred[['db', 'Verify_ALPHAMISSENSE', 'AlphaMissense_score', 'Verify_SIFT', 'SIFT_score', 'Verify_POLYPHEN', 'PolyPhen-2']])

# Drop not necessary columns
prev_pred = prev_pred.drop(['SIFT_prediction', 'POLY_prediction', 'ALPHAMISSENSE_prediction', 'Verify_ALPHAMISSENSE', 'Verify_SIFT', 'Verify_POLYPHEN'], axis=1)
prev_pred.drop_duplicates(inplace=True)
#print(prev_pred)

# Read the file with the predictions in dbNSFP, we want to add the other predictors to the others
other_preds = pd.read_csv('./results/strict_pairs_preds_vs_homolvar.txt')
other_preds.drop_duplicates(inplace=True)
print(other_preds)

# Merge the two DataFrames based on the specified columns
    # This keeps all rows from prev_pred and adds matching columns from other_preds
merged_df = prev_pred.merge(other_preds, 
                      on=['Pfam_code', 'initial_aa', 'pos_align', 'final_aa', 
                          'gene_name', 'pos_prot_change', 'db', 'uniprot',       # MERGE NOT WORKING PROPERLY --> DATA LOOSE
                          'score_substmat', 'seq_identity(%)', 'Pair_type'], 
                      how='left')

#merged_df = merged_df.drop(['score_substmat'], axis=1)
merged_df.drop_duplicates(inplace=True)
#print(merged_df)

# Once we merged the results and we have all the pathogenicity predictors, we obtain an output file
merged_df.to_csv('./results/ALL_predictors_strict_pairs.txt', sep='\t', index=False)