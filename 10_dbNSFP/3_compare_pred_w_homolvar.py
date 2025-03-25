import pandas as pd
from Bio.Data.IUPACData import protein_letters_3to1

# Read file with the variants forming pairs
vars_pairs_strict = pd.read_excel('./data/Vars_from_strict_pairs.xlsx')
# Read the file with the dbNSFP predictions interpreted, with the columns of interest
strict_pred = pd.read_csv('./results/strict_pairs_PREDICTIONS.txt')

new_data = []
# Get the output interpretations for our predictor HomolVar
vars_pairs_strict['HomolVar'] = vars_pairs_strict['db'].map({'gnomAD':'T', 'ClinVar':'D'})

# Iterate over the variants, to add the predictions to the corresponding one
for i, row in vars_pairs_strict.iterrows():
    gene = row['gene_name']
    i_aa = row['initial_aa']
    pos = row['pos_prot_change']
    f_aa = row['final_aa']

    if i_aa in protein_letters_3to1:
            i_aa = protein_letters_3to1[i_aa]
    if f_aa in protein_letters_3to1:
        f_aa = protein_letters_3to1[f_aa]
    try:
        pos = int(pos)
    except:
        print('something went wrong with the pos: ', pos)
        continue
    
    # Filter the conditions to match the specific variant
    gene_condition = (strict_pred['genename'].str.contains(gene, na=False))
    iaa_cond = (strict_pred['aaref'] == i_aa)
    faa_cond = (strict_pred['aaalt'] == f_aa)
    filt_vars = strict_pred[gene_condition & iaa_cond & faa_cond]
    
    #if gene == 'ACTN2':
    #    print(filt_vars[['genename', 'aaref', "aapos", 'aaalt']])

    filt_vars = filt_vars.copy()
    filt_vars["aapos"] = filt_vars["aapos"].astype(str)
    pred_vars = filt_vars[filt_vars["aapos"].str.contains(str(pos), na=False)]

    # Extract only columns with 'output' in the name
    pred_cols = [col for col in pred_vars.columns if 'output' in col]

    # If matching data exists, store it
    if not pred_vars.empty:
        row_data = row.to_dict()
        row_data.update(pred_vars[pred_cols].iloc[0].to_dict())  # Add matching pred columns
        new_data.append(row_data)

# Create a new DataFrame with merged data
vars_pairs_strict = pd.DataFrame(new_data)

# Get the output file
vars_pairs_strict.to_csv('./results/strict_pairs_preds_vs_homolvar.txt', index=False)



###################################################################################
# Let's get a results table comparing the predictors against the annotations in ClinVar and gnomAD

# Identify 'output' columns
pred_columns = [col for col in vars_pairs_strict.columns if 'output' in col]

# Store our previous obtained results for HomolVar
results = [['HomolVar', 0.95, 0]]

for col in pred_columns:
    # Create a valid mask (ignoring empty values and 'A')
    #valid_mask = (vars_pairs_strict[col] != 'A') & (vars_pairs_strict[col] != '') & (~vars_pairs_strict[col].isna())
    # Count matches where values are equal
    #matches = (vars_pairs_strict[col] == vars_pairs_strict['HomolVar']) & valid_mask
    #proportion = matches.sum() / valid_mask.sum() if valid_mask.sum() > 0 else 0 # removing ambiguous

    # Count matches where values are equal
    matches = (vars_pairs_strict[col] == vars_pairs_strict['HomolVar'])

    # Compute proportion
    proportion = matches.sum() / len(vars_pairs_strict[col]) if len(vars_pairs_strict[col]) > 0 else 0   # wo removing ambiguous
    proportion = round(proportion, 2)

    # Compute missing values proportion for this column
    total_values = vars_pairs_strict[col].size  # Total number of rows in the column
    missing_a_count = ((vars_pairs_strict[col] == 'A') | (vars_pairs_strict[col] == '') | (vars_pairs_strict[col].isna())).sum()
    missing_a_proportion = round(missing_a_count / total_values, 2)

    # Store the result
    results.append([col, proportion, missing_a_proportion])

# Create final DataFrame
proportion_df = pd.DataFrame(results, columns=['Column', 'Match Proportion', 'Ambiguous rate'])

# Sort values
proportion_df = proportion_df.sort_values('Match Proportion', ascending=False)

# Display results
#print(proportion_df)

# Get the table of prediction power
proportion_df.to_csv('./results/FINAL_table_comparisons.txt', sep='\t', index=False)