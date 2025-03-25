import pandas as pd


def mutationTaster(value):
    if pd.isna(value) or value == '':  # Handle NaNs and empty values
        return 'A'
    value = value.replace('.','').split(';')
    unique_vals = ''.join(sorted(set(value)))
    if len(unique_vals) > 1 or len(unique_vals) == 0:
        return 'A'
    if unique_vals == 'D':
        return 'D'
    if unique_vals == 'N':
        return 'T'
    return 'A'

def mutationAssessor(value):
    if pd.isna(value) or value == '':  # Handle NaNs and empty values
        return 'A'
    value = value.replace('.','').split(';')
    unique_vals = ''.join(sorted(set(value)))
    if len(unique_vals) > 1 or len(unique_vals) == 0:
        return 'A'
    if unique_vals == 'H':
        return 'D'
    if unique_vals == 'N':
        return 'T'
    else:
        return 'A'

def PROVEAN(value):
    if pd.isna(value) or value == '':  # Handle NaNs and empty values
        return 'A'
    value = value.replace('.','').split(';')
    unique_vals = ''.join(sorted(set(value)))
    if len(unique_vals) > 1 or len(unique_vals) == 0:
        return 'A'
    if unique_vals == 'D':
        return 'D'
    if unique_vals == 'N':
        return 'T'
    else:
        return 'A'

def VEST4(value):
    if pd.isna(value) or value == '':  # Handle NaNs and empty values
        return 'A'
    values = [v for v in value.split(';') if v != '.']  # Remove '.'
    numeric_values = [float(v) for v in values if v]  # Convert to float, ignore empty
    mean_score = sum(numeric_values) / len(numeric_values) if numeric_values else float('nan')

    if mean_score >= 0.819:
        return 'D'
    elif mean_score <= 0.187:
        return 'T'
    else:
        return 'A'

def MetaRNN(value):
    if pd.isna(value) or value == '':  # Handle NaNs and empty values
        return 'A'
    value = value.replace('.','').split(';')
    unique_vals = ''.join(sorted(set(value)))
    if len(unique_vals) > 1 or len(unique_vals) == 0:
        return 'A'
    if unique_vals == 'D':
        return 'D'
    if unique_vals == 'T':
        return 'T'
    return 'A'

def MCAP(value):
    value = value.replace('.','').split(';')
    if pd.isna(value) or value == '':  # Handle NaNs and empty values
        return 'A'
    unique_vals = ''.join(sorted(set(value)))
    if len(unique_vals) > 1 or len(unique_vals) == 0:
        return 'A'
    if unique_vals == 'D':
        return 'D'
    if unique_vals == 'T':
        return 'T'
    return 'A'

def REVEL(value):
    if pd.isna(value) or value == '':  # Handle NaNs and empty values
        return 'A'
    values = [v for v in value.split(';') if v != '.']  # Remove '.'
    numeric_values = [float(v) for v in values if v]  # Convert to float, ignore empty
    mean_score = sum(numeric_values) / len(numeric_values) if numeric_values else float('nan')

    if mean_score >= 0.682:
        return 'D'
    elif mean_score <= 0.086:
        return 'T'
    else:
        return 'A'

def MVP_gMVP(value):
    if pd.isna(value) or value == '':  # Handle NaNs and empty values
        return 'A'
    values = [v for v in value.split(';') if v != '.']  # Remove '.'
    numeric_values = [float(v) for v in values if v]  # Convert to float, ignore empty
    mean_score = sum(numeric_values) / len(numeric_values) if numeric_values else float('nan')

    if mean_score >= 0.75:
        return 'D'
    elif mean_score < 0.75:
        return 'T'
    return 'A'

def DEOGEN2(value):
    if pd.isna(value) or value == '':  # Handle NaNs and empty values
        return 'A'
    value = value.replace('.','').split(';')
    unique_vals = ''.join(sorted(set(value)))
    if len(unique_vals) > 1 or len(unique_vals) == 0:
        return 'A'
    if unique_vals == 'D':
        return 'D'
    if unique_vals == 'T':
        return 'T'
    return 'A'

def LIST_S2(value):
    if pd.isna(value) or value == '':  # Handle NaNs and empty values
        return 'A'
    value = value.replace('.','').split(';')
    unique_vals = ''.join(sorted(set(value)))
    if len(unique_vals) > 1 or len(unique_vals) == 0:
        return 'A'
    if unique_vals == 'D':
        return 'D'
    if unique_vals == 'T':
        return 'T'
    return 'A'

def VARITY(value):
    if pd.isna(value) or value == '':  # Handle NaNs and empty values
        return 'A'
    values = [v for v in value.split(';') if v != '.']  # Remove '.'
    numeric_values = [float(v) for v in values if v]  # Convert to float, ignore empty
    mean_score = sum(numeric_values) / len(numeric_values) if numeric_values else float('nan')

    if mean_score >= 0.5:
        return 'D'
    elif mean_score < 0.5:
        return 'T'
    return 'A'

def ESM1b(value):
    if pd.isna(value) or value == '':  # Handle NaNs and empty values
        return 'A'
    value = value.replace('.','').split(';')
    unique_vals = ''.join(sorted(set(value)))
    if len(unique_vals) > 1 or len(unique_vals) == 0:
        return 'A'
    if unique_vals == 'D':
        return 'D'
    if unique_vals == 'T':
        return 'T'
    return 'A'


def MutScore(value):
    if pd.isna(value) or value == '':  # Handle NaNs and empty values
        return 'A'
    values = [v for v in value.split(';') if v != '.']  # Remove '.'
    numeric_values = [float(v) for v in values if v]  # Convert to float, ignore empty
    mean_score = sum(numeric_values) / len(numeric_values) if numeric_values else float('nan')

    if mean_score >= 0.730:
        return 'D'
    elif mean_score < 0.140:
        return 'T'
    return 'A'

def CADD(value):
    if pd.isna(value) or value == '':  # Handle NaNs and empty values
        return 'A'
    if value >= 20:
        return 'D'
    elif value < 10:
        return 'T'
    return 'A'

def DANN(value):
    if pd.isna(value) or value == '':  # Handle NaNs and empty values
        return 'A'
    if value >= 0.95:
        return 'D'
    elif value < 0.95:
        return 'T'
    return 'A'

def fathmm_XF_coding(value):
    if pd.isna(value) or value == '':  # Handle NaNs and empty values
        return 'A'
    value = value.replace('.','').split(';')
    unique_vals = ''.join(sorted(set(value)))
    if len(unique_vals) > 1 or len(unique_vals) == 0:
        return 'A'
    if unique_vals == 'D':
        return 'D'
    if unique_vals == 'N':
        return 'T'
    return 'A'

def phastCons(value):
    if pd.isna(value) or value == '' or value == '.':  # Handle NaNs and empty values
        return 'A'
    if isinstance(value, str):
        value = float(value)
    if value >= 0.8:
        return 'D'
    elif value < 0.4:
        return 'T'
    return 'A'

def phyloP(value):
    if pd.isna(value) or value == '' or value == '.':  # Handle NaNs and empty values
        return 'A'
    if isinstance(value, str):
        value = float(value)
    if value > 0:
        return 'D'
    elif value < 0:
        return 'T'
    return 'A'