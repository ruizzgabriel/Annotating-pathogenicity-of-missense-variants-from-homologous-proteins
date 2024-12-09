from Bio.Data.IUPACData import protein_letters_3to1, protein_letters_1to3
import plotly.graph_objects as go
import pandas as pd
import numpy as np
import plotly.offline as pyo  # Import Plotly offline mode
import json
import csv
import os
import re

def extract_hom_vars(align_df, fig):

    prot_entries_df = pd.read_csv('/var/www/homologous/data/prot_entry_types_mapped.txt', sep='\t')
    # Filter rows with data in 'hom_vars'
    if align_df is not None:
        hom_vars_df = align_df[align_df['hom_vars'] != '-']

        patho_pos_l = []
        NOpatho_pos_l = []
        reviewd_data_patho = []
        reviewd_data_NOpatho = []

        for index, row in hom_vars_df.iterrows():
            hom_var_data_patho = ''
            hom_var_data_NOpatho = ''
            vars_data = row['hom_vars']
            prot_pos = int(row['prot_pos'])

            if vars_data != '-':
                # Remove parentheses and their content
                print(vars_data)
                vars_data = re.sub(r'\(.*?\)', '', vars_data)
                print(vars_data)
                vars_data = vars_data.replace('[', '').replace(']', '').replace("'", '').split(',')
                for variant in vars_data:
                    print(variant)
                    variant = variant.replace(' ', '').replace('\n', '')
                    hom_prot = variant.split(';')[0]
                    gene = prot_entries_df.loc[prot_entries_df['Uniprot_entry_name'] == hom_prot, 'Gene_name(HGNC)'].values
                    #print(gene)
                    #aa_ref = protein_letters_1to3[aa_ref]
                    #mut_aa = variant.split(';')[1]
                    db = variant.split(';')[2]
                    prot_var = variant.split(';')[3]
                    annot = variant.split(';')[4]
                    if annot != 'Neutral' and annot != 'Pathogenic':
                        if db == 'ClinVar':
                            annot = 'Pathogenic'
                        elif db == 'gnomAD':
                            annot = 'Neutral'

                    if annot == 'Neutral':
                        if prot_pos not in NOpatho_pos_l:
                            NOpatho_pos_l.append(prot_pos)
                            hom_var_data_NOpatho = f"Position: {prot_pos}<br>Annotation: {annot}<br><br>--------- HOMOLOGOUS VARIANTS ---------<br>"
                        hom_var_data_NOpatho += f"{gene[0]} variant: {prot_var}<br>"
                    else:
                        if prot_pos not in patho_pos_l:
                            patho_pos_l.append(prot_pos)
                            hom_var_data_patho = f"Position: {prot_pos}<br>Annotation: {annot}<br><br>--------- HOMOLOGOUS VARIANTS ---------<br>"
                        hom_var_data_patho += f"{gene[0]} variant: {prot_var}<br>"
                
                if hom_var_data_patho:
                    reviewd_data_patho.append(hom_var_data_patho)
                if hom_var_data_NOpatho:
                    reviewd_data_NOpatho.append(hom_var_data_NOpatho)



        """ # Extract x and text data for plotting
        #hom_pos = hom_vars_df['prot_pos']
        pos_l = list(hom_vars_df['prot_pos'].dropna().astype(int).tolist())

        #pos_l = hom_vars_df[hom_vars_df['hom_vars'] != '-', 'prot_pos'].dropna().astype(int).tolist()
        ref_aa = list(hom_vars_df['residue'].dropna())
        #ref_aa = hom_vars_df[hom_vars_df['hom_vars'] != '-', 'residue'].dropna().tolist()
        hom_data = list(hom_vars_df.loc[hom_vars_df['prot_pos'].isin(pos_l), 'hom_vars'])

        reviewd_data_patho = []
        reviewd_data_NOpatho = []

        for prot_pos, aa_ref, position_data in zip(pos_l, ref_aa, hom_data):
            #hom_var_review_data = ['----------------------------']
            hom_var_review_data_patho = ''
            hom_var_review_data_NOpatho = ''
            if position_data != '-':
                position_data = position_data.replace('[', '').replace(']', '').replace("'", '').split(',')
                for hom_var in position_data:
                    hom_var = hom_var.replace(' ', '')
                    hom_prot = hom_var.split(';')[0]
                    prot_code = prot_entries_df.loc[prot_entries_df['Uniprot_entry_name'] == hom_prot, 'Gene_name(HGNC)'].values
                    #print(prot_code)
                    #aa_ref = protein_letters_1to3[aa_ref]
                    mut_aa = hom_var.split(';')[1]
                    db = hom_var.split(';')[2]
                    prot_var = hom_var.split(';')[3]
                    annot = hom_var.split(';')[4]
                    if annot == 'Neutral':
                        annot = 'Non-pathogenic'
                        hom_var_review_data_patho += f"Position: {prot_pos}<br>Gene: {prot_code[0]}<br>Protein variant: {prot_var}<br>Annotation: {annot}<br>---------------------------------------<br>"
                    else:
                        hom_var_review_data_NOpatho += f"Position: {prot_pos}<br>Gene: {prot_code[0]}<br>Protein variant: {prot_var}<br>Annotation: {annot}<br>---------------------------------------<br>"

                reviewd_data_patho.append(hom_var_review_data_patho)
                reviewd_data_NOpatho.append(hom_var_review_data_NOpatho) """
        

        #print(len(patho_pos_l), len(reviewd_data_patho))
        #print(len(NOpatho_pos_l), len(reviewd_data_NOpatho))


        #print(NOpatho_pos_l)
        #print(reviewd_data)
        # Plot triangles at specified positions with hover text
        fig.add_trace(go.Scatter(
            x=NOpatho_pos_l,  # x-positions from 'position' column
            y=[0.8] * len(NOpatho_pos_l),  # y-position set at 5
            mode='markers',
            name='Non-pathogenic homologous variants',
            marker=dict(color='DarkGreen', symbol='triangle-up', size=10),  # Triangles
            hoverinfo='text',
            text=reviewd_data_NOpatho,
        ))

        # Plot triangles at specified positions with hover text
        fig.add_trace(go.Scatter(
            x=patho_pos_l,  # x-positions from 'position' column
            y=[0.9] * len(patho_pos_l),  # y-position set at 5
            mode='markers',
            name='Pathogenic homologous variants',
            marker=dict(color='FireBrick', symbol='triangle-up', size=10),  # Triangles
            hoverinfo='text',
            text=reviewd_data_patho,
        ))

        # Show the plot
        fig.show()
        return
    
    else:
        return None


def pos_mapped_plot_interactive(align_df, vars_df, gene_name, entry, align_type):

    with open('/var/www/homologous/data/entry_prot_len_mapping.json', 'r') as json_file:
        prot_lengths_d = json.load(json_file)

    # Get the length of the protein
    if entry in prot_lengths_d:
        max_pos = prot_lengths_d[entry]
    else:
        return None
    
    if max_pos > 0:
        # Create positions array
        positions = np.arange(1, max_pos + 1)  # Ensure it covers all positions up to max_pos

        #
        # Define the file path
        file_path = "/var/www/homologous/data/prot_entry_types_mapped.txt"

        # Create a dictionary to map Uniprot_entry to Uniprot_entry_name
        entry_to_name = {}

        # Read the file
        with open(file_path, mode='r') as file:
            reader = csv.DictReader(file, delimiter='\t')
            for row in reader:
                # Fill the dictionary with Uniprot_entry -> Uniprot_entry_name mapping
                entry_to_name[row["Uniprot_entry"]] = row["Uniprot_entry_name"]

        entry_name = entry_to_name.get(entry, "Not Found")

        # Define folders
        curated_folder = "/var/www/homologous/data/curated_aligns_by_prots"
        seed_folder = "/var/www/homologous/data/seed_align_by_prots"
        full_folder = "/var/www/homologous/data/full_align_by_prots"

        # Variable to store the protein sequence
        protein_sequence = None

        # Helper function to read the file and extract the residue column as a protein sequence
        def extract_sequence(file_path):
            try:
                # Read the file into a DataFrame
                df = pd.read_csv(file_path, delimiter=',')
                # Extract the 'residue' column and join it into a sequence
                sequence = ''.join(df['residue'].astype(str))
                return sequence
            except Exception as e:
                print(f"Error processing file {file_path}: {e}")
                return None

        # Check for the file in curated_aligns_by_prots
        curated_file = os.path.join(curated_folder, f"{entry_name}_align_CURATED.txt")
        seed_file = os.path.join(seed_folder, f"{entry_name}_align_seed_pfam.txt")
        full_file = os.path.join(full_folder, f"{entry_name}_align_full_pfam.txt")

        if os.path.exists(curated_file):
            # Extract sequence from curated alignment file
            protein_sequence = extract_sequence(curated_file)
        elif os.path.exists(seed_file):
            # Extract sequence from seed alignment file
            protein_sequence = extract_sequence(seed_file)
        elif os.path.exists(full_file):
            # Extract sequence from full alignment file
            protein_sequence = extract_sequence(full_file)
        
        else:
            print(f"Cannot find a sequence for protein {entry_name}.")


        # Create the plot using Plotly
        fig = go.Figure()

        # Assuming that protein_sequence is already defined and has been successfully extracted
        if protein_sequence is None:
            print("Error: protein_sequence is None. Please check the file extraction process.")
            return None  # Exit early if protein_sequence is not found

        # Make sure protein_sequence is not empty
        if len(protein_sequence) < max_pos:
            print(f"Warning: protein_sequence is shorter than the max position ({max_pos}).")
            return None  # Exit early if protein_sequence is too short

        fig.add_trace(go.Scatter(
            x=np.arange(1, max_pos + 1),  # Generate all positions as x-values
            y=[0.15] * max_pos,  # Maintain the y-value for the black line
            mode='lines',
            name='Protein length',
            line=dict(color='Black', width=15),  # Set color to black
            hoverinfo='text',  # Show custom text on hover
            text=[f"Position: {pos}, {protein_sequence[pos-1]}" for pos in range(1, max_pos + 1)],  # Position and AA from protein_sequence
        ))


        # If we have an alignment...
        if align_df is not None:
            align_df = align_df.fillna('-')
            # Map '-' to 0 and others to 1
            try:
                mapped_pos = align_df['curated_align'].apply(lambda x: 0 if x == '-' else 0.5).values
                #align_type = 'CURATED*'
            except KeyError:
                mapped_pos = align_df['pfam_align'].apply(lambda x: 0 if x == '-' else 0.5).values
                #align_type = 'Pfam SEED*'
                if align_type == 'Pfam SEED':
                    align_type = 'Pfam SEED*'
                elif align_type == 'Pfam FULL':
                    align_type = 'Pfam FULL*'

            # Convert to DataFrame for rolling window
            pos_mapped_df = pd.DataFrame({'positions': positions, 'mapped_pos': mapped_pos})

            # Set y-values to 0 or 1 based on the mapped position
            #pos_mapped_df['smoothed_align'] = pos_mapped_df['mapped_pos'].rolling(window=10, min_periods=1).max()  # Use max instead of mean to keep it binary

            # Replace values <= 0.1 with None to hide them
            pos_mapped_df['mapped_pos'] = pos_mapped_df['mapped_pos'].apply(lambda x: x if x > 0.1 else None)


            y_values = np.where(pos_mapped_df['mapped_pos'].notna(), 0.3, np.nan)

            # Extract sequence from alignment (same logic as black line)
            def extract_sequence_from_alignment(file_path):
                try:
                    # Read the file into a DataFrame
                    df = pd.read_csv(file_path, delimiter=',')
                    # Extract the 'align_res' column and join it into a sequence, handling empty values
                    sequence = ''.join(df['align_res'].astype(str).fillna('-').replace('nan', '-'))  # Replace NaN with '-' if needed
                    return sequence
                except Exception as e:
                    print(f"Error processing file {file_path}: {e}")
                    return None
                

            protein_sequence = None
            if os.path.exists(curated_file):
            # Extract sequence from curated alignment file
                protein_sequence = extract_sequence_from_alignment(curated_file)
            elif os.path.exists(seed_file):
                # Extract sequence from seed alignment file
                protein_sequence = extract_sequence_from_alignment(seed_file)
            elif os.path.exists(full_file):
                # Extract sequence from full alignment file
                protein_sequence = extract_sequence_from_alignment(full_file)

            else:
                print(f"Cannot find a sequence for protein {entry_name}.")
                # Add line for mapped regions with only x hoverinfo

            fig.add_trace(go.Scatter(
                x=pos_mapped_df['positions'],
                y=y_values,
                mode='lines',
                name= f'Available alignment ({align_type})',
                line=dict(color='DodgerBlue', width=15),
                hoverinfo='text',  # Show custom text on hover
                text=[f"Position: {pos}, {protein_sequence[pos-1]}" for pos in range(1, max_pos + 1)],  # Position and AA from protein_sequence

            ))
        

        # Extract the positions for pathogenic and non-pathogenic variants
        #patho_vars = vars_df[vars_df['database'] == 'ClinVar']
        #pos_patho_vars = list(patho_vars['position'].dropna().astype(int).tolist())
        #print(len(vars_df))
        #print(vars_df[['gene', 'clinical_significance', 'database']])
        patho_terms = ['Likely pathogenic', 'Likely risk allele', 'Pathogenic', 'Pathogenic/Likely pathogenic', 'Pathogenic/Likely risk allele', 'Pathogenic; risk factor']
        patho_vars = vars_df[vars_df['clinical_significance'].str.strip().isin(patho_terms)]
        pos_patho_vars = list(patho_vars['position'].dropna().astype(int).tolist())
        
        # Extract 'prot_change' for pathogenic variants
        patho_prot_change = list(patho_vars.loc[patho_vars['position'].isin(pos_patho_vars), 'prot_change'])

        #nonpatho_vars = vars_df[vars_df['database'] == 'gnomAD']
        #pos_nonpatho_vars = list(nonpatho_vars['position'].dropna().astype(int).tolist())
        nonpatho_vars = vars_df[vars_df['clinical_significance'].str.strip().isin(['Neutral', 'Benign', 'Benign/Likely benign', 'Likely benign', 'Neutral'])]
        pos_nonpatho_vars = list(nonpatho_vars['position'].dropna().astype(int).tolist())
        
        # Extract 'prot_change' for non-pathogenic variants
        nonpatho_prot_change = list(nonpatho_vars.loc[nonpatho_vars['position'].isin(pos_nonpatho_vars), 'prot_change'])

        #print(f'patho: {len(patho_vars)}   non-patho: {len(nonpatho_vars)}')
        #print(f'patho: {patho_vars}   non-patho: {nonpatho_vars}')

        # Add scatter points for pathogenic variants (red) with hoverinfo
        fig.add_trace(go.Scatter(
            x=pos_patho_vars,
            y=[0.6] * len(pos_patho_vars),
            mode='markers',
            name='Pathogenic variants',
            marker=dict(color='red'),
            hoverinfo='text',
            text=[f"Position: {pos}<br>Protein Change: {prot}" for pos, prot in zip(pos_patho_vars, patho_prot_change)],
            showlegend=False,  # Do not show legend for individual points
            legendgroup='pathogenic'  # Group with other patho variants
        ))

        # Add scatter points for non-pathogenic variants (green) with hoverinfo
        fig.add_trace(go.Scatter(
            x=pos_nonpatho_vars,
            y=[0.5] * len(pos_nonpatho_vars),
            mode='markers',
            name='Non-pathogenic variants',
            marker=dict(color='green'),
            hoverinfo='text',
            text=[f"Position: {pos}<br>Protein Change: {prot}" for pos, prot in zip(pos_nonpatho_vars, nonpatho_prot_change)],
            showlegend=False,  # Do not show legend for individual points
            legendgroup='non-pathogenic'  # Group with other non-patho variants
        ))

        # Add dummy traces to represent the legend entries
        fig.add_trace(go.Scatter(
            x=[None],
            y=[None],
            mode='markers',
            name='Non-pathogenic variants',
            marker=dict(color='green'),
            showlegend=True,  # Show legend for grouped non-patho variants
            hoverinfo='none',  # No hover info for the legend entry
            legendgroup='non-pathogenic'  # Group with other non-patho variants
        ))
        fig.add_trace(go.Scatter(
            x=[None],
            y=[None],
            mode='markers',
            name='Pathogenic variants',
            marker=dict(color='red'),
            showlegend=True,  # Show legend for grouped patho variants
            hoverinfo='none',  # No hover info for the legend entry
            legendgroup='pathogenic'  # Group with other patho variants
        ))


        fig.update_layout(
            #title=f"Positions mapped for {gene_name} in a {align_type} alignment",
            #title=f"Available Pathogenic and Non-pathogenic variants for {gene_name}",
            xaxis_title="Protein position",
            yaxis_title=" ",
            yaxis=dict(showticklabels=False, range=[0, 1], showgrid=False),  # Disable horizontal grid lines
            xaxis=dict(range=[-10, max_pos+10], constrain='domain', gridcolor='lightgray'),  # Enforce x-axis range and set vertical grid color
            legend_title="Features",
            hovermode='closest',
            width=1200,
            height=300,
            margin=dict(l=10, r=10, t=40, b=30),  # Adjust left, right, top, and bottom margins
            legend=dict(
                traceorder='reversed',  # 'normal' for the order they are added
                itemclick='toggle'    # Allows clicking to toggle visibility
            )
        )

        extract_hom_vars(align_df, fig)

        # Generate the plot as a div
        plot_div = pyo.plot(fig, include_plotlyjs=False, output_type='div', config={'responsive': True})

        return plot_div  # Return the div for rendering in your server
    
    return None

"""# Usage
#uniprot_human = 'ANS1A'  # uniprot code with _HUMAN but without it
#gene_name = 'ANKS1A'
uniprot_human = 'NMDZ1'  # uniprot code with _HUMAN but without it
gene_name = 'GRIN1'

align_df = pd.read_csv(f'./{uniprot_human}_HUMAN_align_seed_pfam.txt', sep=',')
vars_df = pd.read_csv('./homologous_NOT_MAPPED_and_grindb.txt', sep=',')
vars_df = vars_df[vars_df['gene'] == gene_name]

# Call the function and get the plot div
output_file = pos_mapped_plot_interactive(align_df, vars_df, gene_name, output_file='output.html')

# In your Flask or other web framework code, you would render this plot_html in your template"""