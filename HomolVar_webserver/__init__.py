import pandas as pd
import os
from markupsafe import escape
from flask import Flask, render_template, request
from flask_sqlalchemy import SQLAlchemy
from sqlalchemy import create_engine
from flask import session, jsonify
from interactive_plot import pos_mapped_plot_interactive

import tools
import checking_tools
from config import *

#APP INIZIALIZATION
app = Flask(__name__)
app.config.from_pyfile('config.py')
db = SQLAlchemy(app)


# Create a SQLAlchemy engine
engine = create_engine(SQLALCHEMY_DATABASE_URI)


@app.route('/' , methods=['GET', 'POST'])
def pathogenicity_predictor():

    ###########
    # Read mapped entries file for displaying the different UniProt entries, NP codes...
    prot_entries = os.path.join(app.root_path, 'data', 'prot_entry_types_mapped.txt')
    prot_entries_df = pd.read_csv(prot_entries, sep='\t')
    prot_entries_df = prot_entries_df.sort_values(by='Gene_name(HGNC)', ascending=True)
    prot_entries_mapped = prot_entries_df[(prot_entries_df['curated_alignment'] == True) | (prot_entries_df['pfam_alignment_seed'] == True) | (prot_entries_df['pfam_alignment_full'] == True)]
    #print(f'Number of proteins that have an alignment: {len(prot_entries_mapped)}')  # 3,519 seed; 18,597 seed + full
    valid_entries_l = pd.concat([prot_entries_mapped['Gene_name(HGNC)'], prot_entries_mapped['Uniprot_entry'], prot_entries_mapped['Uniprot_entry_name'], prot_entries_mapped['RefSeq_NP_code']]).dropna().unique()
    #print('Proteins with alignment displayed', len(valid_entries_l))

    # Get the maximum value of the mutation position from your dataset
    #max_position = variants_df['position'].max()  

    if request.method == "POST":
        gene_or_uniprot = escape(request.form.get("protein"))
        ref_aa = escape(request.form.get("i_aa"))
        mut_aa = escape(request.form.get("f_aa"))
        position = escape(request.form.get("pos"))
        
        if (not position) and (not ref_aa) and (not mut_aa):
            gene_name, entry, entry_name, np_code = tools.get_different_codes(gene_or_uniprot, prot_entries_df)

            if (entry_name is None):
                warn_user = 'Check the entry name (GENE NAME, UNIPROT ENTRY or NCBI REFSEQ), it is not correct.'
                return render_template('/results.html', invalid_aas=warn_user) # I've used the same invalid_aas method in order to avoid adding more conditions to the html file. As this serves.    
            #print(entry_name)
            return entry_align(entry_name)

        ref_aa = str(ref_aa).capitalize()
        mut_aa = str(mut_aa).capitalize()
        pos = checking_tools.check_pos(position)

        ref_aa = checking_tools.get_aa(ref_aa)
        mut_aa = checking_tools.get_aa(mut_aa)

        # If reference or mutated amino acid are None, return invalid input
        if (ref_aa is None) or (mut_aa is None):
            warn_user = 'Please, check the amino acids; at least one of them is incorrect.<br>Alternatively, you can leave the protein variant fields empty to view the positions that can be mapped for this protein.'
            return render_template('/results.html', invalid_aas=warn_user)
        
        elif ref_aa == mut_aa:
            warn_user = 'You have entered a synonymous variant, but this predictor only works for missense variants. Please, try with a missense variant.<br>Alternatively, you can leave the protein variant fields empty to view the positions that can be mapped for this protein.'
            return render_template('/results.html', invalid_aas=warn_user)
        elif (ref_aa == 'fs') or (mut_aa == 'fs'):
            warn_user = 'You have entered a frameshift variant, but this predictor only works for missense variants. Please, try with a missense variant.<br>Alternatively, you can leave the protein variant fields empty to view the positions that can be mapped for this protein.'
            return render_template('/results.html', invalid_aas=warn_user)
        elif (ref_aa == 'Ter') or (mut_aa == 'Ter'):
            warn_user = 'You have entered a truncating variant, but this predictor only works for missense variants. Please, try with a missense variant.<br>Alternatively, you can leave the protein variant fields empty to view the positions that can be mapped for this protein.'
            return render_template('/results.html', invalid_aas=warn_user)

        else:
            prot_var = f'p.{ref_aa}{pos}{mut_aa}'
            gene_name, entry, entry_name, np_code = tools.get_different_codes(gene_or_uniprot, prot_entries_df)
            #print(f'The user wrote: {gene_or_uniprot} and we transformed it into: {entry_name}')

            if (entry_name is None):
                warn_user = 'Check the entry name (GENE NAME, UNIPROT ENTRY or NCBI REFSEQ), it is not correct.'
                return render_template('/results.html', invalid_aas=warn_user) # I've used the same invalid_aas method in order to avoid adding more conditions to the html file. As this serves.
            else: 
                # Check protein type: patho (OMIM) or not patho protein
                form_data = checking_tools.check_protein_type(entry_name, entry, gene_name, prot_var, prot_entries_df)

                # Read variants file
                ## Filtered variants
                variants_wo_mapped_pos_path = os.path.join(app.root_path, 'data', 'homologous_NOT_MAPPED_and_grindb.txt')
                variants_wo_mapped_pos_df = pd.read_csv(variants_wo_mapped_pos_path, sep=',', usecols=['gene', 'clinical_significance', 'prot_change', 'NM', 'initial_aa', 'position', 'final_aa', 'database'])

                
                # Set variables in the session to access them in other functions
                session['gene_name'] = gene_name
                session['entry_name'] = entry_name
                session['ref_aa'] = ref_aa
                session['pos'] = pos
                session['mut_aa'] = mut_aa
                                
                #### DO WE HAVE AN ALIGNMENT FOR THIS PROTEIN?
                #### IS THE VARIANT MATCHING THE POSITION AND REF AMINO ACID FROM THE ALIGNMENT?
                # Check the user's input. We check if we have curated and/or Pfam alignment for this protein and if the reference amino acid
                # given by the user matches with the alignment's one (from UniProt seq).
                warn_user, curated_align_df, pfam_align_df, align_type = checking_tools.check_users_input(entry_name, ref_aa, pos)

                
                #### READ VARIANTS FILE
                # Depending on the align_type, the file of variants will change. That is because the variants that are in 'VARIANTS_seed_homologous_and_grindb.txt' are mixed the variants that can be mapped in seed aligns and CURATED (from GRIN) aligns
                # In the case of FULL alignments, the variants that can be mapped there are in another precomputed file
                
                if align_type == 'CURATED' or 'SEED' in align_type:
                    variants_path = os.path.join(app.root_path, 'data', 'VARIANTS_SEED_homologous_and_grindb.txt')
                    # Read CSV file into a pandas DataFrame
                    variants_df = pd.read_csv(variants_path, sep=',', low_memory=False)
                    variants_df = variants_df[['gene', 'clinical_significance', 'prot_change', 'initial_aa', 'position', 'final_aa', 'Uniprot_entry_name', 'Uniprot_link', 'NM_link', 'Pos_align', 'database', 'Database_link']]
                """ elif 'FULL' in align_type:
                    variants_path = os.path.join(app.root_path, 'data', 'VARIANTS_FULL_align_freq-6_homologous.txt')                     ##### !!!!
                    variants_df = pd.read_csv(variants_path, sep='\t', low_memory=False)
                    variants_df = variants_df[['gene', 'clinical_significance', 'prot_change', 'initial_aa', 'position', 'final_aa', 'Uniprot_entry_name', 'Uniprot_link', 'NM_link', 'Pos_align', 'database', 'Database_link']] """

                #### LET'S GET THE DATA FOR PLOTTING THE ALIGNMENT
                plot_html = False
                # We will prioritize plotting the curated alignment

                if curated_align_df is not None:
                    # Call the function and get the plot div
                    vars_gene_wo_mapped_pos_df = variants_wo_mapped_pos_df[variants_wo_mapped_pos_df['gene'] == gene_name]
                    plot_html = pos_mapped_plot_interactive(curated_align_df, vars_gene_wo_mapped_pos_df, gene_name, entry, align_type)

                elif pfam_align_df is not None:
                    # Call the function and get the plot div
                    vars_gene_wo_mapped_pos_df = variants_wo_mapped_pos_df[variants_wo_mapped_pos_df['gene'] == gene_name]
                    plot_html = pos_mapped_plot_interactive(pfam_align_df, vars_gene_wo_mapped_pos_df, gene_name, entry, align_type)
                else:
                    vars_gene_wo_mapped_pos_df = variants_wo_mapped_pos_df[variants_wo_mapped_pos_df['gene'] == gene_name]
                    plot_html = pos_mapped_plot_interactive(None, vars_gene_wo_mapped_pos_df, gene_name, entry, align_type)


                #### 1.- THE VARIANT IS NOT ANNOTATED YET
                # If the variant is not annotated in ClinVar or gnomAD:
                #if is_annotated is None:
                ###############################################################
                #### 1.1.- THE VARIANT IS OK
                # If all is OK and we don't have to notify anything to the user, we keep on the analysis:
                if warn_user == False:

                    #warn_user_strict_curated = False
                    #warn_user_hom_curated = False
                    #warn_user_strict = False
                    #warn_user_hom = False
                    strict_pairs_df_curated = None
                    homolog_pairs_df_curated = None
                    strict_pairs_df = None
                    homolog_pairs_df = None
                    pairs_df = False

                    if curated_align_df is not None:
                        strict_pairs_df_curated, warn_user_strict_curated = tools.search_for_pairs(entry_name, curated_align_df, variants_df, ref_aa, pos, mut_aa, strict_vs_homologous = 'strict', curated_vs_pfam = 'CURATED')
                        homolog_pairs_df_curated, warn_user_hom_curated = tools.search_for_pairs(entry_name, curated_align_df, variants_df, ref_aa, pos, mut_aa, strict_vs_homologous = 'homologous', curated_vs_pfam = 'CURATED')
                    
                    elif pfam_align_df is not None:
                        strict_pairs_df, warn_user_strict = tools.search_for_pairs(entry_name, pfam_align_df, variants_df, ref_aa, pos, mut_aa, strict_vs_homologous = 'strict', curated_vs_pfam = align_type)
                        homolog_pairs_df, warn_user_hom = tools.search_for_pairs(entry_name, pfam_align_df, variants_df, ref_aa, pos, mut_aa, strict_vs_homologous = 'homologous', curated_vs_pfam = align_type)
                    
                    #strict_pairs_df_joined = None
                    #hom_strict_df_joined = None
                    warn_user = False

                    # If there is a curated df for strict
                    if (strict_pairs_df_curated is not None): # and (strict_pairs_df is not None):
                        #strict_pairs_df_joined = pd.concat([strict_pairs_df_curated, strict_pairs_df], axis = 0)
                        pairs_df = strict_pairs_df_curated
                        
                    # If there is a curated df for homologous
                    if (homolog_pairs_df_curated is not None): # and (homolog_pairs_df is not None):
                        #hom_strict_df_joined = pd.concat([homolog_pairs_df_curated, homolog_pairs_df], axis = 0)
                        pairs_df = homolog_pairs_df_curated
                        
                    # But if they exist strict and homologous, pairs_df is the concat of both of them
                    #if (strict_pairs_df_joined is not None) and (hom_strict_df_joined is not None):
                        #pairs_df = pd.concat([strict_pairs_df_joined, hom_strict_df_joined], axis = 0)
                    if (strict_pairs_df_curated is not None) and (homolog_pairs_df_curated is not None):
                        pairs_df = pd.concat([strict_pairs_df_curated, homolog_pairs_df_curated], axis = 0)
                    
                    if type(pairs_df) == bool:
                        # If there is no a curated align...
                        #if (warn_user_strict == False) and (warn_user_hom == False):
                        if (strict_pairs_df is not None) and (homolog_pairs_df is not None):
                            pairs_df = pd.concat([strict_pairs_df, homolog_pairs_df], axis = 0)
                        #elif warn_user_strict == False:
                        elif strict_pairs_df is not None:
                            pairs_df = strict_pairs_df
                        elif homolog_pairs_df is not None:
                            pairs_df = homolog_pairs_df
                        else:
                            warn_user = 'No strict nor homologous pairs found, i.e., no variants in equivalent position in homologous variants with exact or similar mutated amino acid were found.'

                    # Check if we have pairs. If so, we retrieve the clinical consequence for the variant
                    if warn_user == False:

                        pairs_df = checking_tools.format_pairs_df(pairs_df, prot_var, entry)

                        
                        
                        clin_conseq, final_pairs_df = checking_tools.check_consequence_of_pairs(pairs_df)
                        # Convert the DataFrame to a dictionary
                        pairs_dict = final_pairs_df.to_dict(orient='records')

                        return render_template('/results.html', form_data=form_data, clin_conseq=clin_conseq, pairs_found=pairs_dict, plot_div=plot_html, align_type=align_type)
                    else:
                        return render_template('/results.html', no_pairs_found=warn_user, form_data=form_data, plot_div=plot_html, align_type=align_type)

                #### 1.2.- THE VARIANT IS NOT OK: NOT MATCHING REFERENCE AMINO ACID AND POSITION
                else:
                    if warn_user == 'no_alignment':
                        return render_template('/results.html', form_data=form_data, plot_div=plot_html, align_type=None)
                    return render_template('/results.html', error=warn_user, form_data=form_data, plot_div=plot_html, align_type=None)
                

    # Initialize form_data as empty if it's a GET request
    form_data = {
        'output_user': '',
        'protein': '',
        'prot_var': ''
        }

    return render_template('/user_parameters.html', error=None, form_data=form_data, valid_entries_l=valid_entries_l) # , max_position=max_position)


@app.route('/get-annotation-data')
def get_annotation_data():
    ## Read filtered variants
    variants_wo_mapped_pos_path = os.path.join(app.root_path, 'data', 'homologous_NOT_MAPPED_and_grindb.txt')
    variants_wo_mapped_pos_df = pd.read_csv(variants_wo_mapped_pos_path, sep=',', usecols=['gene', 'clinical_significance', 'prot_change', 'NM', 'initial_aa', 'position', 'final_aa', 'database', 'rsid', 'VariationID'], low_memory=False)
    
    # Access the variables from flask session
    gene_name = session.get('gene_name')
    entry_name = session.get('entry_name')
    ref_aa = session.get('ref_aa')
    pos = session.get('pos')
    mut_aa = session.get('mut_aa')


    #### IS THE VARIANT ALREADY ANNOTATED?
    # Check if the protein variant is already annotated in ClinVar (pathogenic) or gnomAD (neutral)
    is_annotated, existing_annot_df = tools.search_var_annotation(variants_wo_mapped_pos_df, gene_name, ref_aa, pos, mut_aa)
    
    if is_annotated == False:
        is_annotated, existing_annot_df = tools.search_annot_clinvar_gnomad(gene_name, entry_name, ref_aa, pos, mut_aa)

    ## format the output df of existing annotations
    exist_annot_df_curated = tools.format_exist_annot(existing_annot_df, app)

    annot_vars_dict = None
    if exist_annot_df_curated is not None:
        annot_vars_dict = exist_annot_df_curated.to_dict(orient='records')
        #print(annot_vars_dict)

    if annot_vars_dict:  # if data is available
        response = {
            "data": annot_vars_dict  # List of dictionaries, each dictionary is a row for the table
        }
    else:  # if no data available
        response = {
            "message": "This variant was not found annotated in ClinVar nor gnomAD."
        }
    #print(response)
    return jsonify(response)



@app.route('/about')
def about():
    return render_template('about.html')


""" @app.route('/homologous_vars_db/')
def homologous_vars_db():
    # Use the db object to interact with your database
    func_vars = db.session.execute(text("SELECT * FROM homologous_variants_db"))
    rows_hom_vars = [dict(zip(func_vars.keys(), row)) for row in func_vars]  # Convert rows to dictionaries

    return render_template('database_homologous_vars.html', rows_hom_vars=rows_hom_vars)
 """

@app.route('/<entry>')
def entry_align(entry):
    """
    'entry' and 'code_human' will be two UniProt codes:
    e.g. Entry will be code 'Q05586' and code_human will be the code 'NMDZ1_HUMAN'
    """
    prot_entries = os.path.join(app.root_path, 'data', 'prot_entry_types_mapped.txt')
    prot_entries_df = pd.read_csv(prot_entries, sep='\t')
    prot_entries_df = prot_entries_df.sort_values(by='Gene_name(HGNC)', ascending=True)

    gene_name, entry, entry_name, np_code = tools.get_different_codes(entry, prot_entries_df)
    
    variants_wo_mapped_pos_df = pd.read_csv('/var/www/homologous/data/homologous_NOT_MAPPED_and_grindb.txt', sep=',', usecols=['gene', 'clinical_significance', 'position', 'database', 'prot_change'])
    variants_wo_mapped_pos_df = variants_wo_mapped_pos_df[variants_wo_mapped_pos_df['gene'] == gene_name]
    form_data = checking_tools.check_protein_type(entry_name, entry, gene_name, None, prot_entries_df)

    try:
        # Read the CSV file based on the entry. If a curated align exists, keep the curated. Else, the Pfam
        curated_align_path = f'{app.root_path}/data/curated_aligns_by_prots/{entry_name}_align_CURATED.txt'
        seed_align_path = f'{app.root_path}/data/seed_align_by_prots/{entry_name}_align_seed_pfam.txt'
        #full_align_path = f'{app.root_path}/data/full_align_by_prots/{entry_name}_align_full_pfam.txt'        ##### !!!!
        align_type = None
        prot_align = None
        if os.path.exists(curated_align_path):
            prot_align = pd.read_csv(curated_align_path, sep=',')
            align_type = 'CURATED'
        elif os.path.exists(seed_align_path):
            prot_align = pd.read_csv(f'{app.root_path}/data/seed_align_by_prots/{entry_name}_align_seed_pfam.txt', sep=',')
            align_type = 'Pfam SEED'
        """ elif os.path.exists(full_align_path):                                                                                   ##### !!!!
            prot_align = pd.read_csv(f'{app.root_path}/data/full_align_by_prots/{entry_name}_align_full_pfam.txt', sep=',')
            align_type = 'Pfam FULL' """
        
        if align_type != None:
            prot_align['eq_pos_align'] = prot_align['eq_pos_align'].astype(object)
            pd.set_option('future.no_silent_downcasting', True)
            prot_align['eq_pos_align'] = prot_align['eq_pos_align'].fillna(-99).astype(int)
            prot_align['eq_pos_align'] = prot_align['eq_pos_align'].apply(lambda x: int(float(x)) if x != -99 else x)

        # Convert DataFrame to list of lists
        """ data = prot_align.values.tolist()
        columns = prot_align.columns.tolist() """
        
        ##############

        # Call the function and get the plot div
        plot_html = pos_mapped_plot_interactive(prot_align, variants_wo_mapped_pos_df, gene_name, entry, align_type)

        if plot_html is None:
            error = f'There was no proper data for {entry} protein... We recommend you to check this protein in UniProt: <a href="https://www.uniprot.org/uniprotkb/{entry}/entry" target="_blank"> {entry}</a>'
            return render_template('protein_align.html', form_data=form_data, error=error)

        return render_template('protein_align.html', form_data=form_data, plot_div=plot_html, align_type=align_type)  # columns=columns, data=data
    
    except FileNotFoundError:
        plot_html = pos_mapped_plot_interactive(None, variants_wo_mapped_pos_df, gene_name, entry, None)

        if plot_html is None:
            error = f'There was no proper data for {entry} protein... We recommend you to check this protein in UniProt: <a href="https://www.uniprot.org/uniprotkb/{entry}/entry" target="_blank"> {entry}</a>'
            return render_template('protein_align.html', form_data=form_data, error=error)
        else:
            error = f'No alignment was found for {entry} ({gene_name} gene).'
        return render_template('protein_align.html', form_data=form_data, plot_div=plot_html, error=error, align_type=align_type)

        #return render_template('protein_align.html', form_data=form_data, error=error)

