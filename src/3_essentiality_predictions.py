#!/usr/bin/env python3

import cobra
import pandas as pd
import media
import math
import os
from alive_progress import alive_bar
pd.options.mode.chained_assignment = None

__author__ = "Matthew Onorato"
__version__ = "0.0.4"
__date__ = "June 24, 2020"
__maintainer__ = "Matthew Onorato"
__email__ = "monorato2@gmail.com" > "monorato4398@sdsu.edu"
__status__ = "final"


"""
DESCRIPTION: This script simulates single gene deletions in a metabolic model and compares the results to experimental
             data, i.e. the corresponding transposon mutagenesis dataset. Uses pFBA to calculate optimal growth and
             FBA to simulate knock-outs.
"""


def sim_vs_exp(mtb_model, obj_fx, exp_data, in_vivo_medium, ko_method, file_name, sheetname):
    """
    Simulates single gene deletions for a JSON model and compares simulated vs. experimental gene essentiality calls
    :param mtb_model: JSON model
    :param obj_fx: model objective function
    :param exp_data: column name from essentiality_dataset_comparison.xlsx
    :param in_vivo_medium: from media.py
    :param ko_method: please use either 'fba', 'linear moma', or 'linear room'
    :param file_name: name of file to be written
    :param sheetname: name of Excel sheet to be written
    :return: writes dataframe with simulated, experimental, and comparative results as well as overall accuracy and MCC
    """

    model = cobra.io.load_json_model(mtb_model)
    model.objective = obj_fx
    model.reactions.get_by_id(obj_fx).upper_bound = 1000

    # clears model of active exchange reactions (i.e. removes 'preset' media) and adds provided media
    for reaction in model.reactions:
        if 'EX_' in reaction.id:
            reaction.lower_bound = 0

    if exp_data == 'Griffin_2011':
        mtb_medium = media.griffin  # alternatives: iEK_griffinEssen
    elif exp_data == 'Sassetti_2003' or exp_data == 'Zhang_2012':
        mtb_medium = media.m7H10  # m7H9/m7H10/sMtb_CSM (same result), iEK_m7H10
    elif exp_data == 'DeJesus_2017' or exp_data == 'Xu_2017':
        mtb_medium = media.dejesus  # iEK_deJesusEssen
    elif exp_data == 'Minato_2019':
        mtb_medium = media.ym_rich
    else:
        mtb_medium = in_vivo_medium

    for rxn_name, uptake_rate in mtb_medium.items():
        try:
            model.reactions.get_by_id(rxn_name).lower_bound = uptake_rate
        except KeyError:
            # print(f'{rxn_name} not in {model.id}. Please add reaction to model, if appropriate.')
            continue

    if exp_data in vivo_datasets:
        if in_vivo_medium == media.sMtb_in_vivo_mod or in_vivo_medium == media.sMtb_in_vivo:
            if 'predictions' in file_name:
                model.reactions.get_by_id('EX_o2_e').lower_bound = -1.50  # iEK1008: needed for growth in sMtb_in_vivo

    try:
        pfba_model_flux = cobra.flux_analysis.pfba(model)  # model flux with pFBA constraints added
        opt_pfba_growth = pfba_model_flux.fluxes[obj_fx]  # biomass reaction flux w/ pFBA constraints
    except cobra.util.OptimizationError:
        print('Model cannot grow in the provided media. As a result, gene essentiality predictions cannot be made.')
        return

    if opt_pfba_growth is None or opt_pfba_growth < 0.0001:
        print('Model cannot grow in the provided media. As a result, gene essentiality predictions cannot be made.')
        return

    grw_thrs = 0.05 * opt_pfba_growth  # Agudelo et al. 2020 binary classification used

    # performs simulated single gene deletions in appropriate media
    if ko_method == 'moma':
        df_sim = cobra.flux_analysis.single_gene_deletion(model, method='linear moma', solution=pfba_model_flux)
    elif ko_method == 'room':
        df_sim = cobra.flux_analysis.single_gene_deletion(model, method='linear room', solution=pfba_model_flux)
    else:
        df_sim = cobra.flux_analysis.single_gene_deletion(model, method='fba')

    df_sim.rename(index=lambda x: list(x)[0], inplace=True)

    df_sim['growth'] = df_sim['growth'].fillna(0)  # NA in 'growth' if solution is infeasible; infeasible = no grow = 0

    # joins dataframe with experimental results to dataframe with simulated results
    exp_path = '../results/essentiality_dataset_comparison.xlsx'

    if exp_data in vitro_datasets:
        df_exp = pd.read_excel(exp_path, sheet_name='vitro_all_calls', index_col='ids')
    else:
        df_exp = pd.read_excel(exp_path, sheet_name='vivo_all_calls', index_col='ids')

    df_sim_exp = df_sim.join(df_exp[exp_data], how='left')
    df_sim_exp['classification'] = ''

    # excludes genes only in one, blocked reaction (always called NE due to model limitation and not an accurate call)
    genes_to_exclude = []
    blocked_rxns = cobra.flux_analysis.variability.find_blocked_reactions(model, open_exchanges=True)

    for rxn_id in blocked_rxns:
        rxn_gpr = model.reactions.get_by_id(rxn_id).gene_reaction_rule
        if rxn_gpr.count('Rv') == 1:
            if len(model.genes.get_by_id(rxn_gpr).reactions) == 1:
                genes_to_exclude.append(rxn_gpr)

    for index, row in df_sim_exp.iterrows():
        if index in genes_to_exclude:
            df_sim_exp.loc[index, 'classification'] = 'NA'
        elif row['growth'] > grw_thrs and row[exp_data] == 'NE':  # no >= since it would include 0
            df_sim_exp.loc[index, 'classification'] = 'TP'
        elif row['growth'] < grw_thrs and row[exp_data] == 'ES':
            df_sim_exp.loc[index, 'classification'] = 'TN'
        elif row['growth'] >= grw_thrs and row[exp_data] == 'ES':
            df_sim_exp.loc[index, 'classification'] = 'FP'
        elif row['growth'] <= grw_thrs and row[exp_data] == 'NE':
            df_sim_exp.loc[index, 'classification'] = 'FN'
        else:
            df_sim_exp.loc[index, 'classification'] = 'NA'  # for genes without an essential call in original dataset

    counts = df_sim_exp.groupby(['classification'])['classification'].count()

    df_sim_exp['NA'] = ''  # sum of genes only in one blocked reaction and genes without an essential call
    df_sim_exp['NA'][0] = counts['NA']

    try:
        total = counts['TP'] + counts['TN'] + counts['FP'] + counts['FN']  # + counts['NA']
    except KeyError:
        print(f'{model.id} could not make any predictions. This is likely due to the model not being able to grow.')
        return

    df_sim_exp['accuracy'] = ''  # adding empty column and then filling first row only
    df_sim_exp['accuracy'][0] = (counts['TP'] + counts['TN']) / total * 100

    # MCC does not include 'NA' and accounts for different model sizes --> better way to compare model accuracies
    mcc_root = math.sqrt((counts['TP'] + counts['FP']) * (counts['TP'] + counts['FN']) *
                         (counts['TN'] + counts['FP']) * (counts['TN'] + counts['FN']))
    df_sim_exp['MCC'] = ''
    df_sim_exp['MCC'][0] = ((counts['TP'] * counts['TN']) - (counts['FP'] * counts['FN'])) / mcc_root

    file_path = '../results/'

    if file_name not in os.listdir(file_path):
        filewriter = pd.ExcelWriter(file_path + file_name, mode='w')
    else:
        filewriter = pd.ExcelWriter(file_path + file_name, mode='a')

    df_sim_exp.to_excel(filewriter, sheet_name=sheetname)
    filewriter.save()

    model.reactions.get_by_id(obj_fx).upper_bound = 0  # just in case

    return


model_files = [file for file in os.listdir('../../3_base_model/models/') if file.endswith('.json')]

vitro_datasets = ['Sassetti_2003', 'Griffin_2011', 'Zhang_2012', 'DeJesus_2017', 'Xu_2017', 'Minato_2019']

vivo_datasets = ['Sassetti_vivo_2003', 'Zhang_2013', 'ARTIST_2014',
                 'Mendum_Day3_2015', 'Mendum_Day7_2015', 'Subramaniyam_2019']


# calculates gene essentiality accuracies for different biomass reactions using DeJesus dataset & dejesus media
# accuracies will be used as one of many metrics in which to select a single biomass reaction for further analyses
in_vitro_biomass_rxns = ['biomass_iEK1008_60atp', 'biomass_iEK1011_60atp', 'biomass_iNJ661_60atp',
                         'biomass_sMtb_CSM_57atp', 'biomass_sMtb_REB_57atp', 'biomass_sMtbRECON_vitro_57atp']

with alive_bar(len(in_vitro_biomass_rxns), bar='blocks', spinner='notes_scrolling') as bar:
    for model_file in model_files:
        if 'iMtb' in model_file:
            model_path = f'../../3_base_model/models/{model_file}'
            model_name = model_file.split('.')[0]
            combo_name = f'{model_name}_vitro_combinations.xlsx'

            for bm_rxn in in_vitro_biomass_rxns:
                sim_vs_exp(model_path, bm_rxn, 'DeJesus_2017', None, 'fba', combo_name, f'{bm_rxn}')
                bar()


# calculates gene essentiality accuracies for different biomass reactions using ARTIST dataset & different in vivo media
# accuracies will be used as one of many metrics in which to select a single biomass reaction for further analyses
in_vivo_biomass_rxns = {'BM1': 'biomass_iEK1008_60atp',
                        'BM2': 'biomass_iEK1011_60atp',
                        'BM3': 'biomass_iNJ661_60atp',
                        'BM4': 'biomass_iNJ661v_SBML_xml_60atp',
                        'BM5': 'biomass_sMtb_CSI_57atp',
                        'BM6': 'biomass_sMtb_IVB_57atp',
                        'BM7': 'biomass_sMtb_NRC_57atp',
                        'BM8': 'biomass_sMtbRECON_vivo_57atp'}

in_vivo_media = {'M1': ['iEK_inVivo', media.iEK_inVivo],
                 'M2': ['sMtb_CSI', media.sMtb_CSI],
                 'M3': ['sMtb_in_vivo', media.sMtb_in_vivo],
                 'M4': ['sMtb_in_vivo_mod', media.sMtb_in_vivo_mod],
                 'M5': ['zimmerman_in_vivo', media.zimmerman_in_vivo]}

in_vivo_keys = {'Key': list(in_vivo_biomass_rxns.keys()) + list(in_vivo_media.keys())}
in_vivo_names = {'Name': list(in_vivo_biomass_rxns.values()) + [v[0] for v in in_vivo_media.values()]}

df_vivo_summary = pd.DataFrame.from_dict({**in_vivo_keys, **in_vivo_names})

with alive_bar(len(in_vivo_biomass_rxns) * len(in_vivo_media), bar='blocks', spinner='notes_scrolling') as bar:
    for model_file in model_files:
        if 'iMtb' in model_file:
            model_path = f'../../3_base_model/models/{model_file}'
            model_name = model_file.split('.')[0]
            combo_name = f'{model_name}_vivo_combinations.xlsx'

            with pd.ExcelWriter('../results/' + combo_name, mode='w') as writer:
                df_vivo_summary.to_excel(writer, sheet_name='summary', index=False)
            writer.save()

            for bm_key, bm_rxn in in_vivo_biomass_rxns.items():
                for medium_key, med in in_vivo_media.items():
                    sim_vs_exp(model_path, bm_rxn, 'ARTIST_2014', med[1], 'fba', combo_name, f'{bm_key}_{medium_key}')
                    bar()


all_datasets = vitro_datasets + vivo_datasets

with alive_bar(len(model_files) * len(all_datasets), bar='blocks', spinner='notes_scrolling') as bar:
    med = media.sMtb_in_vivo_mod

    for model_file in model_files:
        model_path = f'../../3_base_model/models/{model_file}'
        mod_name = model_file.rsplit('.', maxsplit=1)[0]

        for ko_mthd in ['fba']:  # moma and room not used --> too long to run --> see 4.py for alternative approach
            for dataset in all_datasets:
                if dataset in vitro_datasets:
                    bm_rxn = 'biomass_iNJ661_60atp'
                else:
                    bm_rxn = 'biomass_iNJ661v_SBML_xml_60atp'

                sim_vs_exp(model_path, bm_rxn, dataset, med, ko_mthd, f'{mod_name}_predictions.xlsx', f'{dataset}')
                bar()
