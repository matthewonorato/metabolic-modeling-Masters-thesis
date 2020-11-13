#!/usr/bin/env python3

import os
import cobra
import pandas as pd

__author__ = 'Matthew Onorato'
__version__ = '0.0.1'
__date__ = 'June 06, 2020'
__status__ = 'final'


"""
DESCRIPTION: This script calculates the accuracy of each subsystem's gene essentiality predictions and the percentage of
             unblocked reactions per subsystem.
"""

model_files = [file for file in os.listdir('../../3_base_model/models/') if file.endswith('.json') if 'iMtb' in file]
prediction_files = [file for file in os.listdir('../results/') if 'predictions' in file if 'iMtb' in file]

for prediction_file in prediction_files:
    if 'iMtb' in prediction_file:
        model_file = [file for file in model_files if 'iMtb' in file][0]
    else:
        model_file = [file for file in model_files if 'iEK' in file][0]

    model = cobra.io.load_json_model(f'../../3_base_model/models/{model_file}')

    prediction_path = '../results/' + prediction_file
    df_preds = pd.read_excel(prediction_path, sheet_name=f'{model.id}_all_predictions', index_col=0).fillna('NA')

    # filters out genes not in the model
    model_gene_ids = [gene.id for gene in model.genes]
    df_preds = df_preds.loc[model_gene_ids]

    # maps model genes (and corresponding subsystems) to gene essentiality prediction, excluding NAs
    model_subsystems = set([rxn.subsystem for rxn in model.reactions])
    ess_dataset_names = list(df_preds.columns)
    df_all_subsystem_predictions = pd.DataFrame(index=model_subsystems, columns=ess_dataset_names)

    # maps model genes to reactions then to corresponding subsystems
    model_genes = [gene for gene in model.genes]
    model_genes_rxns_frz = {gene.id: gene.reactions for gene in model_genes}

    gene_subs, model_genes_subs = [], {}
    for gene, rxn_frz in model_genes_rxns_frz.items():
        for rxn in rxn_frz:
            gene_subs.append(rxn.subsystem)

        model_genes_subs[gene] = set(gene_subs)  # removes subsystem duplicates for each gene --> only counted once/gene
        gene_subs = []

    # triple for loop: for each dataset, finds T/F call for each gene and maps T/F call to each subsystem gene is in
    # model_genes_subs_preds Example: 'Rv1185c': [{'Fatty Acid Metabolism': 'TP'}, {'Membrane Metabolism': 'TP'}]
    subsystem_predictions, model_genes_subs_preds = [], {}
    for dataset in ess_dataset_names:
        model_genes_preds = df_preds[dataset].to_dict()
        for gene, subs in model_genes_subs.items():
            for sub in subs:
                ess_pred = model_genes_preds[gene]
                if 'T' in ess_pred or 'F' in ess_pred:
                    subsystem_prediction = {sub: model_genes_preds[gene]}
                    subsystem_predictions.append(subsystem_prediction)

            if len(subsystem_predictions) > 0:
                model_genes_subs_preds[gene] = subsystem_predictions

            subsystem_predictions = []

        # maps each model subsystem to all gene essentiality predictions then calculates each subsystem's accuracy
        subsystem_predictions, model_subs_preds = [], {}
        for subsystem in model_subsystems:
            for sub_preds in model_genes_subs_preds.values():
                for sub_pred in sub_preds:
                    if list(sub_pred.keys())[0] == subsystem:
                        subsystem_predictions.append(sub_pred[subsystem])

            model_subs_preds[subsystem] = subsystem_predictions
            subsystem_predictions = []

        df_model_subs_preds = pd.DataFrame.from_dict(model_subs_preds, orient='index').T

        for subsystem in model_subsystems:
            df_ss_preds = df_model_subs_preds[subsystem].dropna()
            try:
                sub_ess_acc = df_ss_preds.replace({'TP': 1, 'TN': 1, 'FP': 0, 'FN': 0}).sum()/len(df_ss_preds) * 100
                df_all_subsystem_predictions.loc[subsystem, dataset] = round(sub_ess_acc, 0)
            except ZeroDivisionError:
                df_all_subsystem_predictions.loc[subsystem, dataset] = 0

    df_all_subsystem_predictions.rename_axis('Subsystem', axis=0, inplace=True)

    with pd.ExcelWriter(f'../results/{model.id}_subsystem_accuracy.xlsx', mode='w') as writer:
        df_all_subsystem_predictions.to_excel(writer, sheet_name='gene_essentiality')
    writer.save()

    # maps each subsystem to blocked rxn IDs and to all rxns IDs, separately, then calculates each subsystem's accuracy
    df_subsystems_unblk_rxns = pd.DataFrame(index=model_subsystems, columns=['unblocked_reactions'])

    model.objective = 'biomass_iNJ661v_SBML_xml_60atp'  # fewest blocked reactions of in vivo media
    model.reactions.get_by_id('biomass_iNJ661v_SBML_xml_60atp').upper_bound = 1000.0

    blk_rxns = cobra.flux_analysis.variability.find_blocked_reactions(model, open_exchanges=True)

    subsystem_blk_rxns, model_subsystem_blk_rxns = [], {}
    for subsystem in model_subsystems:
        for blk_rxn in blk_rxns:
            if subsystem == model.reactions.get_by_id(blk_rxn).subsystem:
                subsystem_blk_rxns.append(blk_rxn)

        model_subsystem_blk_rxns[subsystem] = subsystem_blk_rxns
        subsystem_blk_rxns = []

    subsystem_rxns, model_subsystems_rxns = [], {}
    for subsystem in model_subsystems:
        for rxn in model.reactions:
            if subsystem == rxn.subsystem:
                subsystem_rxns.append(rxn.id)

        model_subsystems_rxns[subsystem] = subsystem_rxns
        subsystem_rxns = []

    for subsystem, rxns in model_subsystems_rxns.items():
        ss_blk_rxns = len(model_subsystem_blk_rxns[subsystem])
        ss_unblk_rxns = len(rxns) - ss_blk_rxns
        df_subsystems_unblk_rxns.loc[subsystem, 'unblocked_reactions'] = round((ss_unblk_rxns/len(rxns) * 100), 0)

    df_subsystems_unblk_rxns.rename_axis('Subsystem', axis=0, inplace=True)

    with pd.ExcelWriter(f'../results/{model.id}_subsystem_accuracy.xlsx', mode='a') as writer:
        df_subsystems_unblk_rxns.to_excel(writer, sheet_name='unblocked_reactions')
    writer.save()
