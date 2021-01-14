#!/usr/bin/env python3

import cobra
import pandas as pd
from continuous_testing import run_test_suite
from alive_progress import alive_bar

__author__ = 'Matthew Onorato'
__version__ = '0.0.2'
__date__ = 'March 15, 2020'
__maintainer__ = 'Matthew Onorato'
__email__ = 'monorato2@gmail.com' > 'monorato4398@sdsu.edu'
__status__ = 'final'


"""
DESCRIPTION: This script adds new metabolites and new reactions to metabolic model iEK1008.json. It also tracks model
             statistics after each update, a la Memote with continuous testing.
"""


def create_reaction(json_model, rxn_id, rxn_name, ss, mets, lb, ub, gpr):
    """
    Creates a COBRApy Reaction object for each reaction to be added to iEK1008
    :param json_model: JSON-formatted metabolic model
    :param rxn_id:     user-provided reaction ID
    :param rxn_name:   user-provided reaction name
    :param ss:         user-provided subsystem
    :param mets:       user-provided dictionary of Metabolite objects : coefficients
    :param lb:         user-provided reaction lower bound
    :param ub:         user-provided reaction lower bound
    :param gpr:        user-provided gene-protein-reaction relationship (denoted as gene_reaction_rule in model)
    :return:           Reaction object
    """
    reaction = cobra.Reaction(rxn_id)
    reaction.name = rxn_name
    reaction.subsystem = ss
    json_model.add_reactions([reaction])
    reaction.add_metabolites(mets)
    reaction.lower_bound = lb
    reaction.upper_bound = ub
    try:
        reaction.gene_reaction_rule = gpr
    except AttributeError:
        reaction.gene_reaction_rule = ''
    return reaction


def main():
    """
    Adds new reactions and metabolites to iEK1008.json while performing continuous testing
    (note: main() also allows create_reaction() to be imported without main() being run)

    :return: writes new model (iMtb_H37Rv.json) and testing report (iEK1008_test_suite.tsv)
    """
    run_test_suite('../models/iEK1008.json')  # runs test suite with iEK1008.json

    # rewrites iEK1008.json to iMtb_H37Rv.json so original model is not overwritten
    model_iek = cobra.io.load_json_model('../models/iEK1008.json')
    cobra.io.save_json_model(model_iek, '../models/iMtb_H37Rv.json')
    model = cobra.io.load_json_model('../models/iMtb_H37Rv.json')

    # removes 10 imbalanced reactions from iEK1008; all 10 reactions are added back with balanced formulas during update
    rxns_to_bal = [rxn.id for rxn in model.reactions if len(rxn.check_mass_balance()) > 0
                   if 'EX_' not in rxn.id and 'DM_' not in rxn.id and 'BIOMASS' not in rxn.id]

    for rxn_to_bal in rxns_to_bal:
        model.reactions.get_by_id(rxn_to_bal).remove_from_model()
    cobra.io.save_json_model(model, '../models/iMtb_H37Rv.json')

    run_test_suite('../models/iMtb_H37Rv.json', update='imbalanced_reactions_removed')

    # creates COBRApy Metabolite objects for new metabolites
    df_new_mets = pd.read_excel('../data/iEK1008_updates.xlsx', sheet_name='metabolites_added', usecols='A:C')

    new_mets = {}
    for index, row in df_new_mets.iterrows():
        new_met_id = str(row['Metabolite_ID'])
        new_met_name = row['Metabolite_Name']
        new_met_formula = row['Metabolite_Formula']
        if new_met_id.endswith('c'):
            new_met_comp = 'c'
        elif new_met_id.endswith('e'):
            new_met_comp = 'e'
        else:
            print('Metabolite compartment could not be determined. Please check metabolite id.')
            new_met_comp = ''
        new_met = cobra.Metabolite(new_met_id, name=new_met_name, formula=new_met_formula, compartment=new_met_comp)
        new_mets[new_met_id] = new_met

    df_new_rxns = pd.read_excel('../data/iEK1008_updates.xlsx', sheet_name='reactions_added', usecols='A:G')

    with alive_bar(len(df_new_rxns), bar='blocks', spinner='notes_scrolling') as bar:
        for index, row in df_new_rxns.iterrows():
            new_rxn_mets = {}
            new_rxn_form = row['Reaction_Formula']
            if ' --> ' in new_rxn_form:
                new_rxn_form = new_rxn_form.split(' --> ')
            elif ' <=> ' in new_rxn_form:
                new_rxn_form = new_rxn_form.split(' <=> ')
            else:
                print('Unexpected symbol in ' + row['Reaction_Formula'])

            subs = new_rxn_form[0].split(' + ')
            for sub in subs:
                if '.0' in sub:
                    sub_coeff = -1 * float(sub.split(' ')[0])
                    sub_id = sub.split(' ')[-1]
                    try:
                        new_rxn_sub = new_mets[sub_id]
                    except KeyError:  # metabolite is not new, i.e. already in iEK1008
                        new_rxn_sub = model.metabolites.get_by_id(sub_id)
                else:
                    sub_coeff = -1.0
                    try:
                        new_rxn_sub = new_mets[sub]
                    except KeyError:
                        new_rxn_sub = model.metabolites.get_by_id(sub)
                new_rxn_mets[new_rxn_sub] = sub_coeff

            pros = new_rxn_form[1].split(' + ')
            for pro in pros:
                if '.0' in pro:
                    pro_coeff = float(pro.split(' ')[0])
                    pro_id = pro.split(' ')[-1]
                    try:
                        new_rxn_pro = new_mets[pro_id]
                    except KeyError:
                        new_rxn_pro = model.metabolites.get_by_id(pro_id)
                else:
                    pro_coeff = 1.0
                    try:
                        new_rxn_pro = new_mets[pro]
                    except KeyError:
                        new_rxn_pro = model.metabolites.get_by_id(pro)
                new_rxn_mets[new_rxn_pro] = pro_coeff

            # creates new reactions with new COBRApy Reaction and Metabolite objects
            create_reaction(model, row['Reaction_ID'], row['Reaction_Name'], row['Subsystem'], new_rxn_mets,
                            float(row['Lower_Bound']), float(row['Upper_Bound']), row['Gene_Reaction_Rule'])

            cobra.io.save_json_model(model, '../models/iMtb_H37Rv.json')

            run_test_suite('../models/iMtb_H37Rv.json', update=row['Reaction_ID'])

            bar()

    return


if __name__ == '__main__':
    main()
