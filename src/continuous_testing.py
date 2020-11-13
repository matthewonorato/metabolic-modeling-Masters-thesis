#!/usr/bin/env python3

import cobra
import media
import os

__author__ = 'Matthew Onorato'
__version__ = '0.0.3'
__date__ = 'June 17, 2020'
__maintainer__ = 'Matthew Onorato'
__email__ = 'monorato2@gmail.com' > 'monorato4398@sdsu.edu'
__status__ = 'final'

"""
DESCRIPTION: This script runs a series of model tests after each update to metabolic model iEK1008, similar to Memote.
             A single tsv file is written that tracks total counts for each test.
             
             The model info that is tracked is: Model ID
                                                Total Metabolites, Reactions, Genes, Compartments
                                                Metabolic Coverage
                                                Mass & Charge Balance
                                                Disconnected Metabolites
                                                Unbounded Flux In Default Medium (Griffin)
                                                Blocked Reactions
                                                Single Reaction Metabolites
                                                Energy-generating Cycles for ATP, GTP, CTP, and UTP
                                                Biomass & ATPM Flux: iEK_m7H10 vs. m7H10
"""


def run_test_suite(json_model, update='start'):
    """
    Runs a series of tests to determine if any updates significantly impact the metabolic model
    :param json_model: COBRApy-compliant model in JSON format
    :param update: the type of update that is made
    :return: a tsv file chronicling the model update test results
    """
    model = cobra.io.load_json_model(json_model)  # iEK1008.json: Griffin 2011 media pre-set in model

    model_info = {'Model_ID': model.id,
                  'Model_Update': update,
                  'Total_Metabolites': str(len(model.metabolites)),
                  'Total_Reactions': str(len(model.reactions)),
                  'Total_Genes': str(len(model.genes)),
                  'Total_Compartments': str(len(model.compartments.keys())),
                  'Metabolic_Coverage': str(len(model.reactions)/len(model.genes))}

    try:
        model.objective = 'BIOMASS__2'
    except ValueError:
        try:
            model.objective = 'biomass_iEK1008_60atp'
            model.reactions.get_by_id('biomass_iEK1008_60atp').upper_bound = 1000
        except ValueError:
            biomass_rxns = [rxn.id for rxn in model.reactions if 'biomass' in rxn.id.lower()]
            model.objective = biomass_rxns[0]

    to_balance = [rxn.id for rxn in model.reactions if len(rxn.check_mass_balance()) > 0
                  if 'EX_' not in rxn.id and 'DM_' not in rxn.id and 'biomass' not in rxn.id.lower()]

    mets_no_rxn = [met.id for met in model.metabolites if len(met.reactions) == 0]

    try:
        df_fva = cobra.flux_analysis.flux_variability_analysis(model, fraction_of_optimum=1.0)
    except cobra.exceptions.Infeasible:  # needed for models without pre-set media or those that cannot grow
        model.medium = {}  # clears any pre-set media

        for rxn_id, uptake_rate in media.iEK_griffinEssen.items():  # loads Griffin media
            model.reactions.get_by_id(rxn_id).lower_bound = uptake_rate

        df_fva = cobra.flux_analysis.flux_variability_analysis(model, fraction_of_optimum=1.0)

    min_value = 0.99999 * df_fva['minimum'].min()
    max_value = 0.99999 * df_fva['maximum'].max()
    unbound_rxns = (df_fva[df_fva['minimum'] < min_value] + df_fva[df_fva['maximum'] > max_value]).index.to_list()

    model_consistency = {'Imbalanced_Reactions': str(len(to_balance)),
                         'Disconnected_Metabolites': str(len(mets_no_rxn)),
                         'Unbound_Reactions_in_Medium': str(len(unbound_rxns))}

    blocked_rxns = cobra.flux_analysis.variability.find_blocked_reactions(model, open_exchanges=True)

    # single_rxn_mets not the same as orphan & dead-end metabolites --> see updated model results in https://memote.io/
    single_rxn_mets = [met.id for met in model.metabolites if len(met.reactions) == 1]

    model_network_topology = {'Blocked_Reactions': str(len(blocked_rxns)),
                              'Single_Reaction_Metabolites': str(len(single_rxn_mets))}

    # clears model of active exchange reactions (i.e. removes 'preset' media); same as model.medium = {}
    for rxn in model.reactions:
        if 'EX_' in rxn.id:
            rxn.lower_bound = 0

    xtps = {'ATPM': 'Erroneous ATP Flux',
            'NTP3': 'Erroneous GTP Flux',
            'ctp_c': 'Erroneous CTP Flux',
            'utp_c': 'Erroneous UTP Flux'}

    xtp_yields = 0
    model_energy_cycles = {}

    for rxn_or_met in xtps.keys():
        if '_c' not in rxn_or_met:
            model.reactions.get_by_id(rxn_or_met).lower_bound = 0  # Gsmodutils: = 1
            model.objective = rxn_or_met
        elif '_c' in rxn_or_met:
            drain_rxn = cobra.Reaction(f'DR_{rxn_or_met}')
            model.add_reactions([drain_rxn])
            drain_rxn.add_metabolites({cobra.Metabolite(rxn_or_met): -1.0, cobra.Metabolite('h2o_c'): -1.0,
                                       cobra.Metabolite('h_c'): 1.0, cobra.Metabolite('pi_c'): 1.0,
                                       cobra.Metabolite(rxn_or_met.replace('tp_c', 'dp_c')): 1.0})
            drain_rxn.lower_bound = 0
            drain_rxn.upper_bound = 1000
            model.objective = f'DR_{rxn_or_met}'

        xtp_yield = model.optimize().objective_value  # ATPM: closed bounds = 0 // Griffin = 10.875 mmol ATP gDW-1 hr-1
        xtp_yields += xtp_yield

    model_energy_cycles['Erroneous_xTP_Flux'] = str(xtp_yields)

    model_growth = {}
    mtb_media = {'Biomass_Flux_in_iEK_m7H10': media.iEK_m7H10,
                 'Biomass_Flux_in_m7H10': media.m7H10,
                 'ATPM_Flux_in_iEK_m7H10': media.iEK_m7H10,
                 'ATPM_Flux_in_m7H10': media.m7H10}

    for name, mtb_medium in mtb_media.items():
        model.medium = {}

        for ex_rxn_id, uptake_rate in mtb_medium.items():
            try:
                model.reactions.get_by_id(ex_rxn_id).lower_bound = uptake_rate
            except KeyError:  # bypasses exchange reactions in media that are not present in model (e.g. EX_zn2_e)
                continue

        if 'Biomass' in name:
            model.objective = 'BIOMASS__2'
            model_growth[name] = str(model.optimize().objective_value)
        elif 'ATPM' in name:
            model.objective = 'ATPM'
            model_growth[name] = str(model.optimize().objective_value)

    model_tests = {**model_info, **model_consistency, **model_network_topology, **model_energy_cycles, **model_growth}

    header = list(model_tests.keys())
    model_results = list(model_tests.values())

    file_path = '../results/'
    file_name = 'iEK1008_test_suite.tsv'  # f'{model.id}': will write different files when model ID changes

    if file_name not in os.listdir(file_path):
        update_file = open(file_path + file_name, 'w')
        update_file.write('\t'.join(header) + '\n')
        update_file.write('\t'.join(model_results) + '\n')
    else:
        update_file = open(file_path + file_name, 'a')
        update_file.write('\t'.join(model_results) + '\n')
