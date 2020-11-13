#!/usr/bin/env python3

import pandas as pd
import cobra
from continuous_testing import run_test_suite
import re
from alive_progress import alive_bar
from Bio import Entrez
import time

__author__ = 'Matthew Onorato'
__version__ = '0.0.2'
__date__ = 'June 16, 2020'
__maintainer__ = 'Matthew Onorato'
__email__ = 'monorato2@gmail.com' > 'monorato4398@sdsu.edu'
__status__ = 'final'


"""
DESCRIPTION: This script updates gene-reactions-rules based on UniProt-verified annotations. It also adds gene names
             and gene annotations by mapping them from NCBI and UniProt to iMtb_H37Rv.json.
"""

# maps reaction IDs to new GPRs then updates iMtb_H37Rv model with new GPRs
rxn_id_gpr = pd.read_excel('../data/iEK1008_updates.xlsx', sheet_name='reactions_updated',
                           usecols=['Reaction_ID', 'Gene_Reaction_Rule'],
                           index_col='Reaction_ID')['Gene_Reaction_Rule'].to_dict()

model = cobra.io.load_json_model('../models/iMtb_H37Rv.json')

for rxn in model.reactions:
    if rxn.id in rxn_id_gpr.keys():
        if not rxn.gene_reaction_rule == rxn_id_gpr[rxn.id]:
            try:
                rxn.gene_reaction_rule = rxn_id_gpr[rxn.id]
            except AttributeError:  # 'float' object (nan) has no attribute 'strip'
                rxn.gene_reaction_rule = ''
        else:
            continue

cobra.io.save_json_model(model, '../models/iMtb_H37Rv.json')
run_test_suite('../models/iMtb_H37Rv.json', update='GPRs_updated')


# maps UniProt IDs to locus tags (Rv numbers) and gene names using UniProt proteome file
rv_gene_file = '../data/uniprot-proteome_UP000001584_gene_names.tsv'
uniprot_ids_genes = {line.strip().split('\t')[0]: line.strip().split('\t')[1] for line in open(rv_gene_file, 'r')}

uniprot_ids_locus_tags_gene_names = {}
for uniprot_id, locus_or_gene in uniprot_ids_genes.items():
    locus_or_gene = locus_or_gene.split(' ')
    if len(locus_or_gene) == 1 and locus_or_gene[0].startswith('Rv'):  # Rv number only
        uniprot_ids_locus_tags_gene_names[uniprot_id] = [locus_or_gene[0], 'NA']
    elif len(locus_or_gene) == 1 and not locus_or_gene[0].startswith('Rv'):  # gene name only --> None
        print(locus_or_gene)
    else:
        if locus_or_gene[0].startswith('Rv') and locus_or_gene[1].startswith('Rv'):  # remaining Rv_number-only genes
            if '.' not in locus_or_gene[0]:
                uniprot_ids_locus_tags_gene_names[uniprot_id] = [locus_or_gene[0], 'NA']
            else:
                uniprot_ids_locus_tags_gene_names[uniprot_id] = [locus_or_gene[1], 'NA']
        elif re.match('[a-z|P]', locus_or_gene[0]) is not None:  # gene names with lowercase letter or P (for PE/PPE)
                r = re.compile('Rv*')  # matches all gene names starting with Rv
                locus = list(filter(r.match, locus_or_gene))[0]
                uniprot_ids_locus_tags_gene_names[uniprot_id] = [locus, locus_or_gene[0]]
        else:
            if locus_or_gene[0].startswith('Rv'):
                uniprot_ids_locus_tags_gene_names[uniprot_id] = [locus_or_gene[0], 'NA']  # Rvs with other locus tags
            else:
                uniprot_ids_locus_tags_gene_names[uniprot_id] = [locus_or_gene[1], 'NA']  # some start with TB


# annotates genes with gene name, NCBI gene and protein IDs, RefSeq locus tag and gene names, UniProt IDs, and SBO terms
rvs_genes_uniprots = {v[0]: [v[1], k] for k, v in uniprot_ids_locus_tags_gene_names.items()}

model_genes = set([gene for gene in model.genes])

Entrez.email = 'monorato4398@sdsu.edu'
with alive_bar(len(model_genes), bar='blocks', spinner='notes_scrolling') as bar:
    for model_gene in model_genes:
        if model_gene.id in rvs_genes_uniprots.keys():
            if rvs_genes_uniprots[model_gene.id][0] != 'NA':
                model_gene.name = rvs_genes_uniprots[model_gene.id][0]

        if 'ncbigene' not in model_gene.annotation.keys():
            gene_handle = Entrez.esearch(db='gene', term=f'NC_000962 AND {model_gene.id} AND H37Rv', retmax=1)
            gene_record = Entrez.read(gene_handle)
            ncbigene = gene_record['IdList']
            if len(ncbigene) > 0:
                model_gene.annotation['ncbigene'] = ncbigene

        if 'ncbigi' not in model_gene.annotation.keys():
            protein_handle = Entrez.esearch(db='protein', term=f'BioProject: PRJNA57777 AND {model_gene.id}', retmax=1)
            protein_record = Entrez.read(protein_handle)
            ncbigi = protein_record['IdList']

            if len(ncbigi) > 0:
                model_gene.annotation['ncbiprotein'] = ncbigi

        else:
            model_gene.annotation['ncbiprotein'] = model_gene.annotation.pop('ncbigi')

        if 'refseq_locus_tag' not in model_gene.annotation.keys():
            model_gene.annotation['refseq_locus_tag'] = [model_gene.id]

        if 'refseq_name' not in model_gene.annotation.keys():
            if rvs_genes_uniprots[model_gene.id][0] != 'NA':
                model_gene.annotation['refseq_name'] = [rvs_genes_uniprots[model_gene.id][0]]

        if 'sbo' not in model_gene.annotation.keys():
            model_gene.annotation['sbo'] = 'SBO:0000243'

        if 'uniprot' not in model_gene.annotation.keys():
            try:
                model_gene.annotation['uniprot'] = [rvs_genes_uniprots[model_gene.id][1]]
            except KeyError:
                continue

        time.sleep(0.3)

        bar()

# corrects duplicate ncbigene IDs; same ncbigene ID used for Rv2321c/Rv2322c and Rv1104/Rv1105 due to frameshifts
ncbigene_changes = {'Rv2953': '887772', 'Rv0998': '885385', 'Rv2135c': '887273', 'Rv0989c': '885355',
                    'Rv1558': '886363', 'Rv1700': '885049', 'Rv0265c': '886650'}

for gene_id, ncbi_id in ncbigene_changes.items():
    model.genes.get_by_id(gene_id).annotation['ncbigene'] = [ncbi_id]

# fills in missing Uniprot IDs for Rv0815c (same ID as Rv3117) and Rv2922A
uniprot_id_changes = {'Rv0815c': 'P9WHF9', 'Rv2922A': 'P9WQC9'}

for gene_id, uniprot_id in uniprot_id_changes.items():
    model.genes.get_by_id(gene_id).annotation['uniprot'] = [uniprot_id]

cobra.io.save_json_model(model, '../models/iMtb_H37Rv.json')
run_test_suite('../models/iMtb_H37Rv.json', update='gene_names_and_annotations_updated')
