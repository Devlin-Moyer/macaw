# fig_2_S5_data.py
'''
Make CSV files with internal model IDs and gene IDs for all reactions in 
in Human-GEM, Yeast-GEM, and iML1515. For reactions associated with more than
one gene, make multiple rows with that reaction ID so there is only one ID of
each type in each row
'''

import optlang
import cobra
import pandas as pd

# silence annoying optlang message that prints when you read in a model for the
# first time in a given Python session
optlang.glpk_interface.Configuration()

for model_name in ['Human-GEMv1.15', 'yeast-GEMv9.0.0', 'iML1515']:
    # read in model
    model = cobra.io.read_sbml_model(f'GSMMs/{model_name}.xml')
    rxn_to_gene = {'reaction_id': list(), 'gene_id' : list()}
    for r in model.reactions:
        for g in r.genes:
            if model_name == 'Human-GEMv1.15':
                # Human-GEM uses Ensembl gene IDs as the primary gene IDs, but
                # KEGG uses NCBI gene IDs. All but 2 genes in version 1.15 of
                # Human-GEM have NCBI gene IDs in their annotations
                if 'ncbigene' in g.annotation.keys(): 
                    rxn_to_gene['reaction_id'].append(r.id)
                    rxn_to_gene['gene_id'].append(g.annotation['ncbigene'])
            else:
                # KEGG uses the same gene IDs that both Yeast-GEM and iML1515
                # use as the primary IDs for all of their genes
                rxn_to_gene['reaction_id'].append(r.id)
                rxn_to_gene['gene_id'].append(g.id)
    # save as CSV
    pd.DataFrame(rxn_to_gene).to_csv(
        f'figure_data/{model_name}_reactions-to-genes.csv', index = False
    )
