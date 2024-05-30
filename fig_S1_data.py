# fig_S1_data.py
'''
Find all dead-end pathways that contain an irreversible exchange reaction
'''

import sys
import cobra
import pandas as pd

try:
    model_name = sys.argv[1]
except:
    sys.exit('provide name of GSMM')

test_results = pd.read_csv(f'figure_data/{model_name}_results-with-KEGG.csv')
model = cobra.io.read_sbml_model(f'GSMMs/{model_name}.xml')

# irreversible exchange reactions
irrev_exchs = [
    r.id for r in model.reactions if r.boundary and not r.reversibility
]
# pathways containing irreversible exchange reactions flagged by dead-end test
irrev_exch_dead_paths = test_results[
    (test_results['dead_end_test'] != 'ok') &
    test_results['reaction_id'].isin(irrev_exchs)
]['pathway'].unique().tolist()
# all reactions in those pathways flagged by the dead-end test
irrev_exch_dead_rxns = test_results[
    test_results['pathway'].isin(irrev_exch_dead_paths) &
    (test_results['dead_end_test'] != 'ok')
]
irrev_exch_dead_rxns = irrev_exch_dead_rxns[[
    'reaction_id', 'reaction_equation', 'pathway', 'kegg_group'
]]
print(irrev_exch_dead_rxns.sort_values('pathway'))
print(irrev_exch_dead_rxns['kegg_group'].value_counts())
print(test_results['kegg_group'].value_counts())
