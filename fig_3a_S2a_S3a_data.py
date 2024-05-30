# fig_3a_S2a_S3a_data.py
'''
Make tables indicating which test each reaction was flagged by in Human-GEM,
yeast-GEM, and iML1515 to use to color the networks
'''

import sys
import pandas as pd
from macaw_utils import simplify_test_results

def categorize_rxn(row):
    dupe = row['duplicate_test'] != 'ok'
    dead_end = row['dead_end_test'] != 'ok'
    ppi = row['diphosphate_test'] != 'ok'
    loop = row['loop_test'] != 'ok'
    dil = row['dilution_test'] != 'ok'
    if not any([dupe, dead_end, ppi, loop, dil]):
        return('None')
    elif dupe and not any([dead_end, ppi, loop, dil]):
        return('Duplicate')
    elif dead_end and not any([dupe, ppi, loop, dil]):
        return('Dead-End')
    elif ppi and not any([dupe, dead_end, loop, dil]):
        return('Diphosphate')
    elif loop and not any([dupe, dead_end, ppi, dil]):
        return('Loop')
    elif dil and not any([dupe, dead_end, ppi, loop]):
        return('Dilution')
    else:
        return('Multiple')

# get name of model from command-line
model = sys.argv[1]
# silence Pandas' most annoying least necessary warning message
pd.options.mode.chained_assignment = None

edge_list = pd.read_csv(f'figure_data/{model}_edge-list.csv')
all_test_results = pd.read_csv(f'figure_data/{model}_test-results.csv')

# filter down to reactions that appear in the edge list
rxn_ids = set(edge_list['source'].unique().tolist())
rxn_ids.update(edge_list['target'].unique().tolist())
test_results = all_test_results[all_test_results['reaction_id'].isin(rxn_ids)]
# categorize each reaction by the test(s) it was flagged by
test_results['category'] = simplify_test_results(test_results).apply(
    categorize_rxn, axis = 1
)
# save to appropriately-named file
fname_dict = {
    'Human-GEMv1.15' : '3a',
    'yeast-GEMv9.0.0' : 'S2a',
    'iML1515' : 'S3a'
}
test_results[['reaction_id', 'category']].to_csv(
    f'figure_data/fig_{fname_dict[model]}_data.csv', index = False
)
