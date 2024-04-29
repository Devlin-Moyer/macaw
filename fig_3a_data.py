# fig_3a_data.py
'''
Make a table indicating which test each reaction was flagged by in version 1.15
of Human-GEM to use to color the nodes in the graph of all pathway networks
'''

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

# silence Pandas' most annoying least necessary warning message
pd.options.mode.chained_assignment = None

edge_list = pd.read_csv('figure_data/Human-GEMv1.15_edge-list.csv')
all_test_results = pd.read_csv('figure_data/Human-GEMv1.15_test-results.csv')

# filter down to reactions that appear in the edge list
rxn_ids = set(edge_list['source'].unique().tolist())
rxn_ids.update(edge_list['target'].unique().tolist())
test_results = all_test_results[all_test_results['reaction_id'].isin(rxn_ids)]
# categorize each reaction by the test(s) it was flagged by
test_results['category'] = simplify_test_results(test_results).apply(
    categorize_rxn, axis = 1
)
test_results[['reaction_id', 'category']].to_csv(
    'figure_data/fig_3a_nodes.csv', index = False
)
