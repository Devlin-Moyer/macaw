# fig_3a_S2a_S3a_data.py
'''
Make tables indicating which test each reaction was flagged by in Human-GEM,
yeast-GEM, and iML1515 to use to color the networks
'''

import pandas as pd
from macaw.utils import simplify_test_results

def categorize_rxn(row):
    dupe = row['duplicate_test'] != 'ok'
    dead_end = row['dead_end_test'] != 'ok'
    loop = row['loop_test'] != 'ok'
    dil = row['dilution_test'] != 'ok'
    if not any([dupe, dead_end, loop, dil]):
        return('None')
    elif dupe and not any([dead_end, loop, dil]):
        return('Duplicate')
    elif dead_end and not any([dupe, loop, dil]):
        return('Dead-End')
    elif loop and not any([dupe, dead_end, dil]):
        return('Loop')
    elif dil and not any([dupe, dead_end, loop]):
        return('Dilution')
    else:
        return('Multiple')

# silence Pandas' most annoying least necessary warning message
pd.options.mode.chained_assignment = None

for (figure, model) in [
    ('fig_3a', 'Human-GEMv1.15'),
    ('fig_S2a', 'yeast-GEMv9.0.0'),
    ('fig_S3a', 'iML1515')
]:
    edge_list = pd.read_csv(f'figure_data/{model}_edge-list.csv')
    all_test_results = pd.read_csv(f'figure_data/{model}_test-results.csv')
    # ignore the results of the diphosphate test
    all_test_results = all_test_results.drop('diphosphate_test', axis = 1)
    # filter down to reactions that appear in the edge list
    rxn_ids = set(edge_list['source'].unique().tolist())
    rxn_ids.update(edge_list['target'].unique().tolist())
    test_results = all_test_results[all_test_results['reaction_id'].isin(rxn_ids)]
    # categorize each reaction by the test(s) it was flagged by
    test_results['category'] = simplify_test_results(test_results).apply(
        categorize_rxn, axis = 1
    )
    test_results[['reaction_id', 'category']].to_csv(
        f'figure_data/{figure}_node-list.csv', index = False
    )
    edge_list.to_csv(f'figure_data/{figure}_edge-list.csv', index = False)
