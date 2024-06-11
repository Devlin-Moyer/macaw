# investigate_human_dilution.py
'''
The proportion of reactions flagged by the dilution test in Human-GEM jumped
significantly between versions 1.11 and 1.12, but it's not clear which changes
made between those two versions are responsible, so make a list of all the
reactions that were flagged in 1.12 that weren't in 1.11 to see if that
clarifies things (also save the edge lists just in case those are also helpful)
'''

import sys
import optlang
import pandas as pd
import cobra
from macaw_main import dilution_test

try:
    threads = int(sys.argv[1])
except IndexError:
    sys.exit('specify number of threads to use')

optlang.glpk_interface.Configuration()

# get list of IDs of metabolites that are in DMEM or FBS
media_concs = pd.read_csv('media/DMEM-FBS_ingredients.csv')
media_mets = media_concs[
    (media_concs['DMEM'] != '0') | (media_concs['FBS'] != '0')
]['metabolite_id'].to_list()

results_dfs = list()
for version in ['11', '12']:
    print(f'working on version 1.{version}')
    model = cobra.io.read_sbml_model(f'GSMMs/Human-GEMv1.{version}.xml')
    (results, edges) = dilution_test(
        model, media_mets = media_mets, use_names = True,
        add_suffixes = True, threads = threads
    )
    # don't need dead-end test results and want to be able to merge the
    # results from both versions together into a single dataframe, so rename
    # the dilution test column
    results = results.drop('dead_end_test', axis = 1).rename(columns = {
        'dilution_test' : f'dil_test_in_v{version}'
    })
    results_dfs.append(results)
    # also save results separately just in case
    results.to_csv(
        f'figure_data/Human-GEMv{version}_dil-results.csv', index = False
    )
    # might as well just save edge lists to two separate files now
    edge_list.to_csv(
        f'figure_data/Human-GEMv{version}_dil-edges.csv', index = False
    )

# keep all reactions that were removed or added between versions
both_results = results_dfs[0].merge(results_dfs[1], how = 'outer')
# filter for reactions whose test results were different
both_results[
    both_results['dil_test_in_v11'] != both_results['dil_test_in_v12']
].to_csv('figure_data/Human-GEMv1.11-12_dil-results.csv', index = False)
