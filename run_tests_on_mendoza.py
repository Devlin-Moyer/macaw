# run_tests_on_mendoza.py

import sys
from optlang.glpk_interface import Configuration
import cobra
from macaw_main import run_all_tests

try:
    threads = int(sys.argv[1])
except IndexError:
    sys.exit('provide number of threads to use')

# silence annoying optlang message that prints when you read in a model
Configuration()

# TODO: figure out how to parallelize such that each model is on its own
# thread
d = 'GSMMs/mendoza_2019'
for f in os.listdir(d):
    # annoyingly, some of the models were only available as .mats and others
    # were only available as ".sbml" files, and as far as I know, Cobrapy
    # doesn't have a single input function that automatically figures out the
    # file type from the extension
    if f.endswith('mat'):
        model = cobra.io.load_matlab_model(f'{d}/{f}')
    else:
        model = cobra.io.read_sbml_model(f'{d}/{f}')
    # TODO: figure out if they all have ATP maintenance reactions with
    # predictable IDs cuz ideally we would ensure that their lower bounds = 0
    # don't bother with the IDs for redox metabolites, protons, or phosphates,
    # cuz we can just skip the redox duplicates and the diphosphate test and
    # work with the results from the other tests that don't require potentially
    # tedious input
    (test_results, edge_list) = run_all_tests(
        model, timeout = 360, use_names = True, add_suffixes = True
    )
    test_results.to_csv(
        f'figure_data/fig_mendoza_data/{f}_test-results.csv', index = False
    )
