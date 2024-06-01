# run_tests_on_ecoli.py

import sys
from optlang.glpk_interface import Configuration
import cobra
import pandas as pd
from macaw_main import run_all_tests

try:
    threads = int(sys.argv[1])
except IndexError:
    sys.exit('provide number of threads to use')

# silence annoying optlang message that prints when you read in a model
Configuration()
model = cobra.io.read_sbml_model(f'GSMMs/iML1515.xml')

# get list of IDs of metabolites that are in M9 media
in_media = pd.read_csv('input_data/M9_ingredients.csv')['metabolite_id']

# set lower bound on ATP maintenance reaction to 0
model.reactions.get_by_id('ATPM').lower_bound = 0

redox_pairs = [
    # NAD(H), NADP(H), FAD(H2), and FMN(H2)
    ('nad_c', 'nadh_c'), ('nadp_c', 'nadph_c'), ('fad_c', 'fadh2_c'),
    ('fmn_c', 'fmnh2_c'),
    # ubiquinone-8 and ubiquinol-8
    ('q8_c', 'q8h2_c'),
    # menaquinone-8 and menaquinol-8
    ('mqn8_c', 'mql8_c'),
    # riboflavin
    ('ribflv_c', 'rbflvrd_c'),
    # glutathione
    ('gthox_c', 'gthrd_c'),
    # glutaredoxin
    ('grxox_c', 'grxrd_c'),
    # thioredoxin
    ('trdox_c', 'trdrd_c'),
    # flavodoxin (reduced and semi-oxidized)
    ('flxso_c', 'flxr_c'),
    # oxygen and hydrogen peroxide
    ('o2_c', 'h2o2_c')
]
# protons in all compartments
proton_ids = ['h_c', 'h_e', 'h_p']

(test_results, edge_list) = run_all_tests(
    model, redox_pairs, proton_ids, media_mets = in_media, timeout = 360,
    use_names = True, add_suffixes = True, threads = threads
)

fname = 'figure_data/iML1515'
test_results.to_csv(f'{fname}_test-results.csv', index = False)
edge_list.to_csv(f'{fname}_edge-list.csv', index = False)
