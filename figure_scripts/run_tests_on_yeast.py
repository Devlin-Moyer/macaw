# run_tests_on_yeast.py

import sys
from optlang.glpk_interface import Configuration
import cobra
import pandas as pd
from macaw.main import run_all_tests

try:
    (version, threads) = sys.argv[1:]
    threads = int(threads)
except IndexError:
    sys.exit('provide number of threads to use')

pd.set_option('display.max_columns', None)
# silence annoying optlang message that prints when you read in a model
Configuration()
model = cobra.io.read_sbml_model(f'GSMMs/yeast-GEMv{version}.xml')

# get list of IDs of metabolites in the Verduyn minimal mineral medium
media_df = pd.read_csv('figure_data/Table S3.csv')
media_mets = media_df['metabolite_id'].to_list()

# set lower bound on ATP maintenance reaction to 0 for consistency with
# Human-GEM, which doesn't have a comparable reaction/bound
model.reactions.get_by_id('r_4046').lower_bound = 0

redox_pairs = [
    # NAD+ and NADH in every compartment they both exist in
    ('s_1198', 's_1203'), ('s_1199', 's_1204'), ('s_1200', 's_1205'),
    ('s_1202', 's_1206'), ('s_1201', 's_3753'), ('s_2820', 's_2818'),
    # NADP+ and NADPH in every compartment they both exist in
    ('s_1207', 's_1212'), ('s_1208', 's_1213'), ('s_1210', 's_1214'),
    ('s_1211', 's_1215'), ('s_2800', 's_2799'), ('s_2953', 's_2952'),
    # FAD and FADH2 in every compartment they both exist in
    ('s_0687', 's_0689'), ('s_0688', 's_0690'),
    # FMN and FMNH2 in every compartment they both exist in
    ('s_0714', 's_0717'),
    # ubiquinone-6 and ubiquinol-6
    ('s_1537', 's_1535'),
    # ferricytochrome c and ferrocytochrome c
    ('s_0709', 's_0710'),
    # ferricytochrome b5 and ferrocytochrome b5
    ('s_3826', 's_3827'), ('s_4210', 's_4209'),
    # TRX1 and TRX1 disulphide
    ('s_1617', 's_1621'), ('s_1616', 's_1620'), ('s_1618', 's_1622'),
    ('s_1619', 's_1623'),
    # oxygen and hydrogen peroxide
    ('s_1275', 's_0837'), ('s_1279', 's_0840')
]
# protons in all compartments
proton_ids = [m.id for m in model.metabolites if m.name == 'H+']
# IDs for diphosphate in all compartments
ppi_ids = [
    's_0633', 's_0635', 's_0636', 's_0637', 's_0638', 's_2834', 's_2860',
    's_3095', 's_4018', 's_4157', 's_4290', 's_4291'
]
# IDs for inorganic (mono)phosphate
pi_ids = [
    's_1322', 's_1323', 's_1324', 's_1325', 's_1326', 's_1329', 's_2966',
    's_2977', 's_2995', 's_3228', 's_3536'
]
(test_results, edge_list) = run_all_tests(
    model, redox_pairs, proton_ids, ppi_ids, pi_ids, media_mets, timeout = 360,
    use_names = True, add_suffixes = True, threads = threads, verbose = 2
)

fname = f'figure_data/yeast-GEMv{version}'
test_results.to_csv(f'{fname}_test-results.csv', index = False)
edge_list.to_csv(f'{fname}_edge-list.csv', index = False)
