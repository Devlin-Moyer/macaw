# run_tests_on_alteromonas.py

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
model = cobra.io.read_sbml_model('GSMMs/helen-alteromonas.xml')

# set lower bound on ATP maintenance reaction to 0
model.reactions.get_by_id('ATPM').lower_bound = 0

redox_pairs = [
    # NAD(H), NADP(H), and FAD(H2)
    ('cpd00003_c0', 'cpd00004_c0'), ('cpd00006_c0', 'cpd00005_c0'),
    ('cpd00015_c0', 'cpd00982_c0'),
    # ubiquinone-8 and ubiquinol-8
    ('cpd15560_c0', 'cpd15561_c0'),
    # menaquinone-8 and menaquinol-8
    ('cpd15500_c0', 'cpd15499_c0'),
    # cytochrome-bo
    ('cpd18077_c0', 'cpd18078_c0'),
    # glutathione
    ('cpd00111_c0', 'cpd00042_c0'),
    # ferredoxin
    ('cpd11621_c0', 'cpd11620_c0'),
    # oxygen and hydrogen peroxide
    ('cpd00007_c0', 'cpd00025_c0')
]
# protons in all compartments
proton_ids = ['cpd00067_c0', 'cpd00067_e0']
# diphosphate and inorganic monophosphate
ppi_ids = ['cpd00012_c0']
pi_ids = ['cpd00009_c0']

(test_results, edge_list) = run_all_tests(
    model, redox_pairs, proton_ids, ppi_ids, pi_ids, timeout = 360,
    use_names = True, add_suffixes = True, threads = threads
)

fname = 'figure_data/helen-alteromonas'
test_results.to_csv(f'{fname}_test-results.csv', index = False)
edge_list.to_csv(f'{fname}_edge-list.csv', index = False)
