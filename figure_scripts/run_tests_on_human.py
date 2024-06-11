# run_tests_on_human.py

import sys
from optlang.glpk_interface import Configuration
import cobra
import pandas as pd
from macaw_main import run_all_tests

try:
    (version, threads) = sys.argv[1:]
    threads = int(threads)
except IndexError:
    sys.exit('provide number of threads to use')

# silence annoying optlang message that prints when you read in a model
Configuration()
model = cobra.io.read_sbml_model(f'GSMMs/Human-GEMv{version}.xml')

# get list of IDs of metabolites that are in DMEM or FBS
media_df = pd.read_csv('figure_data/Table S1.csv')
media_mets = media_df['metabolite_id'].to_list()

redox_pairs = [
    # NAD+ and NADH in every compartment they both exist in
    ('MAM02552c', 'MAM02553c'), ('MAM02552e', 'MAM02553e'),
    ('MAM02552m', 'MAM02553m'), ('MAM02552n', 'MAM02553n'),
    ('MAM02552r', 'MAM02553r'), ('MAM02552x', 'MAM02553x'),
    # NADP+ and NADPH in every compartment they both exist in
    ('MAM02554c', 'MAM02555c'), ('MAM02554l', 'MAM02555l'),
    ('MAM02554m', 'MAM02555m'), ('MAM02554n', 'MAM02555n'),
    ('MAM02554r', 'MAM02555r'), ('MAM02554x', 'MAM02555x'),
    # FAD and FADH2 in every compartment they both exist in
    ('MAM01802c', 'MAM01803c'), ('MAM01802m', 'MAM01803m'),
    ('MAM01802r', 'MAM01803r'), ('MAM01802x', 'MAM01803x'),
    # ubiquinone and ubiquinol in every compartment they both exist in
    ('MAM03103c', 'MAM03102c'), ('MAM03103e', 'MAM03102e'),
    ('MAM03103m', 'MAM03102m'),
    # ferricytochrome c and ferrocytochrome c
    ('MAM01824c', 'MAM01826c'), ('MAM01824m', 'MAM01826m'),
    # ferricytochrome b5 and ferrocytochrome b5
    ('MAM01823c', 'MAM01825c'),
    # thioredoxin
    ('MAM02990c', 'MAM02666c'), ('MAM02990m', 'MAM02666m'),
    # mitothioredoxin
    ('MAM02487m', 'MAM02486m'),
    # oxygen and hydrogen peroxide
    ('MAM02630c', 'MAM02041c'), ('MAM02630e', 'MAM02041e'),
    ('MAM02630l', 'MAM02041l'), ('MAM02630m', 'MAM02041m'),
    ('MAM02630n', 'MAM02041n'), ('MAM02630r', 'MAM02041r'),
    ('MAM02630x', 'MAM02041x')
]
if version == '1.15':
    # IDs for FMN and FMNH2 changed between versions
    redox_pairs.extend([('MAM20020r', 'MAM20019r'), ('MAM20023i', 'MAM20022i')])
elif version == '1.18':
    redox_pairs.extend([('MAM01828r', 'MAM20019r'), ('MAM01828i', 'MAM20019i')])

# protons in all compartments
proton_ids = [m.id for m in model.metabolites if m.id.startswith('MAM02039')]
# IDs for diphosphate in all compartments of Human-GEM
ppi_ids = [m.id for m in model.metabolites if m.id.startswith('MAM02759')]
# IDs for inorganic (mono)phosphate
pi_ids = [m.id for m in model.metabolites if m.id.startswith('MAM02751')]

(test_results, edge_list) = run_all_tests(
    model, redox_pairs, proton_ids, ppi_ids, pi_ids, media_mets, timeout = 1800,
    use_names = True, add_suffixes = True, threads = threads
)

fname = f'figure_data/Human-GEMv{version}'
test_results.to_csv(f'{fname}_test-results.csv', index = False)
edge_list.to_csv(f'{fname}_edge-list.csv', index = False)
