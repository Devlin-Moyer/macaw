# fig_S4_data.py

import sys
from optlang.glpk_interface import Configuration
import logging
import pandas as pd
import os
import time
import cobra
from macaw_main import dead_end_test, dilution_test, duplicate_test, loop_test

try:
    threads = int(sys.argv[1])
except IndexError:
    sys.exit('provide number of threads to use')

# silence messages from Cobrapy that are not relevant in this context
Configuration()
logging.getLogger('cobra').setLevel(logging.ERROR)

# get list of IDs of metabolites that are in DMEM or FBS
media_concs = pd.read_csv('input_data/DMEM-FBS_ingredients.csv')
media_mets = media_concs[
    (media_concs['DMEM'] != '0') | (media_concs['FBS'] != '0')
]['metabolite_id'].to_list()

# get paths to all versions of Human-GEM
d = 'GSMMs'
model_paths = [f'{d}/{f}' for f in os.listdir(d) if os.path.isfile(f'{d}/{f}')]
# skip any models we already have results for in the output file
out_fname = 'figure_data/fig_S4_data.csv'
if os.path.exists(out_fname):
    already_done = pd.read_csv(out_fname)['model'].to_list()
    model_paths = [p for p in model_paths if p not in already_done]

for model_path in model_paths:
    start_time = time.time()
    model = cobra.io.read_sbml_model(model_path)
    # start by figuring out which of the redox metabolites are actually in this
    # version
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
        # the numeric parts of the IDs for FMN and FMNH2 changed, so include
        # both versions (we'll remove the ones not in this version later)
        ('MAM20020r', 'MAM20019r'), ('MAM20023i', 'MAM20022i'),
        ('MAM01828r', 'MAM20019r'), ('MAM01828i', 'MAM20019i'),
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
    # we know NAD was in every version, but the earlier ones formatted the IDs
    # slightly differently, so if MAM02552c isn't in the model, reformat the IDs
    try:
        model.metabolites.get_by_id('MAM02552c')
    except KeyError:
        redox_pairs = [
            (m1.replace('MAR', 'm'), m2.replace('MAR', 'm'))
            for (m1, m2) in redox_pairs
        ]
    # two of the compartment suffixes changed in version 1.11
    version = int(model_path.split('.')[1])
    if version < 11:
        redox_pairs = [(
            m1.replace('e', 's').replace('x', 'p'),
            m2.replace('e', 's').replace('x', 'p')
        ) for (m1, m2) in redox_pairs]
    # now drop all pairs with IDs that aren't in this particular version
    all_mets = [m.id for m in model.metabolites]
    redox_pairs = [
        pair for pair in redox_pairs
        if (pair[0] in all_mets) and (pair[1] in all_mets)
    ]
    # we can always find all proton IDs by looking at the numeric bit
    proton_ids = [m.id for m in model.metabolites if '02039' in m.id]
    # now run the tests (individually so we skip the network merging)
    (duplicates, _) = duplicate_test(model, redox_pairs, proton_ids)
    (dead_ends, _) = dead_end_test(model)
    (loops, _) = loop_test(model, threads = threads)
    (dilutions, _) = dilution_test(
        model, dead_ends, media_mets, threads = threads
    )
    # dead-end and dilution test results were already merged
    all_test_results = duplicates.merge(dilutions).merge(loops)
    # just summarize how many reactions were flagged by each test in the output
    simplified_results = simplify_test_results(all_test_results)
    any_test = len(simplified_results[simplified_results.loc[
        :, simplified_results.columns.str.contains('test')
    ].apply(lambda col: col != 'ok', axis = 1).any(axis = 1)])
    deads = (simplified_results['dead_end_test'] != 'ok').sum()
    dils = (simplified_results['dilution_test'] != 'ok').sum()
    dupes = (simplified_results['duplicate_test'] != 'ok').sum()
    loops = (simplified_results['loop_test'] != 'ok').sum()
    pd.DataFrame({
        'model_version' : [f'1.{version}'],
        'all_rxns' : [len(model.reactions)],
        'flagged' : [any_test],
        'dead-ends' : [deads],
        'dilution-blocked' : [dils],
        'duplicates' : [dupes],
        'loops' : [loops],
        'redoxes' : [len(redox_pairs)]
    }).to_csv(
        # append so we get something even if it doesn't finish
        out_fname, mode = 'a', header = not os.path.exists(fname), index = False
    )
    msg = f'Took {time_str(start_time, time.time())} to test version {version}.'
    msg += f' Tested {redoxes} pairs of redox mets.'
    print(msg)
