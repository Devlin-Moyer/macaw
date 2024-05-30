# run_tests_on_many.py
'''
Run tests on all models in a particular directory
'''

from optlang.glpk_interface import Configuration
import sys
import logging
import os
import numpy as np
from pebble import ProcessPool, ProcessExpired
import time
import cobra
from macaw_main import dead_end_test, dilution_test, duplicate_test, loop_test
from macaw_utils import time_str, simplify_test_results
import pandas as pd

def handle_one_model(model_path):
    start_time = time.time()
    # annoyingly, some of the models were only available as .mats and others
    # were only available as ".sbml" files, and as far as I know, Cobrapy
    # doesn't have a single input function that automatically figures out the
    # file type from the extension
    if model_path.endswith('mat'):
        model = cobra.io.load_matlab_model(model_path)
    else:
        model = cobra.io.read_sbml_model(model_path)
    # get the appropriate set of pairs of IDs for redox carrier metabolites
    if 'Human-GEM' in model_path:
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
        # see if this version is from before they reformatted metabolite IDs
        try:
            model.metabolites.get_by_id('MAR02552c')
        except KeyError:
            redox_pairs = [
                (pair[0].replace('MAR', 'm'), pair[1].replace('MAR', 'm'))
                for pair in redox_pairs
            ]
        proton_ids = [m.id for m in model.metabolites if '02039' in m.id]
    else:
        # everything other than Human-GEM uses a variety of different formats
        # for metabolite IDs (BiGG, ModelSEED, KEGG, MetaCyc) but doesn't have
        # any subcellular compartments
        redox_pairs = [
            # NAD(H)
            ('nad', 'nadh'), ('C00003', 'C00004'),
            ('cpd00003', 'cpd00004'), ('NAD', 'NADH'), ('00010', '00013'),
            # NADP(H)
            ('nadp', 'nadph'), ('C00006', 'C00005'),
            ('cpd00006', 'cpd00005'), ('NADP', 'NADPH'), ('00002', '00005'),
            # FAD(H2)
            ('fad', 'fadh2'), ('C00016', 'C01352'),
            ('cpd00015', 'cpd00982'), ('FAD', 'FADH2'), ('00251', '00247'),
            # FMN(H2)
            ('fmn', 'fmnh2'), ('C00061', 'C01847'),
            ('cpd00050', 'cpd01270'), ('FMN', 'FMNH2'), ('00622', '00575'),
            # riboflavin
            ('ribflv', 'rbflvrd'), ('C00255', 'C01007'),
            ('cpd00220', 'cpd00739'), ('RIBOFLAVIN', 'CPD-316'),
            # ubiquinones and ubiquinols
            ('q', 'qh2'), ('C00399', 'C00390'), ('cpd28301', 'cpd28300'),
            ('01220', '01121'), ('cpd25985', 'cpd26399'),
            ('CPD0-1118', 'CPD0-2061'), ('cpd28290', 'cpd28015'),
            ('cpd11669', 'cpd11665'), ('UBIQUINONE-2', 'QH2'),
            ('q6', 'q6h2'), ('cpd15290', 'cpd15291'),
            ('UBIQUINONE-6', 'UBIQUINOL-30'), ('cpd25792', 'cpd25913'),
            ('CPD-9717', 'CPD-9955'), ('q8', 'q8h2'), ('cpd15560', 'cpd15561'),
            ('cpd15560', 'cpd29608'), ('UBIQUINONE-8', 'CPD-9956'),
            ('q9', 'q9h2'), ('cpd01351', 'cpd25914'),
            ('UBIQUINONE-9', 'CPD-9957'), ('q10', 'q10h2'),
            ('cpd08232', 'cpd25915'), ('UBIQUINONE-10', 'CPD-9958'),
            # glutathione
            ('gthox', 'gthrd'), ('C00127', 'C00051'), ('cpd00111', 'cpd00042'),
            ('OXIDIZED-GLUTATHIONE', 'GLUTATHIONE'),
            # oxygen and hydrogen peroxide
            ('o2', 'h2o2'), ('C00007', 'C00027'), ('cpd00007', 'cpd00025'),
            ('OXYGEN-MOLECULE', 'HYDROGEN-PEROXIDE'), ('00249', '00410')
        ]
        proton_ids = ['h', 'C00080', 'cpd00067', 'PROTON']
        # some models added compartment suffixes and some didn't, so try
        # different compartment suffixes until at least one of the proton IDs is
        # a metabolite ID that's actually in the model
        all_ids = [m.id for m in model.metabolites]
        if not any(p in all_ids for p in proton_ids):
            for suffix in ['_c', '_c0', '[c]']:
                new_proton_ids = [p + suffix for p in proton_ids]
                if any(p in all_ids for p in new_proton_ids):
                    # update the lists of proton IDs and redox carrier pairs
                    proton_ids = new_proton_ids
                    redox_pairs = [
                        (m1 + suffix, m2 + suffix) for (m1, m2) in redox_pairs
                    ]
                    break
    # it's okay if the given model only contains metabolites with some of these
    # IDs; just filter them out and report how many were left
    all_mets = [m.id for m in model.metabolites]
    redox_pairs = [
        pair for pair in redox_pairs
        if (pair[0] in all_mets) and (pair[1] in all_mets)
    ]
    proton_ids = [m for m in proton_ids if m in all_mets]
    # report how many of the redox IDs were actually in the model
    redoxes = len(redox_pairs)
    # identify sets of reactions that are potentially duplicates of each other
    (duplicates, dupe_edges) = duplicate_test(
        model, redox_pairs, proton_ids, verbose = 0
    )
    # call tests separately cuz we're skipping the diphosphate test and don't
    # need to bother creating pathways
    (dead_ends, dead_end_edges) = dead_end_test(model, verbose = 0)
    # identify reactions that are capable of sustaining non-zero fluxes when all
    # exchange reactions are blocked
    (loops, loop_edges) = loop_test(model, verbose = 0)
    # identify reactions that become incapable of sustaining non-zero fluxes
    # when dilution constraints are added to the model
    (dilutions, dil_edges) = dilution_test(model, dead_ends, verbose = 0)
    # dead-end and dilution test results were already merged, but get the rest
    all_test_results = duplicates.merge(dilutions).merge(loops)
    # get the number of all reactions in the model, the number flagged by any
    # test, and the number flagged by each individual test
    all_rxns = len(model.reactions)
    simplified_results = simplify_test_results(all_test_results)
    flagged = len(simplified_results[simplified_results.loc[
        :, simplified_results.columns.str.contains('test')
    ].apply(lambda col: col != 'ok', axis = 1).any(axis = 1)])
    deads = (simplified_results['dead_end_test'] != 'ok').sum()
    dils = (simplified_results['dilution_test'] != 'ok').sum()
    dupes = (simplified_results['duplicate_test'] != 'ok').sum()
    loops = (simplified_results['loop_test'] != 'ok').sum()
    time_msg = time_str(start_time, time.time())
    return((all_rxns, flagged, deads, dils, dupes, loops, redoxes, time_msg))

# silence annoying optlang message that prints when you read in a model
Configuration()
# get number of threads to use from command-line
try:
    (direc, out_fname, threads, batch_idx, tot_batches) = sys.argv[1:]
    threads = int(threads)
    batch_idx = int(batch_idx)
    tot_batches = int(tot_batches)
except ValueError:
    msg = 'directory containing models, name of output file, # threads, batch '
    msg += 'index, total number of batches (just do 1 for both of the last two '
    msg += 'for no batching)'
    sys.exit(msg)
# by default, the Cobrapy io functions print a LOT of formatting suggestions for
# several of these models, and we don't need all those messages cuz we are not
# here to edit these models in any way, so make sure Cobrapy only prints errors
logging.getLogger('cobra').setLevel(logging.ERROR)

# get a list of all models in the given directory
model_paths = [
    f'{direc}/{f}' for f in os.listdir(direc)
    if os.path.isfile(f'{direc}/{f}')
]
# see if there's already an output file with the appropriate name and skip any
# models we already have results for in that file
fname = f'figure_data/fig_{out_fname}.csv'
if tot_batches > 1:
    fname = fname.replace('.csv', f'_batch-{batch_idx}-of-{tot_batches}.csv')
if os.path.exists(fname):
    already_done = pd.read_csv(fname)['model'].to_list()
    model_paths = [p for p in model_paths if p not in already_done]
# if tot_batches > 1, split the list of paths into tot_batches (approximately)
# equally large groups and only test the models in the specified batch (index)
model_paths = np.array_split(
    model_paths, tot_batches
)[batch_idx - 1].tolist() # SGE task IDs are 1-indexed, Python is 0-indexed

# set up a Pebble ProcessPool to run tests on all models in parallel
pool = ProcessPool(max_workers = threads)
future = pool.map(handle_one_model, model_paths)
# prepare a dict to track the number of reactions flagged by each test
out_dict = {
    'model' : list(),
    'all_rxns' : list(),
    'flagged' : list(),
    'dead-ends' : list(),
    'dilution-blocked' : list(),
    'duplicates' : list(),
    'loops' : list(),
    'redoxes' : list()
}
# keep track of which model we're on
i = 0
iterator = future.result()
while True:
    try:
        # update dict with numbers from the current model
        (rxns, flagged, deads, dils, dupes, loops, rdx, msg) = next(iterator)
        model_name = model_paths[i].split('/')[-1]
        # if we're doing all versions of Human-GEM, just use the version number
        # as the model name
        if 'Human-GEM' in model_name:
            model_name = model_name.split('v')[-1].split('.xml')[0]
        msg = f'Took {msg} to test {model_name} (model {i+1} out of '
        msg += f'{len(model_paths)}). Tested {rdx} pairs of redox mets.'
        print(msg)
        # append to the output file so that we save output as it becomes
        # available instead of having to wait until the end to get anything
        pd.DataFrame({
            'model' : [model_name],
            'all_rxns' : [rxns],
            'flagged' : [flagged],
            'dead-ends' : [deads],
            'dilution-blocked' : [dils],
            'duplicates' : [dupes],
            'loops' : [loops],
            'redoxes' : [rdx]
        }).to_csv(
            fname, mode = 'a', header = not os.path.exists(fname), index = False
        )
    except StopIteration:
        # should only happen if we've reached the end of the list
        break
    except (Exception, ProcessExpired) as error:
        # if we encounter an error for a particular model, just print which one
        # but continue working on the rest
        msg = f'Error when testing {model_name} (model {i+1} out of '
        msg += f'{len(model_paths)}): {error}'
        # some errors don't come with tracebacks, but if this does, print it
        if hasattr(error, 'traceback'):
            msg += f'\n{error.traceback}'
        print(msg)
    finally:
        # make sure we always increment the iterator
        i += 1
# close the ProcessPool
pool.close()
pool.join()
