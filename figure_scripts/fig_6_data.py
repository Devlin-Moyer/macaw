# fig_6_data.py
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
    # the different models from the Mendoza paper used a variety of formats
    # for their metabolite IDs (BiGG, ModelSeed, KEGG, etc.), so include many
    # IDs for each pair then filter down to the ones actually present in the
    # given model
    redox_pairs = [
        # NAD(H)
        ('nad', 'nadh'), ('C00003', 'C00004'), ('cpd00003', 'cpd00004'),
        ('NAD', 'NADH'),
        # NADP(H)
        ('nadp', 'nadph'), ('C00006', 'C00005'), ('cpd00006', 'cpd00005'),
        ('NADP', 'NADPH'),
        # FAD(H2)
        ('fad', 'fadh2'), ('C00016', 'C01352'), ('cpd00015', 'cpd00982'),
        ('FAD', 'FADH2'),
        # FMN(H2)
        ('fmn', 'fmnh2'), ('C00061', 'C01847'), ('cpd00050', 'cpd01270'),
        ('FMN', 'FMNH2'),
        # riboflavin
        ('ribflv', 'rbflvrd'), ('C00255', 'C01007'), ('cpd00220', 'cpd00739'),
        ('RIBOFLAVIN', 'CPD-316'),
        # ubiquinones and ubiquinols
        ('q', 'qh2'), ('C00399', 'C00390'), ('cpd28301', 'cpd28300'),
        ('cpd25985', 'cpd26399'), ('CPD0-1118', 'CPD0-2061'),
        ('cpd28290', 'cpd28015'), ('cpd11669', 'cpd11665'),
        ('UBIQUINONE-2', 'QH2'), ('q6', 'q6h2'), ('cpd15290', 'cpd15291'),
        ('UBIQUINONE-6', 'UBIQUINOL-30'), ('cpd25792', 'cpd25913'),
        ('CPD-9717', 'CPD-9955'), ('q8', 'q8h2'), ('cpd15560', 'cpd15561'),
        ('cpd15560', 'cpd29608'), ('UBIQUINONE-8', 'CPD-9956'), ('q9', 'q9h2'),
        ('cpd01351', 'cpd25914'), ('UBIQUINONE-9', 'CPD-9957'),
        ('q10', 'q10h2'), ('cpd08232', 'cpd25915'),
        ('UBIQUINONE-10', 'CPD-9958'),
        # glutathione
        ('gthox', 'gthrd'), ('C00127', 'C00051'), ('cpd00111', 'cpd00042'),
        ('OXIDIZED-GLUTATHIONE', 'GLUTATHIONE'),
        # oxygen and hydrogen peroxide
        ('o2', 'h2o2'), ('C00007', 'C00027'), ('cpd00007', 'cpd00025'),
        ('OXYGEN-MOLECULE', 'HYDROGEN-PEROXIDE')
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
    # now do the tests
    (dupes, _) = duplicate_test(model, redox_pairs, proton_ids, verbose = 0)
    (dead_ends, _) = dead_end_test(model, verbose = 0)
    (loops, _) = loop_test(model, verbose = 0)
    (dils, _) = dilution_test(model, dead_ends, verbose = 0)
    # dead-end and dilution test results were already merged, but get the rest
    all_test_results = dupes.merge(dils).merge(loops)
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
# if tot_batches > 1, split the list of paths into tot_batches (approximately)
# equally large groups and only test the models in the specified batch (index)
model_paths = np.array_split(
    model_paths, tot_batches
)[batch_idx - 1].tolist() # SGE task IDs are 1-indexed, Python is 0-indexed
# also update the output file name appropriately
if tot_batches > 1:
    out_fname += f'_batch-{batch_idx}-of-{tot_batches}.csv'
else:
    out_fname += '.csv'
out_dir = 'figure_data'
out_fname = f'{out_dir}/{out_fname}'
# skip any models we already have test results for
already_done = set()
for f in os.listdir(out_dir):
    if f.startswith(out_fname):
        already_done.update(pd.read_csv(f'{out_dir}/{f}')['model'].to_list())
model_paths = [p for p in model_paths if p not in already_done]

# set up a Pebble ProcessPool to run tests on all models in parallel
pool = ProcessPool(max_workers = threads)
future = pool.map(handle_one_model, model_paths)
# keep track of which model we're on
i = 0
iterator = future.result()
while True:
    try:
        # update dict with numbers from the current model
        (rxns, flagged, deads, dils, dupes, loops, rdx, msg) = next(iterator)
        model_name = model_paths[i].split('/')[-1]
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
            out_fname,
            mode = 'a',
            header = not os.path.exists(out_fname),
            index = False
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
