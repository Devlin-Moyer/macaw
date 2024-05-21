# run_tests_on_many.py
'''
Run tests on all models in a particular directory
'''

import time
from optlang.glpk_interface import Configuration
import sys
import logging
import os
from pebble import ProcessPool, ProcessExpired
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
    # call tests separately cuz we're skipping the diphosphate test and don't
    # need to bother creating pathways
    (dead_end_results, dead_end_edges) = dead_end_test(
        model, use_names = True, add_suffixes = True, verbose = 0
    )
    # identify sets of reactions that are potentially duplicates of each other
    (duplicates, dupe_edges) = duplicate_test(
        model, use_names = True, add_suffixes = True, verbose = 0
    )
    # identify reactions that are capable of sustaining non-zero fluxes when all
    # exchange reactions are blocked
    (loops, loop_edges) = loop_test(
        model, use_names = True, add_suffixes = True, verbose = 0
    )
    # identify reactions that become incapable of sustaining non-zero fluxes
    # when dilution constraints are added to the model
    (dilution_results, dilution_edges) = dilution_test(
        model, dead_end_results, use_names = True, add_suffixes = True,
        verbose = 0
    )
    # dead-end and dilution test results were already merged, but get the rest
    all_test_results = duplicates.merge(dilution_results).merge(loops)
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
    return((all_rxns, flagged, deads, dils, dupes, loops, time_msg))

setup_start = time.time()
# silence annoying optlang message that prints when you read in a model
Configuration()
# get number of threads to use from command-line
try:
    (direc, out_fname, threads) = sys.argv[1:]
    threads = int(threads)
except IndexError:
    sys.exit('directory containing models, name of output file, # threads')
# by default, the Cobrapy io functions print a LOT of formatting suggestions for
# several of these models, and we don't need all those messages cuz we are not
# here to edit these models in any way, so make sure Cobrapy only prints errors
logging.getLogger('cobra').setLevel(logging.ERROR)

# get a list of all models in the given directory
model_paths = [
    f'{direc}/{f}' for f in os.listdir(direc)
    if os.path.isfile(f'{direc}/{f}')
]
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
    'loops' : list()
}
# keep track of which model we're on
i = 0
print(f'Took {time_str(setup_start, time.time())} to set up')
iterator = future.result()
while True:
    try:
        # update dict with numbers from the current model
        (all_rxns, flagged, deads, dils, dupes, loops, msg) = next(iterator)
        model_name = model_paths.split('/')[-1]
        # if we're doing all versions of Human-GEM, just use the version number
        # as the model name
        if 'Human-GEM' in model_name:
            model_name = model_name.split('v')[-1].split('.xml')[0]
        out_dict['model'].append(model_name)
        out_dict['all_rxns'].append(all_rxns)
        out_dict['flagged'].append(flagged)
        out_dict['dead-ends'].append(deads)
        out_dict['dilution-blocked'].append(dils)
        out_dict['duplicates'].append(dupes)
        out_dict['loops'].append(loops)
        msg = f'Took {msg} to test {model_name} (model {i+1} out of '
        msg += f'{len(model_paths)})'
        print(msg)
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
# turn the dict into a Pandas DataFrame and write to CSV
pd.DataFrame(out_dict).to_csv(f'figure_data/fig_{out_fname}.csv', index = False)
