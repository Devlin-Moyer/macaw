# run_tests_on_agora.py

import sys
from optlang.glpk_interface import Configuration
import logging
import os
from pebble import ProcessPool, ProcessExpired
import time
import cobra
from macaw_main import dead_end_test, duplicate_test, loop_test
from macaw_utils import time_str, simplify_test_results
import pandas as pd

def handle_one_model(model_path):
    model = cobra.io.read_sbml_model(f'GSMMs/AGORA2/{model_path}')
    # TODO: figure out if they all have ATP maintenance reactions
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
    # skip dilution test cuz it takes much longer than others and there are
    # 7,000 AGORA2 models
    all_test_results = duplicates.merge(dead_end_results).merge(loops)
    # just return number of reactions flagged by each test
    simplified_results = simplify_test_results(all_test_results)
    dupes = (simplified_results['duplicate_test'] != 'ok').sum()
    dead_ends = (simplified_results['dead_end_test'] != 'ok').sum()
    loops = (simplified_results['loop_test'] != 'ok').sum()
    all_rxns = len(model.reactions)
    return((all_rxns, dupes, dead_ends, loops))

setup_start = time.time()
# silence annoying optlang message that prints when you read in a model
Configuration()
# get number of threads to use from command-line
try:
    threads = int(sys.argv[1])
except IndexError:
    sys.exit('provide number of threads to use')
# by default, the Cobrapy io functions print a LOT of formatting suggestions for
# several of these models, and we don't need all those messages cuz we are not
# here to edit these models in any way, so make sure Cobrapy only prints errors
logging.getLogger('cobra').setLevel(logging.ERROR)

# set up a Pebble ProcessPool to run tests on all models in parallel
model_paths = os.listdir('GSMMs/AGORA2')
pool = ProcessPool(max_workers = threads)
future = pool.map(handle_one_model, model_paths)
iterator = future.result()
# prepare a dict to track the number of reactions flagged by each test
out_dict = {
    'model' : list(),
    'all_rxns' : list(),
    'duplicates' : list(),
    'dead-ends' : list(),
    'loops' : list()
}
# keep track of which model we're on
i = 0
print(f'Took {time_str(setup_start, time.time())} to set up')
while True:
    start_time = time.time()
    try:
        # update dict with numbers from the current model
        (all_rxns, dupes, dead_ends, loops) = next(iterator)
        out_dict['model'].append(model_paths[i])
        out_dict['all_rxns'].append(all_rxns)
        out_dict['duplicates'].append(dupes)
        out_dict['dead-ends'].append(dead_ends)
        out_dict['loops'].append(loops)
        msg = f'Took {time_str(start_time, time.time())} to test '
        msg += f'{model_paths[i]} (model {i+1} out of {len(model_paths)})'
        print(msg)
    except (StopIteration, IndexError):
        # should only happen if we've reached the end of the list
        break
    except (Exception, ProcessExpired) as error:
        # if we encounter an error for a particular model, just print which one
        # but continue working on the rest
        msg = f'Error after {time_str(start_time, time.time())} when testing '
        msg += f'{model_paths[i]} (model {i+1} out of {len(model_paths)}): '
        msg += f'{error}'
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
pd.DataFrame(out_dict).to_csv('figure_data/fig_agora.csv', index = False)
