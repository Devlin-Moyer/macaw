# run_tests_on_all_human.py
'''
Run the duplicate, dead-end, and loop tests on all versions of Human-GEM
'''

import time
from optlang.glpk_interface import Configuration
import sys
from pebble import ProcessPool, ProcessExpired
import cobra
from macaw_main import dead_end_test, duplicate_test, loop_test
from macaw_utils import time_str, simplify_test_results
import pandas as pd

def handle_one_model(version):
    model = cobra.io.read_sbml_model(f'GSMMs/Human-GEMv{version}.xml')
    # call tests separately cuz we're skipping some and pathway-making
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

# set up a Pebble ProcessPool to run tests on all models in parallel
versions = [f'1.{v}' for v in range(19)]
pool = ProcessPool(max_workers = threads)
future = pool.map(handle_one_model, versions)
iterator = future.result()
# prepare a dict to track the number of reactions flagged by each test
out_dict = {
    'version' : list(),
    'all_rxns' : list(),
    'duplicates' : list(),
    'dead-ends' : list(),
    'loops' : list()
}
# keep track of which version we're on
i = 0
print(f'Took {time_str(setup_start, time.time())} to set up')
while True:
    start_time = time.time()
    try:
        # update dict with numbers from the current model
        (all_rxns, dupes, dead_ends, loops) = next(iterator)
        out_dict['version'].append(f'1.{i}')
        out_dict['all_rxns'].append(all_rxns)
        out_dict['duplicates'].append(dupes)
        out_dict['dead-ends'].append(dead_ends)
        out_dict['loops'].append(loops)
        msg = f'Took {time_str(start_time, time.time())} to test version 1.{i}'
        print(msg)
    except (StopIteration, IndexError):
        # should only happen if we've reached the end of the list
        break
    except (Exception, ProcessExpired) as error:
        # if we encounter an error for a particular model, just print which one
        # but continue working on the rest
        msg = f'Error after {time_str(start_time, time.time())} when testing '
        msg += f'version 1.{i}: {error}'
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
pd.DataFrame(out_dict).to_csv('figure_data/fig_all-human.csv', index = False)
