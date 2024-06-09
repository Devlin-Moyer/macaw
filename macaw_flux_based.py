# macaw_flux_based.py
'''
Tests that involve predicting fluxes through reactions as opposed to just
looking at the structure of the network
'''

from optlang.symbolics import Zero
from macaw_fva import fva
from macaw_utils import add_reaction_equations, edit_dead_end_bounds
import cobra
from macaw_structural import dead_end_test, _dead_end_test_inner
from macaw_dilution import constrain_reversible_rxns, add_leakage_reactions
from pebble import ProcessPool, ProcessExpired
from concurrent.futures import TimeoutError
from macaw_dilution import add_dilution_constraints
import pandas as pd
from optlang.interface import OPTIMAL, UNBOUNDED
import math

def loop_test(
    given_model, zero_thresh = 10**-8, corr_thresh = 0.9, use_names = False,
    add_suffixes = False, threads = 1, verbose = 1
):
    '''
    Do FVA on the given model after setting both bounds on all exchange
    reactions to 0 to identify all reactions that are involved in fully internal
    "loops" of reactions that can sustain fluxes independently of the rest of
    the network (also known as type III extreme pathways)
    Remove any objective functions and any positive non-zero lower bounds on
    reactions before starting to ensure that the model is feasible with no
    exchange fluxes
    Generate 1,000 possible solutions to the GSMM using OptGP, get pairwise
    correlations between all reactions' distributions of possible fluxes, and
    create an edge list connecting all pairs of reactions with correlations
    higher than corr_thresh that also share at least one metabolite (in an
    attempt to separate out reactions that participate in different loops)
    '''
    if verbose > 0:
        print('Starting loop test...')
    # work with a copy of the given model
    model = given_model.copy()
    # make sure model has no objective function cuz most objective functions
    # are unlikely to be feasible with no exchange fluxes
    model.objective = Zero
    for r in model.reactions:
        # set both bounds on all exchange reactions to 0 so the only reactions
        # that will be capable of fluxes will be in internal loops
        if r.boundary:
            r.lower_bound = 0
            r.upper_bound = 0
        if r.lower_bound > 0:
            # look for reactions with positive non-zero lower bounds cuz those
            # are probably ATP maintenance reactions and will probably make the
            # model infeasible with no exchange fluxes and stop us from finding
            # loops
            r.lower_bound = 0
    # do FVA and highlight all reactions with a non-zero minimum or maximum
    fva_results = fva(
        model, zero_thresh = zero_thresh, threads = threads,
        verbose = verbose - 1
    ).reset_index(names = 'reaction_id')
    fva_results['loop_test'] = fva_results.apply(
        lambda row: 'in loop'
            if (row['minimum'] != 0) or (row['maximum'] != 0) else 'ok',
        axis = 1
    )
    # remove all reactions that couldn't have flux and then get 1,000 possible
    # solutions to the remaining model so we can look at correlations between
    # those possible fluxes to separate out all the flagged reactions into the
    # different loops they comprise
    model.remove_reactions(
        fva_results[fva_results['loop_test'] == 'ok']['reaction_id'].to_list()
    )
    solutions = cobra.sampling.sample(model, 1000, processes = threads)
    # round fluxes that are suspiciously close to zero to make sure we
    # definitely skip all reactions that could never have flux
    solutions = solutions.applymap(lambda x: 0 if abs(x) < zero_thresh else x)
    non_zero_rxn_fluxes = solutions[
        solutions.columns[(solutions != 0).all(axis = 0)]
    ]
    # get correlation matrix
    corr_mat = non_zero_rxn_fluxes.corr()
    # get list of tuples of reaction IDs that had |correlations| > 0.9
    corr_mat.index = corr_mat.columns
    corr_pairs = corr_mat.stack().reset_index(name = 'correlation')
    edge_df = corr_pairs[
        # skip rows for the correlation of a reaction with itself
        (corr_pairs['level_0'] != corr_pairs['level_1']) &
        (corr_pairs['correlation'].abs() > corr_thresh)
    ]
    edge_list_bad = list(zip(
        edge_df['level_0'].to_list(), edge_df['level_1'].to_list()
    ))
    # filter out edges between reactions that don't share any metabolites cuz
    # some loops can involve dozens if not hundreds of reactions and we def
    # don't want to connect every single pair of reactions in those
    edge_list = list()
    for (r1, r2) in edge_list_bad:
        r1_mets = model.reactions.get_by_id(r1).metabolites
        r2_mets = model.reactions.get_by_id(r2).metabolites
        if set(r1_mets).intersection(set(r2_mets)):
            edge_list.append((r1, r2))
    # don't include the minimum and maximum columns in the output and add
    # a column for the reaction equations
    out_df = add_reaction_equations(
        fva_results[['reaction_id', 'loop_test']], given_model,
        use_names = use_names, add_suffixes = add_suffixes
    )
    if verbose > 0:
        count = (fva_results['loop_test'] == 'in loop').sum()
        print(f' - Found {count} reactions involved in infinite loops.')
    return((out_df, edge_list))

def dilution_test(
    given_model, dead_end_results = None, media_mets = None,
    zero_thresh = 10**-8, timeout = 1800, max_attempts = 3, use_names = False,
    add_suffixes = False, verbose = 1, threads = 1
):
    '''
    Use dead-end test results to set both bounds to zero for all reactions that
    were part of dead-ends (sometimes adding dilution reactions allows them to
    carry fluxes, which significantly increases the space of possible solutions
    to the model and thus the time it takes to do FVA, so blocking them first
    speeds things up significantly) and to set the appropriate bounds to zero
    for the reversible reactions that couldn't go in just one of their two
    directions due to dead-ends.
    Then do FVA to determine the range of feasible fluxes for each reaction both
    before and after imposing dilution constraints. Report all reactions whose
    fluxes were only capable of being zero with dilution constraints but could
    be non-zero without, zero without dilution constraints but non-zero with, or
    only zero either way.
    GLPK sometimes hangs when doing FVA on models with dilution constraints, so
    the timeout argument is the number of seconds to wait for the minimum or
    maximum flux for a single reaction before killing that process and starting
    over (each reaction is handled on a separate process through Pebble)
    Any minimum or maximum fluxes within zero_thresh of zero are rounded to zero
    to handle rounding errors
    '''
    model = given_model.copy()
    if dead_end_results is None:
        if verbose > 0:
            msg = 'Output of dead_end_test was not provided to dilution_test, '
            msg += 'so dead_end_test is being run on the given model in '
            msg += 'preparation for running dilution_test.'
            print(msg)
        (dead_end_results, dead_end_edges) = dead_end_test(
            model, use_names, add_suffixes, verbose
        )
    # set both bounds to zero for all reactions found to be dead-ends
    # also set appropriate bound to zero for reversible reactions found to be
    # structurally prevented from carrying flux in one direction
    model = edit_dead_end_bounds(model, dead_end_results)
    if verbose > 0:
        print('Starting dilution test...')
    # if given a list of metabolites we should allow uptake of, find the
    # exchange reactions for those metabolites and set their lower bounds to
    # -1000 and set the lower bounds on all other exchange reactions to 0
    if media_mets is not None:
        i = 0
        for r in model.boundary:
            met = list(r.metabolites)[0]
            id_match = any(m == met.id for m in media_mets)
            name_match = any(m == met.name for m in media_mets)
            if id_match or name_match:
                r.lower_bound = -1000
                i += 1
            else:
                r.lower_bound = 0
        if verbose > 0:
            msg = ' - Adjusted lower bounds on exchange reactions so that only '
            msg += f'{i} metabolites could be consumed.'
            if i < len(media_mets):
                msg += f' {len(media_mets) - i} of the metabolite IDs provided '
                msg += 'in the media_mets argument were not in the GSMM.'
            print(msg)
    # ensure the model has no objective function
    model.objective = Zero
    if verbose > 1:
        msg = ' - Finding ranges of possible fluxes for all reactions in model '
        msg += 'without dilution constraints...'
        print(msg)
    fva_before = fva(
        model, zero_thresh = zero_thresh, threads = threads,
        verbose = verbose - 1
    )
    # see docstrings for what these functions do
    model = constrain_reversible_rxns(
        model, fva_before, zero_thresh, verbose = verbose - 1
    )
    model = add_leakage_reactions(model, verbose = verbose - 1)
    # prepare dicts with all metabolite and reaction IDs as keys and strings as
    # values indicating whether or not those metabolites and reactions were
    # blocked by dilution constraints
    met_dict = {m.id : '' for m in given_model.metabolites}
    rxn_dict = {r.id : '' for r in given_model.reactions}
    mets_to_test = list(met_dict.keys())
    # set up a Pebble ProcessPool to add a dilution reaction and constraint for
    # each metabolite in parallel so we can see which dilution constraints block
    # fluxes through which reactions
    pool = ProcessPool(
        max_workers = threads,
        initializer = _pool_init,
        initargs = (model, zero_thresh)
    )
    future = pool.map(dilution_test_inner, mets_to_test, timeout = timeout)
    iterator = future.result()
    # do a nested while loop so that we can retry metabolites we encounter
    # errors for on the first try cuz optimizing models with dilution
    # constraints sometimes leads to rare intermittent errors
    attempts = 0
    fva_before = fva_before.rename(
        columns = {'minimum' : 'min_before', 'maximum' : 'max_before'}
    ).reset_index(names = 'reaction_id')
    while mets_to_test:
        attempts += 1
        to_redo = list()
        # we'll get output in the same order as the metabolite IDs we mapped to 
        # the pool, so keep track of which metabolite we're on
        i = 0
        while True:
            try:
                fva_after = next(iterator)
                # make sure all mins and maxes in fva_after weren't NA
                all_mins_na = fva_after['min_after'].isna().all()
                all_maxs_na = fva_after['max_after'].isna().all()
                # compare FVA results for the reactions involving this
                # metabolite from before and after adding its dilution
                # constraint
                fva_both = fva_after.merge(fva_before)
                min_zero_before = (fva_both['min_before'] == 0).all()
                max_zero_before = (fva_both['max_before'] == 0).all()
                both_zero_before = min_zero_before and max_zero_before
                min_zero_after = (fva_both['min_after'] == 0).all()
                max_zero_after = (fva_both['max_after'] == 0).all()
                both_zero_after = min_zero_after and max_zero_after
                if all_mins_na or all_maxs_na:
                    result = 'error'
                elif both_zero_after and not both_zero_before:
                    result = 'blocked by dilution'
                elif both_zero_before and not both_zero_after:
                    result = 'unblocked by dilution'
                elif both_zero_before and both_zero_after:
                    result = 'always blocked'
                else:
                    result = 'ok'
                # add this result for this metabolite and all the reactions it
                # participates in
                met_dict[mets_to_test[i]] = result
                met_obj = given_model.metabolites.get_by_id(mets_to_test[i])
                for r in met_obj.reactions:
                    # if we've already labeled this reaction as "blocked by
                    # dilution", then don't overwrite it with an "ok" because a
                    # different metabolite that participates in this reaction
                    # isn't blocked by its dilution constraint
                    if rxn_dict[r.id] != 'blocked by dilution':
                        # then see if this reaction was always blocked, cuz it
                        # is possible for some but not all reactions that a
                        # particular metabolite participates in to always be
                        # blocked
                        one_fva = fva_both[fva_both['reaction_id'] == r.id]
                        if (
                            (one_fva['min_before'] == 0) &
                            (one_fva['max_before'] == 0) &
                            (one_fva['min_after'] == 0) &
                            (one_fva['max_after'] == 0)
                        ).all():
                            rxn_dict[r.id] = 'always blocked'
                        else:
                            rxn_dict[r.id] = result
            except (StopIteration, IndexError):
                # should only happen if we've reached the end of the list
                break
            except TimeoutError as error:
                # doing FBA on models with dilution constraints sometimes makes
                # Python/Cobrapy/optlang/GLPK/who knows what hang; add it to the
                # list of metabolites to redo
                to_redo.append(mets_to_test[i])
                if verbose > 1:
                    met_name = model.metabolites.get_by_id(mets_to_test[i]).name
                    if timeout >= 60:
                        t = f'{int(timeout/60)} minutes'
                    else:
                        t = f'{timeout} seconds'
                    msg = f' - Took longer than {t} to test dilution '
                    msg += 'constraint for {met_name} ({mets_to_test[i]}).'
                    print(msg)
            except ProcessExpired as error:
                # same thing for other kinds of errors; just a different message
                to_redo.append(mets_to_test[i])
                if verbose > 1:
                    met_name = model.metabolites.get_by_id(mets_to_test[i]).name
                    msg = ' - Encountered an error when testing dilution '
                    msg += f'constraint for {met_name} ({mets_to_test[i]}): '
                    msg += f'{error}'
                    # some errors don't come with tracebacks, but if they do,
                    # print them cuz those are useful
                    if hasattr(error, 'traceback'):
                        msg += f'\n{error.traceback}'
                    print(msg)
            finally:
                # increment iterator regardless of success or failure
                i += 1
        # make sure every single metabolite's result wasn't N/A; that always
        # means something is very wrong
        if all(v == 'error' for v in met_dict.values()):
            msg = 'Every single dilution constraint rendered the model '
            msg += 'infeasible; something is seriously wrong with either the '
            msg += 'code or the model.'
            raise ValueError(msg)
        elif all(v == '' for v in met_dict.values()):
            msg = 'There was an error for every single metabolite; increase '
            msg += 'verbosity if you see no error messages'
            raise ValueError(msg)
        # now see if we have any metabolites to redo
        if to_redo:
            mets_to_test = to_redo
            if verbose > 0:
                msg = ' - Encountered errors when testing dilution constraints '
                msg += f'for {len(to_redo)} metabolites'
                if attempts < max_attempts:
                    msg += '; trying to test them again.'
                else:
                    msg += '; result will be "error" for those reactions'
                    for m in met_dict.keys():
                        if met_dict[m] == '':
                            met_dict[m] = 'error'
                        for r in model.metabolites.get_by_id(m).reactions:
                            if rxn_dict[r] == '':
                                rxn_dict[r] = 'error'
                print(msg)
            # avoid infinite loops if some errors are not in fact intermittent
            if attempts >= max_attempts:
                break
        else:
            if (attempts > 1) and (verbose > 0):
                msg = ' - Successfully tested all dilution constraints after '
                msg += f'{attempts} tries.'
                print(msg)
            break
    pool.close()
    pool.join()
    # make lists of blocked metabolites and reactions
    blocked_mets = [
        given_model.metabolites.get_by_id(m) for m in met_dict.keys()
        if met_dict[m] == 'blocked by dilution'
    ]
    blocked_rxns = [
        given_model.reactions.get_by_id(r) for r in rxn_dict.keys()
        if (rxn_dict[r] == 'blocked by dilution') and ('leakage' not in r)
    ]
    # count number of metabolites that are blocked by their own dilution
    # constraint and not by an upstream or downstream metabolite's constraint
    direct_mets = len(blocked_mets)
    # use dead-end test algorithm to figure out if any of these blocked_mets
    # also block reactions beyond those that directly involve them
    for met in given_model.metabolites:
        _dead_end_test_inner(met, blocked_mets, blocked_rxns, list(), list())
    # remove all metabolites and reactions from these lists that were already
    # already flagged as dead-ends cuz _dead_end_test_inner will find those in
    # addition to other reactions indirectly blocked by dilution constraints
    # for metabolites that do not participate in those reactions
    dead_end_mets = set()
    for result in dead_end_results['dead_end_test']:
        if not result.startswith('only') and (result != 'ok'):
            # dead-end reactions have semicolon-delimited lists of all dead-end
            # metabolites that participate in them
            dead_end_mets.update(result.split(';'))
    blocked_mets = [m.id for m in blocked_mets if m.id not in dead_end_mets]
    blocked_rxns = [
        r.id for r in blocked_rxns
        if all(m.id not in dead_end_mets for m in r.metabolites)
    ]
    # just as we did with dead-end test, make a list of the blocked metabolites
    # for each blocked reaction to make it clear why each reaction was blocked
    results_dict = dict()
    edge_list = list()
    for r in given_model.reactions:
        if r.id in blocked_rxns:
            mets = [m.id for m in r.metabolites if m.id in blocked_mets]
            results_dict[r.id] = ';'.join(mets)
            edge_list.extend([(m, r.id) for m in mets])
        else:
            results_dict[r.id] = rxn_dict[r.id]
    dilution_test_results = pd.DataFrame({
        'reaction_id' : [r.id for r in model.reactions]
    })
    dilution_test_results['dilution_test'] = dilution_test_results[
        'reaction_id'
    ].map(results_dict)
    # merge with the dead-end test results before returning
    out_df = dilution_test_results.merge(dead_end_results)
    # reorder columns for consistency with other tests
    out_df = out_df[[
        'reaction_id', 'reaction_equation', 'dead_end_test', 'dilution_test'
    ]]
    # print a summary of the findings
    if verbose > 0:
        msg = f' - Found {direct_mets} metabolites for which adding a '
        msg += 'dilution constraint prevented fluxes through all associated '
        msg += f'reactions.\n - Found {len(blocked_rxns)} reactions that '
        msg += 'are blocked by one or more dilution constraints.'
        print(msg)
    # don't need to add a column for reaction equations cuz it's already in the
    # dead-end test results
    return((out_df, edge_list))

def _pool_init(m, z):
    global model, zero_thresh
    (model, zero_thresh) = (m, z)

def dilution_test_inner(met_id):
    '''
    Create a dilution reaction and constraint for a single metabolite in the
    given model (expected as a global variable)
    Do FVA for all reactions the given metabolite participates in
    Round all FVA results within zero_thresh (expected as a global variable)
    to zero to handle rounding errors
    Don't use macaw_fva.fva() cuz that creates a ProcessPool and nesting
    those seems to dramatically increase the amount of memory used
    '''
    # create a dilution reaction and constraint for this metabolite
    dil_model = add_dilution_constraints(
        # don't need to do FVA and set bounds on reversible reactions or add
        # leakage reactions cuz we already did both of those things
        model, [met_id], preprocess = False, leak_flux = 0, verbose = 0
    )
    dil_model.objective = Zero
    # now do FVA for all reactions this metabolite participates in (except the
    # leakage reactions cuz those were just added and not in the original model)
    rxns = [
        r for r in dil_model.metabolites.get_by_id(met_id).reactions
        if ('leakage' not in r.id) and ('dilution' not in r.id)
    ]
    dil_model.solver.objective.direction = 'min'
    mins = [handle_one_reaction(dil_model, r) for r in rxns]
    dil_model.solver.objective.direction = 'max'
    maxes = [handle_one_reaction(dil_model, r) for r in rxns]
    fva_results = pd.DataFrame({
        'reaction_id' : [r.id for r in rxns],
        'min_after' : mins,
        'max_after' : maxes
    })
    fva_results['min_after'] = fva_results['min_after'].apply(
        lambda x: 0 if abs(x) < zero_thresh else x
    )
    fva_results['max_after'] = fva_results['max_after'].apply(
        lambda x: 0 if abs(x) < zero_thresh else x
    )
    return(fva_results)

def handle_one_reaction(model, reaction):
    model.solver.objective.set_linear_coefficients({
        reaction.forward_variable : 1, reaction.reverse_variable : -1
    })
    model.slim_optimize()
    # handle non-optimal solver statuses
    if model.solver.status == OPTIMAL:
        obj_val = model.solver.objective.value
    elif model.solver.status == UNBOUNDED:
        # this means the objective value was infinite
        if model.solver.objective.direction == 'min':
            obj_val = -float('inf')
        else:
            obj_val = float('inf')
    else:
        # if it wasn't optimal or unbounded, use NaN
        obj_val = float('nan')
    # reset objective coefficients to 0
    model.solver.objective.set_linear_coefficients({
        reaction.forward_variable: 0, reaction.reverse_variable : 0
    })
    return(obj_val)
