# macaw_structural.py
'''
Functions for finding (and optionally, fixing) sets of reactions with problems
that are detectable just by examining their definitions, bounds, and/or
immediate context in the network (i.e. as opposed to detecting problems by
predicting fluxes through reactions)
'''

import pandas as pd
from macaw_utils import add_reaction_equations
import itertools as it
from macaw_utils import flip_reaction

def dead_end_test(
    given_model, use_names = False, add_suffixes = False, verbose = 1
):
    '''
    Find all metabolites that have 0 or 1 reactions associated with them, as
    well as any reactions involving those metabolites. Do so recursively to
    identify whole dead-end pathways. Also identify all reversible reactions
    that involve metabolites that can only be produced or only be consumed by
    all other reactions they participate in, make those reversible reactions
    irreversible in the appropriate direction, and see if that creates any new
    dead-ends.
    '''
    if verbose > 0:
        print('Starting dead-end test...')
    # work with a copy of the given model so the original object is left as-is
    model = given_model.copy()
    dead_end_rxns = list()
    dead_end_mets = list()
    fwd_only = list()
    rev_only = list()
    # loop over metabolites and assess their dead-ended-ness and the dead-ended
    # ness of any reactions that involve them
    for met in model.metabolites:
        _dead_end_test_inner(
            met, dead_end_mets, dead_end_rxns, fwd_only, rev_only
        )
    # make a dict where the keys are reaction IDs and the values are either
    # "ok" if the reaction wasn't involved in any dead ends or a list of the
    # metabolites that participate in that reaction that only participate in
    # other dead-end reactions so you can trace out pathways of connected
    # dead-end reactions to y'know maybe find the root cause. also indicate if
    # reversible reactions had just one direction blocked by a dead-end
    test_results = dict()
    edge_list = list()
    for r in model.reactions:
        if r in dead_end_rxns:
            mets = [m.id for m in r.metabolites if m in dead_end_mets]
            test_results[r.id] = ';'.join(mets)
            edge_list.extend([(m, r.id) for m in mets])
        elif r.id in fwd_only:
            test_results[r.id] = 'only when going backwards'
        elif r.id in rev_only:
            test_results[r.id] = 'only when going forwards'
        else:
            test_results[r.id] = 'ok'
    out_df = pd.DataFrame({'reaction_id' : [r.id for r in model.reactions]})
    out_df['dead_end_test'] = out_df['reaction_id'].map(test_results)
    # print a summary of the findings unless verbose = 0
    if verbose > 0:
        irrevs = out_df['dead_end_test'].str.contains('only').sum()
        msg = f' - Found {len(dead_end_mets)} dead-end metabolites.\n - Found '
        msg += f'{len(dead_end_rxns)} reactions incapable of sustaining steady-'
        msg += 'state fluxes in either direction due to these dead-ends.\n - '
        msg += f'Found {irrevs} reversible reactions that can only carry '
        msg += 'steady-state fluxes in a single direction due to dead-ends.'
        print(msg)
    # add a column for the reaction equations, e.g.
    # glucose + ATP + H2O -> glucose-6-phosphate + ADP
    out_df = add_reaction_equations(
        out_df, given_model, use_names = use_names, add_suffixes = add_suffixes
    )
    return((out_df, edge_list))

def _dead_end_test_inner(
    met_to_check, dead_end_mets, dead_end_rxns, fwd_only, rev_only
):
    '''
    Check if a given metabolite is only associated with a single reaction, then
    check all the other metabolites in that reaction but avoid checking the
    same metabolite twice so it doesn't recurse infinitely
    '''
    # first make sure this isn't already on our list of dead-end metabolites
    if met_to_check not in dead_end_mets:
        # if any reactions this metabolite participates in are already on our
        # list of dead-ends, don't count them when deciding if this metabolite
        # is also a dead-end
        rxns_to_check = [
            r for r in met_to_check.reactions if r not in dead_end_rxns
        ]
        # now see if the metabolite still has at least 2 non-dead-end reactions
        if len(rxns_to_check) < 2:
            # this metabolite is only associated with 0 or 1 reactions that we
            # haven't already marked as dead-ends, so it's a dead-end metabolite
            dead_end_mets.append(met_to_check)
            # if there's only one reaction associated with this metabolite that
            # isn't already on our list of dead-end reactions, add it
            if len(rxns_to_check) == 1:
                dead_end_rxns.append(rxns_to_check[0])
                # check the other metabolites in this reaction
                for m in rxns_to_check[0].metabolites:
                    if m not in dead_end_mets:
                        _dead_end_test_inner(
                            m, dead_end_mets, dead_end_rxns, fwd_only, rev_only
                        )
        else:
            # this metabolite participates in at least 2 reactions that we
            # haven't already flagged as being dead-ends, but see if it
            # participates in 0 or 1 reversible reactions and all the
            # irreversible reactions produce it or they all consume it
            rev_rxns = [r for r in rxns_to_check if r.reversibility]
            irrev_rxns = [r for r in rxns_to_check if r not in rev_rxns]
            all_prod = all(met_to_check in r.products for r in irrev_rxns)
            all_cons = all(met_to_check in r.reactants for r in irrev_rxns)
            if (len(rev_rxns) == 0) and (all_prod or all_cons):
                # this metabolite can only ever be produced or only ever be
                # consumed, so it and all reactions it participates in are
                # dead-ends
                dead_end_mets.append(met_to_check)
                dead_end_rxns.extend(rxns_to_check)
                # see if the other metabolites that participate in these
                # reactions are also dead-ends now that we know that all of
                # these reactions are
                mets_to_check = set([
                    m for r in rxns_to_check
                    for m in r.metabolites
                    if m not in dead_end_mets
                ])
                for met in mets_to_check:
                    _dead_end_test_inner(
                        met, dead_end_mets, dead_end_rxns, fwd_only, rev_only
                    )
            elif (len(rev_rxns) == 1) and all_prod:
                # this reaction participates in exactly one reversible reaction
                # but all other reactions can only produce it, so make that
                # reversible reaction irreversible in the direction that
                # consumes this metabolite
                if met_to_check in rev_rxns[0].reactants:
                    # now it'll be an irreversible reaction if we see it again
                    rev_rxns[0].lower_bound = 0
                    fwd_only.append(rev_rxns[0].id)
                elif met_to_check in rev_rxns[0].products:
                    # just setting the upper bound to 0 would mean the products
                    # could only be consumed and the reactants could only be
                    # produced, which is confusing, so instead flip the reaction
                    # around (i.e. switch products and reactants and bounds),
                    # then set the lower bound to zero
                    flip_reaction(rev_rxns[0])
                    rev_rxns[0].lower_bound = 0
                    rev_only.append(rev_rxns[0].id)
                # then make sure the other metabolites that participate in the
                # formerly reversible reaction aren't dead-ends now
                for m in rev_rxns[0].metabolites:
                    if (m not in dead_end_mets) and (m != met_to_check):
                        _dead_end_test_inner(
                            m, dead_end_mets, dead_end_rxns, fwd_only, rev_only
                        )
            elif (len(rev_rxns) == 1) and all_cons:
                # make the one reversible reaction irreversible in the direction
                # that produces this metabolite
                if met_to_check in rev_rxns[0].products:
                    rev_rxns[0].lower_bound = 0
                    fwd_only.append(rev_rxns[0].id)
                elif met_to_check in rev_rxns[0].reactants:
                    flip_reaction(rev_rxns[0])
                    rev_rxns[0].lower_bound = 0
                    rev_only.append(rev_rxns[0].id)
                for m in rev_rxns[0].metabolites:
                    if (m not in dead_end_mets) and (m != met_to_check):
                        _dead_end_test_inner(
                            m, dead_end_mets, dead_end_rxns, fwd_only, rev_only
                        )

def duplicate_test(
    model, redox_pairs = list(), proton_ids = list(), use_names = False,
    add_suffixes = False, verbose = 1
):
    '''
    Finds sets of reactions that may be duplicates of each other because they:
    - Involve exactly the same metabolites with exactly the same stoichiometric
      coefficients (but potentially different associated genes)
    - Involve exactly the same metabolites, but go in different directions and/
      or some are reversible and some are not
    - Involve exactly the same metabolites, but with different stoichiometric
      coefficients
    - Represent the oxidation and/or reduction of the same metabolite, but use
      different electron acceptors/donors from the given list of pairs of
      oxidized and reduced forms of various electron carriers (e.g. NAD(H),
      NADP(H), FAD(H2), ubiquinone/ubiquinol, cytochromes).
    There are sometimes cases where having sets of reactions that fall into one
    of the above categories is a legitimate representation of real biochemistry
    (e.g. different transport proteins for import vs export of a particular
    metabolite across a membrane, enzymes that can use NAD(H) and NADP(H)
    interchangeably to catalyze a particular redox reaction), but in my
    experience, there are quite a few cases where at least one reaction in a
    set of reactions flagged by this function should not exist.
    '''
    if verbose > 0:
        print('Starting duplicate test...')
    # make sure the redox_pairs argument is a list or tuple of at least two
    # lists or tuples that each contain exactly two strings
    if not (isinstance(redox_pairs, list) or isinstance(redox_pairs, tuple)):
        msg = 'The redox_pairs argument to duplicate_test must be either a '
        msg += f'list or a tuple, not a {type(redox_pairs).__name__}.'
        raise TypeError(msg)
    elif not all(
        isinstance(p, list) or isinstance(p, tuple) for p in redox_pairs
    ):
        bad_types = 's or '.join(set([
            type(p).__name__ for p in redox_pairs
            if not isinstance(p, list) or not isinstance(p, tuple)
        ]))
        msg = 'The elements of the redox_pairs argument to duplicate_test must '
        msg += f'all be lists or tuples, not {bad_types}s.'
        raise TypeError(msg)
    elif not all(all(isinstance(m, str) for m in p) for p in redox_pairs):
        bad_types = set()
        for pair in redox_pairs:
            for thing in pair:
                if not isinstance(thing, str):
                    bad_types.add(type(thing).__name__)
        bad_types = 's or '.join(bad_types)
        msg = 'The elements of the redox_pairs argument to duplicate_test must '
        msg += 'be lists or tuples of *strings* of metabolite IDs in the given '
        msg += 'Cobra.Model object, not {bad_types}s.'
        raise TypeError(msg)
    elif ((len(redox_pairs) < 2) or (proton_ids == list())) and (verbose > 0):
        msg = 'Can\'t look for pairs of redox reactions that use different '
        msg += 'redox carriers (e.g. NAD(P), FAD, ubiquinone) but are otherwise'
        msg += ' identical unless the redox_pairs argument contains at least '
        msg += 'two pairs of metabolite IDs corresponding to the oxidized and '
        msg += 'reduced forms of a redox carrier (e.g. [(NAD+, NADH), (FAD, '
        msg += 'FADH2)]) and the proton_ids argument contains the ID of at '
        msg += 'least one metabolite ID that represents a proton. Will still '
        msg += 'look for other kinds of potentially duplicate reactions.'
        print(msg)
        redox_mets = list()
    else:
        # make a list of all the redox carrier and proton metabolite objects
        redox_mets = list()
        for pair in redox_pairs:
            for met in pair:
                try:
                    redox_mets.append(model.metabolites.get_by_id(met))
                except KeyError:
                    msg = f'The metabolite ID {met} in the pair {pair} passed '
                    msg += 'to the redox_pairs argument of duplicate_test did '
                    msg += 'not match the ID of any metabolite in the given '
                    msg += 'GSMM. Verify that there are no typos in the ID and '
                    msg += 'that you provided the correct GSMM and try again.'
                    raise ValueError(msg)
        for met in proton_ids:
            try:
                redox_mets.append(model.metabolites.get_by_id(met))
            except KeyError:
                msg = f'The metabolite ID {met} in the list of metabolite IDs '
                msg += 'passed to the proton_ids argument to duplicate_test did'
                msg += 'not match the ID of any metabolite in the given GSMM. '
                msg += 'Verify that there are no typos in teh ID and that you '
                msg += 'provided the correct GSMM and try again.'
                raise ValueError(msg)
        # also replace the metabolite IDs in redox_pairs with the corresponding
        # metabolite objects
        redox_pairs = [
            (model.metabolites.get_by_id(m1), model.metabolites.get_by_id(m2))
            for (m1, m2) in redox_pairs
        ]
    # make dicts with one key per reaction and lists of the IDs of other
    # reactions that are (potentially) the specified kind of duplicate
    exact_dupes = {r.id : list() for r in model.reactions}
    direction_dupes = {r.id : list() for r in model.reactions}
    stoich_dupes = {r.id : list() for r in model.reactions}
    redox_dupes = {r.id: list() for r in model.reactions}
    # loop over all pairs of reactions
    edge_list = list()
    for (r1, r2) in it.combinations(model.reactions, 2):
        # check if they involve exactly the same metabolites
        if set(r1.metabolites.keys()) == set(r2.metabolites.keys()):
            # see if they're exactly they same or if they have different
            # coefficients, directions, and/or reversibilities
            same_stoich = all(
                r1.metabolites[met] == r2.metabolites[met]
                for met in r1.metabolites.keys()
            )
            opp_dir = all(m in r2.products for m in r1.reactants)
            same_rev = r1.reversibility == r2.reversibility
            if same_stoich and not opp_dir and same_rev:
                exact_dupes[r1.id].append(r2.id)
                exact_dupes[r2.id].append(r1.id)
                edge_list.append((r1.id, r2.id))
            if not same_stoich:
                stoich_dupes[r1.id].append(r2.id)
                stoich_dupes[r2.id].append(r1.id)
                edge_list.append((r1.id, r2.id))
            if opp_dir or not same_rev:
                direction_dupes[r1.id].append(r2.id)
                direction_dupes[r2.id].append(r1.id)
                edge_list.append((r1.id, r2.id))
        else:
            # will be an empty list if no pairs of metabolite IDs were given
            if redox_mets:
                # see if both reactions involve one of the given pairs
                both_redox = any([
                    (m1 in r1.metabolites) and (m2 in r1.metabolites)
                    for (m1, m2) in redox_pairs
                ]) and any([
                    (m1 in r2.metabolites) and (m2 in r2.metabolites)
                    for (m1, m2) in redox_pairs
                ])
                # see if these two reactions involve exactly the same
                # metabolites if you ignore the redox carriers and protons
                mets_1 = {m for m in r1.metabolites if m not in redox_mets}
                mets_2 = {m for m in r2.metabolites if m not in redox_mets}
                # if both sets are empty, we're looking at reactions where two
                # of the redox pairs react with each other or reactions that
                # transport them between compartments, which we don't want to
                # report as potential duplicates, cuz those are 100% realistic
                if both_redox and (mets_1 == mets_2) and mets_1 and mets_2:
                    redox_dupes[r1.id].append(r2.id)
                    redox_dupes[r2.id].append(r1.id)
                    edge_list.append((r1.id, r2.id))
    # turn those dicts into a dataframe that has the list of IDs of potential
    # duplicate reactions separated by semicolons in appropriately-named columns
    # for each category of duplicate reactions (or "ok" for reactions that we
    # found no duplicates for)
    out_df = pd.DataFrame({'reaction_id' : [r.id for r in model.reactions]})
    out_df['duplicate_test_exact'] = out_df['reaction_id'].map(
        lambda r: ';'.join(exact_dupes[r]) if exact_dupes[r] else 'ok'
    )
    out_df['duplicate_test_directions'] = out_df['reaction_id'].map(
        lambda r: ';'.join(direction_dupes[r]) if direction_dupes[r] else 'ok'
    )
    out_df['duplicate_test_coefficients'] = out_df['reaction_id'].map(
        lambda r: ';'.join(stoich_dupes[r]) if stoich_dupes[r] else 'ok'
    )
    if redox_mets:
        out_df['duplicate_test_redox'] = out_df['reaction_id'].map(
            lambda r: ';'.join(redox_dupes[r]) if redox_dupes[r] else 'ok'
        )
    else:
        out_df['duplicate_test_redox'] = 'N/A'
    # some reactions may be more than one kind of duplicate, so count the total
    # number of reactions that were at least one kind of duplicate
    any_dupe = 0
    for r in model.reactions:
        if (
            exact_dupes[r.id] or direction_dupes[r.id] or
            stoich_dupes[r.id] or redox_dupes[r.id]
        ):
            any_dupe += 1
    if verbose > 0:
        exacts = (out_df['duplicate_test_exact'] != 'ok').sum()
        direcs = (out_df['duplicate_test_directions'] != 'ok').sum()
        stoichs = (out_df['duplicate_test_coefficients'] != 'ok').sum()
        redoxes = (~out_df['duplicate_test_redox'].isin(['ok', 'N/A'])).sum()
        msg = f' - Found {any_dupe} reactions that were some type of duplicate:'
        msg += f'\n   - {exacts} were completely identical to at least one '
        msg += f'other reaction.\n   - {direcs} involve the same metabolites '
        msg += 'but go in the opposite direction or have the opposite '
        msg += f'reversibility as at least one other reaction.\n   - {stoichs} '
        msg += f'involve the same metabolites but with different coefficients '
        msg += f'as at least one other reaction.\n   - {redoxes} redox '
        msg += 'reactions oxidize and/or reduce the same metabolites as another'
        msg += ' redox reaction but use a different pair of the given redox '
        msg += 'carriers.'
        print(msg)
    # add a column for the reaction equations, e.g.
    # glucose + ATP + H2O -> glucose-6-phosphate + ADP
    out_df = add_reaction_equations(
        out_df, model, use_names = use_names, add_suffixes = add_suffixes
    )
    return((out_df, edge_list))

def diphosphate_test(
    given_model, ppi_ids, pi_ids, use_names = False, add_suffixes = False,
    verbose = 1
):
    '''
    Given lists of the metabolite IDs that represent diphosphate (ppi_ids) and
    inorganic (mono)phosphate ions (pi_ids) in given_model, identifies all
    reversible reactions that involve one of the diphosphate metabolites (to
    exclude reactions that just transport diphosphate between compartments).
    
    Returns a dataframe with one row per reaction and a "diphosphate_test"
    column that says:
    - "should be irreversible" if diphosphate is a product
    - "should be flipped and made irreversible" if diphosphate is a reactant
    - "ok" if the reaction is already irreversible or doesn't involve
      diphosphate

    Nominally, hydrolysis of a molecule with a diphosphate group is almost
    perfectly reversible at physiological pHs, but real cells express lots of
    highly active diphosphatases specifically to drive this reaction in the
    direction of hydrolysis. This is crucial in real life because DNA and RNA
    polymerases catalyze the reaction
        DNA/RNA + (d)NTP + H2O <-> DNA/RNA-NMP + PPi
    and crucial for steady-state modeling, since, unless you're explicitly
    including thermodynamic constraints, things like FBA will abuse the
    reversibility of these reactions to generate ATP from dozens of utterly
    implausible sources.
    '''
    # first see if we got empty lists for the metabolite IDs
    if (not ppi_ids or not pi_ids) and (verbose > 0):
        msg = 'Can\'t look for reversible reactions that involve diphosphate '
        msg += 'without a list of the IDs of the metabolite objects in the '
        msg += 'given GSMM that represent diphosphate ions and a list of the '
        msg += 'IDs of the metabolites that represent phosphate ions.'
        print(msg)
        # create a dataframe that has reaction IDs in one column and NAs in the
        # other (i.e. the same format as the one we'd return if we actually did
        # the test)
        out_df = pd.DataFrame({
            'reaction_id' : [r.id for r in given_model.reactions]
        })
        out_df['diphosphate_test'] = 'N/A'
        return((out_df, given_model))
    if verbose > 0:
        print('Starting diphosphate test...')
    # copy the given model so the original object remains as it was
    model = given_model.copy()
    # make a dict with with a list of reaction IDs and a list of test results
    # for each reaction (to be made into a Pandas DataFrame later)
    result_dict = {'reaction_id' : [], 'diphosphate_test' : []}
    for r in model.reactions:
        # reactions that involve both diphosphate and phosphate are probably
        # diphosphatase or antiport reactions, so those are fine
        pi_in_rxn = any(m.id in pi_ids for m in r.metabolites)
        # also skip exchange and irreversible reactions
        if (r not in model.boundary) and not pi_in_rxn and r.reversibility:
            # figure out if diphosphate is a product or reactant (or both)
            ppi_prod = any(m.id in ppi_ids for m in r.products)
            ppi_reac = any(m.id in ppi_ids for m in r.reactants)
            if ppi_prod and not ppi_reac:
                test_result = 'should be irreversible'
            elif ppi_reac and not ppi_prod:
                test_result = 'should be flipped and made irreversible'
            else:
                # if there was a diphosphate metabolite in both the products and
                # the reactants for this reaction, assume this is a diphosphate
                # transport reaction (i.e. not problematic)
                test_result = 'ok'
        else:
            test_result = 'ok'
        result_dict['reaction_id'].append(r.id)
        result_dict['diphosphate_test'].append(test_result)
    # turn the dict of results into a dataframe
    out_df = pd.DataFrame(result_dict)
    # print a summary of the findings unless verbose = 0
    if verbose > 0:
        i = out_df['diphosphate_test'].str.contains('irreversible').sum()
        msg = f' - Found {i} suspiciously reversible reactions that involve '
        msg += 'diphosphate.'
        print(msg)
    # add a column for the reaction equations, e.g.
    # glucose + ATP + H2O -> glucose-6-phosphate + ADP
    out_df = add_reaction_equations(
        out_df, model, use_names = use_names, add_suffixes = add_suffixes
    )
    return(out_df)
