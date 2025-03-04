# utils.py
'''
Miscellaneous minor utility functions
'''

import math
import pandas as pd

def time_str(start, end):
    '''
    Given two time.time objects, return a string with the hours, minutes and
    seconds between the two times
    '''
    (hrs, rem) = divmod(end - start, 3600)
    (mins, secs) = divmod(rem, 60)
    secs = round(secs)
    if hrs == 0:
        if mins == 0:
            msg = f'{secs} seconds'
        else:
            msg = f'{mins} minutes and {secs} seconds'
    else:
        msg = f'{hrs} hours, {mins} minutes, and {secs} seconds'
    return(msg)

def flip_reaction(reaction):
    '''
    Switch the products and reactants and upper and lower bounds for the given
    Cobra.Reaction object
    '''
    orig_ub = reaction.upper_bound
    orig_lb = reaction.lower_bound
    # make a dict with 2 * the current stoichiometric coefficient of each
    # metabolite, then "subtract" this dict from the current reaction to
    # negate the existing coefficients on all metabolites
    new_met_dict = {m : 2 * s for (m, s) in reaction.metabolites.items()}
    reaction.subtract_metabolites(new_met_dict)
    reaction.lower_bound = -1 * orig_ub
    reaction.upper_bound = -1 * orig_lb
    # everything should be modified in-place

def sigfig_round(x, sfs):
    '''
    Round the given number to the given number of significant figures and also
    take care of things that are suspiciously close to zero
    '''
    if abs(x) < 10**-8:
        out = 0.0
    elif abs(x) == float('Inf'):
        # don't try to do actual math with infinities
        out = x
    else:
        out = round(x, sfs-int(math.floor(math.log10(abs(x))))-1)
    return(out)

def add_reaction_equations(
    df, given_model, id_col = 'reaction_id', use_names = False,
    add_suffixes = False
):
    '''
    Given a Pandas DataFrame that contains a column called "reaction_id" and a
    Cobra.Model object that contains reactions with those IDs, add a new column
    to the dataframe called "reaction_equation" like "A + B <-> C + D" where A,
    B, C, and D are either metabolite IDs or metabolite names and can optionally
    have suffixes indicating the compartment they're in added to them (e.g.
    "pyruvate [c]" vs "pyruvate [m]" for cytosolic and mitochondrial pyruvate)
    '''
    # make a copy of the given model so that, if we edit the metabolite names or
    # IDs, the original model object remains unchanged
    model = given_model.copy()
    if add_suffixes:
        for m in model.metabolites:
            if use_names:
                # in case this gets called on the same model multiple times in
                # a row or the metabolite names were already like that
                if not m.name.endswith(f' [{m.compartment}]'):
                    m.name += f' [{m.compartment}]'
            else:
                # cobrapy (reasonably) disapproves of metabolite IDs that
                # contain whitespace characters, so use an underscore instead
                if not m.id.endswith(f'_[{m.compartment}]'):
                    m.id += f'_[{m.compartment}]'
    rxn_string_dict = {
        r.id : r.build_reaction_string(use_names)
        for r in model.reactions
        if r.id in df[id_col].to_list()
    }
    # df.insert raises an exception if you try to insert a column with the
    # same name as an existing column, so catch that and skip this without
    # raising an exception
    if 'reaction_equation' in df.columns:
        msg = 'The dataframe passed to add_reaction_equations already had a '
        msg += 'column called "reaction_equation"; returning it as-is.'
        print(msg)
    else:
        # insert the new column immediately after the reaction ID column
        df.insert(
            loc = df.columns.get_loc(id_col) + 1,
            column = 'reaction_equation', 
            value = df[id_col].map(rxn_string_dict)
        )
    # technically I think the dataframe was modified in-place but whatever
    return(df)

def edit_dead_end_bounds(given_model, results):
    '''
    Given a Cobrapy Model and a Pandas DataFrame containing the results of
    running the dead-end test on that Model:
        - Set both bounds to zero for all reactions found to be dead-ends
        - Set the appropriate bound to zero for reversible reactions found to
          be incapable of carrying steady-state fluxes in only one of their
          directions due to dead-ends. If a reversible reaction could only
          sustain fluxes in the reverse direction, switch its products and
          reactants and then set its lower bound to 0 rather than setting its
          upper bound to 0 and only allow it to run backwards
    '''
    # make sure the dataframe has a column called dead_end_test
    if 'dead_end_test' not in results.columns:
        msg = 'The DataFrame passed to the results argument of macaw.utils.'
        msg += 'edit_dead_end_bounds did not have a column called "dead_end_'
        msg += 'test".'
        raise ValueError(msg)
    model = given_model.copy()
    for (idx, row) in results.iterrows():
        rxn = model.reactions.get_by_id(row['reaction_id'])
        if row['dead_end_test'] != 'ok':
            if row['dead_end_test'] == 'only when going backwards':
                rxn.lower_bound = 0
            elif row['dead_end_test'] == 'only when going forwards':
                rxn.upper_bound = 0
                flip_reaction(rxn)
            else:
                # if it wasn't "ok" and wasn't "only when going forwards/
                # backwards", it must be a list of metabolite IDs and thus a
                # reaction that's incapable of going in either direction
                rxn.lower_bound = 0
                rxn.upper_bound = 0
    return(model)

def simplify_test_results(given_df):
    '''
    Simplify the columns of test results to just say "bad" or "ok"
    '''
    df = given_df.copy()
    if 'dead_end_test' in df.columns:
        # don't count reversible reactions the dead-end test found to only be
        # capable of carrying flux in one direction as "bad"
        df['dead_end_test'] = df['dead_end_test'].apply(lambda x:
            'bad' if (x != 'ok') and not x.startswith('only') else 'ok'
        )
    if 'dilution_test' in df.columns:
        # if we read the test results from a file, there could be NAs that won't
        # have startswith methods, so get those out of the way first
        df['dilution_test'] = df['dilution_test'].apply(
            lambda x: 'bad' if pd.isna(x) else x
        )
        df['dilution_test'] = df['dilution_test'].apply(lambda x:
            'bad' if (x != 'ok') and not x.startswith('always') else 'ok'
        )
    if 'diphosphate_test' in df.columns:
        # leave it alone if it's all NAs
        if not df['diphosphate_test'].isna().all():
            df['diphosphate_test'] = df['diphosphate_test'].apply(
                lambda x: 'bad' if x.startswith('should') else 'ok'
            )
    if 'loop_test' in df.columns:
        df['loop_test'] = df['loop_test'].apply(
            lambda x: 'bad' if x != 'ok' else x
        )
    # simplifying the duplicate test results is more involved cuz they're split
    # across multiple columns
    if any(c.startswith('duplicate_test') for c in df.columns):
        ok = ['ok', float('nan')]
        df['duplicate_test'] = df.apply(
            lambda row: 'bad' if (
                (
                    row['duplicate_test_exact'] not in ok
                ) or (
                    row['duplicate_test_directions'] not in ok
                ) or (
                    row['duplicate_test_coefficients'] not in ok
                ) or (
                    row['duplicate_test_redox'] not in ok
                )
            ) else 'ok',
            axis = 1
        )
        # now we can drop the other duplicate test columns
        df = df.drop(columns = [
            'duplicate_test_exact', 'duplicate_test_directions',
            'duplicate_test_coefficients', 'duplicate_test_redox'
        ])
    return(df)
