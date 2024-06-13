# macaw_dilution.py
'''
Adds a "dilution reaction" for each metabolite in the model (or a given list of
metabolites in the given model) and constrains the flux through that reaction to
be 1/dil_factor * the sum of the absolute values of the fluxes through all
reactions that involve that metabolite weighted by that metabolite's
stoichiometric coefficients in each reaction.

Requires the model to be capable of net production of all metabolites that
participate in at least one reaction with nonzero predicted flux; this prevents
"perfect" recycling of cofactors by requiring at least some net production of
them. This can improve predictions of knockout phenotypes by requiring fluxes
through biosynthetic pathways that might otherwise not be required to have flux.

Mitigates but does not completely solve the problem of unbounded loop fluxes.

This concept was introduced in https://doi.org/10.1186/gb-2010-11-4-r43 and this
implementation of the idea is based on the implementation presented in
https://doi.org/10.1371/journal.pcbi.1003126
'''

import cobra
from macaw.macaw_fva import fva
import re
import pandas as pd
import itertools as it
from optlang.symbolics import Zero

def add_dilution_constraints(
    given_model, mets_to_dilute = None, leak_flux = 1, dil_factor = 1000,
    preprocess = True, fva_results = None, zero_thresh = 10**-8, debug = False,
    debug_rxn = '', threads = 1, verbose = 1
):
    '''
    Creates a "dilution" reaction for each metabolite in mets_to_dilute that
    consumes it and produces nothing whose flux is constrained to be the sum
    of the absolute values of all fluxes through reactions that involve that
    metabolite divided by dil_factor.

    If no mets_to_dilute are given, will use the list of all metabolites in the
    model except for ones that have "tRNA" or "cytochrome" (case-insensitive)
    in their IDs and/or names. See docstring for add_dilution_reactions for an
    explanation.

    Before doing any of that, does FVA for all internal (i.e. not exchange)
    reversible reactions and sets their bounds to their minimum and maximum
    possible fluxes; see docstring for constrain_reversible_rxns for details.
    Note that skipping this preprocessing step can make solutions to the model
    with dilution constraints take about twice as long to compute, so only skip
    it if you have a compelling reason.
    '''
    # make sure dil_factor is > 0, cuz there's no point in making dilution
    # constraints with dilution factors of 0 (they'd divide by zero)
    if dil_factor <= 0:
        if verbose > 0:
            msg = ' - The dil_factor passed to add_dilution_constraints was '
            msg += '<= 0; dilution constraints can only be imposed if '
            msg += 'dil_factor > 0. Returning given_model exactly as it was.'
            print(msg)
        return(given_model)
    model = given_model.copy()
    if preprocess:
        # do FVA on all reversible reactions (see docstring for details)
        model = constrain_reversible_rxns(
            model, fva_results, zero_thresh, threads, verbose
        )
    # if no list of metabolites to add dilution constraints for was given, use
    # the list of all metabolites in the model except for ones with "tRNA" in
    # their names, since those usually lack a biosynthetic pathway (since it's
    # extremely complicated and frequently different for each tRNA) and can thus
    # only participate in a loop between being charged with an amino acid and
    # being returned to their original state by a "translation" reaction, which
    # tends to be required for biomass flux
    if not mets_to_dilute:
        to_skip = ['trna', 'cytochrome']
        mets_to_dilute = [
            m.id for m in model.metabolites
            # be case-insensitive
            if all(s not in m.name.lower() for s in to_skip)
            # check metabolite names and IDs, just to be thorough
            and all(s not in m.id.lower() for s in to_skip)
        ]
        if verbose > 0:
            msg = ' - Since no mets_to_dilute were passed to '
            msg += 'add_dilution_constraints, adding dilution constraints for '
            msg += f'all {len(mets_to_dilute)} metabolites that don\'t have '
            msg += '"tRNA" or "cytochrome" in their names or IDs.'
            print(msg)
    elif not all(isinstance(m, str) for m in mets_to_dilute):
        # if this was a list of Metabolite objects from given_model, replace it
        # with a list of their IDs, since those objects are technically not
        # the same as the metabolite objects in model, and that'll cause drama
        # later
        try:
            mets_to_dilute = [m.id for m in mets_to_dilute]
        except:
            msg = 'Something was wrong with the contents of mets_to_dilute in '
            msg += 'add_dilution_reactions(); make sure it is either a list of '
            msg += 'cobra.Metabolite objects from given_model or a list of '
            msg += 'strings corresponding to the IDs of cobra.Metabolite '
            msg += 'objects that are present in given_model.'
            raise ValueError(msg)
    # now add dilution reactions, leakage reactions, and dilution constraints
    model = add_dilution_reactions(model, mets_to_dilute, verbose)
    model = add_leakage_reactions(model, leak_flux, verbose)
    if debug:
        # if debug is True, add each dilution reaction one at a time,
        # check to see if the given reporter reaction can still
        # sustain flux, and stop if/when a particular dilution
        # constraint blocks flux through the given reporter reaction
        model = debug_dilution(
            model, debug_rxn, mets_to_dilute, dil_factor
        )
    else:
        # otherwise, just make all the dilution constraints in a list comp
        dil_consts = [
            make_dilution_constraint(model, met_id, dil_factor)
            for met_id in mets_to_dilute
        ]
        model.add_cons_vars(dil_consts)
    return(model)

def constrain_reversible_rxns(
    given_model, fva_results = None, zero = 10**-8, threads = 1, verbose = 1
):
    '''
    Do FVA on all internal (i.e. not exchange) reversible reactions and set
    their bounds to their minimum and maximum possible fluxes to mitigate the
    extent to which loop fluxes betweeen their forward and reverse halves
    "cheat" dilution constraints.
    
    We want dilution fluxes to scale with the sum of the absolute values of
    all fluxes that involve each metabolite without imposing any non-linear
    constraints on the model (because solving models with non-linear constraints
    is substantially more computationally expensive, especially for larger
    GSMMs), so we don't want to use an absolute value operator. Cobrapy
    internally represents all reversible reactions as pairs of irreversible ones
    that do exactly the opposite reaction (e.g. A <-> B is represented as A -> B
    and B -> A). Since both of these half-reactions are irreversible, they can
    only have non-negative fluxes, so adding their fluxes will always give us a
    positive number to use for the overall reversible reaction in any dilution
    constraints that include it.

    However, as far as I know, there's no way to prevent both "halves" of a
    single reversible reaction from carrying a non-zero flux at the same time
    without imposing non-linear constraints on the model, so this sum of the
    forward and reverse fluxes for a reversible reaction frequently exceeds the
    net flux through the reversible reaction by a substantial margin. Setting
    the bounds on each reversible reaction to the maximum possible *net* flux in
    each direction significantly limits the extent to which these sum of the
    forward and reverse flux can exceed the net flux through a reversible
    reaction without involving any binary or other non-linear constraints.

    Notably, this significantly restricts the solution space of the model with
    dilution constraints, which, at least sometimes, seems to halve the time it
    takes to compute solutions to the model with dilution constraints, which can
    be quite significant if you're doing, say, FVA on a model with several
    thousand reactions.
    '''
    # work with a copy of the given model to y'know be polite or whatever
    model = given_model.copy()
    rev_rxns = [
        r.id for r in model.reactions
        # leave exchange reactions out of this, since we also ignore them when
        # imposing dilution constraints
        if r.reversibility and not r.boundary
    ]
    # do FVA on all of those reactions unless FVA results were provided
    if fva_results is None:
        if verbose > 0:
            print(' - Doing FVA for all reversible non-exchange reactions.')
        fva_results = fva(
            model, rev_rxns, zero_thresh = zero, threads = threads,
            verbose = verbose - 1
        )
    else:
        # if we got FVA results for more than just the non-exchange reversible
        # reactions, filter down to those
        fva_results = fva_results[fva_results.index.isin(rev_rxns)].copy()
    if verbose > 0:
        msg = ' - Setting bounds on reversible reactions to minimum and maximum'
        msg += ' possible fluxes.'
    # now set the lower bounds to the minima and upper bounds to maxima
    # assume NAs mean the minimum or maximum was infinite
    fva_results['minimum'] = fva_results['minimum'].fillna(-float('Inf'))
    fva_results['maximum'] = fva_results['maximum'].fillna(float('Inf'))
    # count how many reactions fall into each category
    fwd_only = 0
    rev_only = 0
    neither = 0
    for (rxn_id, row) in fva_results.iterrows():
        if (abs(row['minimum']) <= zero) and (abs(row['maximum']) <= zero):
            neither += 1
            # set both bounds on this reaction to 0
            model.reactions.get_by_id(rxn_id).lower_bound = 0
            model.reactions.get_by_id(rxn_id).upper_bound = 0
            if verbose > 1:
                msg = f' - Reaction {rxn_id} can\'t sustain flux in either '
                msg += 'direction; setting both bounds to zero.'
                print(msg)
        elif (row['minimum'] >= -zero) and (row['maximum'] > zero):
            fwd_only += 1
            # set the lower bound to zero
            model.reactions.get_by_id(rxn_id).lower_bound = 0
            if verbose > 1:
                msg = f' - Reaction {rxn_id} could never go backwards; setting '
                msg += 'lower bound to zero.'
                print(msg)
        elif (row['minimum'] < -zero) and (row['maximum'] <= zero):
            rev_only += 1
            # setting the upper bound to zero would mean this reaction could
            # only have negative fluxes, which can cause all sorts of bizarre
            # problems when you least expect it, so switch the products and
            # reactants and then set the lower bound to zero
            rxn = model.reactions.get_by_id(rxn_id)
            new_met_dict = {
                met : -1 * coef for (met, coef) in rxn.metabolites.items()
            }
            # "add" this dict twice to wind up with the opposite of the original
            rxn.add_metabolites(new_met_dict)
            rxn.add_metabolites(new_met_dict)
            # now we can set the lower bound to 0
            rxn.lower_bound = 0
            if verbose > 1:
                msg = f' - Reaction {rxn_id} could only sustain negative '
                msg += 'fluxes; switching products and reactants and setting '
                msg += 'lower bound to zero.'
                print(msg)
        else:
            # if this reaction could sustain flux in either direction, set both
            # bounds to those numbers to prevent the separate forward and
            # reverse halves of this reaction we're going to make from
            # sustaining arbitrarily large loop fluxes between each other
            rxn = model.reactions.get_by_id(rxn_id)
            # round to avoid setting bounds that are only non-zero due to
            # rounding errors (e.g. -1.572935*10^-8)
            rxn.lower_bound = round(row['minimum'], 3)
            rxn.upper_bound = round(row['maximum'], 3)
    if verbose > 0:
        start_rev_rxns = len([
            r for r in given_model.reactions
            if r.reversibility and not r.boundary
        ])
        msg = f' - According to FVA results, out of {start_rev_rxns} nominally '
        msg += 'reversible reactions (excluding exchange reactions), '
        msg += f'{fwd_only} could only sustain forward fluxes, {rev_only} could'
        msg += f' only sustain reverse fluxes, and {neither} couldn\'t sustain '
        msg += 'fluxes in either direction.'
        print(msg)
    return(model)

def add_dilution_reactions(given_model, mets_to_dilute = None, verbose = 1):
    '''
    Given a Cobrapy Model object and, optionally, a list of IDs of Metabolite
    objects in that Model, add a "dilution" reaction for each metabolite in
    mets_to_dilute that irreversibly consumes that metabolite and produces
    nothing to represent "dilution" of that metabolite out of/away from the
    cell, or loss to a side reaction (e.g. oxidative damage).

    If no mets_to_dilute are given, make a dilution reaction for all metabolites
    in the model except for those with "tRNA" or "cytochrome" in their names.
    Those metabolites are frequently present in models but only capable
    of participating in short loops of reactions (e.g. redox loops for
    cytochromes, charging with an amino acid then losing the amino acid in
    a "translation" reaction for tRNAs) and lack biosynthesis pathways. Creating
    dilution constraints for metabolites that can only participate in loops of
    reactions and never be net produced (e.g. by a biosynthetic pathway or
    exchange reaction) prevents flux through all reactions such metabolites
    participate in. Since tRNA and cytochrome metabolites are frequently
    included in GSMMs, rarely have accompanying biosynthetic pathways, and are
    especially likely to participate in particularly important reactions (i.e.
    the electron transport chain leading up to ATP synthase for cytochromes and
    "translation" reactions associated with the biomass reaction for tRNAs),
    it's generally more trouble than it's worth to try and include dilution
    constraints for cytochromes and tRNAs.
    '''
    # start by making sure given_model doesn't already have dilution reactions
    probably_dilution_rxns = [
        r for r in given_model.reactions
        if r.id.endswith('_dilution') and (len(r.metabolites) == 1)
    ]
    if any(probably_dilution_rxns):
        if verbose > 0:
            msg = ' - The given_model passed to add_dilution_reactions appears '
            msg += 'to already have at least one dilution reaction; returning '
            msg += 'the model as-is.'
            print(msg)
        return(given_model)
    # otherwise, make a copy of given_model to add dilution reactions to
    model = given_model.copy()
    # now make an irreversible reaction that consumes each metabolite in
    # mets_to_dilute whose ID is <metabolite ID>_dilution so they're easy to
    # identify later
    dil_rxns = [make_dilution_reaction(model, m) for m in mets_to_dilute]
    # adding a list of many reactions is much faster than separately adding
    # many individual reactions due to the way Cobrapy Models work
    model.add_reactions(dil_rxns)
    return(model)

def make_dilution_reaction(model, met_id):
    '''
    Make a reaction that consumes the given metabolite and produces nothing
    '''
    met_obj = model.metabolites.get_by_id(met_id)
    dil_rxn = cobra.Reaction(f'{met_id}_dilution')
    dil_rxn.name = f'{met_obj.name} Dilution'
    dil_rxn.lower_bound = 0
    dil_rxn.upper_bound = float('Inf')
    dil_rxn.add_metabolites({met_obj : -1.0})
    return(dil_rxn)

def add_leakage_reactions(given_model, bound = 1, verbose = 1):
    '''
    Finds all pairs of Metabolite objects in given_model that represent the
    same real-world compound in different subcellular compartments and creates
    a new reversible reaction that reversibly interconverts them with a flux of
    up to +/- bound (which should be a relatively small non-negative non-zero 
    number). This represents "leakage" of this metabolite across the boundary
    separating the two compartments.

    This is important if the only way for a particular pair of metabolites to
    to move between a particular pair of compartments is via antiport, i.e.
    A[out] + B [in] <-> A[in] + B[out]
    and the only way to produce B[in] is to derive it from A[out], which might
    be the case for an NTP and its corresponding NDP.

    Dilution constraints would require a small amount of A[in] to be "wasted"
    via its dilution reaction, thus ensuring that there is always less B[in]
    than the amount of A[out] that initially entered, so there is never enough
    B[in] to exchange for more A[out], preventing the antiport reaction from 
    ever having fluxes greater than 0.

    Leakage reactions resolve this conundrum this by allowing an arbitrary
    small flux of each metabolite to move between compartments independently
    of any possible antiport scheme to make up for these dilution fluxes. They
    are also hypothetically plausible representations of reality, since no
    biological membrane is always completely impemetrable to any particular
    metabolite, and we're usually trying to model steady states of metabolic
    networks.
    '''
    # first, make sure given_model doesn't already have leakage reactions
    if any(
        ('--' in r.id) and r.id.endswith('_leakage')
        for r in given_model.reactions
    ):
        if verbose > 0:
            msg = ' - The given_model passed to add_leakage_reactions appears '
            msg += 'to already contain leakage reactions; returning the model '
            msg += 'as-is.'
            print(msg)
        return(given_model)
    # then make sure bound isn't zero
    if bound == 0:
        if verbose > 0:
            msg = '- The bound argument passed to add_leakage_reactions was 0, '
            msg += 'so leakage reactions will not be added.'
            print(msg)
        return(given_model)
    # otherwise, as usual, modify a copy of the given model
    model = given_model.copy()
    leakage_rxns = list()
    # in a perfect world, metabolites in different compartments that represent
    # the same real-world metabolite would contain some reference to each other
    # in their IDs, but we do not live in a perfect world, so hope that these
    # metabolites can be linked with their names. start by making a list of all
    # metabolites that have the same name as at least one other metabolite
    all_names = list()
    for m in model.metabolites:
        # if metabolite names appear to end with their compartments, drop the
        # compartment names so we still find matches
        regexp = re.compile(r' ?[\(\[\{]?' + m.compartment + r'[\)\]\}]?$')
        all_names.append(regexp.sub('', m.name))
    name_counts = pd.Series(all_names).value_counts()
    shared_names = name_counts[name_counts > 1].index.unique().tolist()
    # now loop over this list of names and get the corresponding metabolites
    for met_name in shared_names:
        met_objs = [m for m in model.metabolites if met_name == m.name]
        # loop over each possible pair of these metabolites
        for (m1, m2) in it.combinations(met_objs, 2):
            # a Cobrapy Metabolite object has a "reactions" attribute that is a
            # frozenset of all the Cobrapy Reaction objects that involve that
            # Metabolite object
            shared_rxns = m1.reactions.intersection(m2.reactions)
            # see if any reaction involves both of these metabolites
            if shared_rxns:
                # create a reversible reaction that interconverts these
                leak_rxn = cobra.Reaction(
                    f'{m1.id}--{m2.id}_leakage', name = f'{met_name} Leakage'
                )
                leak_rxn.add_metabolites({m1 : -1, m2 : 1})
                # only let it have a flux of up to 1/dil_factor
                leak_rxn.lower_bound = -bound
                leak_rxn.upper_bound = bound
                leakage_rxns.append(leak_rxn)
    if verbose > 0:
        msg = f' - Added a "leakage" reaction for all {len(leakage_rxns)} pairs'
        msg += ' of metabolites in different compartments that seem to '
        msg += 'represent the same real-world metabolite and set their bounds '
        msg += f'to +/- {bound}.'
        print(msg)
    # add all leakage reactions in one go cuz that's much faster than adding
    # them one at a time
    model.add_reactions(leakage_rxns)
    return(model)

def make_dilution_constraint(model, met_id, dil_factor):
    # start out with an optlang/SymPy expression that's just zero, and add in
    # the optlang/SymPy variables associated with the reactions that this
    # metabolite participates in
    expression = Zero
    # passing around Metabolite objects directly has frequently given us weird
    # errors, so we're passing around the metabolite IDs instead
    met_obj = model.metabolites.get_by_id(met_id)
    for r in met_obj.reactions:
        if 'dilution' not in r.id:
            # add both the forward and reverse variable so this is always
            # positive, regardless of which direction the reaction is going in
            expression += r.forward_variable + r.reverse_variable
        else:
            # if this is the metabolite's dilution reaction, subtract its
            # flux times the dilution factor (since it's irreversible, reverse
            # variable should always be 0, but subtracted just in case)
            expression -= dil_factor * (r.forward_variable - r.reverse_variable)
    # set the upper and lower bounds on this constraint to 0 so that the
    # dilution flux (scaled by the dilution factor) must equal the sum of fluxes
    # through all other reactions involving this metabolite
    dilution_constraint = model.problem.Constraint(
        expression, lb = 0, ub = 0, name = f'{met_id}_dilution_constraint'
    )
    return(dilution_constraint)
