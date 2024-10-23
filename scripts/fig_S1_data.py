# fig_S1_data.py
'''
Make toy models to demonstrate how dilution test works
'''

import optlang
import cobra

# silence annoying optlang message
optlang.glpk_interface.Configuration()

# make four metabolites
a = cobra.Metabolite('a', name = 'A', compartment = 'x')
b = cobra.Metabolite('b', name = 'B', compartment = 'x')
c = cobra.Metabolite('c', name = 'C', compartment = 'x')
d = cobra.Metabolite('d', name = 'D', compartment = 'x')
# reaction that produces A from nothing
v1 = cobra.Reaction('v1', name = 'v1', lower_bound = 0, upper_bound = 100)
v1.add_metabolites({a : 1})
# reaction that turns A and C into B and D
v2 = cobra.Reaction('v2', name = 'v2', lower_bound = 0, upper_bound = 100)
v2.add_metabolites({a : -1, b : 1, c : -1, d : 1})
# reaction that consumes B and produces nothing
v3 = cobra.Reaction('v3', name = 'v3', lower_bound = 0, upper_bound = 100)
v3.add_metabolites({b : -1})
# reaction that turns D into C
v4 = cobra.Reaction('v4', name = 'v4', lower_bound = 0, upper_bound = 100)
v4.add_metabolites({c : 1, d : -1})
# add the reactions to a model (also implicitly adds metabolites)
model = cobra.Model('cofactor_loop')
model.add_reactions([v1, v2, v3, v4])
# save this as the model to use for panel A
cobra.io.save_json_model(model, 'data/fig_S1a_model.json')

# for each metabolite, add a dilution reaction and constraint and find the
# maximum possible dilution flux
for (met_id, panel) in [('a', 'b'), ('b', 'c'), ('c', 'd'), ('d', 'e')]:
    # make a fresh copy of the model each time just to be sure everything stays
    # separate
    dil_model = model.copy()
    met_obj = dil_model.metabolites.get_by_id(met_id)
    dil_rxn = cobra.Reaction('v_dil', name = 'Dilution')
    dil_rxn.add_metabolites({met_obj : -1})
    dil_model.add_reactions([dil_rxn])
    dil_exp = optlang.symbolics.Zero
    for r in met_obj.reactions:
        if r.id != 'v_dil':
            dil_exp += r.forward_variable + r.reverse_variable
        else:
            # if this is the metabolite's dilution reaction, subtract its
            # flux times the dilution factor (since it's irreversible, reverse
            # variable should always be 0, but subtracted just in case)
            dil_exp -= 1000 * (r.forward_variable - r.reverse_variable)
    # set the upper and lower bounds on this constraint to 0 so that the
    # dilution flux must equal the sum of fluxes (scaled by the dilution factor)
    # through all other reactions involving this metabolite
    dil_const = model.problem.Constraint(
        dil_exp, lb = 0, ub = 0, name = 'dilution_constraint'
    )
    dil_model.add_cons_vars([dil_const])
    # find a solution that maximizes the dilution flux
    dil_model.objective = dil_rxn
    fluxes = dil_model.optimize().fluxes.reset_index()
    fluxes.columns = ['reaction_id', 'flux']
    # save the model and fluxes to use for the appropriate panel
    cobra.io.save_json_model(dil_model, f'data/fig_S1{panel}_model.json')
    fluxes.to_csv(f'data/fig_S1{panel}_data.csv', index = False)
