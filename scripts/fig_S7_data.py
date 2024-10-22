# fig_S7_data.py
'''
Make a toy model with two metabolites that can exist in two different
compartments and only move between those compartments via an antiport reaction
to show how dilution constraints can block fluxes through some antiport
reactions but adding "leakage" reactions prevents this
'''

import optlang
import cobra
from macaw.dilution import add_dilution_constraints

# silence annoying optlang message that prints the first time you initialize a
# Cobrapy Model in a given Python session
optlang.glpk_interface.Configuration()

# make metabolites
a_o = cobra.Metabolite('a_out', name = 'A (out)', compartment = 'out')
b_o = cobra.Metabolite('b_out', name = 'B (out)', compartment = 'out')
a_i = cobra.Metabolite('a_in', name = 'A (in)', compartment = 'in')
b_i = cobra.Metabolite('b_in', name = 'B (in)', compartment = 'in')
# make reactions
v1 = cobra.Reaction(
    'v1', name = 'v1 (0, 100)', lower_bound = 0, upper_bound = 100
)
v1.add_metabolites({a_o : 1})
v2 = cobra.Reaction(
    'v2', name = 'v2 (0, inf)', lower_bound = 0, upper_bound = float('Inf')
)
v2.add_metabolites({a_o : -1, a_i : 1, b_i : -1, b_o : 1})
v3 = cobra.Reaction(
    'v3', name = 'v3 (0, inf)', lower_bound = 0, upper_bound = float('Inf')
)
v3.add_metabolites({a_i : -1, b_i : 1})
v4 = cobra.Reaction(
    'v4', name = 'v4 (0, inf)', lower_bound = 0, upper_bound = float('Inf')
)
v4.add_metabolites({b_o : -1})
# add reactions to model (will also add metabolites as needed)
normal_model = cobra.Model('leakage_demo')
normal_model.add_reactions([v1, v2, v3, v4])
# set objective to maximizing flux through v4
normal_model.objective = 'v4'

# make one version with dilution reactions but not leakage reactions and one
# with both (in addition to the existing one with neither)
dil_only_model = add_dilution_constraints(
    normal_model, mets_to_dilute = ['a_out', 'a_in', 'b_in', 'b_out'],
    preprocess = False, leak_flux = 0, dil_factor = 1000
)
dil_leak_model = add_dilution_constraints(
    normal_model, mets_to_dilute = ['a_out', 'a_in', 'b_in', 'b_out'],
    preprocess = False, leak_flux = 1, dil_factor = 1000
)

# solve each model and save the fluxes to a file
panels = ['a', 'b', 'c']
i = 0
for model in [normal_model, dil_only_model, dil_leak_model]:
    fluxes = model.optimize().fluxes.reset_index()
    fluxes.columns = ['reaction_id', 'flux']
    fluxes.to_csv(f'data/fig_S7{panels[i]}_data.csv', index = False)
    i += 1

# only need to save the model with both leakage and dilution reactions to
# make Escher maps cuz it contains all reactions present in the other two
cobra.io.save_json_model(dil_leak_model, 'data/fig_S7_model.json')
