# fig_S7_data.py
'''
Make a copy of version 9.0.0 of yeast-GEM where you fix all the errors in lipoic
acid metabolism
'''

import optlang
import cobra

# silence annoying optlang message
optlang.glpk_interface.Configuration()

model = cobra.io.read_sbml_model('GSMMs/yeast-GEMv9.0.0.xml')

# add compartment suffixes to metabolite names
for m in model.metabolites:
    m.name += f' [{m.compartment}]'
# replace reaction names with IDs so Escher map looks less terrible
for r in model.reactions:
    r.name = r.id
# save this as the model for figure S5A
cobra.io.save_json_model(model, 'data/fig_S7a_model.json')

# replace Protein N6-(lipoyl)lysine with dihydrolipoylprotein in the
# biosynthesis reaction that we're keeping
model.reactions.get_by_id('r_4323').add_metabolites({
    # add coefficient of -1 to remove a product
    model.metabolites.get_by_id('s_3947') : -1,
    # then add a coefficient of +1 to add a new product
    model.metabolites.get_by_id('s_0628') : 1
})

# add a new reaction to convert lipoylprotein into lipoamide
new_rxn = cobra.Reaction('r_new')
new_rxn.add_metabolites({
    model.metabolites.get_by_id('s_1098') : -1,
    model.metabolites.get_by_id('s_1097') : 1
})
model.add_reactions([new_rxn])

# remove inaccurate/redundant reactions and metabolites
model.remove_reactions([
    # duplicate biosynthesis reactions
    'r_4324', 'r_4281',
    # duplicate GCS reactions
    'r_0504', 'r_0509',
    # scavenging reactions
    'r_4256', 'r_4258'
])
# Protein N6-(lipoyl)lysine in [c] and [m], lipoyl-AMP, and lipoyl-ACP
to_remove = ['s_3947', 's_3887', 's_3916', 's_3946']
model.remove_metabolites([model.metabolites.get_by_id(m) for m in to_remove])

# save this as the model for figure S5B
cobra.io.save_json_model(model, 'data/fig_S7b_model.json')
