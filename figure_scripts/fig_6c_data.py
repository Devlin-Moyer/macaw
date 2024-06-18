# fig_6c_data.py
'''
1. Predict impacts of knocking out LIPT1 and LIAS on maximum possible fluxes
through Pyruvate Dehydrogenase (PDH) and Glycine Cleavage System (GCS)
reactions in version 1.15 of Human-GEM.
2. Make all the changes we proposed to lipoic acid metabolism that weren't
incorporated into version 1.18 of Human-GEM.
3. Predict impacts of knocking out LIPT1, LIAS, GLRX5, and IBA57 on maximum
possible fluxes through PDH and GCS reactions in this version of Human-GEM.
'''

import sys
import optlang
import cobra
from macaw.fva import fva
from macaw.dilution import add_dilution_constraints
import pandas as pd

def fix_lipoate_biosynth(old_model):
    new_model = old_model.copy()
    # create all new reaction objects in one go
    rxns = [cobra.Reaction(f'new_rxn_{i}') for i in range(18)]
    # get all the existing metabolite objects and create all the new ones we're
    # gonna need to make and edit all these reactions
    oct_gcsh = new_model.metabolites.get_by_id('MAM00210c')
    lipoyl_gcsh = new_model.metabolites.get_by_id('MAM00209c')
    dhlipoyl_gcsh = new_model.metabolites.get_by_id('MAM00208c')
    cyto_acp = new_model.metabolites.get_by_id('MAM00184m')
    mito_acp = new_model.metabolites.get_by_id('MAM02484m')
    ac_coa = new_model.metabolites.get_by_id('MAM01261m')
    coa = new_model.metabolites.get_by_id('MAM01597m')
    ac_acp = cobra.Metabolite(
        'MAM01258m', name = 'acetyl-[ACP]', compartment = 'm'
    )
    mal_acp = new_model.metabolites.get_by_id('MAM02442m')
    proton_m = new_model.metabolites.get_by_id('MAM02039m')
    co2 = new_model.metabolites.get_by_id('MAM01596m')
    nadph = new_model.metabolites.get_by_id('MAM02555m')
    nadp = new_model.metabolites.get_by_id('MAM02554m')
    water = new_model.metabolites.get_by_id('MAM02040m')
    acac_acp = cobra.Metabolite(
        'MAM01254m', name = 'acetoacetyl-[ACP]', compartment = 'm'
    )
    ho_but = cobra.Metabolite(
        'MAM00158m', name = '3-hydroxybutanoyl-[ACP]', compartment = 'm'
    )
    bute_acp = cobra.Metabolite(
        'MAM01409m', name = '(2E)-butenoyl-[ACP]', compartment = 'm'
    )
    buty_acp = cobra.Metabolite(
        'MAM01411m', name = 'butyryl-[ACP]', compartment = 'm'
    )
    ox_hex = cobra.Metabolite(
        'MAM00881m', name = '3-oxohexanoyl-[ACP]', compartment = 'm'
    )
    ho_hex = cobra.Metabolite(
        'MAM01635m', name = '3-hydroxyhexanoyl-[ACP]', compartment = 'm'
    )
    hexe_acp = cobra.Metabolite(
        'MAM00052m', name = '(2E)-hexenoyl-[ACP]', compartment = 'm'
    )
    hexa_acp = cobra.Metabolite(
        'MAM02121m', name = 'hexanoyl-[ACP]', compartment = 'm'
    )
    ox_oct = cobra.Metabolite(
        'MAM00891m', name = '3-oxooctanoyl-[ACP]', compartment = 'm'
    )
    ho_oct = cobra.Metabolite(
        'MAM00161m', name = '3-hydroxyoctanoyl-[ACP]', compartment = 'm'
    )
    octe_acp = cobra.Metabolite(
        'MAM00058m', name = '(2E)-octenoyl-[ACP]', compartment = 'm'
    )
    octa_acp_m = cobra.Metabolite(
        'MAM02643m', name = 'octanoyl-[ACP]', compartment = 'm'
    )
    cys = new_model.metabolites.get_by_id('MAM01628m')
    fe_2 = new_model.metabolites.get_by_id('MAM01821m')
    ala = new_model.metabolites.get_by_id('MAM01307m')
    two_isc = cobra.Metabolite('2Fe2S', compartment = 'm')
    atp = new_model.metabolites.get_by_id('MAM01371m')
    four_isc = cobra.Metabolite('4Fe4S', compartment = 'm')
    adp = new_model.metabolites.get_by_id('MAM01285m')
    p_i = new_model.metabolites.get_by_id('MAM02751m')
    gsh = new_model.metabolites.get_by_id('MAM02026m')
    fad = new_model.metabolites.get_by_id('MAM01802m')
    hs = cobra.Metabolite('MAM02042m', name = 'HS-', compartment = 'm')
    gssh = cobra.Metabolite('GSSH', compartment = 'm')
    fadh2 = new_model.metabolites.get_by_id('MAM01803m')
    o2 = new_model.metabolites.get_by_id('MAM02630m')
    sulfite = new_model.metabolites.get_by_id('MAM02949m')
    deoxy_a_m = cobra.Metabolite(
        'MAM01098m', name = '5-deoxyadenosine', compartment = 'm'
    )
    deoxy_a_c = new_model.metabolites.get_by_id('MAM01098c')
    sam = new_model.metabolites.get_by_id('MAM02877m')
    met = new_model.metabolites.get_by_id('MAM02471m')
    apo_a1 = new_model.metabolites.get_by_id('MAM01350c')
    octa_acp_c = new_model.metabolites.get_by_id('MAM02643c')
    nad_c = new_model.metabolites.get_by_id('MAM02552c')
    nad_m = new_model.metabolites.get_by_id('MAM02552m')
    nadh_c = new_model.metabolites.get_by_id('MAM02553c')
    nadh_m = new_model.metabolites.get_by_id('MAM02553m')
    proton_c = new_model.metabolites.get_by_id('MAM02039m')
    # move octanoyl-GCSH, lipoyl-GCSH, and dihydrolipoyl-GCSH from cyto to mito
    oct_gcsh.id = 'MAM00210m'
    oct_gcsh.compartment = 'm'
    lipoyl_gcsh.id = 'MAM00209m'
    lipoyl_gcsh.compartment = 'm'
    dhlipoyl_gcsh.id = 'MAM00208m'
    dhlipoyl_gcsh.compartment = 'm'
    # add in the currently absent mitochondrial fatty acid synthesis pathway
    # the first reaction in the pathway exists as a dead-end, so edit it, but
    # the rest are gonna have to be entirely new
    new_model.reactions.get_by_id('MAR01517').add_metabolites({
        # want to replace [ACP] with mito_acp, but it's a reactant, so "add"
        # one [ACP] and "subtract" one mito_acp
        cyto_acp : 1, mito_acp : -1
    })
    rxns[0].add_metabolites({ac_coa : -1, mito_acp : -1, ac_acp : 1, coa : 1})
    rxns[1].add_metabolites({
        mal_acp : -1, ac_acp : -1, proton_m : -1,
        acac_acp : 1, co2 : 1, mito_acp : 1
    })
    rxns[2].add_metabolites({
        acac_acp : -1, nadph : -1, proton_m : -1, ho_but : 1, nadp : 1
    })
    rxns[3].add_metabolites({ho_but : -1, bute_acp : 1, water : 1})
    rxns[4].add_metabolites({
        bute_acp : -1, nadph : -1, proton_m : -1, buty_acp : 1, nadp : 1
    })
    rxns[5].add_metabolites({
        buty_acp : -1, mal_acp : -1, proton_m : -1,
        ox_hex : 1, co2 : 1, mito_acp : 1
    })
    rxns[6].add_metabolites({
        ox_hex : -1, nadph : -1, proton_m : -1, ho_hex : 1, nadp : 1
    })
    rxns[7].add_metabolites({ho_hex : -1, hexe_acp : 1, water : 1})
    rxns[8].add_metabolites({
        hexe_acp : -1, nadph : -1, proton_m : -1, hexa_acp : 1, nadp : 1
    })
    rxns[9].add_metabolites({
        hexa_acp : -1, mal_acp : -1, proton_m : -1,
        ox_oct : 1, co2 : 1, mito_acp : 1
    })
    rxns[10].add_metabolites({
        ox_oct : -1, nadph : -1, proton_m : -1, ho_oct : 1, nadp : 1
    })
    rxns[11].add_metabolites({ho_oct : -1, octe_acp : 1, water : 1})
    rxns[12].add_metabolites({
        octe_acp : -1, nadph : -1, proton_m : -1, octa_acp_m : 1, nadp : 1
    })
    # add in missing iron-sulfur cluster biosynthesis reactions & metabolites
    rxns[13].add_metabolites({
        cys : -2, fe_2 : -2, nadph : -1,
        ala : 2, two_isc : 1, nadp : 1, proton_m : 1
    })
    two_isc_gpr = 'ENSG00000004779 and ENSG00000136003 and ENSG00000161513 and '
    two_isc_gpr += 'ENSG00000165060 and ENSG00000214113 and ENSG00000244005 '
    two_isc_gpr += 'and ENSG00000267673'
    rxns[13].gene_reaction_rule = two_isc_gpr
    rxns[14].add_metabolites({
        atp : -2, water : -2, nadph : -1, two_isc : -2,
        adp : 2, proton_m : 3, nadp : 1, four_isc : 1, p_i : 2
    })
    four_isc_gpr = 'ENSG00000100209 and ENSG00000109519 and ENSG00000113013 '
    four_isc_gpr += 'and ENSG00000135070 and ENSG00000161513 and '
    four_isc_gpr += 'ENSG00000165898 and ENSG00000181873 and ENSG00000182512 '
    four_isc_gpr += 'and ENSG00000267673'
    rxns[14].gene_reaction_rule = four_isc_gpr
    # add in missing sulfur metabolism reactions & metabolites
    rxns[15].add_metabolites({gsh : -1, fad : -1, hs : -1, gssh : 1, fadh2 : 1})
    rxns[15].gene_reaction_rule = 'ENSG00000137767'
    rxns[16].add_metabolites({
        gssh : -1, o2 : -1, water : -1, gsh : 1, sulfite : 1, proton_m : 2
    })
    rxns[16].gene_reaction_rule = 'ENSG00000105755'
    # add in missing 5'-deoxyadenosine metabolite & reaction
    rxns[17].add_metabolites({deoxy_a_m : -1, deoxy_a_c : 1})
    # just completely redo the LIAS reaction
    lias_rxn = new_model.reactions.get_by_id('MAR06403')
    lias_rxn.subtract_metabolites(lias_rxn.metabolites)
    lias_rxn.add_metabolites({
        oct_gcsh : -1, four_isc : -1, sam : -2, nadph : -1, proton_m : -3,
        dhlipoyl_gcsh : 1, fe_2 : 4, deoxy_a_m : 2, met : 2, hs : 2, nadp : 1
    })
    lias_gpr = 'ENSG00000121897 and ENSG00000137714 and ENSG00000140905 and '
    lias_gpr += 'ENSG00000161513 and (ENSG00000135070 or ENSG00000165898 or '
    lias_gpr += 'ENSG00000169599)'
    lias_rxn.gene_reaction_rule = lias_gpr
    # move the rest of the LIPT2 reaction from cyto to mito and fix GPR
    lipt2_rxn = new_model.reactions.get_by_id('MAR06402')
    # also drop apoA1 from the reactants and add mitoACP as a product
    lipt2_rxn.add_metabolites({
        apo_a1 : 1, mito_acp : 1, octa_acp_c : 1, octa_acp_m : -1
    })
    lipt2_gpr = 'ENSG00000004779 and ENSG00000140905 and ENSG00000175536'
    lipt2_rxn.gene_reaction_rule = lipt2_gpr
    # move rest of DLD reaction from cytosol to mitochondria
    new_model.reactions.get_by_id('MAR06404').add_metabolites({
        # NADH and proton are reactants, so coeffs are currently -1, so adding
        # with coeffs of one brings them to 0, i.e. removes them
        nadh_c : 1, nadh_m : -1, proton_c : 1, proton_m : -1,
        # and switch signs for NAD+ cuz it's a product
        nadh_c : -1, nadh_m : 1
    })
    new_model.add_reactions(rxns)
    return(new_model)

try:
    threads = int(sys.argv[1])
except:
    sys.exit('provide number of threads for doing FVA')

# read in version 1.18 of Human-GEM and make a copy with all the rest of my
# proposed changes to lipoic acid metabolism
print('setting up')
optlang.glpk_interface.Configuration()
old_model = cobra.io.read_sbml_model('GSMMs/Human-GEMv1.15.xml')
newer_model = cobra.io.read_sbml_model('GSMMs/Human-GEMv1.18.xml')
new_model = fix_lipoate_biosynth(newer_model)
# set lower bounds on exchange reactions to -1000 for everything that's in
# DMEM or FBS or 0 for everything else
media_concs = pd.read_csv('input_data/DMEM-FBS_ingredients.csv')
in_media = media_concs[
    (media_concs['DMEM'] != '0') | (media_concs['FBS'] != '0')
]['metabolite_id'].to_list()
for model in [old_model, new_model]:
    for r in model.boundary:
        if list(r.metabolites)[0].id in in_media:
            r.lower_bound = -1000
        else:
            r.lower_bound = 0
# ensure both models have no objective function
old_model.objective = optlang.symbolics.Zero
new_model.objective = optlang.symbolics.Zero
# do FVA once on each model so we don't have to redo it each time we impose
# dilution constraints
old_fva = fva(old_model, threads = threads)
new_fva = fva(new_model, threads = threads)
fva_dict = {'Human-GEM 1.15' : old_fva, 'Human-GEM 1.18+' : new_fva}

# get % reduction in maximum possible flux through GCS and PDH reactions in
# both models after knocking out LIPT1 with and without also imposing dilution
# constraints
model_dict = {'Human-GEM 1.15' : old_model, 'Human-GEM 1.18+' : new_model}
gene_dict = {
    'GLRX5': 'ENSG00000182512', 'IBA57' : 'ENSG00000181873',
    'LIAS' : 'ENSG00000121897', 'LIPT1' : 'ENSG00000144182'
}
# start with observed consequences of mutations in each of these genes in real
# patients to compare the predicted knockout results to
ko_results = {
    'condition' : ['Patients'] * 4,
    'gene' : ['GLRX5', 'IBA57', 'LIAS', 'LIPT1'],
    'impact' : ['Both blocked'] * 3 + ['PDH blocked, no impact on GCS'],
}
for gene in gene_dict.keys():
    for model_name in ['Human-GEM 1.15', 'Human-GEM 1.18+']:
        for dilution in ['without', 'with']:
            # get appropriate model and make a copy
            model = model_dict[model_name].copy()
            # impose dilution constraints if we're supposed to
            if dilution == 'with':
                model = add_dilution_constraints(
                    model, fva_results = fva_dict[model_name],
                    dil_factor = 1000, verbose = 0
                )
            # get max possible fluxes through GCS and PDH reactions before
            # knocking anything out
            with model as wt_model:
                wt_model.objective = 'MAR08433'
                wt_gcs = round(wt_model.slim_optimize(), 2)
                wt_model.objective = 'MAR08746'
                wt_pdh = round(wt_model.slim_optimize(), 2)
            # knock out the appropriate gene and reassess max fluxes
            with model as ko_model:
                # the old model doesn't have most of these genes in it
                try:
                    ko_model.genes.get_by_id(gene_dict[gene]).knock_out()
                    ko_model.objective = 'MAR08433'
                    ko_gcs = round(ko_model.slim_optimize(), 2)
                    ko_model.objective = 'MAR08746'
                    ko_pdh = round(ko_model.slim_optimize())
                except KeyError:
                    ko_gcs = wt_gcs
                    ko_pdh = wt_pdh
            # get the % differences
            if wt_gcs != 0:
                gcs_diff = round(100 * (wt_gcs - ko_gcs) / wt_gcs)
            else:
                gcs_diff = 0
            if wt_pdh != 0:
                pdh_diff = round(100 * (wt_pdh - ko_pdh) / wt_pdh)
            else:
                pdh_diff = 0
            msg = f'Knocking out {gene} in {model_name} model {dilution} '
            msg += 'dilution constraints reduced maximum possible fluxes '
            msg += f'through GCS by {gcs_diff}% and PDH by {pdh_diff}%.'
            print(msg)
            # update output dict
            ko_results['condition'].append(f'{model_name} {dilution} Dilution')
            ko_results['gene'].append(gene)
            if (wt_gcs == ko_gcs) and (wt_pdh == ko_pdh):
                ko_results['impact'].append('No impact on either')
            elif (wt_gcs > ko_gcs) and (wt_pdh == ko_pdh):
                ko_results['impact'].append('GCS blocked, no impact on PDH')
            elif (wt_gcs == ko_gcs) and (wt_pdh > ko_pdh):
                ko_results['impact'].append('PDH blocked, no impact on GCS')
            elif (wt_gcs > ko_gcs) and (wt_pdh > ko_pdh):
                ko_results['impact'].append('Both blocked')
            else:
                ko_results['impact'].append(f'GCS: {gcs_diff}; PDH: {pdh_diff}')

# save knockout results as DataFrame and edited version of Human-GEM 1.18 model
pd.DataFrame(ko_results).to_csv('figure_data/fig_6c_data.csv', index = False)
cobra.io.save_json_model(new_model, 'figure_data/fig_6bc_model.json')
