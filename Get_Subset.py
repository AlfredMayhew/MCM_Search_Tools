"""Creates a mechanism containing reaction pathways involving specified 
compounds from a larger mechanism"""
#imports
import sys
import re

if len(sys.argv) == 4:
    thisfname, in_mech, out_mech, target_specs = sys.argv
else:
    raise Exception(f"""Input invalid: '{" ".join(sys.argv)}' .
                    3 arguments required (in this order):
                    - Path to the input mechanism
                    - Path to save the new mechanism
                    - List of species to select subset for (comma separated)""")

#list of inorganic species whose chemistry will automatically be carried forward
inorgs = ['CO', 'H2', 'H2O2', 'HNO3', 'HO2', 'HO2NO2', 'HONO', 'HSO3', 'N2O5',
          'NA', 'NO', 'NO2', 'NO3', 'O', 'O1D', 'O3', 'OH', 'SA', 'SO2', 'SO3']
rxn_header = "* Reaction definitions. ;\n"
ro2_header = "* Peroxy radicals. ;\n"
###############################################################################
#read the input mechanism
with open(in_mech) as file:
    in_mech_lines = file.readlines()
    
#select only the reaction lines from the mechainsm
rxn_pat = "%.*:(.*)=(.*);\n"
rxn_lines = [x for x in in_mech_lines if re.match(rxn_pat, x)]

#go through the reactions of the selected species and store any reactions 
#involved in the reaction cascade
def find_reactions(spec):
    loss_rxns = []
    prod_rxns = []
    out_reacts = set()
    out_prods = set()
    
    for l in rxn_lines:
        #select the reactants and products for this reaction
        reacts = re.split(" *\+ *", re.match(rxn_pat, l)[1])
        reacts = [x.strip() for x in reacts]
        prods = re.split(" *\+ *", re.match(rxn_pat, l)[2])
        prods = [x.strip() for x in prods]

        if spec in reacts:
            loss_rxns.append(l)
            out_prods.update(prods)
        if spec in prods:
            prod_rxns.append(l)
            out_reacts.update(reacts)

    return loss_rxns, prod_rxns, out_prods, out_reacts

def return_all_specs(in_rxn_lines):
    specs = set()
    for l in in_rxn_lines:
        #select the reactants and products for this reaction
        reacts = re.split(" *\+ *", re.match(rxn_pat, l)[1])
        reacts = [x.strip() for x in reacts if x.strip()]
        prods = re.split(" *\+ *", re.match(rxn_pat, l)[2])
        prods = [x.strip() for x in prods if x.strip()]

        specs.update(reacts)
        specs.update(prods)
    return specs

#list to store all reactions to keep
sel_rxns = []

#select any inorganic reactions that should be carried forward
for l in rxn_lines:
    #select the reactants and products for this reaction
    reacts = re.split(" *\+ *", re.match(rxn_pat, l)[1])
    reacts = [x.strip() for x in reacts if x.strip()]
    prods = re.split(" *\+ *", re.match(rxn_pat, l)[2])
    prods = [x.strip() for x in prods if x.strip()]

    if all([x in inorgs for x in reacts]) and all([x in inorgs for x in prods]):
        sel_rxns.append(l)

#initialise a 'products' set that starts with our species of interest
next_prods = set(target_specs.split(","))
compl_prods = []

#iterate through the 'next prods' list and find all of the loss processes 
#removing the products as we go, until we have no new products left
while next_prods:
    sel_p = list(next_prods)[0]
    
    sel_loss_rxns, sel_prod_rxns, sel_prods, sel_reacts = find_reactions(sel_p)
    
    sel_rxns += sel_loss_rxns
    compl_prods.append(sel_p)
    
    next_prods.remove(sel_p)
    next_prods.update([x for x in sel_prods if ((x not in compl_prods) and 
                                                 (x not in inorgs))])

rxn_header_idx = in_mech_lines.index(rxn_header)

#edit the RO2 list to only include species in the new mechanism
ro2_header_idx = in_mech_lines.index(ro2_header)

ro2_lines = [x.strip() for x in in_mech_lines[ro2_header_idx:rxn_header_idx] if not x.startswith("*")]

old_ro2s = [x.strip() for x in "".join(ro2_lines).replace("RO2 = ","").replace(";","").split("+")]

new_ro2s = [x for x in old_ro2s if x in return_all_specs(sel_rxns)]

new_mech_lines = in_mech_lines[:ro2_header_idx+1]
new_mech_lines.append(f"RO2 = {' + '.join(new_ro2s)} ;\n")
new_mech_lines.append(rxn_header)
new_mech_lines += sel_rxns

with open(out_mech, "w") as file:
    file.writelines(new_mech_lines)

