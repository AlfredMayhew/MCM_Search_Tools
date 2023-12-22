"""Searches all species in the MCM for one that matches a provided SMILES 
string (ignoring stereochemistry)"""
#imports
import os
import search_tools
from rdkit import Chem

args = os.sys.argv
print(f"Supplied arguments: {args[1:]}")
if len(args) != 2:
    raise Exception("""One argument must be provided: A SMILES string to search for.""")
else:
    this_fname, in_smiles = args
    
###############################################################################
#read the MCM masses file
mass_df = search_tools.MCM_Masses()

#convert inputted SMILES to cannonical SMILES (and removing stereochemistry) to 
#allow comparison to the MCM df
in_canon_smiles = search_tools.MolToSmilesNoStereo(Chem.MolFromSmiles(in_smiles))

###############################################################################
#Search for any matching MCM species
match_smiles = list(mass_df[mass_df["Canon_SMILES"] == in_canon_smiles].index)

if len(match_smiles) > 0:
    print(f"Found the following matches: {', '.join(match_smiles)}")
else:
    print("No matching MCM species found.")


