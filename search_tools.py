"""Provides information about species from the MCM"""
from rdkit import Chem
import re
import pandas as pd

masses_path = "mcm_data/Whole_MCM_Masses.txt"
fac_path = "mcm_data/Whole_MCM.fac"

def MolToSmilesNoStereo(mol):
    Chem.RemoveStereochemistry(mol)
    return Chem.MolToSmiles(mol, canonical = True)
                                 
def MCM_Masses():
    spec_df = pd.read_csv(masses_path, skiprows=18, sep="\s+",
                      header=None, index_col=0, engine="python")
    spec_df.columns = ["SMILES", "ID", "mass"]
    spec_df["Mols"] = spec_df["SMILES"].apply(Chem.MolFromSmiles)
    spec_df["Canon_SMILES"] = spec_df["Mols"].apply(MolToSmilesNoStereo)
    
    return spec_df
        
def MCM_Reactions():
    """Produces a list of MCM reactions, with the rate, products, and reactants"""
    #read in compounds and masses from mcm output file
    with open(fac_path,"r") as file: 
        mcm_lines=file.read()
        mcm_lines = mcm_lines.split("* Reaction definitions. ;")[1].split("\n")
        
        
    reactions = []
    for line in mcm_lines:
        correct_format = re.compile("%.*:.*\=.*;")
        
        if re.match(correct_format, line):
            rate = line.split(":")[0].strip(" ;%")
            reaction = line.split(":")[1].strip(" ;%")
            reactants = [x.strip(" ;") for x in reaction.split("=")[0].split("+")]
            try:
                products = [x.strip(" ;") for x in reaction.split("=")[1].split("+")]
            except IndexError: #if the reaction has no products then an error will be thrown
                products = [""]
            
            reactions.append({"Rate" : rate,
                              "Reactants" : reactants,
                              "Products" : products})
    
    
    return reactions

def Get_Formation_Reactions(species):
    """Returns all formation reactions for a species"""
    all_reactions = MCM_Reactions()
    output_reactions = []
    for r in all_reactions:
        if species in r["Products"]:
            output_reactions.append(r)
            
    return output_reactions
