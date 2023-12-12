"""Provides information about species from the MCM"""
from rdkit import Chem
import re

masses_path = "mcm_data/Whole_MCM_Masses.txt"
fac_path = "mcm_data/Whole_MCM.fac"

def SMILESDict(subset="whole"):
    """Produces a dictionary of MCM compounds and their SMILES strings"""
    #read in compounds and masses from mcm output file
    if subset == "whole":
        with open(masses_path,"r") as file: 
            mcm=file.readlines() 
            for i,s in enumerate(mcm):
                mcm[i]=s.split()
    elif subset == "isoprene":
        with open("Isoprene_MCM_Masses.txt","r") as file: 
            mcm=file.readlines() 
            for i,s in enumerate(mcm):
                mcm[i]=s.split()
    else:
        raise Exception(f"Unrecognised subset {subset}")
        
    #make dictionary of MCM compounds and SMILES
    mcmdict={}
    for l in mcm:
        mcmdict[l[0]]=l[1]
        
    return mcmdict

def StructuralDict(subset = "whole"):
    """Produces a dictionary of MCM compounds and a selection of their 
    structural features/ functional groups"""
    #read in compounds and masses from mcm output file
    if subset == "whole":
        with open(masses_path,"r") as file: 
            mcm=file.readlines() 
            for i,s in enumerate(mcm):
                mcm[i]=s.split()
    elif subset == "isoprene":
        with open("Isoprene_MCM_Masses.txt","r") as file: 
            mcm=file.readlines() 
            for i,s in enumerate(mcm):
                mcm[i]=s.split()
    else:
        raise Exception(f"Unrecognised subset {subset}")
            
    #list of functional groups to search for
    #number of different atoms (C,H,N,O)
    car=Chem.MolFromSmarts("[C,c]")
    hyd=Chem.MolFromSmarts("[H]")
    oxy=Chem.MolFromSmarts("[O,o]")
    nit=Chem.MolFromSmarts("[N,n]")
    
    #nitrate
    nitrate = Chem.MolFromSmarts("[$([NX3](=[OX1])(=[OX1])O),$([NX3+]([OX1-])(=[OX1])O)]")
    #alcohol
    alcohol = Chem.MolFromSmarts("[#6][OX2H]")
    #hydroperoxide
    hydroperox= Chem.MolFromSmarts("[OX2][OX2H]")
    #Carbonyl
    carbonyl=Chem.MolFromSmarts("[CX3]=[OX1]")
    #carbon-carbon double bond
    ccdoub=Chem.MolFromSmarts("[$([CX3]=[CX3]),$([CX2](=C)=C)]")
    #peroxy acetyl nitrate (pan)
    pan=Chem.MolFromSmarts("[$(C(=O)OO[NX3](=[OX1])(=[OX1])),$(C(=O)OO[NX3+]([OX1-])(=[OX1]))]")
    
    #peroxy radical
    peroxyrad=Chem.MolFromSmarts("[OX2][OX1+0]")
    #alcoxy radical
    alcoxyrad=Chem.MolFromSmarts("[#6][OX1+0]")
    
    #compile into list of functional groups
    groups=[car,hyd,oxy,nit,nitrate,alcohol,hydroperox,carbonyl,ccdoub,pan,
            peroxyrad,alcoxyrad]
    groupnames=["C","H","O","N","NO3","OH","OOH","CO",
                "CCDoubleBond","PAN","OO.","O."]

    #make dictionary of MCM compounds and structural desciptors
    mcmdict={}
    for l in mcm:
        mcmdict[l[0]]={}
        #Get molecular formula in rdkit
        mol=Chem.rdmolops.AddHs(Chem.MolFromSmiles(l[1]))
        
        #Go through functional groups and count occurences
        for i,g in enumerate(groups): 
            if g != ccdoub:
                mcmdict[l[0]][groupnames[i]]=len(mol.GetSubstructMatches(g))
            else: #need to divide the number of CCdouble bonds by two (two atoms in each bond)
                mcmdict[l[0]][groupnames[i]]=len(mol.GetSubstructMatches(g))/2
        
    return mcmdict
    
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
