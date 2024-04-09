#imports
import os
import search_tools
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image
from PIL import ImageFont
import textwrap
from PIL import ImageDraw
from rdkit.Chem.Draw import IPythonConsole 

#Structural features accepted
func_smarts = {"ONO2" : "[$([NX3](=[OX1])(=[OX1])O),$([NX3+]([OX1-])(=[OX1])O)]", #nitrate
               "OH" : "[#6][OX2H]", #alcohol
               "OOH" : "[OX2][OX2H]", #hydroperoxide
               "C=O" : "[CX3]=[OX1]", #Carbonyl
               "C=C" : "[C]=[C]", #carbon-carbon double bond
               "C(=O)OONO2" : "[$(C(=O)OO[NX3](=[OX1])(=[OX1])),$(C(=O)OO[NX3+]([OX1-])(=[OX1]))]", #peroxy acetyl nitrate (pan) 
               "OO." : "[OX2][OX1+0]", #peroxy radical   
               "O." : "[#6][OX1+0]", #alcoxy radical
               }
synonyms = {"ONO2" : ["nitrate", "organonitrate", "NO3"],
            "OH" : ["alcohol", "hydroxy", "hydroxyl"],
            "OOH" : ["hydroperoxide", "hydroperoxy", "hydroperoxyl"],
            "C=O" : ["CO", "carbonyl"],
            "C=C" : ["CCDoubleBond", "carbon-carbon double bond", "double bond",
                     "carbon double bond"], 
            "C(=O)OONO2" : ["PAN"], 
            "OO." : ["peroxy radical", "peroxyl radical", "peroxyl", "peroxy"],  
            "O." : ["alkoxy radical", "alkoxyl radical", "alkoxy", "alkoxyl"], 
            }

args = os.sys.argv
print(f"Supplied arguments: {args[1:]}")
if len(args) <= 1:
    raise Exception(f"""No arguments provided.
List the number of structural features of interest. E.g. 'alcohol=1 C=O=2 C=6',

Accepted features are: {', '.join(func_smarts.keys())}

Any other arguments will be treated as chemical element symbols (e.g. C or H).""")

#filter arguments and add them to a dictionary
argdict = {}
all_syns = [x for y in synonyms.values() for x in y]
for i,a in enumerate(args):
    if i == 0: #the first argument will be this file's name
        pass
    else:
        kwd = a.rsplit("=", 1)[0]
        val = a.rsplit("=", 1)[1]
        
        #check that the value can be converted to a float
        try:
            float(val)
        except ValueError:
            raise Exception(f"Invalid value of {val} for provided keyword argument of {kwd}.")
        
        if (kwd.casefold() in [x.casefold() for x in func_smarts.keys()]):
            argdict[kwd] = float(val)
        elif (kwd.casefold() in [x.casefold() for x in all_syns]):
            #check which synonym this is
            for k,v in synonyms.items():
                if kwd.casefold() in [x.casefold() for x in v]:
                    argdict[k] = float(val)
                    break
        else:
            #try to parse the argument as an elemental symbol.
            pt = Chem.GetPeriodicTable()
            try:
                atom_num = pt.GetAtomicNumber(kwd[0].upper() + kwd[1:].lower())
            except:
                raise Exception(f"""Invalid Input: {kwd} is not a recognised species property
                                Valid inputs are: {', '.join(func_smarts.keys())} or elemental symbols.""")
            
            argdict[atom_num] = float(val)

###############################################################################
#get the mcm compound dataframe
mass_df = search_tools.MCM_Masses()

#function to check whether a given dict entry matches for a given molecule
def check_single_match(keyval, mol):
    k = keyval[0]
    v = keyval[1]
    
    if type(k) == str:
        matches = mol.GetSubstructMatches(Chem.MolFromSmarts(func_smarts[k]))
        return len(matches) == v
    elif type(k) == int:
        matches = mol.GetSubstructMatches(Chem.MolFromSmarts(f"[#{k}]"))
        return len(matches) == v
    else:
        raise Exception(f"Unrecognised type for argdict: {type(k)}")
#function to check all of the key-value pairs in the arg dict
def check_all_matches(mol):
    h_mol = Chem.AddHs(mol)
    for kv in argdict.items():
        if not check_single_match(kv, h_mol):
            return False
    return True

#apply this function to the mass df to find MCM species that fit the criteria
applied_series = mass_df["Mols"].apply(check_all_matches)
matches = applied_series[applied_series].index

if len(matches) > 0:
    print(f"""{len(matches)} Matching MCM Species:
      
          {", ".join(matches)}""")
else:
    raise Exception("There are no MCM species with these requirements. Please consider broadening your search.")
      
if len(matches) > 100:
    raise Exception("List of MCM matches is longer than 100 compounds. Please consider refining your search in order to draw the output.")
else: 
    print("Drawing structures to output file...")

###############################################################################
#save structures to output file
molecules = mass_df.loc[matches, "Mols"]
img=Draw.MolsToGridImage(molecules.values,molsPerRow=6,subImgSize=(200,200),
                         legends=[x for x in molecules.index],maxMols=100,
                         returnPNG=False) 


# add title to image stating functional groups selected
draw = ImageDraw.Draw(img)
titletext=f"Functional Groups: {', '.join(args[1:])}"
font = ImageFont.load_default().font

#wrap text
margin = offset = 5
for line in textwrap.wrap(titletext, width=int(img.size[0]/font.getsize("a")[0][0])-margin):
    draw.text((margin, offset), line,font=ImageFont.load_default(),fill="#0000")
    offset += font.getsize(line)[0][1]
img.save('MCM_Search_Result.png')

print("...Complete! Look for `MCM_Search_Result.png`.")

