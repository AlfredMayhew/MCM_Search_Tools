#imports
import os
import MCMSpeciesInfo
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image
from PIL import ImageFont
import textwrap
from PIL import ImageDraw
from rdkit.Chem.Draw import IPythonConsole 

#Structural features accepted
acceptedkwds = ["subset","C","H","O","N","NO3","OH","OOH","CO",
              "CCDoubleBond","PAN","OO.","O."]

args = os.sys.argv
print(f"Supplied arguments: {args[1:]}")
if len(args) == 1:
    raise Exception(f"""No arguments provided.
List the number of structural features of interest. E.g. 'C=6 Alcohol=1 Carbonyl=2 ',
An MCM subset can also be specified.
Accepted features are: {', '.join(acceptedkwds)}""")

#filter arguments and add them to a dictionary
argdict = {}

for i,a in enumerate(args):
    if i == 0:
        pass
    else:
        kwd = a.split("=")[0]
        val = a.split("=")[1]
        if (kwd in acceptedkwds) and (kwd != "subset"):
            argdict[kwd] = float(val)
        elif (kwd in acceptedkwds) and (kwd == "subset"):
            argdict[kwd] = val
        else:
            raise Exception(f"""Invalid Input: {kwd} is not a recognised species property
                            Valid inputs are: {', '.join(acceptedkwds)}""")


#get dictionary of mcm compound properties
if "subset" in argdict.keys():
    propDict = MCMSpeciesInfo.StructuralDict(subset=argdict["subset"])
else:
    propDict = MCMSpeciesInfo.StructuralDict()
argdict.pop('subset', None)

#Go through dictionary of properties and match compounds with matching properties
matchlist = []

for s in propDict:
    match = True #match will be set to false on disagreement for this compound
    for a in argdict:
        if propDict[s][a] != argdict[a]:
            match = False
            break
    
    if match == True:
        matchlist.append(s)

if matchlist:
    print(f"""{len(matchlist)} Matching MCM Species:
      
          {", ".join(matchlist)}""")
else:
    raise Exception("There are no MCM species with these requirements. Please consider broadening your search")
      
if len(matchlist) > 999:
    raise Exception("List of MCM matches is longer than 999 compounds. Please consider refining your search")
else: 
    print("Drawing structures to output file...")

#save structures to output file
mcmdict = MCMSpeciesInfo.SMILESDict()
molecules = [Chem.MolFromSmiles(mcmdict[x]) for x in matchlist]
img=Draw.MolsToGridImage(molecules,molsPerRow=6,subImgSize=(200,200),
                         legends=[x for x in matchlist],maxMols=999,
                         returnPNG=False) 


# add title to image stating functional groups selected
draw = ImageDraw.Draw(img)
titletext=f"Functional Groups: {', '.join(args[1:])}"
font = ImageFont.load_default().font

#wrap text
margin = offset = 5
for line in textwrap.wrap(titletext, width=int(img.size[0]/font.getsize("a")[0])-margin):
    draw.text((margin, offset), line,font=font,fill="#0000")
    offset += font.getsize(line)[1]
img.save('MCM_Search_Result.png')

print("...Complete! Look for `MCM_Search_Result.png`.")

