"""Convert a stoichiometric representation to BondGraphTools"""
import BondGraphTools as bgt
import numpy as np
import datetime

def model(s,filename=None):
    """ Converts a stoichiometric matrix  to the BondGraphTools bond graph representation
    """

    ## Extract info
    N = s["N"]                  # Stoichiometric matrix

    ## Copy or create forward and reverse stoichiometric matrices
    if "Nf" in s:
        #print("Using Nf")
        Nf = s["Nf"]
    else:
        Nf = -N*(N<0)

    if "Nr" in s:
        Nr = s["Nr"]
    else:
        Nr = N*(N>0)
        
    shape = N.shape
    n_X = shape[0]              # Number of species (states)
    n_V = shape[1]              # Number of reactions (flows)
 
    species = s["species"]      # List of species names
    reaction = s["reaction"]    # List of reaction names

    if "name" in s.keys():
        name = s["name"]
    else:
        name = "none"

    ## Sanity check
    if not ((n_X == len(species))):
        print("Number of species (",len(species),") not equal to number of rows (", n_X,").")

    if not ((n_V == len(reaction))):
        print("Number of reactions (",len(reaction),") not equal to number of columns (", n_V,").")


    ## Create model header
    indent = "    "
    delim = '"""'
    headStr = ("import BondGraphTools as bgt\n"
              "import sympy as sp\n\n"
              "def model():\n"
              "{0}{3} Acausal bond graph {1}.py\n"
              "{0}Created by stoichBondGraph at {2}\n\n"
              "{0}Usage:\n"
              "{0}import {1}; model = {1}.model()\n"
              "{0}{3}\n\n"
              '{0}model = bgt.new(name="{1}")\n'
     )

    head = headStr.format(indent,name, datetime.datetime.now().ctime(),delim)

    ## print(head)

    ## Formatting
    CeStr = ("\n{0}### Ce:{1}\n"
             "\n{0}## Component Ce:{2}{1}\n"
             "{0}K_{1} =  sp.symbols('K_{1}')\n"
             "{0}RT = sp.symbols('RT')\n"
             "{0}{2}{1} = bgt.new('Ce',name='{1}',value={{'k':K_{1},'R':RT,'T':1}},library='BioChem')\n"
             "{0}model.add({2}{1})\n"
    )

    junStr = ("\n{0}## Junction {1}:{2}\n"
              "{0}{2} = bgt.new('{1}')\n"
              "{0}model.add({2})\n"
    )

    bondStr = ("\n{0}## Bond from {1} to {2}\n"
               "{0}bgt.connect({1},{2})\n"
               )

    conCR = ("\n{0}### Ce:{1} to Re:{2}\n")
    conRC = ("\n{0}### Re:{1} to Ce:{2}\n")

    prepend = "comp_"
    junName = "{1}{0}_jun"  # Junction name
    junF = "{1}{0}_junF"  # Forward junction name
    junR = "{1}{0}_junR"  # Reverse junction name

    specStr = indent+"### Species\n"

    ## Create species
    for spec in species:
        specStr += CeStr.format(indent,spec,prepend)
        specStr += junStr.format(indent,"0",junName.format(spec,prepend))
        specStr += bondStr.format(indent,junName.format(spec,prepend),prepend+spec)

        ## print(specStr)
        
    ## Create reactions
    ReStr = ("\n{0}### Re:{1}\n"
             "\n{0}## Component Re:{2}{1}\n"
             "{0}kappa_{1} =  sp.symbols('kappa_{1}')\n"
             "{0}RT = sp.symbols('RT')\n"
             "{0}{2}{1} = bgt.new('Re',name='{1}',value={{'r':kappa_{1},'R':RT,'T':1}},library='BioChem')\n"
             "{0}model.add({2}{1})\n"
    )

    reacStr = indent+"### Reactions\n"
    for reac in reaction:
        reacStr += ReStr.format(indent,reac,prepend)
        reacStr += junStr.format(indent,"1",junF.format(reac,prepend))
        portF = "("+prepend+reac+",0)"
        reacStr += bondStr.format(indent,junF.format(reac,prepend),portF)
        reacStr += junStr.format(indent,"1",junR.format(reac,prepend))
        portR = "("+prepend+reac+",1)"
        reacStr += bondStr.format(indent,portR,junR.format(reac,prepend))
    ## print(reacStr)

                
    ## Create the connections
    connectStr = indent+"### Connections\n"
    for i,row in enumerate(Nf):
        for j,coeff in enumerate(row):
            if (coeff>0):
                for k in range(coeff):
                    #print(species[i],reaction[j])
                    connectStr += conCR.format(indent,species[i],reaction[j])
                    connectStr += bondStr.format(indent,junName.format(species[i],prepend),junF.format(reaction[j],prepend))

    for i,row in enumerate(Nr):
        for j,coeff in enumerate(row):
            if (coeff>0):
                for k in range(coeff):
                    connectStr += conRC.format(indent,reaction[j],species[i])
                    connectStr += bondStr.format(indent,junR.format(reaction[j],prepend),junName.format(species[i],prepend))



            # if (coeff>0):
            #     for k in range(abs(coeff)):
            #         connectStr += conRC.format(indent,reaction[j],species[i])
            #         connectStr += bondStr.format(indent,junR.format(reaction[j],prepend),junName.format(species[i],prepend))

    ## print(connectStr)

    ## Print to file
    if filename is None:
        Filename = name+".py"
    else:
        Filename = filename+".py"
        
    f = open(Filename,'w')
    f.write(head)
    f.write(specStr)
    f.write(reacStr)
    f.write(connectStr)
    f.write("\n"+indent+"return model\n")
    f.close

