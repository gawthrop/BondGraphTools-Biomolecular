import numpy as np
import pandas as pd
import copy

def Integer(s):
    """ Make a stoichiometric matrix with integer elements
    """

    N = s["N"]                  # Stoichiometric matrix
    reaction = s["reaction"]
    species = s["species"]
    
    NT = N.T                    # Transpose
    for j,row in enumerate(NT):
        # I = np.nonzero(row)
        # Coeff = row[I]
        I_nonInt = np.nonzero(row-row.astype(int))
        if len(I_nonInt[0])>0:
            if len(I_nonInt[0])>1:
                print("Extract.Integer only handles one non-integer per reaction")

            i = I_nonInt[0][0]
            # print(i,reaction[j],species[i],row[i])
            #print(row[np.nonzero(row)])
            factor = 1/abs(row[i])
            row *= factor
            print("Multiplying reaction",reaction[j], "(",j,")",
                  "by",factor,"to avoid non-integer species",
                  species[i],"(",i,")")
            #print(row[np.nonzero(row)])

    N = NT.T
    s['N'] = N.astype(int)

    return s

def extract(filename='ecoli_core_model_no_biomass.xlsx',quiet=False):

    if not quiet:
        print('Extracting stoichiometric matrix from:',filename)
    
    df = pd.read_excel(filename,index_col=0)
    #print(df)
    species = [x.upper() for x in list(df.index)]
    reaction = [x.upper() for x in list(df)]
    N = df.to_numpy()
    nX = len(species)
    nV = len(reaction)
    if not quiet:
        print("nX:", nX)
        print("nV:", nV)


    ## Remove () and [] and - from the  names
    for i,spec in enumerate(species):
        species[i] = spec.replace('(','_').replace(')','').replace('[','').replace(']','').replace('-','')
    for i,reac in enumerate(reaction):
        reaction[i] = reac.replace('(','').replace(')','').replace('[','').replace(']','').replace('-','')

    ## Look for reactions with no RHS and replace by chemostats
    Nf = -((N<0)*N)
    Nr = (N>0)*N
    NfT = Nf.T
    NrT = Nr.T
    J = []
    chemostats = []
    reac = []
    #zap_reac = ['ATPM','ADK1']
    zap_reac = []
    for j in range(nV-1,-1,-1):
        if np.max(NrT[j]) == 0:
            i = np.flatnonzero(NfT[j])[0]
            if not quiet:
                print('Deleting reaction:',reaction[j],'and adding chemostat:', species[i])
            chemostats.append(species[i])
            N = np.delete(N,j,1)
        elif reaction[j] in zap_reac:
            print('Deleting reaction:',reaction[j])
            N = np.delete(N,j,1)
        else:
            reac.append(reaction[j])

    ## Reactions now in wrong order: reverse them
    reac = list(reversed(reac))
    nX = len(species)
    nV = len(reac)
    if not quiet:
        print("nX:", nX)
        print("nV:", nV)


    if not quiet:
        print(len(chemostats),'chemostats')
    #print(len(chemostats),'chemostats:',chemostats)
    s = {}
    s['N'] = N
    s['species'] = species
    s['reaction'] = reac
    s['chemostats'] = chemostats
    return s

def getSpecies(s,reaction):
    """ Extract the species associated with a list of reactions
    """

    N = s["N"]
    NT = N.T
    Species = s["species"]
    Reaction = s["reaction"]
    
    species = []
    for reac in reaction:
        Nj = NT[Reaction.index(reac)]
        I = np.nonzero(Nj)[0]
        for i in I:
            species.append(Species[i])
    species = list(set(species)) # Make unique
    species.sort()               # Sort
    return species
        
    
def choose(s,reaction=[]):
    """ Choose subset reactions.
    """

    NN = s["N"]
    Species = s["species"]
    Reaction = s["reaction"]

    species = getSpecies(s,reaction)
    n_X = len(species)
    n_V = len(reaction)
    
    N = np.zeros((n_X,n_V))
    for i,spec in enumerate(species):
        for j,reac in enumerate(reaction):
            N[i,j] = NN[Species.index(spec),Reaction.index(reac)]

    N = N.astype(int)
    
    ## Forward and reverse
    Nf = -((N<0)*N)
    Nr = (N>0)*N

    sub = {}
    sub["N"] = N
    sub["Nf"] = Nf
    sub["Nr"] = Nr
    sub["species"] = species
    sub["reaction"] = reaction
    
    return sub
