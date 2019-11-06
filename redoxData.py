import numpy as np
import stoich as st
import scipy.constants as const

def mV(E):
    return(str(int(round(1000*E))))

def data():
    """Redox midpoint (pH 7) data
    """

    redoxData = {}
    
    ## Data from Bla14

    ## NADP/NADPH
    dat= {}
    dat['oxidised'] = 'NADP'
    dat['reduced'] = 'NADPH'
    dat['E7'] = -0.324
    dat['electrons'] = 2
    dat['protons'] = 1
    dat['source'] = 'Bla14:TableA1.2'
    redoxData['NADP/NADPH'] = dat

    ## NAD/NADH
    dat= {}
    dat['oxidised'] = 'NAD'
    dat['reduced'] = 'NADH'
    dat['E7'] = -0.320
    dat['electrons'] = 2
    dat['protons'] = 1
    dat['source'] = 'FalRav07:Table4.1'
    redoxData['NAD/NADH'] = dat
    
    ## O2/2H2O
    dat = {}
    dat['oxidised'] = 'O2'
    dat['reduced'] = 'H2O'
    dat['E7'] = 0.816
    dat['electrons'] = 4
    dat['protons'] = 4
    dat['source'] = 'Bla14:TableA1.2'
    redoxData['O2/2H2O'] = dat

    ## P700/P700+
    dat = {}
    dat['photon'] = True
    dat['oxidised'] = 'P700'
    dat['reduced'] = 'P700+'
    dat['E7'] = 0.49
    dat['electrons'] = 1
    dat['protons'] = 0
    dat['source'] = 'Bla14:TableA1.2'
    redoxData['P700/P700+'] = dat

    ## P870/P870+
    dat = {}
    dat['photon'] = True
    dat['oxidised'] = 'P870'
    dat['reduced'] = 'P870+'
    dat['E7'] = 0.45
    dat['electrons'] = 1
    dat['protons'] = 0
    dat['source'] = 'Bla14:TableA1.2'
    redoxData['P870/P870+'] = dat

    ## P680/P680+
    dat = {}
    dat['photon'] = True
    dat['oxidised'] = 'P680'
    dat['reduced'] = 'P680+'
    dat['E7'] = 1.1
    dat['electrons'] = 1
    dat['protons'] = 0
    dat['source'] = 'Bla14:TableA1.2'
    redoxData['P680/P680+'] = dat

    ## P700+/P700*
    dat = {}
    dat['photon'] = True
    dat['oxidised'] = 'P700+'
    dat['reduced'] = 'P700*'
    dat['E7'] = 1.26
    dat['electrons'] = 1
    dat['protons'] = 0
    dat['source'] = 'Bla14:TableA1.3'
    redoxData['P700+/P700*'] = dat

    ## P870+/P870*
    dat = {}
    dat['photon'] = True
    dat['oxidised'] = 'P870+'
    dat['reduced'] = 'P870*'
    dat['E7'] = 0.94
    dat['electrons'] = 1
    dat['protons'] = 0
    dat['source'] = 'Bla14:TableA1.3'
    redoxData['P870+/P870*'] = dat

    ## P680+/P680*
    dat = {}
    dat['photon'] = True
    dat['oxidised'] = 'P680+'
    dat['reduced'] = 'P680*'
    dat['E7'] = 0.8
    dat['electrons'] = 1
    dat['protons'] = 0
    dat['source'] = 'Bla14:TableA1.3'
    redoxData['P680+/P680*'] = dat

    ## PQ/PQH2 (AKA Plastoquinone/Plastoquinol, AKA UQ/UQH2)
    dat = {}
    dat['oxidised'] = 'PQ'
    dat['reduced'] = 'PQH2'
    dat['E7'] = 0
    dat['electrons'] = 2
    dat['protons'] = 2
    dat['source'] = 'Bla14:TableA1.2'
    redoxData['PQ/PQH2'] = dat

    
    ## Plastocyanin (ox/red), PcOx/PcRed)
    dat = {}
    dat['oxidised'] = 'PcOx'
    dat['reduced'] = 'PcRed'
    dat['E7'] = 0.380
    dat['electrons'] = 1
    dat['protons'] = 0
    dat['source'] = 'FalRav07:Table4.1'
    redoxData['PcOx/PcRed'] = dat

    ## Feredoxin (ox/red), FdOx/FdRed)
    dat = {}
    dat['oxidised'] = 'FdOx'
    dat['reduced'] = 'FdRed'
    dat['E7'] = -0.430
    dat['electrons'] = 1
    dat['protons'] = 0
    dat['source'] = 'FalRav07:Table4.1'
    redoxData['FdOx/FdRed'] = dat


    return redoxData

def Phi(couple):
    """ Net potential for couple without protons.
    """

    dat = data()[couple]
    phi_ph7 = st.V_N()*np.log(1e-7) # phi of protons at pH 7
    Phi7 = dat['electrons']*dat['E7']   # Phi with protons
    Phi = Phi7 - dat['protons']*phi_ph7 # Phi without protons
    #E = Phi/dat['electrons']
    return Phi

def phi():
    """ 
    Returns dict of phi deduced from redox poentials
    phis are relative to the couple
    """
    phi = {}
    for couple,dat in data().items():
 
        ## Extract names of oxidised and reduced species
        ox = dat["oxidised"]
        re = dat["reduced"]
        # print(ox,re)

        ## Extract stoichiometry
        oxre = couple.split('/')
        # print(oxre[0],oxre[1])
        if ox in [oxre[0]]:
            n_ox = 1
        else:
            n_ox = int(oxre[0][0])
        if re in [oxre[1]]:
            n_re = 1
        else:
            n_re = int(oxre[1][0])
        # print(n_ox,n_re)

        Phi_oxre = Phi(couple)
        if "photon" in dat.keys():
            phi_ox = Phi_oxre
            phi_re = 0
        else:
            phi_ox = (Phi_oxre/2)/n_ox
            phi_re = -(Phi_oxre/2)/n_re

        ## Put phi into dict
        phi[ox] = phi_ox
        phi[re] = phi_re
        

    ## Put in zero potential for the E components
    phi['E1'] = 0
    phi['E2'] = 0

    return phi
    
def E(couple):
    """ Redox potential for couple without protons.
    """
    dat = data()[couple]
    E = Phi(couple)/dat['electrons']
    return E

def E7(couple):
    """ Redox potential for couple with protons at pH 7
    """
    dat = data()[couple]
    E7 = dat['E7']
    return E7

def VpH(pH=7):
    """
    Voltage corresponding to pH
    """
    conc = np.exp(-pH*np.log(10)) # Concentration
    return st.V_N()*np.log(conc)

def EpH(couple,pH=7):
    """ Redox potential for couple with protons at given pH
    """
    E0 = E(couple)
    dat = data()[couple]
    EpH = E0 + VpH(pH)*(dat['protons']/dat['electrons'])
    return EpH

def V_photon(wavelength=680):
    """Voltage corresponding to the energy of a coulomb of photons
    wavelength in nm
    """
    F = const.physical_constants['Faraday constant'][0]
    N = const.physical_constants['Avogadro constant'][0]
    h = const.physical_constants['Planck constant'][0]
    c = const.physical_constants['speed of light in vacuum'][0]

    v = c/(wavelength*1e-9)
    
    return (N*h*v)/F

def print_raw():
    """ Print raw redox potential data
    """

    for couple, dat in data().items():
        print("\nRedox couple:",couple)
        print("Phi = ", int(1000*Phi(couple)),"mV")
        print("Phi/n = ", int(1000*Phi(couple)/dat["electrons"]),"mV")
        for key, val in dat.items():
            if key is "E7":
                print("\t",key+":",int(1000*val),"mV")
            else:
                print("\t",key+":",val)

def table():
    """ LaTeX table of redox potentials
    """
    sep = " & "
    head = ["E_7", "n", "m", "source"]
    header = "Couple"
    for h in head:
        header = header + sep + h
    print("\\begin{array}{|l|l|l|l|l|}")
    print("\\hline")
    print(header,"\\\\")
    print("\\hline")
    print_keys = ["electrons","protons","source"]
    for couple, dat in data().items():
        row = couple+sep+mV(dat['E7'])
        for key in print_keys:
            row = row+sep+str(dat[key])
        print(row,"\\\\")
    print("\\hline")
    print("\\end{array}")
 
