import numpy as np
#import scipy.constants as con

def thermoConstants(T=273.15):
    """Some useful thermodynamic constants."""
    
    R = 8.3144621;		# Gas constant J K^-1 mol^-1
    F = 9.64853399e4;           # Faraday constant C mol^-1
    RT = R*T;                   # Product of R and T
    V_N = RT/F;                 # Faraday voltage
    A = 6.02214129e23; 		# Avogadro's constant mol^-1

    const = {}                  # Initialise dictionary
    const["R"] = R;
    const["RT"] = RT;
    const["F"] = F;
    const["V_N"] = V_N;
    const["A"] = A;

    return const

def Gibbs_Std():
    """Gibbs free energy at standard conditions."""
    
    ## Gibbs free energy. Table A3 WuYanVin07
    ## (298.15 K, 1 M reactants, I = 0.17 M, P = 1 atm)
    ## Units are J/mol
    Gibbs = {}                  # Dictionary initialise.
    Gibbs["H2O"] = -235.74e3
    Gibbs["O2"] = 16.40e3
    Gibbs["NADH"] = 39.31e3
    Gibbs["NAD"] = 18.10e3
    Gibbs["QH2"] = -23.30e3
    Gibbs["Q"] = 65.17e3        # COQ
    Gibbs["ATP"] = -2771.00e3
    Gibbs["ADP"] = -1903.96e3
    Gibbs["AMP"] = -1034.66e3
    Gibbs["Fe2"] = -27.41e3	# Cred
    Gibbs["Fe3"] = -6.52e3	# Cox
    Gibbs["PI"] = -1098.27e3
    Gibbs["Cr"] = -252.68e3
    Gibbs["FADH2"] = -67.60e3
    Gibbs["FAD"] = 19.55e3
    Gibbs["COASH"] = -0.72e3
    Gibbs["ACCOA"] = -178.19e3
    Gibbs["OAA"] = -794.41e3
    Gibbs["CIT"] = -1165.59e3
    Gibbs["ICIT"] = -1158.94e3
    Gibbs["AKG"] = -793.41e3
    Gibbs["SCOA"] = -507.55e3
    Gibbs["SUC"] = -690.44e3
    Gibbs["FUM"] = -603.32e3
    Gibbs["MAL"] = -842.66e3
    Gibbs["ASP"] = -692.26e3
    Gibbs["GLU"] = -692.40e3
    Gibbs["PYR"] = -470.82e3
    Gibbs["GLC"] = -907.21e3
    Gibbs["G6P"] = -1758.87e3
    Gibbs["CO2"] = -530.71e3


    ## From LiWuQi11 Table 5. - (April 2019 translation)
    Gibbs["GLC"] = -916.39e3
    ##Gibbs["ATP"] = -2770.03e3
    Gibbs["F6P"] = -1563.96e3
    Gibbs["F16P"] = -2401.08e3
    Gibbs["DHAP"] = -1194.65e3
    Gibbs["GAP"] = -1187.64e3
    Gibbs["BPG"] = -2276.28e3
    ##Gibbs["NADH"] = 43.91e3
    Gibbs["PG3"] = -1429.38e3
    Gibbs["PG2"] = -1423.49e3
    Gibbs["PEP"] = -1190.84e3
    Gibbs["PYR"] = -393.85e3
    Gibbs["OAA"] = -714.06e3
    Gibbs["CIT"] = -1022.44e3
    Gibbs["ISCIT"] = -1016.69e3
    Gibbs["NADP"] = -834.79e3
    Gibbs["AKG"] = -676.45e3
    Gibbs["SUC"] = -589.56e3
    Gibbs["GDP"] = -1904.22e3
    Gibbs["FUM"] = -500.99e3
    Gibbs["MAL"] = -741.56e3
    Gibbs["G6P"] = -1567.08e3
    Gibbs["COAS"] = -57.17e3
    Gibbs["NADPH"] = -787.35e3
    Gibbs["SUCCoA"] = -471.06e3
    Gibbs["E4P"] = -1306.62e3
    Gibbs["PGLT"] = -1575.83e3
    Gibbs["PGN"] = -1782.55e3
    Gibbs["R5P"] = -1435.72e3
    Gibbs["RU5P"] = -1434.72e3
    Gibbs["S7P"] = -1685.66e3
    Gibbs["X5P"] = -1436.00e3

    # ## From LiWuQi11 Table 5. (Old version)
    # Gibbs["GLC"] = -917.44e3;
    # Gibbs["ATP"] = -2768.78e3;
    # Gibbs["ADP"] = -1905.45e3;
    # Gibbs["F6P"] = -1760.7995e3;
    # Gibbs["F16P"] = -2598.81e3;
    # Gibbs["DHAP"] = -1293.57e3;
    # Gibbs["GAP"] = -1286.54e3;
    # Gibbs["BPG"] = -2354.55e3;
    # Gibbs["PG3"] = -1501.32e3;
    # Gibbs["PG2"] = -1495.07e3;
    # Gibbs["PEP"] = -1262.42e3;
    # Gibbs["PYR"] = -465.99e3;
    # Gibbs["OAA"] = -791.77e3;
    # Gibbs["CIT"] = -1157.36e3;
    # Gibbs["ISCIT"] = -1151.6e3;
    # Gibbs["NADP"] = -836.21e3;
    # Gibbs["AKG"] = -786.12e3;
    # Gibbs["SUC"] = -685.82e3;
    # Gibbs["GDP"] = -1904.79e3;
    # Gibbs["FUM"] = -600.03e3;
    # Gibbs["MAL"] = -840.59e3;


    ## Synonyms
    Gibbs["FDP"] = Gibbs["F16P"]
    Gibbs["G3P"] = Gibbs["GAP"]
    Gibbs["GLCN"] = Gibbs["GLC"]
    Gibbs["RU5PD"] = Gibbs["RU5P"]
    Gibbs["XU5PD"] = Gibbs["X5P"]
    Gibbs["6PGL"] = Gibbs["PGLT"]
    Gibbs["6PGC"] = Gibbs["PGN"]
    
    
    ## Inorganic phosphate
    Gibbs["HPO4"] = Gibbs["PI"]		# Assume that they are the same thing
    
    ## H has zero Gibbs free energy
    Gibbs["Hx"] = 0
    Gibbs["Hi"] = 0
    Gibbs["H"] = 0

    ## Generic redox stuff
    redox=["E1","E2"]
    for sp in redox:
        Gibbs[sp] = 0

    return Gibbs

def conc_nom():
    """Concentrations at nominal conditions."""

    conc = {}
    ## Concentrations from ParRubXu16: nchembio.2077-S4
    conc["ADP"] = 5.69E-4		# ParRubXu16
    conc["ATP"] = 4.67E-3		# ParRubXu16
    conc["CO2c"] = 7.63E-3		# ParRubXu16
    conc["CO2m"] = 6.53E-3		# ParRubXu16
    conc["HPO4"] = 5.83E-3		# ParRubXu16
    conc["NAD"] = 5.02E-4		# ParRubXu16
    conc["NADH"] = 7.50E-5		# ParRubXu16
    conc["PYR"] = 5.88E-3              # ParRubXu16
    
    conc["CO2"] = conc["CO2m"]
    
    ## Concentrations from Mur09
    conc["O2_buffer"] = 200e-6	# Mur 09
    conc["O2_min   "] = 3e-6		# Mur 09
    conc["O2_max   "] = 30e-6		# Mur 09
    conc["O2"] = 25e-6		# Mur 09
    conc["sO2"] = 100e-12		# Mur 09

    ## Protons
    pHx = 7.78			# PorGheZan05
    pHi = 6.88			# PorGheZan05
    conc["Hx"] = 10**-pHx
    conc["Hi"] = 10**-pHi
    conc["H"] = 10**-7

    ## Concentrations from BazBeaVin16
    conc["Fe2"] = (2.7e-3)/2;
    conc["Fe3"] = (2.7e-3)/2;
    conc["QH2"] = (20e-3)/2;
    conc["Q"] = (20e-3)/2;
    conc["NADH"] = (3e-3)/2;
    conc["NAD"] = (3e-3)/2;

    ## Concentration of water!
    conc["H2O"] = 1
    
    return conc

def phi_Std():
    """Faraday-equivalent potential at standard conditions."""

    Gibbs = Gibbs_Std()         # Gibbs free energies
    const = thermoConstants()   # Useful constants
    F = const["F"]              # Faraday constant

    phi = {}
    for species, mu in Gibbs.items():
        phi[species] = mu/F

    return phi

def phi_nom(conc_nom):
    """Faraday-equivalent potential at nominal conditions."""

    phi = phi_Std()             # Faraday-equivalent potential
    const = thermoConstants()   # Useful constants
    V_N = const["V_N"]          # Faraday voltage
    
    for species, conc in conc_nom.items():
        if species in phi:
            phi[species] += V_N*np.log(conc)

    return phi

def phi_species(phi_nom,species):
    """Faraday-equivalent potentials in order of species"""

    phi = np.zeros((len(species),1))
    for i,sp in enumerate(species):
        if sp in phi_nom.keys():
            phi[i] = phi_nom[sp]
        else:
            print("phi not available for",sp+'.',"Set to zero")
            phi[i] = 0

    return phi

def known_species():
    """List the known species"""

    Gibbs = Gibbs_Std()
    species = list(Gibbs.keys())
    species.sort()
    return species

def ParRubXu16_conc():
    """Table 5 of ParRubXu16"""

    conc = {}
    conc['13DPG'] = 2.24e-06
    conc['BPG'] = 0.000237
    conc['23DHB'] = np.nan
    conc['2DHGULN'] = 2.77e-06
    conc['PG2'] = 9.49e-06
    conc['PG3'] = 0.000375
    conc['PSERL'] = 0.00044
    conc['4HBZ'] = np.nan
    conc['6PGC'] = 1.65e-05
    conc['AACOA'] = np.nan
    conc['ACCOA'] = 2.88e-05
    conc['ACTP'] = np.nan
    conc['ACONC'] = 1.1e-05
    conc['ACSER'] = 5.72e-05
    conc['ADE'] = np.nan
    conc['ADN'] = np.nan
    conc['APS'] = np.nan
    conc['ADP'] = 0.000569
    conc['ADPGLC'] = np.nan
    conc['AKG'] = 0.000797
    conc['ALAL'] = 0.00698
    conc['AMP'] = 4.23e-05
    conc['ANTH'] = np.nan
    conc['ARGL'] = 0.000255
    conc['ASNL'] = 0.000215
    conc['ASPL'] = 0.0149
    conc['ATP'] = 0.00467
    conc['CBASP'] = np.nan
    conc['CO2'] = 0.00763
    conc['CO2_M'] = 0.00653
    conc['CIT'] = 0.000584
    conc['CITRL'] = np.nan
    conc['CMP'] = 1.18e-05
    conc['COA'] = np.nan
    conc['COA_M'] = 0.00404
    conc['CTP'] = 0.000897
    conc['CAMP'] = 1.3e-07
    conc['CYSL'] = 8.4e-05
    conc['CYTD'] = np.nan
    conc['CSN'] = np.nan
    conc['DAMP'] = 1.68e-05
    conc['DATP'] = 9.74e-07
    conc['DCDP'] = 1.82e-06
    conc['DCMP'] = 3.71e-05
    conc['DCTP'] = np.nan
    conc['DAD2'] = np.nan
    conc['DGSN'] = np.nan
    conc['2DR5P'] = np.nan
    conc['DGMP'] = np.nan
    conc['DHORS'] = 0.000735
    conc['DHAP'] = 0.00163
    conc['DTDP'] = np.nan
    conc['DTMP'] = 1.18e-05
    conc['DTTP'] = np.nan
    conc['E4P'] = 1.03e-05
    conc['FAD'] = 5.6e-06
    conc['FMN'] = np.nan
    conc['FDP'] = 0.00152
    conc['F6P'] = 9.69e-05
    conc['FUM'] = np.nan
    conc['FUM_M'] = 0.000485
    conc['GDP'] = 3.02e-05
    conc['GLCN'] = 0.000211
    conc['ID66'] = np.nan
    conc['GAM6P'] = np.nan
    conc['G6P'] = 0.000675
    conc['GLUL'] = 0.0638
    conc['GLNL'] = 0.0172
    conc['GTHRD'] = 0.00309
    conc['GTHOX'] = 1.8e-05
    conc['G3P'] = 0.000141
    conc['GLYCR'] = np.nan
    conc['GLY'] = 0.00371
    conc['GMP'] = 1.81e-05
    conc['GTP'] = 0.000677
    conc['GUA'] = np.nan
    conc['GSN'] = 1.35e-06
    conc['ID80'] = 0.00107
    conc['HISL'] = 0.00041
    conc['HISTD'] = np.nan
    conc['HCYSL'] = np.nan
    conc['ID84'] = np.nan
    conc['IDP'] = np.nan
    conc['IMP'] = 1.23e-05
    conc['INS'] = 1.33e-06
    conc['ICIT'] = np.nan
    conc['ICIT_M'] = 3.21e-05
    conc['ILEL'] = 0.00176
    conc['ID91'] = 0.00352
    conc['ITP'] = np.nan
    conc['LEUL'] = 0.00176
    conc['LYSL'] = 0.000506
    conc['MALL'] = 0.00139
    conc['MALCOA'] = 4.95e-06
    conc['METL'] = 0.000639
    conc['ID98'] = 7.26e-05
    conc['INOST'] = np.nan
    conc['ACGAM1P'] = 7.47e-06
    conc['ID101'] = 6.24e-06
    conc['ID102'] = 1.03e-05
    conc['NACASP'] = 0.0029
    conc['ACORN'] = np.nan
    conc['NAD'] = 0.000502
    conc['NADH'] = 7.5e-05
    conc['NADP'] = 2.84e-05
    conc['NADPH'] = 6.54e-05
    conc['ORN'] = np.nan
    conc['OROT'] = 8.41e-06
    conc['OAA'] = np.nan
    conc['OAA_M'] = 2.01e-06
    conc['ID113'] = np.nan
    conc['PHEL'] = 0.00084
    conc['PHPYR'] = 0.00177
    conc['PI'] = 0.00583
    conc['PEP'] = 1.16e-05
    conc['PROL'] = 0.00123
    conc['PPCOA'] = np.nan
    conc['PRPP'] = np.nan
    conc['PYR'] = 0.00588
    conc['QULN'] = np.nan
    conc['RIBFLV'] = np.nan
    conc['R5P'] = 2.84e-05
    conc['RU5PD'] = 5.27e-06
    conc['AHCYS'] = 5.71e-07
    conc['AMET'] = np.nan
    conc['S7P'] = 1.81e-05
    conc['SERL'] = 0.00486
    conc['SKM'] = np.nan
    conc['GLYC3P'] = np.nan
    conc['SUCC'] = 0.000352
    conc['SUCCOA'] = np.nan
    conc['SUCCOA_M'] = 6.8e-06
    conc['TAUR'] = np.nan
    conc['THRL'] = 0.00669
    conc['THYMD'] = 2.64e-06
    conc['TRE'] = np.nan
    conc['TRPL'] = 0.00018
    conc['TYRL'] = 0.000938
    conc['UDP'] = 0.000133
    conc['UDPG'] = 0.00153
    conc['UDPGLCUR'] = 9.75e-05
    conc['UACGAM'] = 0.00897
    conc['UMP'] = 1.45e-05
    conc['URI'] = np.nan
    conc['UTP'] = 0.00176
    conc['VALL'] = 0.00151
    conc['XU5PD'] = 2.99e-05

    # ## Synonyms
    # conc["FDP"] = conc["F16P"]
    # conc["G3P"] = conc["GAP"]
    # conc["GLCN"] = conc["GLC"]
    # conc["RU5PD"] = conc["RU5P"]
    # conc["XU5PD"] = conc["X5P"]
    # conc["6PGL"] = conc["PGLT"]
    # conc["6PGC"] = conc["PGN"]
    

    return conc

    
def GarGri16_conc():
        """Table 18.2 of GarGri16
        """

        conc = {}
        conc['GLC'] = 5
        conc['GLCN'] = conc['GLC']
        conc['G6P'] = 0.083
        conc['F6P'] = 0.014
        conc['FDP'] = 0.031
        conc['DHAP'] = 0.14
        conc['G3P'] = 0.019
        conc['BPG'] = 0.001
        conc['PG3'] = 0.12
        conc['PG2'] = 0.030
        conc['PEP'] = 0.023
        conc['PYR'] = 0.051
        conc['LAC'] = 2.9

        conc['ATP'] = 1.85
        conc['ADP'] = 0.14
        conc['PI'] = 1.0


        for spec in conc.keys():
            conc[spec] /= 1000

        ## These are not in the table
        CONC = ParRubXu16_conc()
        conc['NAD'] = CONC['NAD']
        conc['NADH'] = CONC['NADH']

        return conc
    

    
def Phi_ParRubXu16():
    PHI = {}
    ##############################
    ### Mammalian iBMK (Mammalian) ###
    ##############################
    Phi = {}
    ## PGI (G6P => F6P)
    Phi['PGI'] = 0.00549306

    ## PFK (F6P + ATP => FBP + ADP + H)
    Phi['PFK'] = 0.13857

    ## FBA (FBP -=> DHAP + GAP)
    Phi['FBA'] = 0.050474

    ## TPI (DHAP => GAP)
    Phi['TPI'] = 0.00818777

    ## GAPD (GAP + NAD + Pi -> 13BPG + NADH + H)
    Phi['GAPD'] = 0.0210395

    ## PGK (13BPG + ADP => 3PG + ATP)
    Phi['PGK'] = 0.0228014

    ## PGM (3PG => 2PG)
    Phi['PGM'] = 0.0549306

    ## ENO (2PG => PEP + H2O)
    Phi['ENO'] = 0.0561743

    ## PYK (PEP + ADP + H => Pyr + ATP)
    Phi['PYK'] = 0.0389697

    ## PDH (Pyr + NAD + CoA => AcCoA + NADH + CO2)
    Phi['PDH'] = 0.293827

    ## GND (6PG + NADP => Ru5P + NADPH + CO2)
    Phi['GND'] = 0.101984

    ## RPE (Ru5P => Xu5P)
    Phi['RPE'] = 0.000207285

    ## RPI (Ru5P => R5P)
    Phi['RPI'] = 0.0197958

    ## TKT1 (Xu5P + R5P => GAP + S7P)
    Phi['TKT1'] = 0.000621856

    ## TALA (S7P + GAP => E4P + F6P)
    Phi['TALA'] = 0.0135772

    ## TKT2 (Xu5P + E4P => GAP + F6P)
    Phi['TKT2'] = 0.00404206

    ## CS (AcCoA + OAA + H2O => Cit + CoA + H)
    Phi['CS'] = 0.426904

    ## ACONT (Cit => Icit)
    Phi['ACONT'] = 0.000207285

    ## ICDYR (Icit + NAD => aKG + NADH + CO2)
    Phi['ICDYR'] = 0.0590763

    ## ICDHYR (Icit + NADP => aKG + NADPH + CO2)
    Phi['ICDHYR'] = np.nan

    ## AKGD (aKG + NAD + CoA => SuccCoA + NADH + CO2)
    Phi['AKGD'] = 0.353422

    ## SUCOAS (SuccCoA + ADP + Pi => Succ + CoA + ATP)
    Phi['SUCOAS'] = np.nan

    ## FUM (Fum + H2O => Mal)
    Phi['FUM'] = 0.00155464

    ## MDH (Mal + NAD => OAA + NADH + H)
    Phi['MDH'] = 0.00538942

    ## ME1 (Mal + NAD => Pyr + NADH + CO2)
    Phi['ME1'] = np.nan

    ## ? (Ser + thf => Gly + mlthf + H2O)
    Phi['?'] = 0.0178265

    PHI['Mammalian'] = Phi

    ##############################
    ### Yeast (Yeast) ###
    ##############################
    Phi = {}
    ## PGI (G6P => F6P)
    Phi['PGI'] = 0.00186557

    ## PFK (F6P + ATP => FBP + ADP + H)
    Phi['PFK'] = 0.220033

    ## FBA (FBP -=> DHAP + GAP)
    Phi['FBA'] = 0.00870599

    ## TPI (DHAP => GAP)
    Phi['TPI'] = 0.00673677

    ## GAPD (GAP + NAD + Pi -> 13BPG + NADH + H)
    Phi['GAPD'] = 0.0691297

    ## PGK (13BPG + ADP => 3PG + ATP)
    Phi['PGK'] = 0.0712025

    ## PGM (3PG => 2PG)
    Phi['PGM'] = 0.0469501

    ## ENO (2PG => PEP + H2O)
    Phi['ENO'] = 0.0484011

    ## PYK (PEP + ADP + H => Pyr + ATP)
    Phi['PYK'] = 0.099497

    ## PDH (Pyr + NAD + CoA => AcCoA + NADH + CO2)
    Phi['PDH'] = 0.590349

    ## GND (6PG + NADP => Ru5P + NADPH + CO2)
    Phi['GND'] = 0.113489

    ## RPE (Ru5P => Xu5P)
    Phi['RPE'] = 0.00207285

    ## RPI (Ru5P => R5P)
    Phi['RPI'] = 0.0189666

    ## TKT1 (Xu5P + R5P => GAP + S7P)
    Phi['TKT1'] = 0.0021765

    ## TALA (S7P + GAP => E4P + F6P)
    Phi['TALA'] = 0.0636366

    ## TKT2 (Xu5P + E4P => GAP + F6P)
    Phi['TKT2'] = 0.00943148

    ## CS (AcCoA + OAA + H2O => Cit + CoA + H)
    Phi['CS'] = 0.136083

    ## ACONT (Cit => Icit)
    Phi['ACONT'] = 0.0363786

    ## ICDYR (Icit + NAD => aKG + NADH + CO2)
    Phi['ICDYR'] = 0.0609419

    ## ICDHYR (Icit + NADP => aKG + NADPH + CO2)
    Phi['ICDHYR'] = np.nan

    ## AKGD (aKG + NAD + CoA => SuccCoA + NADH + CO2)
    Phi['AKGD'] = 0.623514

    ## SUCOAS (SuccCoA + ADP + Pi => Succ + CoA + ATP)
    Phi['SUCOAS'] = np.nan

    ## FUM (Fum + H2O => Mal)
    Phi['FUM'] = -9.32784e-06

    ## MDH (Mal + NAD => OAA + NADH + H)
    Phi['MDH'] = -0.0234232

    ## ME1 (Mal + NAD => Pyr + NADH + CO2)
    Phi['ME1'] = 0.0950403

    ## ? (Ser + thf => Gly + mlthf + H2O)
    Phi['?'] = 0.0342021

    PHI['Yeast'] = Phi

    ##############################
    ### E. coli (Ecoli) ###
    ##############################
    Phi = {}
    ## PGI (G6P => F6P)
    Phi['PGI'] = 0.0165828

    ## PFK (F6P + ATP => FBP + ADP + H)
    Phi['PFK'] = 0.256101

    ## FBA (FBP -=> DHAP + GAP)
    Phi['FBA'] = 0.0205213

    ## TPI (DHAP => GAP)
    Phi['TPI'] = 0.00818777

    ## GAPD (GAP + NAD + Pi -> 13BPG + NADH + H)
    Phi['GAPD'] = 0.0136808

    ## PGK (13BPG + ADP => 3PG + ATP)
    Phi['PGK'] = 0.0147173

    ## PGM (3PG => 2PG)
    Phi['PGM'] = 0.0328547

    ## ENO (2PG => PEP + H2O)
    Phi['ENO'] = 0.0285017

    ## PYK (PEP + ADP + H => Pyr + ATP)
    Phi['PYK'] = 0.0734827

    ## PDH (Pyr + NAD + CoA => AcCoA + NADH + CO2)
    Phi['PDH'] = 0.290614

    ## GND (6PG + NADP => Ru5P + NADPH + CO2)
    Phi['GND'] = 0.156293

    ## RPE (Ru5P => Xu5P)
    Phi['RPE'] = 0.000829142

    ## RPI (Ru5P => R5P)
    Phi['RPI'] = 4.14571e-05

    ## TKT1 (Xu5P + R5P => GAP + S7P)
    Phi['TKT1'] = 0.00414571

    ## TALA (S7P + GAP => E4P + F6P)
    Phi['TALA'] = 0.056278

    ## TKT2 (Xu5P + E4P => GAP + F6P)
    Phi['TKT2'] = 0.0166865

    ## CS (AcCoA + OAA + H2O => Cit + CoA + H)
    Phi['CS'] = 0.37954

    ## ACONT (Cit => Icit)
    Phi['ACONT'] = 0.0219722

    ## ICDYR (Icit + NAD => aKG + NADH + CO2)
    Phi['ICDYR'] = np.nan

    ## ICDHYR (Icit + NADP => aKG + NADPH + CO2)
    Phi['ICDHYR'] = 0.0615638

    ## AKGD (aKG + NAD + CoA => SuccCoA + NADH + CO2)
    Phi['AKGD'] = 0.122091

    ## SUCOAS (SuccCoA + ADP + Pi => Succ + CoA + ATP)
    Phi['SUCOAS'] = 0.222832

    ## FUM (Fum + H2O => Mal)
    Phi['FUM'] = 0.000103643

    ## MDH (Mal + NAD => OAA + NADH + H)
    Phi['MDH'] = 0.00538942

    ## ME1 (Mal + NAD => Pyr + NADH + CO2)
    Phi['ME1'] = 0.100948

    ## ? (Ser + thf => Gly + mlthf + H2O)
    Phi['?'] = 0.0899619

    PHI['Ecoli'] = Phi

    return PHI
    
