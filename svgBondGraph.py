"""Convert bond graph diagram in SVG format to BondGraphTools format"""
import BondGraphTools as bgt
import numpy as np
import sympy as sp
import svgpathtools as svgp
from lxml import etree
import datetime

def indices(list,element):
    """ Find all indices of matching elements in list
    """

    indexList = [i for i, x in enumerate(list) if x == element]
    return indexList
    
def getBonds(svg,allowedColour=["#000000","#ff0000"]):
    """ Extract the bonds from the SVG file.
    """


    paths, attributes = svgp.svg2paths(svg)
    bondList = []
    for i_path,path in enumerate(paths):
        colour = attributes[i_path]["stroke"]
        isBond = (colour in allowedColour) and (len(path)>1)

        if isBond:
            segLen = []
            for i,segment in enumerate(path):
                #print(i,segment,segment.length())
                segLen.append(segment.length())
                #print(segLen,segLen[0],segLen[-1])
                if segLen[0]>segLen[-1]:
                    ## Tail at start
                    bondTail = path[0][0]
                    bondHead = path[-2][1]
                else:
                    bondTail = path[-1][1]
                    bondHead = path[1][0]
                    #print("Tail:",bondTail, "Head:", bondHead )
            bondList.append([bondTail,bondHead])
    return bondList

def getComponents(svg,allowedColour=["#000000","#ff0000"],convertR=False,convertCe=False):
    """ Extract BG components from the svg file
    """
    prefix = 'BGT'
    convertList = ['Re','re']
    if convertR:
        convertList.append('R')
        
    convertListCe = []
    if convertCe:
        convertListCe.append('Ce')

    compList = []
    compLocation = [];
    ReList = []
    CeList = []
    portList = []
    portLocation = [];
    junIndex = 0
    tree = etree.parse(svg)
    root = tree.getroot()

    for i,item in enumerate(root[-1]):
        tag = item.tag
        if isinstance(tag, str):

            if '}g' in tag:
                ## Rotated component
                ## So extract the translated coordinates
                # print("Found g",item.tag)
                # print(item.attrib)
                transform = item.attrib['transform']
                ## Extract x,y coords
                s = transform.split('(')
                s = s[1].split(',')
                x_trans = int(s[0])
                s = s[1].split(')')
                y_trans = int(s[0])

                Item = item.getchildren()[0]
                tag = Item.tag
            else:
                Item = item
                x_trans = 0
                y_trans = 0
                
            if '}text' in tag:
                comp = Item.text
                colour = Item.attrib["fill"]
                isComponent = ( ((":" in comp)
                                or (comp in ["0","1"] ))
                               and (colour in allowedColour )
                )
                isPort = (comp[0] in ['[']) and (comp[-1] in [']'])

                if isPort:
                    portList.append(comp[1:-1])
                    x = int(Item.attrib['x'])
                    y = int(Item.attrib['y'])
                    portLocation.append(complex(x+x_trans,y+y_trans))
                    
                if isComponent:
                    if ":" not in comp:
                        comp = comp+":"+prefix+str(junIndex)
                        junIndex += 1

                    compList.append(comp)
                    s = comp.split(':')
                    if s[0] in convertList:
                        ReList.append(s[1])
                    if s[0] in convertListCe:
                        CeList.append(s[1])

                    x = int(Item.attrib['x'])
                    y = int(Item.attrib['y'])
                    compLocation.append(complex(x+x_trans,y+y_trans))
                    
    return compList,compLocation,portList,portLocation,ReList,CeList


def makeConnections(bondList,
                    compList,compLocation,
                    portList,portLocation,
                    convertR=False,convertCe=False,quiet=False):
    """Create connections in BondGraphTools format
    """
    
    indent = "    "
    bondList = np.array(bondList)
    compLocation = np.array(compLocation)
    twoPortList = ['Re','re']
    if convertR:
        twoPortList.append('R')
    headList = []
    tailList = []
    tailPort = {}
    headPort = {}
    compPort = {}
    tails = bondList.T[0]
    heads = bondList.T[1]
    tailhead = np.hstack((tails,heads))
    n_bond = len(tails)
    for i,port in enumerate(portList):
        ## Locate nearest bond end to port
        distance = np.abs(tailhead-portLocation[i])
        nearest = np.argmin(distance)
        if nearest<n_bond:
            tailPort[nearest] = port
        else:
            headPort[nearest-n_bond] = port

        ## Locate nearest component to port
        distance = np.abs(compLocation-portLocation[i])
        nearest = np.argmin(distance)
        if port in compPort.keys():
            compPort[port].append(compList[nearest])
        else:
            compPort[port] = [compList[nearest]]
        
    # print("tailPort:", tailPort)
    # print("headPort:", headPort)

    for i_bond,bond in enumerate(bondList):
        comps = []; types = []; names = []
        for tailHead in bond:
            distance = np.abs(compLocation-tailHead)
            # print('tailHead:',tailHead)
            # print('compLocation:',compLocation)
            # print('distance:',distance)
            nearest = np.argmin(distance)
            # print(nearest) 
            comp = compList[nearest]
            s = comp.split(":")
            comps.append(comp)
            types.append(s[0])
            names.append(s[1])
        # print(comps)


        ## Tail
        if types[0] in twoPortList:
            tail = "("+names[0]+",1)"
        else:
            if i_bond in tailPort.keys():
                tail = "("+names[0]+",'"+tailPort[i_bond]+"')"
            else:
                tail = names[0]
            
        ##Head
        if types[1] in twoPortList:
            head = "("+names[1]+",0)"
        else:
            if i_bond in headPort.keys():
                head = "("+names[1]+",'"+headPort[i_bond]+"')"
            else:
                head = names[1]
 
        headList.append(head)
        tailList.append(tail)

    # print("compPort:",compPort)
    return headList,tailList,compPort

def makeBond(headList,tailList):
    """Create bonds in BondGraphTools format
    """
    indent = "    "
    strBond = indent+"## Bonds\n"
    bond = "{0}bgt.connect({1},{2})\n"
    for i,tail in enumerate(tailList):
        strBond += bond.format(indent,tail,headList[i])

    return strBond

def makeSubsystem(compType,compName,compPort):
    """ Create subsystem in BondGraphTools format
    """

    indent = "    "
    print("Creating subsystem:",compType+":"+compName)
    subStr = ("\n{0}## Subsystem {1}:{2}\n"
              "{0}import {1}_abg\n"
              "{0}{2} = {1}_abg.model()\n"
              "{0}{2}.name = '{2}'\n"
              "{0}model.add({2})\n"
         )
    compStr = "{0}:{1}"
    exposeStr = "{0}bgt.expose({1} / '{2}','{2}')\n"

    strSub = subStr.format(indent,compType,compName)

    for port,sysName in compPort.items():
        if compStr.format(compType,compName) in sysName:
            strSub += exposeStr.format(indent,compName,port)

    
    #print(strSub)
    return strSub

def makeComponents(compList,compPort,convertR=False,convertCe=False):
    """Create components in BondGraphTools format
    """
    indent = "    "
    strComp = ""

    SeStr = ("\n{0}## Component Se:{1}\n"
             "{0}{1} = bgt.new('Se',name='{1}')\n"
             "{0}model.add({1})\n"
    )

    SfStr = ("\n{0}## Component Sf:{1}\n"
             "{0}{1} = bgt.new('Sf',name='{1}')\n"
             "{0}model.add({1})\n"
    )
    
    CStr = ("\n{0}## Component C:{1}\n"
            "{0}{1} =  sp.symbols('{1}')\n"
            "{0}{1} = bgt.new('C',name='{1}',value={{'C':{1}}})\n"
            "{0}model.add({1})\n"
    )

    RStr = ("\n{0}## Component R:{1}\n"
            "{0}{1} =  sp.symbols('{1}')\n"
            "{0}{1} = bgt.new('R',name='{1}',value={{'r':{1}}})\n"
            "{0}model.add({1})\n"
    )

    RStrTwoPort = ("\n{0}## Component R:{1} (Two port version: convertR=True)\n"
                   "{0}kappa_{1} =  sp.symbols('kappa_{1}')\n"
                   "{0}RT = sp.symbols('RT')\n"
                   "{0}{1} = bgt.new('Re',name='{1}',value={{'r':kappa_{1},'R':RT,'T':1}},library='BioChem')\n"
                   "{0}model.add({1})\n"

    )

    
    CeStr = ("\n{0}## Component Ce:{1}\n"
             "{0}K_{1} =  sp.symbols('K_{1}')\n"
             "{0}RT = sp.symbols('RT')\n"
             "{0}{1} = bgt.new('Ce',name='{1}',value={{'k':K_{1},'R':RT,'T':1}},library='BioChem')\n"
             "{0}model.add({1})\n"
    )

    ## If Ce has an argument (eg C:A:2) create additional TF and 0 to give appropriate stoichiometry
    CeArgStr = ("\n{0}## Component Ce:{1}\n"
                "{0}K_{1} =  sp.symbols('K_{1}')\n"
                "{0}RT = sp.symbols('RT')\n"
                "{0}{1}_Ce = bgt.new('Ce',name='{1}',value={{'k':K_{1},'R':RT,'T':1}},library='BioChem')\n"
                "{0}model.add({1}_Ce)\n"
                "{0}{1}_TF = bgt.new('TF',value={2},name='{1}_TF')\n"
                "{0}model.add({1}_TF)\n"
                "{0}{1} = bgt.new('0')\n"
                "{0}model.add({1})\n"
                "{0}bgt.connect(({1}_TF,0),{1}_Ce)\n"
                "{0}bgt.connect({1},({1}_TF,1))\n"
    )



    ReStr = ("\n{0}## Component Re:{1}\n"
             "{0}kappa_{1} =  sp.symbols('kappa_{1}')\n"
             "{0}RT = sp.symbols('RT')\n"
             "{0}{1} = bgt.new('Re',name='{1}',value={{'r':kappa_{1},'R':RT,'T':1}},library='BioChem')\n"
             "{0}model.add({1})\n"
    )

    junStr = ("\n{0}## Junction {1}:{2}\n"
              "{0}{2} = bgt.new('{1}')\n"
              "{0}model.add({2})\n"

    )

    compNames = []
    for comp in compList:
        s = comp.split(":")
        compType = s[0]
        compName = s[1]
        if len(s)>2:
            compArg = s[2]
        else:
            compArg = None
        #print(comp,compArg)

        compNames.append(compName)
        if compType in ["C"]:
            strComp = strComp+CStr.format(indent,compName)
        elif (compType in ["R"]):
            if convertR:
                strComp = strComp+RStrTwoPort.format(indent,compName)
            else:
                strComp = strComp+RStr.format(indent,compName)
        elif compType in ["Ce"]:
            if compArg is None:
                strComp = strComp+CeStr.format(indent,compName)
            else:
                strComp = strComp+CeArgStr.format(indent,compName,compArg)
        elif compType in ["Re"]:
            strComp = strComp+ReStr.format(indent,compName)
        elif compType in ["0","1"]:
            strComp = strComp+junStr.format(indent,compType,compName)
        elif compType in ["Se"]:
            strComp = strComp+SeStr.format(indent,compName)
        elif compType in ["Sf"]:
            strComp = strComp+SfStr.format(indent,compName)
        else:
            #print("UNKNOWN:",compType+":"+compName)
            strComp += makeSubsystem(compType,compName,compPort)
            
    # strComp = (strComp+"\n"+indent+"## Component list\n"
    #            +indent+"components = (\n"
    # )

    # for compName in compNames:
    #     strComp = strComp+"      "+compName+",\n"

    # strComp = (strComp[:-2]+indent+"\n"+indent+")\n\n"
    #                +indent+"bgt.add(model, *components)\n"
    #)

    return strComp

def convertRe(ReList,headList,tailList,quiet=False):
    """ Convert one-port Re to two-port 
    """
    compList = []              # Extra junctions
    for re in ReList:
        port0 = "({0},0)".format(re)
        port1 = "({0},1)".format(re)
        if port1 not in tailList:
            if not quiet:
                print("Converting one-port "+re, "to two-port")
                
            ## Find junction at tail of bond
            iRe = headList.index(port0)
            jun = tailList[iRe]

            ## New reverse junction
            junR = jun+"r"
            compList.append("1:"+junR)

            ## Connect the tails of  bonds to new junction
            for i,tail in enumerate(tailList):
                if (tail in [jun]) and (headList[i] not in [port0]):
                    #print("    Swap", jun, "for", junR, "("+headList[i]+")")
                    tailList[i] = junR
                    
            ## Connect the Re to the new junction
            tailList.append(port1)
            headList.append(junR)

    return headList,tailList,compList

def ConvertCe(CeList,headList,tailList,quiet=False):
    """ Convert Ce to zero junction and Ce
    """
    compList = []              # Extra junctions
    for Ce in CeList:
        if not quiet:
            print("Appending zero junction to Ce:"+Ce)
            
        junCe = Ce+"z"
        compList.append("0:"+junCe)
        if Ce in headList:
            iCe = indices(headList,Ce)
            #print('iCe:',iCe)
            for i in iCe:
                headList[i] = junCe
            headList.append(Ce)
            tailList.append(junCe)
            
        if Ce in tailList:
            iCe = indices(tailList,Ce)
            #print('iCe:',iCe)
            for i in iCe:
                tailList[i] = junCe
            tailList.append(Ce)
            headList.append(junCe)
             
    return headList,tailList,compList

def model(svg,convertR=False,convertCe=False,quiet=False):
    """ Converts the SVG graphical BG to the BondGraphTools computational BG
    """
    
    indent = "    "
    ## File handling
    s = svg.split('_abg')
    name = s[0]
    filename = name+"_abg.py"
    f = open(filename,'w')

    header = ("import BondGraphTools as bgt\n"
              "import sympy as sp\n\n"
              "def model():\n"
              "{0}{3} Acausal bond graph {1}_abg.py\n"
              "{0}Created by svgBondGraph at {2} from {1}_abg.svg\n\n"
              "{0}Usage:\n"
              "{0}import {1}_abg; model = {1}_abg.model()\n"
              "{0}{3}\n\n"
              '{0}model = bgt.new(name="{1}")\n'
     )
    
    ## get Bonds and the location of head and tail
    bondList = getBonds(svg)

    ## Get components and their location
    compList,compLocation,portList,portLocation,ReList,CeList = getComponents(svg,convertR=convertR,convertCe=convertCe)
        
    # ## Create the components - return as string
    # strComp = makeComponents(compList,convertR=convertR)

    ## Create the bond connections
    if len(bondList)>0:
        headList,tailList,compPort = makeConnections(bondList,
                                                     compList,compLocation,
                                                     portList,portLocation,
                                                     convertR=convertR,convertCe=convertCe)
    else:
        headList=[]
        tailList=[]
        compPort={}

    ## Convert one port Re to two port
    headList,tailList,compListNew = convertRe(ReList,headList,tailList,quiet=quiet)
    compList.extend(compListNew)
    
    ## Convert one port Ce to zero + Ce
    headList,tailList,compListNew = ConvertCe(CeList,headList,tailList,quiet=quiet)
    compList.extend(compListNew)
    
    ## Make bond connection str
    strBond = makeBond(headList,tailList) 
    # print(compList)
    # print(sorted(compList))
    
    ## Create the components - return as string
    strComp = makeComponents(sorted(compList),compPort,convertR=convertR,convertCe=convertCe)
    
    ## Write out the _abg.py file
    ## Header
    f.write(header.format(indent,name,datetime.datetime.now().ctime(),'"""'))
    
    ## Body
    f.write(strComp+"\n"+strBond+"\n")
    
    ## End
    f.write(indent+"return model\n")
    
    f.close
