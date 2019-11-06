## File ETC_abg.py

import BondGraphTools as bgt
import PII_abg
import Cyt_abg
import PI_abg
import Fer_abg

def model():
    """ 
    Model of chloroplast electron transport chain
    """

    ETC = bgt.new(name='ETC')   # Create system
    PII = PII_abg.model()
    Cyt = Cyt_abg.model()
    PI = PI_abg.model()
    Fer = Fer_abg.model()
    bgt.add(ETC,PII,Cyt,PI,Fer)
    
    return ETC
