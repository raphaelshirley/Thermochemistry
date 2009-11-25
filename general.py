# This contains all the required packages for the thermochemistry scripts as well
# as constants that are required. It is not currently used by any of the scripts

# Created by ras81 25NOV2009 


############### PACKAGES ###################
#IMPORT ALL THE PACKAGES USED IN ANY OF THE THERMO SCRIPTS
#(THESE MUST ALL BE INSTALLED)
import MySQLdb
import openbabel, pybel
from Cantera import *
from math import *
import random
from pylab import *
from scipy import *
from numpy import *
import private  # Database info
############### PACKAGES ###################

############### CONSTANTS ##################
#e.g. to get boltzmann's constant write constants.kb
class constants:
    """A simple class to specify some standard physical constants"""
    kb = 1.3806503E-23         #Boltzman's constant SI
    h  = 6.626068E-34          #Planck's constant   SI
    c  = 299792458             #Speed of light      SI
    Na = 6.0221415E23          #Avogadro's Number 
    amu = 1.66053886E-27       #mass of carbon12/12 SI  used to be called c12over12
    a0 = 5.29177E-11           #bohr radius         SI
    atm = 101325               #1atm = 101325Pa     SI
    hbar = 6.626068E-34 /(2*pi)#h over 2 pi         SI

constants = constants()   # initialize object (instance of the class constants -
                          # doesn't matter that the names are the same 
############### CONSTANTS ##################

################ UNIT CONVERTERS ############
#e.g. to get 4.0 cm^-1 in to Hz write unitconverter.cmtohz(4.0)
class unitconverter:
    """Class defining a number of functions to convert between units"""
    def cmtohz(self,cm):
        result = cm*100*constants.h*constants.c/constants.kb
        return result
    
unitconverter = unitconverter() # initialize object (instance of the class unitconverter -
                                # doesn't matter that the names are the same 

############### UNIT CONVERTERS ############
