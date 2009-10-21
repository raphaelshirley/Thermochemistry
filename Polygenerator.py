# Program written by Raph to fit polynomials to Cp and other
# thermodynamic properties for use in chemkin and cantera
# from MySQL database. It will run for all species in database: thermonames
# it expects all values in KiloJoule/mol except As they are from the
# thermochemistry script

#Started 23 April 2008.
#Last modified 11 March 2009

############################################################################
##NOTE THIS REQUIRES THE VERSION OF NUMPY BEFORE THE CURRENT VERSION########
##################(AS AT APRIL 2008)########################################
############################################################################

# 15 May 2008 - enforced second derivative to be conserved accross Tswitch.
# I'd originally neglected this constraint required to make derivatives
# smooth as well as continuous - 15 May 2008

# 20 june 2008 - removed constraint on second derivative of Cp. I'm fairly
# sure it is unnecessary 20 june 2008

# 23 june 2008 - changed program to generally read in all species with
# a given project name - 23 june 2008

# 11 March 2009 - ras81 - changed units to SI. molar gas constant was in
# calories before. All Iactually did was change R to kJ/{mol.K} and
# multiplied S and Cp by 1000 to get it back to how it was- ras81 - 11 March 2009

# 20 JULY 2009 - ras81 - added database names as user variables to simplify changing them
# - ras81 20 JULY 2009

# 21 OCT 2009 - took sensitiv data out and put it in separate file that isn't made public - ras81

#############################################################################
##################################Input Section##############################
mmmyy = 'JUN08'                              
note = 'From Raphs program that reads sql database - 23jun2008 at 19:44'
project = 'aluminium'   # silver or silica 23rd june
thermonamedb = 'RAS81thermoname_B3LYP'
thermodatadb = 'RAS81thermodata_B3LYP'
frequencydb = 'RAS81frequency_B3LYP'
knownentdb = 'RAS81knownent'
comoentdb = 'RAS81comoent_B3LYP'
###############################END input Section############################
############################################################################

# constants
kb = 1.3806503E-23 #Boltzman's constant SI
hp = 6.626068E-34 #Planck's constant SI
c = 299792458 #Speed of light SI
#R = 1.987242 #Calories per mol per kelvin
R = 8.314472 #Joules per mol per kelvin
#R = 8.314472E-3 #KiloJoules per mol per kelvin



import MySQLdb
from math import *
from pylab import *
from scipy import *
from scipy.optimize import *
import openbabel, pybel
from numpy import *
import private  # Database info

#######################OPEN LINK TO SQL#####################################
conn = MySQLdb.connect (host = private.defaulthost,
                        user = private.defaultuser,
                        passwd = private.defaultpasswd,
                        db = private.defaultdb);
cursor = conn.cursor ()
#Get the list of species in all formats from chemkinnames on sql
sql1="select Name, ID, canname, project from %s" %(thermonamedb)  
cursor.execute(sql1)
specieslist = cursor.fetchall()
numspecies = len(specieslist)
allfinaldata = []
#############################END OPEN LINK TO SQL#########################

##############################LOOP OVER ALL SPECIES#######################
for i in range(numspecies):
    
    if not specieslist[i][3] == project:
        continue
    
    #print
    #print "species",i
    smile = specieslist[i][0]
    mol = pybel.readstring('smi', smile)
    #print smile
    #print mol.OBMol.GetSpacedFormula(1,'')
            
    #### this removes implicit hydrogens.
    #### ..which we don't want to do when reading in a SMILE
    ##    totalspin=0
    ##    for atom in mol.atoms:
    ##        valenceunsatisfied=atom.OBAtom.ImplicitHydrogenCount()
    ##        if valenceunsatisfied:
    ##            atom.OBAtom.SetSpinMultiplicity(valenceunsatisfied+1)
    ##            # multiplicity of 2 if 1 missing H atom
    ##            totalspin+=valenceunsatisfied
    ##    print "  sum of spins: %d" % totalspin
    ##    mol.OBMol.SetFormula( mol.OBMol.GetSpacedFormula(1,'') )

    formula = mol.OBMol.GetSpacedFormula()
    elements123 = formula.split( )
    numelements = len(elements123)/2

   

    elemental = [' ', ' ', ' ', ' ']
    El = [0.0,0.0,0.0,0.0]
    for run16 in range(0,numelements):
        El[run16] = int(elements123[2*run16 +1])
        elemental[run16] = elements123[2*run16]

    element1 = elemental[0]
    element2 = elemental[1]
    element3 = elemental[2]
    element4 = elemental[3]

    El1 = int(El[0])
    El2 = int(El[1])
    El3 = int(El[2])
    El4 = int(El[3])



    speciesname = specieslist[i][2]
    filename = specieslist[i][0]
    print speciesname, filename
    Hstandard = [0.0, 0.0] 
    sql2="select temp, entropy, heatcapacity, enthalpy, freeenergy from %s where ID='%s'" %(thermodatadb ,specieslist[i][1])
    cursor.execute(sql2)
    row2=cursor.fetchall()
    Ntemps = len(row2)
    #If it's an atom then only one temp and no need to do polyfitting
    if Ntemps == 1:
        allfinaldata = allfinaldata + [[filename, speciesname, element1, '1', element2, '0', element3, '0', element4, '0',
                       25, 100, 4000, row2[0][2]/R, 0, 0, 0, 0, row2[0][3]/R, row2[0][1]/R, row2[0][2]/R, 0, 0, 0, 0, row2[0][3]/R,
                                        row2[0][1]/R, 'This is an atom and therm does not change with time']]
        continue
    
    else:
        T = range(0,Ntemps)
        S = range(0,Ntemps)
        Cp = range(0,Ntemps)
        H = range(0,Ntemps)
        G = range(0,Ntemps)
        for run in range(0,Ntemps):
            T[run] = row2[run][0]
            S[run] = row2[run][1]*1000.0     # just multiplied by 1000 to get it back to before ras81-11 MAR 2009
            Cp[run] = row2[run][2]*1000.0    # same here thought this was easiest way to fix bug ras81-11 MAR 2009
            H[run] = row2[run][3]
            G[run] = row2[run][4]


    sql3="select Enthalpy from %s where id='%s'" %(knownentdb, specieslist[i][1]) 
    cursor.execute(sql3)
    enthalpy = cursor.fetchone()

    whodidthis = 'Enthalpy is from the literature'
    if project == 'silica':
        whodidthis = 'Enthalpy is from the literature (Ho and Melius)'
    if enthalpy == None:
        sql4="select Enthalpy from %s where id='%s'" %(comoentdb, specieslist[i][1]) 
        cursor.execute(sql4)
        enthalpy = cursor.fetchone()
        whodidthis = 'Enthalpy calculated by CoMo using isodesmic/isogyric reactions where posible'

    if enthalpy == None:
        enthalpy = [0.0]
        whodidthis = 'We have no enthalpy for this species'
            
    Hstandard[0] = enthalpy[0]
   
       
    #############################Define Cp function to be minimized##################
    #Define polynomial function for Cp

    def Polytest(Tswitch, a, b, x):
        if x < Tswitch:
            result = polyval(a, x)
        if x >= Tswitch:
            result = polyval(b, x)
        return result  

       
    #sum of squares
    def Minfunc(params):
        #define a1 and a2 so as to ensure continuity and smoothness
        Tswitch2 = abs(params[0]) +500# abs(params[0]) +500 or params[0]  or  1000 ####################################################################################################################
        a3 = params[1]  #removed may 15 2008 to enforce second derivative continuity
        a4 = params[2]
        a5 = params[3]
        b1 = params[4]
        b2 = params[5]
        b3 = params[6]
        b4 = params[7]
        b5 = params[8]
        #a3 = b3 + 3*b4*Tswitch2 + 6*b5*Tswitch2**2 - 3*a4*Tswitch2 + 6*a5*Tswitch2**2   # added 15th may to enforce second derivative continuity
        a2 = (b2 + 2*b3*Tswitch2 + 3*b4*Tswitch2**2 + 4*b5*Tswitch2**3 -
              2*a3*Tswitch2 - 3*a4*Tswitch2**2 - 4*a5*Tswitch2**3)        
        a1 = (b1 + b2*Tswitch2 + b3*Tswitch2**2 + b4*Tswitch2**3 + b5*Tswitch2**4 -
              a2*Tswitch2 - a3*Tswitch2**2 - a4*Tswitch2**3 - a5*Tswitch2**4)
        atest = [a5, a4, a3, a2, a1]
        btest = [b5, b4, b3, b2, b1]
        squares = 0.0
        for run6 in range(0,Ntemps):
            deltaCp = Cp[run6] - Polytest(Tswitch2, atest, btest, T[run6])
            squares = squares + deltaCp**2
        return squares

    ###########################End function minimized#############################

    #################################Produce Guess###############################
    #At the moment this splits at 1000 but adds variable split later
    #This section is to produce first guess which should be very good
    #define vector of values up to 1000K

    Ta = range(0,41)
    Sa = range(0,41)
    Cpa = range(0,41)
    Ha = range(0,41)
    Ga = range(0,41)
    for run4 in range(0,41):
        Ta[run4] = T[run4]
        Sa[run4] = S[run4]
        Cpa[run4] = Cp[run4]
        Ha[run4] = H[run4]
        Ga[run4] = G[run4]
    # the best fit line of Cp vs T of order N for the first half (up to 1000K)
    # must set constraint to make curves line up.
    N = 4
    coeffs1 = polyfit(Ta,Cpa,N)
    bestCpa = polyval(coeffs1, T)

    #define vector of values over 1000K
    #31 so overlaps
    Tb = range(0,(Ntemps - 31))
    Sb = range(0,(Ntemps - 31))
    Cpb = range(0,(Ntemps - 31))
    Hb = range(0,(Ntemps - 31))
    Gb = range(0,(Ntemps - 31))
    for run5 in range (31,Ntemps):   
        Tb[run5-31] = T[run5]
        Sb[run5-31] = S[run5]
        Cpb[run5-31] = Cp[run5]
        Hb[run5-31] = H[run5]
        Gb[run5-31] = G[run5]

    #Cpb[119] is last value

    # the best fit line of Cp vs T of order N for the second half (1000-4000K)
    #must set constraint to make curves line up.
    N = 4
    coeffs2 = polyfit(Tb,Cpb,N)
    bestCpb = polyval(coeffs2, T)

##    coeffs1[2]= (coeffs2[2] + 3*coeffs2[1]*1000 + 6*coeffs2[0]*1000             
##                        - 3*coeffs1[1]*1000 - 6*coeffs1[0]*1000)    # added 15may to enforce second deriv cont
##    coeffs1[3]=(coeffs2[3] + 2*coeffs2[2]*1000 +
##                3*coeffs2[1]*1000**2 + 4*coeffs2[0]*1000**3 -
##                2*coeffs1[2]*1000 - 3*coeffs1[1]*1000**2 -
##                4*coeffs1[0]*1000**3)
##    coeffs1[4]= (coeffs2[4] + coeffs2[3]*1000 + coeffs2[2]*1000**2 +
##                coeffs2[1]*1000**3 + coeffs2[0]*1000**4 -
##                coeffs1[3]*1000 - coeffs1[2]*1000**2 -
##                coeffs1[1]*1000**3 - coeffs1[0]*1000**4)
##    #initialCpa = polyval(coeffs1, T)
##    #initialCpb hasn't changes (same as bestCpb)
##    guess = [500,coeffs1[1],coeffs1[0],
##         coeffs2[4],coeffs2[3], coeffs2[2],coeffs2[1],coeffs2[0]] # changed 15 May2008
    #coeffs1[2]= (coeffs2[2] + 3*coeffs2[1]*1000 + 6*coeffs2[0]*1000             
    #                        - 3*coeffs1[1]*1000 - 6*coeffs1[0]*1000)    # added 15may to enforce second deriv cont
    coeffs1[3]=(coeffs2[3] + 2*coeffs2[2]*1000 +
                3*coeffs2[1]*1000**2 + 4*coeffs2[0]*1000**3 -
                2*coeffs1[2]*1000 - 3*coeffs1[1]*1000**2 -
                4*coeffs1[0]*1000**3)
    coeffs1[4]= (coeffs2[4] + coeffs2[3]*1000 + coeffs2[2]*1000**2 +
                 coeffs2[1]*1000**3 + coeffs2[0]*1000**4 -
                 coeffs1[3]*1000 - coeffs1[2]*1000**2 -
                 coeffs1[1]*1000**3 - coeffs1[0]*1000**4)
    #initialCpa = polyval(coeffs1, T)
    #initialCpb hasn't changes (same as bestCpb)

    ########################################################Guess two fitted polynomials that aren't smooth etc#############################################500/1000
    guess = [500, coeffs1[2],coeffs1[1],coeffs1[0],
             coeffs2[4],coeffs2[3], coeffs2[2],coeffs2[1],coeffs2[0]]
    #print 'this is the first sum of squares', Minfunc(guess)
    #########################END Produce guess####################################

    ##############################Minimise  Cp###################################
    #fmin(func, x0, args=(), xtol=0.0001, ftol=0.0001, maxiter=None, maxfun=None,
    #full_output=0, disp=1, retall=0, callback=None)
    final = fmin(Minfunc, guess, args=(), xtol=0.00000001, ftol=0.00000001, maxiter=None,
                    maxfun=None, full_output=0, disp=1, retall=0, callback=None)

    final[0] =abs(final[0]) + 500
##
##    aaa3 = (final[5] + 3*final[6]*final[0] + 6*final[7]*final[0]             
##                        - 3*final[1]*final[0] - 6*final[2]*final[0])  # added 15may to enforce second deriv cont
##    aaa2 =(final[4] + 2*final[5]*final[0] +
##                3*final[6]*final[0]**2 + 4*final[7]*final[0]**3 -
##                2*aaa3*final[0] - 3*final[1]*final[0]**2 -
##                4*final[2]*final[0]**3)
##    aaa1 = (final[3] + final[4]*final[0] + final[5]*final[0]**2 +
##                 final[6]*final[0]**3 + final[7]*final[0]**4 -
##                 aaa2*final[0] - aaa3*final[0]**2 -
##                 final[1]*final[0]**3 - final[2]*final[0]**4)
##    aaa = [final[2],final[1],aaa3,aaa2,aaa1]
##    bbb = [final[7],final[6],final[5],final[4],final[3]]
    #aaa3 = (final[5] + 3*final[6]*final[0] + 6*final[7]*final[0]             
    #                        - 3*final[1]*final[0] - 6*final[2]*final[0])  # added 15may to enforce second deriv cont
    aaa2 =(final[5] + 2*final[6]*final[0] +
                3*final[7]*final[0]**2 + 4*final[8]*final[0]**3 -
                2*final[1]*final[0] - 3*final[2]*final[0]**2 -
                4*final[3]*final[0]**3)
    aaa1 = (final[4] + final[5]*final[0] + final[6]*final[0]**2 +
                 final[7]*final[0]**3 + final[8]*final[0]**4 -
                 aaa2*final[0] - final[1]*final[0]**2 -
                 final[2]*final[0]**3 - final[3]*final[0]**4)
    aaa = [final[3],final[2],final[1],aaa2,aaa1]
    bbb = [final[8],final[7],final[6],final[5],final[4]]
    Cpfinal = range(0,Ntemps)
    for run7 in range(0,Ntemps):
        Cpfinal[run7] = Polytest(final[0],aaa,bbb,T[run7])

    ###############################END minimize Cp##############################

    ################################Find a6/b6##################################
    # This bit uses the input H at the start at a given temp (from isodesmic etc)
    # to find a6 and b6
    Hcoeffs1 =[aaa[0]/5,aaa[1]/4,aaa[2]/3,aaa[3]/2,aaa[4],0]
    Hcoeffs2 =[bbb[0]/5,bbb[1]/4,bbb[2]/3,bbb[3]/2,bbb[4],0]
    b6ini = polyval(Hcoeffs1, final[0]) - polyval(Hcoeffs2, final[0])
    Hcoeffs2 = [bbb[0]/5,bbb[1]/4,bbb[2]/3,bbb[3]/2,bbb[4],b6ini]
    ##UNITS BELOW COULD BE WRONG
    a6 = (Hstandard[0]-Polytest(final[0], Hcoeffs1, Hcoeffs2, Hstandard[1])/1000)*1000   # Think/Hope this expects H in kJ/mol ras81-11MAR2009
    if Hstandard[1] < final[0]:
        b6 = a6 + b6ini
    if Hstandard[1] >= final[0]:
        b6 = a6
        a6 = a6 - b6ini   #not much thought went into this
            
    Hcoeffs1 = [aaa[0]/5,aaa[1]/4,aaa[2]/3,aaa[3]/2,aaa[4],a6]
    Hcoeffs2 = [bbb[0]/5,bbb[1]/4,bbb[2]/3,bbb[3]/2,bbb[4],b6]



    Hbestfit = range(0,Ntemps)
    for run10 in range(0,Ntemps):
        Hbestfit[run10] =  Polytest(final[0], Hcoeffs1, Hcoeffs2, T[run10])/1000.0

    #################################END find a6/b6#############################

    ##################################Find a7###################################
    # This uses S to get a7 and make it continuous

    Scoeffs1 =[aaa[0]/4,aaa[1]/3,aaa[2]/2,aaa[3],aaa[4],0]
    Scoeffs2 =[bbb[0]/4,bbb[1]/3,bbb[2]/2,bbb[3],bbb[4],0]
    Spolycoeffs1 = [aaa[0]/4,aaa[1]/3,aaa[2]/2,aaa[3],0]
    Spolycoeffs2 = [bbb[0]/4,bbb[1]/3,bbb[2]/2,bbb[3],0]
    b7ini = ((polyval(Spolycoeffs1, final[0]) + aaa[4]*log(final[0])) -
                     (polyval(Spolycoeffs2, final[0]) + bbb[4]*log(final[0])))
    Scoeffs2 =[bbb[0]/4,bbb[1]/3,bbb[2]/2,bbb[3],bbb[4],b7ini]
    Spolycoeffs2 = [bbb[0]/4,bbb[1]/3,bbb[2]/2,bbb[3],b7ini]

    def Ent(temperature):
        if temperature < final[0]:
            result2 = (aaa[4]*log(temperature) + polyval(Spolycoeffs1,temperature))                  
        if temperature >= final[0]:
            result2 = (bbb[4]*log(temperature) + polyval(Spolycoeffs2,temperature))
        return result2


    def Minent(a7b7):
        Ssquares = 0
        for run8 in range(0,Ntemps):
            deltaS = S[run8] - (Ent(T[run8])+a7b7)
            Ssquares = Ssquares + deltaS**2
        return Ssquares

    a7b7final = fmin(Minent, 0, args=(), xtol=0.000001, ftol=0.000001, maxiter=None,
                            maxfun=None, full_output=0, disp=1, retall=0, callback=None)

    #print 'a7 and b7', a7b7final
    Scoeffs1 =[aaa[0]/4,aaa[1]/3,aaa[2]/2,aaa[3],aaa[4],a7b7final]
    Scoeffs2 =[bbb[0]/4,bbb[1]/3,bbb[2]/2,bbb[3],bbb[4],(a7b7final+b7ini)]
    Spolycoeffs1 = [aaa[0]/4,aaa[1]/3,aaa[2]/2,aaa[3],a7b7final]
    Spolycoeffs2 = [bbb[0]/4,bbb[1]/3,bbb[2]/2,bbb[3],(a7b7final+b7ini)]

    Sbestfit = range(0,Ntemps)
    for run9 in range(0,Ntemps):
        Sbestfit[run9] =  Ent(T[run9])

    ##################################END find a7/b7############################

    #############################FINAL COEFFICIENTS#############################
    a1 = aaa[4]/R
    a2 = aaa[3]/R
    a3 = aaa[2]/R
    a4 = aaa[1]/R
    a5 = aaa[0]/R
    a6 = a6/R
    a7 = a7b7final[0]/R
    b1 = bbb[4]/R
    b2 = bbb[3]/R
    b3 = bbb[2]/R
    b4 = bbb[1]/R
    b5 = bbb[0]/R
    b6 = b6/R
    b7 = (a7b7final[0]+b7ini)/R
    ###########################END FINAL COEFFICIENTS###########################

    ############################Print coefficients##############################
    print 'polynomial coefficients:'
    print '(a = below switchover, b = above switchover)'
    print 'Switchover T', final[0]
    print 'a1 =', "%1.9E"%(a1)
    print 'a2 =', "%1.9E"%(a2)
    print 'a3 =', "%1.9E"%(a3)
    print 'a4 =', "%1.9E"%(a4)
    print 'a5 =', "%1.9E"%(a5)
    print 'a6 =', "%1.9E"%(a6)
    print 'a7 =', "%1.9E"%(a7)
    print 'b1 =', "%1.9E"%(b1)
    print 'b2 =', "%1.9E"%(b2)
    print 'b3 =', "%1.9E"%(b3)
    print 'b4 =', "%1.9E"%(b4)
    print 'b5 =', "%1.9E"%(b5)
    print 'b6 =', "%+1.6E"%(b6)
    print 'b7 =', "%+1.6E"%(b7)

    numelement1 = str(El1)
    numelement2 = str(El2)
    numelement3 = str(El3)
    numelement4 = str(El4)

    #define element of allfinaldata containing all info needed for Cantera/ChemKin
                   #                  0             1         2          3           4            5         6          7           8           9           
    allfinaldata = allfinaldata + [[filename, speciesname, element1, numelement1, element2, numelement2, element3, numelement3, element4, numelement4,
                       min(T), final[0], T[Ntemps-1], a1, a2, a3, a4, a5, a6, a7, b1, b2, b3, b4, b5, b6, b7, whodidthis]]
                   #    10      11          12      13  14  15  16  17  18  19  20  21  22  23  24  25  26     27

###############################END LOOP OVER ALL SPECIES#####################################



###################################PRINT IN CANTERA FORM#####################################

print '#'
print '# Generated by raphs program that reads sql database'
print '# ', mmmyy
print '#'
print 'units(length = "cm", time = "s", quantity = "mol", act_energy = "cal/mol")'
print ''
print ''
print 'ideal_gas(name = "gas",'
print '      elements = " ',element1, element2, element3, element4, 'Ar",'
print '      species = """ '
numsilverspecies = len(allfinaldata)
for run44 in range(numsilverspecies):
    print allfinaldata[run44][1]
print '                   """,'
print '      reactions = "all",'
print '     initial_state = state(temperature = 300.0,'
print '                        pressure = OneAtm)    )'


print '#####################################################################################'
print '############################## CANTERA SPECIES ######################################'
print '#####################################################################################'
for run17 in range(len(allfinaldata)):
    filename = allfinaldata[run17][0]
    speciesname = allfinaldata[run17][1]
    element1 = allfinaldata[run17][2]
    numelement1 = allfinaldata[run17][3]
    element2 = allfinaldata[run17][4]
    numelement2 = allfinaldata[run17][5]
    element3 = allfinaldata[run17][6]
    numelement3 = allfinaldata[run17][7]
    element4 = allfinaldata[run17][8]
    numelement4 = allfinaldata[run17][9]
    T[0] = allfinaldata[run17][10]
    final[0] = allfinaldata[run17][11]
    T[Ntemps-1] = allfinaldata[run17][12]
    a1 = allfinaldata[run17][13]
    a2 = allfinaldata[run17][14]
    a3 = allfinaldata[run17][15]
    a4 = allfinaldata[run17][16]
    a5 = allfinaldata[run17][17]
    a6 = allfinaldata[run17][18]
    a7 = allfinaldata[run17][19]
    b1 = allfinaldata[run17][20]
    b2 = allfinaldata[run17][21]
    b3 = allfinaldata[run17][22]
    b4 = allfinaldata[run17][23]
    b5 = allfinaldata[run17][24]
    b6 = allfinaldata[run17][25]
    b7 = allfinaldata[run17][26]
    whodidthis = allfinaldata[run17][27]
    print '#'
    print '#Species number ', run17
    print '#', "%s"%(note)
    print '#',"%s %s"%(filename, whodidthis), ' From Raphs program'
    if speciesname == 'H':
        print '#HIGHLY DUBIOUS - BEWARE###############################################################################################LOOK HERE BEWARE'
    print 'species(name = "'+speciesname+'",'
    if numelement2 == '0':
        print '    atoms = " ',element1+':'+ numelement1,  '",'
    elif numelement3 == '0':
        print '    atoms = " ',element1+':'+ numelement1, element2+':'+ numelement2,  '",'
    elif numelement4 == '0':
        print '    atoms = " ',element1+':'+ numelement1, element2+':'+ numelement2,element3+':'+ numelement3,  '",'
    else:
        print '    atoms = " ',element1+':'+ numelement1, element2+':'+ numelement2,element3+':'+ numelement3,element4+':'+ numelement4,  '",'
    print '    thermo = ('
    print '       NASA( [  ',T[0],',  ',final[0],'], [',"%1.9E"%(a1),',  ',"%1.9E"%(a2),',' 
    print '              ',"%1.9E"%(a3),',  ',"%1.9E"%(a4),', ',"%1.9E"%(a5),','
    print '              ',"%1.9E"%(a6),', ',"%1.9E"%(a7),'] ),'
    print '       NASA( [ ',final[0],',  ',T[Ntemps-1],'], [ ',"%1.9E"%(b1),',  ',"%1.9E"%(b2),', '
    print '              ',"%1.9E"%(b3),',  ',"%1.9E"%(b4),', ',"%1.9E"%(b5),','
    print '              ',"%1.9E"%(b6),', ',"%1.9E"%(b7),'] )'
    print '             ),'
    print '    note = "',filename, mmmyy,'"'
    print '       )'

#######################################END PRINT IN CANTERA FORM######################################

##########################################PRINT IN CHEMKIN FORM########################################

print '!#####################################################################################'
print '!############################## CHEMKIN SPECIES ######################################'
print '!#####################################################################################'
for run17 in range(len(allfinaldata)):
    filename = allfinaldata[run17][0]
    speciesname = allfinaldata[run17][1]
    element1 = allfinaldata[run17][2]
    numelement1 = allfinaldata[run17][3]
    element2 = allfinaldata[run17][4]
    numelement2 = allfinaldata[run17][5]
    element3 = allfinaldata[run17][6]
    numelement3 = allfinaldata[run17][7]
    element4 = allfinaldata[run17][8]
    numelement4 = allfinaldata[run17][9]
    T[0] = allfinaldata[run17][10]
    final[0] = allfinaldata[run17][11]
    T[Ntemps-1] = allfinaldata[run17][12]
    a1 = allfinaldata[run17][13]
    a2 = allfinaldata[run17][14]
    a3 = allfinaldata[run17][15]
    a4 = allfinaldata[run17][16]
    a5 = allfinaldata[run17][17]
    a6 = allfinaldata[run17][18]
    a7 = allfinaldata[run17][19]
    b1 = allfinaldata[run17][20]
    b2 = allfinaldata[run17][21]
    b3 = allfinaldata[run17][22]
    b4 = allfinaldata[run17][23]
    b5 = allfinaldata[run17][24]
    b6 = allfinaldata[run17][25]
    b7 = allfinaldata[run17][26]
    whodidthis = allfinaldata[run17][27]
    
    if numelement2 == '0':
        chemkinlist = "%-2s%3s"%(element1,numelement1)
    elif numelement3 == '0':
        chemkinlist = "%-2s%3s%-2s%3s"%(element1,numelement1,element2,numelement2)
    elif numelement4 == '0':
        chemkinlist = "%-2s%3s%-2s%3s%-2s%3s"%(element1,numelement1,element2,numelement2,element3,numelement3)
    else:
        chemkinlist = "%-2s%3s%-2s%3s%-2s%3s%-2s%3s"%(element1,numelement1,element2,numelement2,element3,numelement3,element4,numelement4)
        
    #chemkinlist = "%-2s%3s%-2s%3s%-2s%3s%-2s%3s"%(element1,numelement1,element2,numelement2,element3,numelement3,element4,numelement4)
    print '!'
    print '!', "%s"%(note)
    print '!Species number ', run17
    print '!',"%s %s"%(filename, whodidthis), ' From Raphs program'
    if len(speciesname) >= 17:
        print '!NAME TOO LONG FOR CHEMKIN - CHANGE NAME TO LESS THAN 17 CHARS ##########################################LOOK HERE BEWARE'            
    if filename == 'H':
        print '!HIGHLY DUBIOUS - BEWARE##################################################################################LOOK HERE BEWARE'
    print "%-18s%-6s%-20sG%10g%10g%8g%5s 1"%(speciesname,mmmyy,chemkinlist,T[0],T[Ntemps-1],final[0],'')
    print "%+1.6E %+1.6E %+1.6E %+1.6E %+1.6E"%(b1,b2,b3,b4,b5),'    2'
    print "%+1.6E %+1.6E %+1.6E %+1.6E %+1.6E"%(b6,b7,a1,a2,a3),'    3'
    print "%+1.6E %+1.6E %+1.6E %+1.6E"%(a4,a5,a6,a7),'                   4'


##########################################END PRINT IN CHEMKIN FORM########################################




