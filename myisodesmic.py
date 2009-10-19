#################################################################################################
####                                                                                         ####
####                   ISODESMIC REACTION FINDER AND ENTHALPY CALCULATION                    ####
####                                                                                         ####
#################################################################################################

# An isodesmic reaction finder
# designed to find the smallest possible isodesmic reaction that links the given species, 
# or as close as possible
# Richard West - April 2008 
# R.H.West.00@cantab.net or rwest@alum.mit.edu


# # for parsing gaussian log files
# import cclib
# # http://cclib.sourceforge.net/wiki/index.php/Using_cclib


# PuLP : A Linear Programming modeler in Python
# http://www.jeannot.org/~js/code/index.en.html#PuLP
# Modified by Markus Sander
# Modified by Tim Totton (13/3/2009) to give correct reading of dH(0->298K) from thermodata
# database and to write the known enthalpyies (converted to 0K) to the comoent database

#################################################################################################
# Import libraries
from pulp import *

from math import *
import re
import glob
import os
import MySQLdb;
#################################################################################################
# MySQL tables specification

thermonameSQL = 'blank'
thermodataSQL = 'blank'
frequencySQL = 'blank'
knownentSQL = 'blank'
comoentSQL = 'blank'

#################################################################################################
# Local directories

## change these strings to indcate directories which contain known and unknown species mol files

known = 'enthalpyofformation/project/known/'
unknown = 'enthalpyofformation/project/unknown/'

#################################################################################################
# Use SQL database have usesql = 1 other wise usesql = 0

usesql = 1

#################################################################################################






#################################################################################################
## Defining SQL database use
#################################################################################################

if usesql:    
    conn = MySQLdb.connect (host = "127.0.0.1",
                           user = "prime",
                           passwd = "m4sterpl4n",
                           db = "prime");
    cursor = conn.cursor ()

listoffilesknown=glob.glob(known+'*.mol')

listoffilesunknown=glob.glob(unknown+'*.mol')
#listoffiles.extend(glob.glob(unknown+'*.mol'))
#listoffiles.extend(glob.glob('TEOS_geometries/unknownmp2/*.mol'))
#listoffiles.extend(glob.glob('TEOS_geometries/knownmp2/*.mol'))

print '-------------------------------------------------------------------------------------'
print
print 'ISODESMIC ONE REACTION OUTPUT'
print 'Written by R. West, M. Sander and T. Totton'
print
# print 'List of files'
print
# print listoffiles
print
print '-------------------------------------------------------------------------------------'
print

periodicTable={ 1: 'H',  2: 'He',  3: 'Li',  4: 'Be',  5: 'B',  6: 'C',  7: 'N',  8: 'O',  9: 'F',  10: 'Ne',  11: 'Na',  12: 'Mg',  13: 'Al',  14: 'Si',  15: 'P',  16: 'S',  17: 'Cl',  18: 'Ar',  19: 'K',  20: 'Ca',  21: 'Sc',  22: 'Ti',  23: 'V',  24: 'Cr',  25: 'Mn',  26: 'Fe',  27: 'Co',  28: 'Ni',  29: 'Cu',  30: 'Zn',  31: 'Ga',  32: 'Ge',  33: 'As',  34: 'Se',  35: 'Br',  36: 'Kr',  37: 'Rb',  38: 'Sr',  39: 'Y',  40: 'Zr',  41: 'Nb',  42: 'Mo',  43: 'Tc',  44: 'Ru',  45: 'Rh',  46: 'Pd',  47: 'Ag',  48: 'Cd',  49: 'In',  50: 'Sn',  51: 'Sb',  52: 'Te',  53: 'I',  54: 'Xe',  55: 'Cs',  56: 'Ba',  57: 'La',  58: 'Ce',  59: 'Pr',  60: 'Nd',  61: 'Pm',  62: 'Sm',  63: 'Eu',  64: 'Gd',  65: 'Tb',  66: 'Dy',  67: 'Ho',  68: 'Er',  69: 'Tm',  70: 'Yb',  71: 'Lu',  72: 'Hf',  73: 'Ta',  74: 'W',  75: 'Re',  76: 'Os',  77: 'Ir',  78: 'Pt',  79: 'Au',  80: 'Hg',  81: 'Tl',  82: 'Pb',  83: 'Bi',  84: 'Po',  85: 'At',  86: 'Rn',  87: 'Fr',  88: 'Ra',  89: 'Ac',  90: 'Th',  91: 'Pa',  92: 'U',  93: 'Np',  94: 'Pu',  95: 'Am',  96: 'Cm',  97: 'Bk',  98: 'Cf',  99: 'Es',  100: 'Fm',  101: 'Md',  102: 'No',  103: 'Lr',  104: 'Rf',  105: 'Db',  106: 'Sg',  107: 'Bh',  108: 'Hs',  109: 'Mt',  110: 'Ds',  111: 'Rg',  112: 'Uub',  113: 'Uut',  114: 'Uuq',  115: 'Uup',  116: 'Uuh',  117: 'Uus',  118: 'Uuo'}

import openbabel, pybel
# please cite:
# Pybel: a Python wrapper for the OpenBabel cheminformatics toolkit
# Noel M O'Boyle, Chris Morley and Geoffrey R Hutchison
# Chemistry Central Journal 2008, 2:5
# doi:10.1186/1752-153X-2-5



for unknownspecies in listoffilesunknown:

    speciestable={}
    gaussiantable={}
    enthalpy={}
    energy={}
    enthalpyroom={}
    molpath='SDFfiles'
    os.path.isdir(molpath) or os.mkdir(molpath)
    goodpath=os.path.join(molpath,'correct')
    os.path.isdir(goodpath) or os.mkdir(goodpath)
    badpath=os.path.join(molpath,'dodgy')
    os.path.isdir(badpath) or os.mkdir(badpath)
    okpath=os.path.join(molpath,'unchecked')
    os.path.isdir(okpath) or os.mkdir(okpath)
    listoffiles=[]

    for item in listoffilesknown:
            listoffiles.append(item)
    listoffiles.extend(glob.glob(unknownspecies))

    for filename in listoffiles:
        basename=os.path.basename(filename)
        sdffilename=os.path.join(goodpath, basename+'.sdf')
        if os.path.isfile(sdffilename):
            mol = pybel.readfile('sdf', sdffilename).next()
            speciestable[mol.OBMol.GetFormula()]=mol
            listoffiles.remove(filename)

    print 'Reading in files: Sum of spins'
    print

    for filename in listoffiles:
        basename=os.path.basename(filename)
        # read in the first molecule in a mol file
        mol = pybel.readfile("mol", filename).next()
        print os.path.basename(filename),

        totalspin=0
        for atom in mol.atoms:
            valenceunsatisfied=atom.OBAtom.ImplicitHydrogenCount()
            if valenceunsatisfied:
                atom.OBAtom.SetSpinMultiplicity(valenceunsatisfied+1)
                # multiplicity of 2 if 1 missing H atom
                totalspin+=valenceunsatisfied

        print "  sum of spins: %d" % totalspin

        #reset formula
        mol.OBMol.SetFormula( mol.OBMol.GetSpacedFormula(1,'') )
        
        smile=mol.write('smi').strip()    
        molname=basename.replace('.mol','')
        gaussianname=basename.replace('.mol','.g03%')
        mol.OBMol.SetTitle(molname)
        
    #    if mol.OBMol.GetFormula().find('H')>=0 :
    #        print "Strange Hydrogen stuff"
    #        break 

    #    if filename.find('unknown')>=0:
        if 'unknown' in filename:
            mol.known=0
        else:
            mol.known=1   

        if totalspin and mol.OBMol.GetTotalSpinMultiplicity() != totalspin+1:
            print "Sum of atomic spin multiplicities (%d) != molecule spin multiplicity (%d)" % (totalspin+1,mol.OBMol.GetTotalSpinMultiplicity())
            print  mol.write('molreport')
            mol.write('sdf',filename=os.path.join(badpath,basename+'.sdf'),overwrite=True)
            speciestable[molname]=mol ####################????????????????????
        else:
            #store molecule
            speciestable[molname]=mol  # or smile or mol.OBMol.GetFormula()
            gaussiantable[molname]=gaussianname
            mol.write('sdf',filename=os.path.join(okpath,basename+'.sdf'),overwrite=True)

                
                    

    #    if len(speciestable)>20: break # stop prematurely for debugging
       # if os.path.basename(filename)=='ti2o3cl2-b971-6311+gdp.freq.log': break
        
        # tio2cl2trig2-ub3lyp-6311+gdp.log
     
    # # read in logfiles using cclib   
    # for filename in listoffiles:
    #     myfile=cclib.parser.ccopen(filename)
    #     import logging
    #     myfile.logger.setLevel(logging.ERROR) # reduce amount of logging messages
    #     data=myfile.parse()
    #         
    #     # try to calculate bond orders
    #     analysis=cclib.method.MBO(myfile)
    #     if analysis.calculate():
    #         print "HURRAH!!!"
    #    # del data
    #    # del myfile

    print
    print "Have finished reading in log files. Setting up optmization"
    print
    print '-------------------------------------------------------------------------------------'
    print


    # print stuff out
    bondtable={}
    elementtable={}
    for smile in speciestable:
        mol=speciestable[smile]
        for atom in mol.atoms:
            # print periodicTable[atom.OBAtom.GetAtomicNum()]
            elementtable[ periodicTable[atom.OBAtom.GetAtomicNum()] ]=1

        if mol.OBMol.NumBonds():
            for bond in openbabel.OBMolBondIter(mol.OBMol):
                # standardbond=sorted(( bond.GetBeginAtom().GetType(),bond.GetEndAtom().GetType() ))
                standardbond=sorted(( periodicTable[bond.GetBeginAtom().GetAtomicNum()],periodicTable[bond.GetEndAtom().GetAtomicNum() ] ))
                bondname= "%s -%d- %s"% (standardbond[0], bond.GetBondOrder(), standardbond[1])
                if not bondtable.has_key(bondname): bondtable[bondname]={}
                if not bondtable[bondname].has_key(smile): bondtable[bondname][smile]=0
                bondtable[bondname][smile]+=1
    bondnames=bondtable.keys()         
    elements=elementtable.keys()
        
    species = speciestable.keys()        
                
    # species=['Ti2O2Cl4','TiOCl2','TiCl4','TiO2', 'OClO'] 
    #,'TiOCl', 'TiOCl2', 'TiOCl3', 'TiO2Cl2', 'TiO2Cl3', 'Ti2O2Cl3', 'Ti2O2Cl4','Ti2O3Cl2', 'Ti2O3Cl3', 'Ti3O4Cl4', 'Ti5O6Cl8']
    spi=range(len(species))
    #elements=['Ti','O','Cl']
    eli=range(len(elements))

    composition=[[0 for s in spi] for e in eli]

    for e in eli:
        for s in spi:
            count=0
            mymatch=re.compile(elements[e]+'(\d*)')
            for m in mymatch.findall(speciestable[species[s]].OBMol.GetFormula()): #species[s]
                count+= int(m or 1)
            print "No. of %s in %s is %d" % (elements[e], species[s], count)
            composition[e][s]=count 



    #################################################################################################
    ## Finding Isodesmic reactions
    #################################################################################################


    print
    print '-------------------------------------------------------------------------------------'
    print
    print 'Finding Isodesmic Reactions'
    print


    #################################################################################################
    ## Defining getinvspecies
    #################################################################################################

        
    #for target in spi:
    def getinvspecies(target):
        # don't bother with known targets
    #    if speciestable[species[target]].known:
    #        continue
        prob = LpProblem("Isodesmic Reaction Search", LpMinimize)        
        boi=range(len(bondnames))
        bonds=[[0 for s in spi] for b in boi]

        for b in boi:
            for s in spi:
                bondname=bondnames[b]
                speciesname=species[s]
                if bondtable[bondname].has_key(speciesname):
                    bonds[b][s]=bondtable[bondname][speciesname]
                else:
                    bonds[b][s]=0
    #    e= elements.index('O')
    #    b= bondnames.index('Ti=O')
    #    bonds[b][s]= composition[e][s]  # assume all O atoms are double bonded to Ti
    #    e=elements.index('Cl')
    #    b=bondnames.index('Ti-Cl')
    #    bonds[b][s]= composition[e][s]  # assume all Cl atoms are bonded to Ti

    # A vector of integer variables - the stoichiometries of each species
        stoich = LpVariable.matrix("stoich", spi, None, None) #, LpInteger) 

    # A vector of positive integer variables - the absolute (positive) value of the stoichiometries 
        modstoich = LpVariable.matrix("modstoich", spi, 0, None) #, LpInteger) 
        for i in range(len(stoich)):
            prob+= modstoich[i]>=stoich[i]
            prob+= modstoich[i]>= -1*stoich[i]


    # A vector of positive integer variables - the absolute (positive) value of the net bonds made 
        modbondsmade = LpVariable.matrix("modbondsmade", boi, 0, None) #, LpInteger) 
        for b in boi:
            prob+= modbondsmade[b]>= lpDot( bonds[b], stoich)
            prob+= modbondsmade[b]>= -1 * lpDot( bonds[b], stoich)


    # Constraints:
    # elemental creation must be 0 for each element
        creation=[0 for e in eli]
        for e in eli:
            creation[e] = lpDot( composition[e] , stoich)
            prob += creation[e] == 0

    # objective function
    #prob += lpSum(modstoich)+stoich[target]
        prob += lpSum(modbondsmade) + 0.001*lpSum([lpDot( bonds[b], modstoich) for b in boi])
    #prob += lpSum( [lpDot( bonds[b], stoich) for b in boi ] )



        # remove all old constraints
        for s in spi: 
            stoich[s].lowBound=None
            stoich[s].upBound=None
        # don't use unknown species
        for s in spi:
            if speciestable[species[s]].known==0:
                    stoich[s].upBound = 0
                    stoich[s].lowBound = 0
        # reaction must make some of target!
        stoich[target].upBound = None
        stoich[target].lowBound = 1 


        # solve the problem
        prob.solve(COIN(msg=0)) # msg = 0
        if prob.status<0: prob.solve() # with messages
        involvedspecies={}
        for s in spi:
            if value(stoich[s])!=0:
                    involvedspecies[species[s]]=1
    #    print species[target]
    #    for s in involvedspecies:
    #          print s
        return involvedspecies


    #################################################################################################
    ## Defining Isodesmic Enthalpy Calculation
    #################################################################################################



    def getenthalpy(target,notuse,bondnames):
    #bondnames=['Ti=O','Ti-Cl','O-Cl']
        prob = LpProblem("Isodesmic Reaction Search", LpMinimize)
        boi=range(len(bondnames))
        bonds=[[0 for s in spi] for b in boi]
            
        for b in boi:
            for s in spi:
                bondname=bondnames[b]
                speciesname=species[s]
                if bondtable[bondname].has_key(speciesname):
                    bonds[b][s]=bondtable[bondname][speciesname]
                else:
                    bonds[b][s]=0
    #    e= elements.index('O')
    #    b= bondnames.index('Ti=O')
    #    bonds[b][s]= composition[e][s]  # assume all O atoms are double bonded to Ti
    #    e=elements.index('Cl')
    #    b=bondnames.index('Ti-Cl')
    #    bonds[b][s]= composition[e][s]  # assume all Cl atoms are bonded to Ti

    # A vector of integer variables - the stoichiometries of each species
        stoich = LpVariable.matrix("stoich", spi, None, None) #, LpInteger) 

    # A vector of positive integer variables - the absolute (positive) value of the stoichiometries 
        modstoich = LpVariable.matrix("modstoich", spi, 0, None) #, LpInteger) 
        for i in range(len(stoich)):
            prob+= modstoich[i]>=stoich[i]
            prob+= modstoich[i]>= -1*stoich[i]


    # A vector of positive integer variables - the absolute (positive) value of the net bonds made 
        modbondsmade = LpVariable.matrix("modbondsmade", boi, 0, None) #, LpInteger) 
        for b in boi:
            prob+= modbondsmade[b]>= lpDot( bonds[b], stoich)
            prob+= modbondsmade[b]>= -1 * lpDot( bonds[b], stoich)


    # Constraints:
    # elemental creation must be 0 for each element
        creation=[0 for e in eli]
        for e in eli:
            creation[e] = lpDot( composition[e] , stoich)
            prob += creation[e] == 0

    # objective function
    #prob += lpSum(modstoich)+stoich[target]
        prob += lpSum(modbondsmade) + 0.001*lpSum([lpDot( bonds[b], modstoich) for b in boi])
    #prob += lpSum( [lpDot( bonds[b], stoich) for b in boi ] )


    #target species index
    #    target=1

    # print prob

     #   print notuse
        # remove all old constraints
        for s in spi:
            stoich[s].lowBound=None
            stoich[s].upBound=None
        # don't use unknown species
        for s in spi:
            if speciestable[species[s]].known==0 or species[s] in notuse:
    #           if speciestable[species[s]].known==0:
                stoich[s].upBound = 0
                stoich[s].lowBound = 0
        # reaction must make some of target!
        stoich[target].upBound = None
        stoich[target].lowBound = 1


        # solve the problem
        prob.solve(COIN(msg=0)) # msg = 0
        if prob.status<0: prob.solve() # with messages
                     


        print 
        print "To find %s use: \n"%species[target],
    #    print speciestable[species[1]].write('smi').split('\t')[0]
        for s in spi:
            if usesql: 
                if speciestable[species[s]].known==1:
                    sql="select Enthalpy from %s where name='%s'" %(knownentSQL,speciestable[species[s]].write('smi').split('\t')[0])
                    cursor.execute(sql)
                    if cursor.rowcount==0:
                        print "No entry for %s in %s found" %(speciestable[species[s]].write('smi'),knownentSQL)
            
                    row=cursor.fetchone()
                    enthalpy[species[s]]=row[0]
            else:
                enthalpy[species[s]]=0
            if usesql:   
                sql="select Energy,ID from %s where molefile like '%s'" %(thermonameSQL,gaussiantable[species[s]])
                cursor.execute(sql)
                row=cursor.fetchone()
                if cursor.rowcount==0:
                    print "No entry for %s (%s) in %s found" %(gaussiantable[species[s]],species[s],thermonameSQL)
                energy[species[s]]=float(row[0])
                sql="select enthalpy from %s where ID='%s' and temp=298.15 " %(thermodataSQL,row[1])
                cursor.execute(sql)
                row=cursor.fetchone()
                enthalpyroom[species[s]]=row[0]
            else:
                enthalpyroom[species[s]]=0.0
                energy[species[s]]=0.0
        targetent=0.0

        print
        rxnspecies=0.0
        elembal=[]
        for s in spi:
            n=value(stoich[s]) 
            if n<0:
                if speciestable[species[s]].known==1:
                    print "%+.2f %s E=%s kJ/mol H(298K)=%s kJ/mol dH(0->298K)=%s kJ/mol \n"%(-n,species[s],energy[species[s]],enthalpy[species[s]],enthalpyroom[species[s]]),
                    targetent=targetent+(enthalpy[species[s]]-energy[species[s]]-enthalpyroom[species[s]])*(-n)
                    # N.B. Gaussian species energies in kJ/mol (known enthalpies in kJ/mol)
                    rxnspecies=rxnspecies+1.0
                else:
                    print "%+.2f %s E=%s kJ/mol H=unknown kJ/mol dH(0->298K)=%s kJ/mol \n"%(-n,species[s],energy[species[s]],enthalpyroom[species[s]]),
                    targetent=targetent-energy[species[s]]*(-n)
                    search=species[s]
                    rxnspecies=rxnspecies+1.0
                for e in eli:
                    count=0
                    mymatch=re.compile(elements[e]+'(\d*)')
                    for m in mymatch.findall(speciestable[species[s]].OBMol.GetFormula()): #species[s]
                        count+= int(m or 1)
                    elembal=elembal+[elements[e],count,n]
                        
        print "<=>\n", 
        for s in spi:
            n=value(stoich[s])
            if n>0: 
                if speciestable[species[s]].known==1:
                    print "%+.2f %s E=%s kJ/mol H(298K)=%s kJ/mol dH(0->298K)=%s kJ/mol \n"%(-n,species[s],energy[species[s]],enthalpy[species[s]],enthalpyroom[species[s]]),
                    # It is assumed that the enthalpy of the known species is given at 298K in kJ/mol
                    targetent=targetent+(-enthalpy[species[s]]+energy[species[s]]+enthalpyroom[species[s]])*n
                    rxnspecies=rxnspecies+1.0
                else:
                    print "%+.2f %s E=%s kJ/mol H=unknown kJ/mol dH(0->298K)=%s kJ/mol \n"%(-n,species[s],energy[species[s]],enthalpyroom[species[s]]),
                    targetent=targetent+energy[species[s]]*n
                    search=species[s]
                    rxnspecies=rxnspecies+1.0
                for e in eli:
                    count=0
                    mymatch=re.compile(elements[e]+'(\d*)')
                    for m in mymatch.findall(speciestable[species[s]].OBMol.GetFormula()): #species[s]
                        count+= int(m or 1)
                    elembal=elembal+[elements[e],count,n]

    ## work out if the number of elements balance with error criterion

        print
        print 'Check element balance:'
        
        rxnelement=len(elembal)/(3*rxnspecies)

        error=0.0
        
        for t in range(int(rxnelement)):
            sum1=0.000000000000000
            for i in range(int(rxnspecies)):
                sum1=sum1+round((elembal[1+3*(t+i*int(rxnelement))]*elembal[2+3*(t+i*int(rxnelement))]),5)
            print str(elembal[0+t*3])+': \t'+str(sum1)
            if abs(round(sum1,5))>1e-3:
                error=error+1.0
        print
        if error>0:
            print '!! There is an error in element balance !! \n'
            targetent='no enthalpy calculated'
            targetent298='no enthalpy calculated'
            

    ## Find enthalpy at 298K
        n=value(stoich[target])
        sql="select ID from %s where molefile like '%s'" %(thermonameSQL,gaussiantable[species[target]])
        cursor.execute(sql)
        row=cursor.fetchone()
        if cursor.rowcount==0:
            print "No entry for %s (%s) in %s found" %(gaussiantable[species[target]],species[s],thermonameSQL)
        sql="select enthalpy from %s where ID='%s' and temp=298.15 " %(thermodataSQL,row[0])
        cursor.execute(sql)
        row=cursor.fetchone()
        enthalpyroom[species[target]]=row[0]
        if targetent!='no enthalpy calculated':
            targetent298=targetent+enthalpyroom[species[target]]*n
        else:
            targetent298='no enthalpy calculated'
        

    ## Find Bonds
        #    print "bond names",
        #    print bondnames
        print "Bond type \tReact's\t->Products \tChange"
        for b in boi:
            bd_total = value(lpDot( bonds[b], modstoich))
            bd_net = value(lpDot( bonds[b], stoich))
            if bd_total:
                print "%s \t %.2f \t-> %.2f \t %.2f"% (bondnames[b], (bd_total-bd_net)/2, (bd_total+bd_net)/2, bd_net)
        bd_total = value(lpSum([lpDot( bonds[b], modstoich) for b in boi]))
        bd_net = value(lpSum(modbondsmade))
        print "Totals  \t %.2f \t-> %.2f \t %.2f \n"% ( (bd_total-bd_net)/2, (bd_total+bd_net)/2, bd_net) 
    #    H[search]=targetent
        print 'Enthalpy at 0K = '+str(targetent)
        print 'Enthalpy at 298K = '+str(targetent298)
        return  [targetent,targetent298,bd_net]




    #################################################################################################
    ## Run Isodesmic reactions
    #################################################################################################


    #invspecies=getinvspecies(1)
    #for s in invspecies:
    #   print [s]
    #   getenthalpy(1,[s])
    H={}
    numreactions={} 
    ##for s in spi: 
    ##    H[species[s]]=0.0 

    for s in spi:
        
        j=0
        print '-------------------------------------------------------------------------------------'
        print 
        print 'Species: \t'+str(species[s])
        print 
    #    print speciestable[species[s]].known
        if speciestable[species[s]].known==1:
            print'Species enthalpy of formation already known'
            print
        if speciestable[species[s]].known==0:
            H[species[s],j]=getenthalpy(s,[''],bondnames)
            invspecies=getinvspecies(s)
    #        print invspecies
            for i in invspecies:
                if speciestable[i].known==1:
                    j=j+1
                    H[species[s],j]=getenthalpy(s,[i],bondnames)
            numreactions[species[s]]=j+1                              



    #################################################################################################
    ## Run Isodesmic reactions statistics
    #################################################################################################


    print 
    print '-------------------------------------------------------------------------------------'
    print '-------------------------------------------------------------------------------------'
    print
    print 'Isodesmic reaction statistics'
    print 
    print '-------------------------------------------------------------------------------------'
    print '-------------------------------------------------------------------------------------'
    print 

    Haverage={}
    Hdev={}
    for s in spi:
      #  print 
       # print '-------------------------------------------------------------------------------------'
        average=0.0
        average298=0.0
        dev=0
       # if speciestable[species[s]].known==1:
       #     print 
       #     print 'Species: \t'+str(species[s])
       #     print 'Enthalpy of formation already known from literature' 
        if speciestable[species[s]].known==0:
    #           print  species[s] +"    =     " +str(H[species[s],0])
            print 
            print 'Species: \t'+str(species[s])
            print 'No. Rxns.: \t'+str(numreactions[species[s]])
            energies=[]
            energies298=[]
            counter=0
            for j in range(0,numreactions[species[s]]):
                if H[species[s],j][2]==0.0:                         ##Only use isodesmic reactions
                 if H[species[s],j][0]!='no enthalpy calculated':
                    if (energies.count(H[species[s],j][0])==0):
                        print '%s' %(str(H[species[s],j]))
                        average=average+H[species[s],j][0]
                        average298=average298+H[species[s],j][1]
                        counter=counter+1
                    energies.append(H[species[s],j][0])
                    energies298.append(H[species[s],j][1])
            energies=[]
            average=average/counter
            average298=average298/counter
            for j in range(0,numreactions[species[s]]):
                if H[species[s],j][0]!='no enthalpy calculated':
                 if H[species[s],j][2]==0.0:                         ##Only use isodesmic reactions
                    if (energies.count(H[species[s],j][0])==0):
                        dev=dev+(1.0/counter)*(H[species[s],j][0]-average)*(H[species[s],j][0]-average)
                    energies.append(H[species[s],j][0])
            print 'Average at 0K = \t'+str(average)
            print 'Average at 298K = \t'+str(average298)
            
            dev=sqrt(dev)
            Hdev[species[s]]=dev
            Haverage[species[s]]=average
            print 'Deviation = \t \t'+str(dev)
            sql="Insert into %s (Name,Enthalpy0K,error,Enthalpy298K) VALUES ('%s','%f','%f','%f')" %(comoentSQL,species[s],average,dev,average298)
            cursor.execute(sql)

    print
    print '-------------------------------------------------------------------------------------'
    print
    print 'Enthalpies of formation at 0K (kJ/mol):'
    print
    for s in spi:
        if speciestable[species[s]].known==0:
            print "H("+ species[s] +")="+str(Haverage[species[s]])
    #        sql="Insert into %s (Name,Enthalpy) VALUES ('%s','%f')" %(comoentSQL,gaussiantable[species[s]].replace('%',''),H[species[s]])
    #        cursor.execute(sql)

    ##########################################################################################################################
    ## Find known enthalpies at 0K
            
    for s in spi:
        if speciestable[species[s]].known==1:
            sql1="select Enthalpy from %s where name='%s'" %(knownentSQL,speciestable[species[s]].write('smi').split('\t')[0])
            cursor.execute(sql1)
            
            if cursor.rowcount==0:
                print "No entry for %s in %s found" %(speciestable[species[s]].write('smi'),knownentSQL)
            
            row=cursor.fetchone()
            enthalpy[species[s]]=row[0]
            
            sql2="select ID from %s where molefile like '%s'" %(thermonameSQL,gaussiantable[species[s]])
            cursor.execute(sql2)
            row=cursor.fetchone()
            if cursor.rowcount==0:
                print "No entry for %s (%s) in %s found" %(gaussiantable[species[s]],species[s],thermonameSQL)

    #        sql3="select enthalpy from %s where ID='%s' and temp='298.15'" %(thermodataSQL,row[0])
    #        cursor.execute(sql3)
    #        row=cursor.fetchone()
    #        enthalpyroom[species[s]]=row[0]
    #        H[species[s]]=enthalpy[species[s]]-enthalpyroom[species[s]]
    #        print "H("+ species[s] +")="+str(H[species[s]])
    #        sql="Insert into %s (Name,Enthalpy0K,Enthalpy298K) VALUES ('%s','%f','%f')" %(comoentSQL,gaussiantable[species[s]].replace('%',''),H[species[s]],enthalpy[species[s]])
    #        cursor.execute(sql)












    #for s in spi:
    #       smile= speciestable[species[s]].write('smi').split('\t')[0]
    #        herzler= speciestable[species[s]].write('smi').split('\t')[1].split('\n')[0]
    #       print smile
    #       print herzler    
    #       print "\n"
    #       sql="select * from chemkinnames where Herzler='%s'" % (herzler)
    #       cursor.execute(sql)
    #       if (cursor.rowcount==0):
    #               print herzler
    #               sql="insert into chemkinnames (Herzler,Chemkin,Smile) VALUES ('%s','%s','%s') " % (herzler,herzler,smile)
    #               print sql   
    #               cursor.execute(sql) 

    ##for s in spi:
    ##    if speciestable[species[s]].known==0:
    ##        smile= speciestable[species[s]].write('smi').split('\t')[0]
    ##        sql="select latex from chemkinnames where Smile='%s'" % (smile)
    ##        cursor.execute(sql)
    ##        row=cursor.fetchone()
    ##        print "%s && $%s \pm %s$&\\\\" % (row[0],str(round(Haverage[species[s]])).split('.')[0],str(round(Hdev[species[s]])).split('.')[0])
    #       print row[0]
    #       print "\n"

    # print H
    #    print "total bonds\t",
    #    print [value(lpDot( bonds[b], modstoich)) for b in boi],
    #    
    #    print "sum total bonds\t",
    #    print value(lpSum([lpDot( bonds[b], modstoich) for b in boi]))
    #    
    #    print "net bonds made\t",
    #    print [value(lpDot( bonds[b], stoich)) for b in boi],
    #    
    #    print "sum mod net bonds made",
    #    print value(lpSum(modbondsmade))
            

