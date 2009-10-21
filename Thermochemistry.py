# Program written by Raph to calculate partition function and thermo properties
# from the data in the SQL database.

#started March 08

# 09 JULY 08 - extended to calculate thermo properties
# as well as partition functions - 09 JULY 08

# 10 JULY 08 - bugfixing + made to produce CSV file - 10 JULY 08

# 11 JULY 2008 - started making it read SQL/write to SQL +
# take account of internal rotors - 11 JULY 2008

# 5 SEPT 2008 - finally figured out what was wrong with entropy
# translational partition function is multiplied by volume and
# therefore introduces pressure term using basic ideal gas assumption
# this introduces power of T - 5 SEPT 2008

# 25 FEB 2009 - modified slightly by ss663 to go through all the species in
# referred SQL database to write thermodata into SQL database

# 3 MARCH 2009 - modified by tst25 to read symmetry number from SQL
# database and to be able to handle atoms. Also made inputs for SQL table names.

# 26 MARCH 2009 - fixed bug with multiplicity - now giving correct entropy

# 21 OCT 2009 - took sensitiv data out and put it in separate file that isn't made public - ras81

##########################INPUTS###################################################
Ireduced = [0]     #0 if no internal rotors, [1,2,3] if 3 of them.
Tmax = 4000.0
Tstep = 20.0

# SQL information
project = 'blank'
thermonameSQL = 'blank'
thermodataSQL = 'blank'
frequencySQL = 'blank'

#s = 2              #Symmetry factor
CSV = True        #writes CSV file results.csv
SQL = True        #writes to SQL database
plotgraphs = False  #plots all thermo properties against T
#############################END INPUTS############################################

import MySQLdb
import openbabel, pybel
from Cantera import *
from math import *
import random
from pylab import *
from scipy import *
from numpy import *
import private  # Database info

# define constants
kb = 1.3806503E-23         #Boltzman's constant
h  = 6.626068E-34          #Planck's constant
c  = 299792458             #Speed of light
Na = 6.0221415E23          #Avogadro's Constant
c12over12 = 1.66053886E-27 #mass of carbon12/12
a0 = 5.29177E-11           #bohr radius
P = 101325                 #1atm = 101325Pa

convertcm = 100*h*c/kb # converts the wavenumbers for powers ie convert to joules and divide by k.
hbar = h/(2*pi)

########################DEFINE ARRAY OF TEMPS######################################
numTs = Tmax/Tstep + 1
Tlist=[]
Tlow=[]
Thigh=[]
for n in range(1,numTs):
    if Tstep*n < 298.15:
        Tlow = Tlow + [Tstep*n]
    else:
        Thigh = Thigh + [Tstep*n]

Tlist=Tlow + [298.15]+ Thigh
Hlist = range(len(Tlist))
Slist = range(len(Tlist))
Cplist = range(len(Tlist))
Glist = range(len(Tlist))
qviblist = range(len(Tlist))
qrotlist = range(len(Tlist))
qtranslist = range(len(Tlist))
dq_vibBYdToverq_viblist = range(len(Tlist))
##########################END DEFINE ARRAY OF TEMPS######################################


#######################OPEN LINK TO SQL#####################################
conn = MySQLdb.connect (host = private.defaulthost,
                        user = private.defaultuser,
                        passwd = private.defaultpasswd,
                        db = private.defaultdb);
cursor = conn.cursor ()
#Get the list of species in all formats from chemkinnames on sql
sql1="select Name, ID, spin, Ia, Ib, Ic, mass, symmetrynumber from %s where project='%s'" %(thermonameSQL, project)  
cursor.execute(sql1)
speciesdata1 = cursor.fetchall()
num_species = len(speciesdata1)
#print num_species

for i in range(0,num_species):
    M = speciesdata1[i][6]                  #Mass
    MOI = range(0,3)
    MOI[0] = speciesdata1[i][3]             #Ix
    MOI[1] = speciesdata1[i][4]             #Iy
    MOI[2] = speciesdata1[i][5]             #Iz
    s = speciesdata1[i][7]                  #Symmetry number
    g = float(speciesdata1[i][2])                  # multiplicity
    #print g

    sql2="select frequency from %s where ID='%i'" %(frequencySQL, speciesdata1[i][1])
    cursor.execute(sql2)
    Vini = cursor.fetchall()
    #print Vini
    V = range(len(Vini))
    for run8 in range(len(Vini)):
        V[run8] = Vini[run8][0]    
    Ireduced = [0]        # internal rotor reduced MOIs

    if MOI[0] == 0 or MOI[1] == 0 or MOI[2] == 0:
        print 'linear molecule'
    else:
        print 'non-linear molecule'

    for run in range(numTs):
        T=Tlist[run]

        #####################################PARTITION FUNCTION#########################################    

        # Calculate translational part of partition function
        # Factor E6 converts from m to cm (industry standard) dimension (L^-3) comes from emitted volume 
        # Removed factor E6 it has no effect on Thermo properties - 10 July 2008
        # does affect the entropy and free energy - 14 july 2008
        q_trans = (( 2 * pi * M * c12over12 * kb * T )**(1.5))/((h**3.0))#* 1000000) 
        #print q_trans

        # Calculate rotational part of partition function
        # accounting for whether or not the molecule is linear
        if MOI[0] == 0.0 or MOI[1] == 0.0 or MOI[2] == 0.0:
            
            # If just an atom all MOI=0 and q_rot=1 with m=1.5
            if MOI[0] == 0.0 and MOI[1] == 0.0 and MOI[2] == 0.0:
                q_rot = 1.0
                m = 1.5
            else:
                q_rot = ((8 * pi * pi * ((MOI[0]+MOI[1]+MOI[2])/2)* (c12over12 * a0**2.0) * kb * T)/(s * h**2.0))
                m = 2.5
            #print q_rot
            
        else:
            q_rot = 8 * pi**2.0 * (( (8 * (pi**3.0) * MOI[0]*MOI[1]*MOI[2]* ((c12over12 * (a0**2))**3))**0.5) *( (kb*T)**1.5) / (s* (h**3.0) )) #(c12over12 * (a0**2))**3)
            m = 3.0
            #print q_rot

        # Calculate vibrational part of partition function
        q_vib = 1.0
        for mode in V:
            if mode <= 0.0:
                print 'negative/zero mode rejected'
            else:
                q_vib = q_vib  / (1.0 - exp(-convertcm*mode/T)) # Taking E0 = 0
                #q_vib = q_vib * exp(-1.43877*mode/(2*T)) / (1 - exp(-1.43877*mode/T)) # Taking E0 = hf/2 (including ZPVE)
        
        # Calculate total partition function.
        # q factorises because assume that modes are independent (fairly good approx)
        q =  g * q_trans * q_rot * q_vib
        #print 'translational partition function', q_trans
        #print 'rotational partition function', q_rot
        #print 'vibrational partition function', q_vib
        #print 'partition function', q

        #############################INTERNAL ROTOR####################
        if not Ireduced == [0]:
            n = len(Ireduced)
            m = m + n
            q_IR = 1.0
            for run6 in range(n):
                q_IR = q_IR*(8*(pi**3)*Ireduced[run6]*kb*T)**(0.5)/h
            q = q*q_IR        
        #########################END INTERNAL ROTOR####################

        #################################END OF PARTITION FUNCTION#######################################

        #############################FIRST DERIV OF PARTITION FUNCTION################################### 
        #first calculate first derivative of vibrational part
        dq_vibBYdToverq_vib = 0.0
        for mode2 in V:
            if mode2 <= 0.0:
                print 'negative mode rejected'
            else:
                dq_vibBYdToverq_vib = dq_vibBYdToverq_vib + (convertcm*mode2/(T**2))*exp(-convertcm*mode2/T) / (1.0 - exp(-convertcm*mode2/T))


        #dqBYdT=g * q_trans * q_rot * (4/T) * q_vib + g * q_trans * q_rot * dq_vibBYdT #Basic product rule

        dqBYdToverq = m/T + dq_vibBYdToverq_vib

        
        #########################END OF FIRST DERIV OF PARTITION FUNCTION################################### 
    

        #############################SECOND DERIV OF PARTITION FUNCTION################################### 
        #first calculate second derivative of vibrational part
        # THIS ALGEBRA IS HORRIFIC!!!!
        # the second derivative of q_vib is a number of terms relating to the first deriv
        # plus a sum over the square of the sum in the first deriv.

        #this is sum over first sum by omega*consts
        Horridsum1 = 0.0
        for mode3 in V:
            if mode3 <= 0.0:
                print 'negative mode rejected'
            else:
                Horridsum1 = Horridsum1 + ((-convertcm*mode3/(T**2))*exp(-convertcm*mode3/T) / (1.0 - exp(-convertcm*mode3/T)))*(-convertcm*mode3/(T**2))

        # this is sum over square of first sum
        Horridsum2 = 0.0
        for mode4 in V:
            if mode4 < 0.0:
                print 'negative mode rejected'
            else:
                Horridsum2 = Horridsum2 + ((-convertcm*mode4/(T**2))*exp(-convertcm*mode4/T) / (1.0 - exp(-convertcm*mode4/T)))**2


    
        d2qBYdTsquaredoverq = ((m**2 -m)/(T**2)
                               + ((2*m)/T)*dq_vibBYdToverq_vib
                               + dq_vibBYdToverq_vib**2
                               - (2.0/T)*dq_vibBYdToverq_vib
                               + Horridsum1
                               + Horridsum2)

        
        #########################END OF SECOND DERIV OF PARTITION FUNCTION################################### 
    
        #############################THERMODYNAMIC PROPERTIES################################### dqBYdToverq  d2qBYdTsquaredoverq
        H = 0
        S = 0
        Cp = 0
        G = 0
        H = (Na * kb * T**2)*dqBYdToverq + Na*kb*T

        S = Na * kb * (log(q) + T*dqBYdToverq + log(kb*T) - log(P) + 1.0)  

        Cp = Na*kb + Na*kb*T*(2.0*dqBYdToverq - T*((dqBYdToverq)**2) + T*d2qBYdTsquaredoverq)
    
        G = H - T*S

        #These are in SI units
        Hlist[run] = H/1000.0        #kJ/mol
        Slist[run] = S/1000.0        #kJ/mol.K               # S        #J/mol.K
        Cplist[run] = Cp/1000.0      #kJ/mol.K             # Cp       #J/mol.K
        Glist[run] = G/1000.0        #kJ/mol
    
        qviblist[run] =  q_vib
        qrotlist[run] =  q_rot 
        qtranslist[run] = q_trans
        dq_vibBYdToverq_viblist[run] = dq_vibBYdToverq_vib
        #########################END OF THERMODYNAMIC PROPERTIES###################################

    if CSV:
        print 'Writing CSV file in kJ/mol.K and kJ/mol'
        file2 = open('thermodata.csv','w')
        writeCSV(file2, ['Temp (K)', 'H (kJ/mol)', 'S (J/mol.K)', 'Cp (J/mol.K)', 'G (kJ/mol)', 'q/qvib' ] )
        for run4 in range(numTs):
            writeCSV(file2, [Tlist[run4], Hlist[run4], Slist[run4], Cplist[run4], Glist[run4], qviblist[run4], qrotlist[run4], qtranslist[run4], dq_vibBYdToverq_viblist[run4] ] )
        file2.close()



    if SQL:
        print 'writing to SQL in kJ/mol.K and kJ/mol'
        for run5 in range(numTs):
            writesql="INSERT into %s (ID, name, temp, entropy, heatcapacity, enthalpy, freeenergy) VALUES ('%i', '%s', '%f', '%f', '%f', '%f', '%f')" %(thermodataSQL, speciesdata1[i][1], speciesdata1[i][0], Tlist[run5], Slist[run5], Cplist[run5], Hlist[run5], Glist[run5])
            cursor.execute(writesql)
            
    

    if plotgraphs:
        print 'Plotting graphs in kJ/mol.K and kJ/mol'
        fig = figure()
        plot(Tlist, Hlist)
        xlabel('T (K)')
        ylabel('H (kJ/mol)')
        fig = figure()
        plot(Tlist, Slist)
        xlabel('T (K)')
        ylabel('S (kJ/mol.K)')
        fig = figure()
        plot(Tlist, Cplist)
        xlabel('T (K)')
        ylabel('Cp (kJ/mol.K)')
        fig = figure()
        plot(Tlist, Glist)
        xlabel('T (K)')
        ylabel('G (kJ/mol)')
    
#sql2="select frequency from frequency where ID='%i'" %(speciesdata1[1])
#masssql="Update Thermoname set mass=%f where ID='%i'" %(float(mass),id)
#cursor.execute(masssql) 
