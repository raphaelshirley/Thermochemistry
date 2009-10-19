# Code to calculate the equilibrium composition
## MODIFIED BY ss663 on 04/03/2009 to avoid the tedious job of
## individually entering species data. This version is more general.


from Cantera import *
from Cantera.Reactor import *
from Cantera.Func import *
from math import *
import random
from pylab import *
from scipy import *
import numpy


##### USER VARIABLES ##########################################################

## PROCESS PARAMETERS
P1 = 1.01E5 #Pressure
T1 = 500    #Lower-bound temp
T2 = 4000   #Upper-bound temp
dT = 25     #Step size
T = range(T1, (T2+dT), dT)
Ntemps = (T2 - T1)/dT + 1

## DATA FROM CTI FILES
cti_file = 'fullchemnewTEST.cti' # cti file
initialcomposition = 'AlCl3:0.05,Cl4Ti:0.475, O2:0.475'

plotgraphs = False  #plots concentrations against T

##### END OF USER VARIABLES ##########################################################

## Get all species names out of cti file - following method gets from thermo if it is formatted as in Raphs
# polygenerator script
chemfile = open('./' + cti_file, 'r')
speciesname = []
species = []
for line in chemfile:
    
    if line[0:14] == 'species(name =':
        speciesname = line.split()[2]
        speciesname = speciesname.strip('",')
        species = species + [speciesname]
    else:
        continue
chemfile.close()

print 'species = ', species



###############################

#### DEFINING VARIABLES #####

## Define number of species
num_species = len(species)
print 'there are', num_species, 'species'

## Import cti file with species data
gas = importPhase(cti_file, 'gas')

## Define array to hold species index
ispecies = range(num_species)

## Define array to hold mole-fractions [this is actually a 2D-array]
mf_species = range(Ntemps)

################################

#### Get values of species index
for i in range(0, num_species):
    ispecies[i] = gas.speciesIndex(species[i])
    #print ispecies[i]

#### WRITE CSV FILE ####
    
file2 = open('results.csv','w')
speciesandtemp = range(num_species+1)
speciesandtemp = ['Temp (K)'] + species 
writeCSV(file2, speciesandtemp )

## Define a two-dimensional array to hold mole-fractions of each species at all T
num_speciesplusone = num_species + 1
for i in range(Ntemps):
    mf_species[i] = range(0, num_speciesplusone)
    

## AND EQUILIBRATE!
gas.set(T = T1, P = P1, X = "%s" %(initialcomposition) )
gas.equilibrate('TP', maxsteps=10000)





for run in range(Ntemps):
    gas.setTemperature(T[run])
    gas.equilibrate('TP', maxsteps=1000)
    #writeCSV(file2, [T[run]] )
    mf_species[run][0] = T[run]
    for i in range(1, num_speciesplusone):
        mf_species[run][i] = gas.moleFraction(ispecies[i-1])
    writeCSV(file2,   mf_species[run] )
  
file2.close()

###############################

##### PLOT MOLE_FRACTION Vs T #######

## Define an array with line_types for the graph [worst job ever!!!]
line_type = ('-b', '-g', '-r', '-c', '-m', '-y', '-k',
             '--b', '--g', '--r', '--c', '--m', '--y', '--k',
             '-.b', '-.g', '-.r', '-.c', '-.m', '-.y', '-.k',
             ':b', ':g', ':r', ':c', ':m', ':y', ':k',
             'ob', 'og', 'or', 'oc', 'om', 'oy', 'ok',
             '+b', '+g', '+r', '+c', '+m', '+y', '+k',
             'xb', 'xg', 'xr', 'xc', 'xm', 'xy', 'xk',
             'vb', 'vg', 'vr', 'vc', 'vm', 'vy', 'vk',
             'db', 'dg', 'dr', 'dc', 'dm', 'dy', 'dk',
             '1b', '1g', '1r', '1c', '1m', '1y', '1k',
             '2b', '2g', '2r', '2c', '2m', '2y', '2k',
             '3b', '3g', '3r', '3c', '3m', '3y', '3k',
             '4b', '4g', '4r', '4c', '4m', '4y', '4k',
             'pb', 'pg', 'pr', 'pc', 'pm', 'py', 'pk',
             'hb', 'hg', 'hr', 'hc', 'hm', 'hy', 'hk',
             '.b', '.g', '.r', '.c', '.m', '.y', '.k',
             '^b', '^g', '^r', '^c', '^m', '^y', '^k',
             '<b', '<g', '<r', '<c', '<m', '<y', '<k',
             '>b', '>g', '>r', '>c', '>m', '>y', '>k',
             'sb', 'sg', 'sr', 'sc', 'sm', 'sy', 'sk',
             'Db', 'Dg', 'Dr', 'Dc', 'Dm', 'Dy', 'Dk',
             'Hb', 'Hg', 'Hr', 'Hc', 'Hm', 'Hy', 'Hk',
             ',b', ',g', ',r', ',c', ',m', ',y', ',k',
             '_b', '_g', '_r', '_c', '_m', '_y', '_k',
             '|b', '|g', '|r', '|c', '|m', '|y', '|k'
             )


if plotgraphs == True:
    speciesconcentrations = range(Ntemps)
    for i in range(num_species):
        for run in range(Ntemps):
            speciesconcentrations[run] = mf_species[run][i+1]
        semilogy(T, speciesconcentrations, line_type[i], label=species[i])
    axis([500, 4000, 1E-8, 1])
    xlabel('T (K)')
    ylabel('mol fraction')

    figure()

    for i in range(num_species):
        for run in range(Ntemps):
            speciesconcentrations[run] = mf_species[run][i+1]
        semilogy(T, speciesconcentrations, line_type[i], label=species[i])
    axis([500, 4000, 1E-8, 1])
    legend(loc='best')
    xlabel('T (K)')
    ylabel('mol fraction')


##    -     # solid line
##    --    # dashed line
##    -.    # dash-dot line
##    :     # dotted line

##    b  # blue
##    g  # green
##    r  # red
##    c  # cyan
##    m  # magenta
##    y  # yellow
##    k  # black
##    w  # white

##    .     # points
##    ,     # pixels
##    o     # circle symbols
##    ^     # triangle up symbols
##    v     # triangle down symbols
##    <     # triangle left symbols
##    >     # triangle right symbols
##    s     # square symbols
##    +     # plus symbols
##    x     # cross symbols
##    D     # diamond symbols
##    d     # thin diamond symbols
##    1     # tripod down symbols
##    2     # tripod up symbols
##    3     # tripod left symbols
##    4     # tripod right symbols
##    h     # hexagon symbols
##    H     # rotated hexagon symbols
##    p     # pentagon symbols
##    |     # vertical line symbols
##    _     # horizontal line symbols
##    steps # use gnuplot style 'steps' # kwarg only
## 
##The following color abbreviations are supported::
## 



