# Script to simulate constant pressure reactor for a given time it spits out
# temps and species concentrations at all the time steps.

# 17 NOV 2009 Started reactor simulator based on Equilibrium.py and an
# old reactor simulator I wrote for the TEOS work - ras81 


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
P1 = 2.36E5 #Pressure Pa
Tini = 1561.0    #Temperature K
time = 0.3     # total time to run s
dt = 0.001   # step size s

Ntimesteps = time/dt + 1   #Num timesteps + 1 for 0s


## DATA FROM CTI FILES
cti_file = 'prf33.cti' # cti file

# SET INITIAL CONCENTRATIONS
#  the following expects mole fractions:
initialcomposition = 'C7H16:1.35003802023096E-20, C7H15:4.25011032373343E-23, C7H14:1.76004520822144E-23, C7H15O2:2.15005606963173E-24, C7H14OOH:9.19022590343797E-26, O2C7H14OOH:9.62023491754841E-26, C7KET:1.7700488710513E-23, C5H11CO:1.0200271376408E-20, C5H11:8.76021688932753E-21, C8H18:8.90023517170702E-19, C8H17:3.22008502280254E-22, C8H16:1.07002895317081E-21, C8H17O2:3.89010495127297E-23, C8H16OOH:3.99010308279325E-25, O2C8H16OOH:2.39006148445196E-25, C8KET:3.14008321786246E-23, C6H13CO:1.81004702375146E-20, C6H13:1.66004157716141E-18, C3H7:7.27017158579661E-10, C3H6:2.42006147386201E-06, C2H4:1.09003077929077E-08, C2H3:2.74007419316208E-10, CH2O:1.03002530093089E-10, HCO:8.71022607287704E-10, CO:0.00148004164070106, CO2:0.111002710587097, H2O2:2.90007780304224E-06, HO2:0.000289007414021237, H:0.000151004163011112, OH:0.00125003438917093, H2O:0.126003255246103, O2:0.0198005429646148, N2:0.740018070580647 '

plotgraphs = True  #plots concentrations against T
LowC = 1E-60   #Lowest concentraion on log plots

molefractions = False    # If False writes mass fractions

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
mf_species = range(Ntimesteps)

################################

#### Get values of species index
for i in range(0, num_species):
    ispecies[i] = gas.speciesIndex(species[i])
    #print ispecies[i]

#### WRITE CSV FILE ####
    
file2 = open('results.csv','w')
speciesandtemp = range(num_species+1)
speciesandtemp = ['Time (s)'] + ['Temp (K)'] + species 
writeCSV(file2, speciesandtemp )

## Define a two-dimensional array to hold mole-fractions and temps of each species at all time steps
num_speciesplustwo = num_species + 2   # + 1 for time and  + 1 for temp
for i in range(Ntimesteps):
    mf_species[i] = range(0, num_speciesplustwo)
    

## set initial conditions!
gas.set(T = Tini, P = P1, X = "%s" %(initialcomposition) )


######################### SET UP THE REACTOR
r   = Reactor(gas)

gas_b = Air()
gas_b.set(P = P1)
env = Reservoir(gas_b)


# Define a wall between the reactor and the environment, and
# make it flexible, so that the pressure in the reactor is held
# at the environment pressure.
w = Wall(r,env)
w.set(K = 1.0e6)   # set expansion parameter. dV/dt = KA(P_1 - P_2)
w.set(A = 1.0)
######################### SET UP THE REACTOR

sim = ReactorNet([r])
time = 0.0
tim = zeros(Ntimesteps,'d')
#Temp = range(Ntimesteps)


for run in range(Ntimesteps):
    tim[run] = time
    mf_species[run][0] = time
    mf_species[run][1] = gas.temperature()
    for i in range(2, num_speciesplustwo):
        if molefractions == True:
            mf_species[run][i] = gas.moleFraction(ispecies[i-2])  # mole fractions
        else:
            mf_species[run][i] = gas.massFraction(ispecies[i-2])# mass fractions
    writeCSV(file2,   mf_species[run] )
    time += dt
    sim.advance(time)

 
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
    timeslist = range(Ntimesteps)
    for run in range(Ntimesteps):
        timeslist[run] = mf_species[run][0]

    speciesconcentrations = range(Ntimesteps)
    for i in range(num_species):
        for run in range(Ntimesteps):
            speciesconcentrations[run] = mf_species[run][i+1]
        semilogy(timeslist, speciesconcentrations, line_type[i], label=species[i])
    axis([0, time, LowC, 1])
    xlabel('time (s)')
    ylabel('mol fraction')

    figure()
    
    for i in range(num_species):
        for run in range(Ntimesteps):
            speciesconcentrations[run] = mf_species[run][i+1]
        semilogy(timeslist, speciesconcentrations, line_type[i], label=species[i])
    axis([0, time, LowC, 1])
    legend(loc='best')
    xlabel('time (s)')
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



