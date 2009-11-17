#! /usr/bin/env python

# give the list of log files as arguments on the command line, eg.
# python extract.py *.log

# it identifies transition states according to filename: eg. ts_reactant__product.log or ts_reactant1_reactant2__product1_product2.log
# command line option -r to reset the sqlite database


import getopt, sys, os, subprocess, re

#from BackwardsReader import BackwardsReader

from QuantumDatabase import *


def readlogfile(logfilename):
    fin=open(logfilename,'r')
    resultblocks=list()
    for line in fin:
        resultblock=''
        if line[0:5]==u' 1\\1\\':
            while line.strip()[-1]!=u'@':
                resultblock=resultblock+line.strip()
                line=fin.next()
            resultblock=resultblock+line.strip() # don't forget the last line!
            resultblocks.append(resultblock)
    fin.close()
    
    results=list()
    for i in resultblocks:
        temp=i.split('\\\\')
        res=list()
        for j in temp:
            res.append(j.split('\\'))
        results.append(res)
    return results


def makeReaction():
    return 1

#
#def main():

try:
    opts, args = getopt.getopt(sys.argv[1:], "rho:", ["reset" "help", "output="])
except getopt.GetoptError, err:
    # print help information and exit:
    print str(err) # will print something like "option -a not recognized"
    usage()
    sys.exit(2)
output = None
reset=True #### this should probably normally be false
for o, a in opts:
    if o == "-v":
        verbose = True
    elif o in ("-h", "--help"):
        usage()
        sys.exit()
    elif o in ("-o", "--output"):
        output = a
    elif o in ("-r", "--reset"):
        reset=True
    else:
        assert False, "unhandled option"

if reset:
    print "Resetting database"
else:
    print "Not resetting database"
connectToDatabase(reset=reset)


if len(args):
    print "Importing files ",args
else:
    args=['CH4.log','CH3J.log',
          'CH3NHCH3.log','CH3NJCH3.log', 
          'ts_CH3NHCH3_CH3J__CH3NJCH3_CH4.log']
    print "No files specified. Testing with ",args
for filename in args:

    results=readlogfile(filename)
    jobname=os.path.basename(filename)
    (jobname,extension)=os.path.splitext(jobname)
    
    if re.search('\Wts\W',results[0][1][0]):
        print "I think %s is a transition state"%jobname
        species=None
        reactionname=jobname
        ### look up reaction by name. not 100% safe
        try: 
            reaction=Reaction.byName(reactionname)
            print "reaction %s already found. Assume its reactants are correct!"%reactionname
        except SQLObjectNotFound:
            print "reaction %s doesn't exist yet. Creating now."%reactionname
            reaction=Reaction(name=reactionname)
            (reactants,products)= reactionname.split('__')
            reactants=reactants.split('_')
            products=products.split('_')
            reaction.makeFromReactantsProducts(reactants,products)
            species=None ### currently essential!!! should make better structured
        
    else:   
        print "I think %s is a stable species"%jobname
        reaction=None
        # it's a species
        speciesname=jobname
        try: 
            species=Species.byName(speciesname)
        except SQLObjectNotFound:
            print "species %s doesn't exist yet. Creating now."%speciesname
            species=Species(name=speciesname)
        
    for res in results:
        #assert jobname==res[2][0], "job name has changed"
        if species and species.formula:
            assert species.formula==res[0][6], "species formula has changed"
        elif species:
            species.formula=res[0][6]
            
        c=Calculation(species=species, # species may be "None" or 
                      reaction=reaction, # reaction may be "None"
                      jobType=res[0][3],
                      method=res[0][4],
                      basisSet=res[0][5],
                      commandLine=res[1][0]
                      )
        # geometry is stored in res[3]
        for item in res[4]:
            (name,value)=item.split('=')
            r=Result(calculation=c,
                     name=name,
                     value=value
                     )
    
#if __name__ == "__main__":
#    main()
#    