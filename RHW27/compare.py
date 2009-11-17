#! /usr/bin/env python

# used to compare species and find the energy change and barrier height
# of a list of reactions in reactionlist.txt

import getopt, sys, os, subprocess, re

from QuantumDatabase import *
from UnitUtilities import *
dreload(Calculation)

# def main():

try:
    opts, args = getopt.getopt(sys.argv[1:], "ho:", ["help", "output="])
except getopt.GetoptError, err:
    # print help information and exit:
    print str(err) # will print something like "option -a not recognized"
    usage()
    sys.exit(2)
output = None
verbose = False
saveonly = False
template = ''
for o, a in opts:
    if o == "-v":
        verbose = True
    elif o in ("-h", "--help"):
        usage()
        sys.exit()
    elif o in ("-o", "--output"):
        output = a
    else:
        assert False, "unhandled option"


print "Connecting to database"
connectToDatabase(reset=False, debug=False)

# if len(args):
#     print "Reaction",args
# else:
#     args=['CH4 = CH3J + HJ']
#     print "No reaction specified. Testing with ",args
# for reaction in args:
#     (reactants,products)= reaction.split('=')
#     reactants=reactants.split('+')
#     products=products.split('+')
#     
    
##### TESTING

print "Testing the arithmetic of calculations: ",
a=Calculation.get(1)
assert float(a['HF'].value)<0
b=a+a
c=b-3*a # = -a
d=a+c # = 0
assert round(d['HF'].value,13)==0 # let it vary by 10^-13 
print "passed OK"

print "Testing the arithmetic of species:"
a=Species.get(1)
print "HF energy of %s = %f"%(a.name,float(a.calculations[0]['HF'].value))
assert float(a.calculations[0]['HF'].value)<0
b=a+a
c=b-3*a # = -a
d=a+c # = 0
print "HF energy of %s = %f"%(d.name,float(d.calculations[0]['HF'].value))
assert round(d.calculations[0]['HF'].value,13)==0, "Calculation problem"
print "passed OK"

### end of testing



### interpreting results

hline = r'\hline'
hline = ''

latexfile=open('latex/scheme.tex','w')
latexfile.write("%!TEX root =  Document.tex\n")

## list of reactions

### reading prime formatted xml
# import xml.dom.minidom
# doc=xml.dom.minidom.parse("Reactions.xml")
# rxlist=doc.getElementsByTagNameNS("http://purl.org/NET/prime/","reaction")
# for reaction in rxlist[]:
#     reactants=reaction.getElementsByTagName("reactants")[0].getElementsByTagName("speciesLink")
#     for reactant in reactants:
#         stoichiometry=reactant.childNodes[0]
#         assert stoichiometry.nodeType==stoichiometry.TEXT_NODE, "Oops. Was expecting a text node containing the stoichiometry"
#         stoichiometry=float(stoichiometry.nodeValue)
#         reactantname=reactant.attributes['preferredKey'].value
    
reactionsfile=open('reactionlist.txt','r')
for reactionstring in reactionsfile:
    #reactionstring="CH3NHCH3 CH3J = CH3NJCH3 CH4"
    
    if reactionstring.find(' = ')==-1: # couldn't find search term
        latexcomment=r'\multicolumn{6}{p{\textwidth}}{'+reactionstring+r'} \\ '+hline+"\n"
        latexfile.write(latexcomment)
        continue

    (reactants,products)= reactionstring.split(' = ')
    reactants=reactants.split() # split on whitespace
    products=products.split() # split on whitespace
                
    reactionname='ts'
    for s in reactants:
        reactionname+='_'+s
    reactionname+='_'
    for s in products:
        reactionname+='_'+s
    
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


    r = reaction
#for r in Reaction.select():
    print
    print "Reaction:",r.name
    rs=Reactant.select(Reactant.q.reaction==r)
    answer=None
    latexreactant=None
    latexproduct=None
    
    desiredSettings=DerivedCalculation(jobType='Mixed', method='CBS-QB3', basisSet='CBS-QB3')
    
    for reactant in rs:
        print "  Reactant: %s, Calculations: %d"%(reactant.species.name,reactant.species.calculations.count())
        # list(reactant.species.calculations)
        if answer:
            answer += reactant.species * reactant.stoichiometry
        else:
            answer = reactant.species * reactant.stoichiometry
        
        #### LATEX STUFF
        latexthis=reactant.species.latex()
                
        try:
            reactant.species.findCalculation(desiredSettings)
        except IndexError:
            latexthis="\colorbox{yellow}{%s}"%latexthis
        
        if abs(reactant.stoichiometry)!=1: # add stoichiometry if not 1
            latexthis="%g %s \n"%(abs(reactant.stoichiometry),latexthis) 
        if reactant.stoichiometry<0:
            if latexreactant:
                latexreactant += " + " + latexthis
            else:
                latexreactant = latexthis
        else:
            if latexproduct:
                latexproduct += " + " + latexthis
            else:
                latexproduct = latexthis     
            
    print "Sum:",answer.name
    print "Common Calculations:",len(answer.calculations)
    reactionChange=answer
    
    
    try:
        reactionChange.findCalculation(desiredSettings)
    except IndexError:
        print "No CBS-QB3 calculation available"
        CBSQB3energychange=None
    else:
        CBSQB3energychange=convertHartreeToKcalmol(reactionChange.findCalculation(desiredSettings)['CBSQB3'].value)
        print "CBS-QB3 Energy change: %f kCal/Mol" % CBSQB3energychange

 
    
    print "Transition state:"
    answer=r*1
    for reactant in rs.filter(Reactant.q.stoichiometry<0):
        answer += reactant.species * reactant.stoichiometry
    print "Sum:",answer.name
    print "Common Calculations:",len(answer.calculations)
    transitionStateChange=answer
    

    desiredSettings=DerivedCalculation(jobType='Mixed', method='CBS-QB3', basisSet='CBS-QB3')
    try:
        transitionStateChange.findCalculation(desiredSettings)
        CBSQB3barrier=convertHartreeToKcalmol(transitionStateChange.findCalculation(desiredSettings)['CBSQB3'].value)
        print "CBS-QB3 Barrier height: %f kCal/Mol" % CBSQB3barrier
    except IndexError:
        print "No CBS-QB3 calculation available"
        CBSQB3barrier = None;
    
    ##### LATEX STUFF 
    latexreaction= "%d &\n"%r.id  # print the sqlobject reaction id
    # latexreaction = " & %% blank column\n" # don't print the reaction id
    latexreaction+= " %s &\n %s &\n"%(latexreactant,latexproduct) 
    if CBSQB3energychange:
        latexreaction+= "$%.1f$ & %% CBS-QB3 energy change \n"%CBSQB3energychange
    else:
        latexreaction+= "--- & %% insufficient CBS-QB3 calcs to find this\n"
    if CBSQB3barrier:
        latexreaction+= "$%.1f$ & %% CBS-QB3 barrier height \n"%CBSQB3barrier 
        latexreaction+= "$%.1f$  %% CBS-QB3 barrier height in exothermic direction \n"%( CBSQB3barrier-CBSQB3energychange if (CBSQB3energychange>0) else CBSQB3barrier )
    else:
        latexreaction+= "--- %% insufficient CBS-QB3 calcs to find this\n"
    latexreaction+= r"\\"+ hline + "\n"
    latexfile.write(latexreaction)

latexfile.close()

print "Species not used in these reactions:"
for sp in Species.select("""species.id not in (select species_id from reactant)"""):
    string= "Name: %20s  Formula: %s , Calculations: %d"%(sp.name, sp.formula, sp.calculations.count())
    print string


print "typesetting..."
try:
    output = subprocess.Popen(["pdflatex", 'document.tex'], cwd='latex', stdout=subprocess.PIPE).communicate()[0]
    print "I have the output '%s' to report"%output[-100:-1].strip()
    print "Now trying to open a preview with quicklook"
    output=subprocess.Popen(["qlmanage", "-p", 'document.pdf'],cwd='latex', stderr=subprocess.PIPE)
except OSError, e:
    print >>sys.stderr, "Execution failed:", e


# if __name__ == "__main__":
#     main()
#     