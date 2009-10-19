#! /usr/bin/env python

## set up classes for species, reactions, calculations, results, etc.
# these are based on SQLobjects, so their properties are all stored in a 
# relational database (SQLite). 
# Also define operators such as Add and Subtract,
# eg. to find energy changes subtract one species from another, and query the energy of the resulting "DerivedSpecies"

import getopt, sys, os, re

from sqlobject import *

""" set up the database """

class DerivedSpecies():
    def __init__(self,name=None):
        self.name=name
        self.calculations=list()
        self.formula=None
    def __add__(self,other):
        answer=DerivedSpecies(name="%s + %s"%(self.name,other.name))
        for c in self.calculations:
            for d in other.calculations:
                if c.sameAs(d):
                    answer.calculations.append(c+d)
        return answer
    def __sub__(self,other):
        answer=self+(-other)
        return answer
    def __neg__(self):
        answer=DerivedSpecies(name="- %s"%self.name)
        for c in self.calculations:
            answer.calculations.append(-c)
        return answer
    def __mul__(self,other):
        try:    
            factor=float(other)
        except (ValueError,TypeError):
            return NotImplemented
        answer=DerivedSpecies(name="%.4g * %s"%(other,self.name))
        for c in self.calculations:
            answer.calculations.append(c*factor)
        return answer
    def __rmul__(self,other):
        return self*other 
    def findCalculation(self,copySettingsFrom):
        for c in self.calculations:
            if c.sameAs(copySettingsFrom):
                return c
        raise IndexError, "No calculation found with those settings"
        
    def latex(self):
        structurefile=os.path.join(os.curdir,'structures', self.name+'.pdf')
        if os.path.exists(structurefile):
            latexthis="\includegraphics[scale=1.5]{%s.pdf} "%self.name
        else:
            latexthis=re.sub(r'(\d)',r'$_\1$',self.name) # subscript numbers
        return latexthis

class Species(SQLObject,DerivedSpecies):
    name=StringCol(alternateID=True,
                   unique=True,
                   notNone=True)
    reactions=SQLRelatedJoin(otherClass='Reaction',
                          joinColumn='species_id', 
                          otherColumn='reaction_id', 
                          intermediateTable='reactant', 
                          createRelatedTable=False)
    calculations=SQLMultipleJoin('Calculation')
    formula=StringCol(length=20,
                      default='')

class Reaction(SQLObject,DerivedSpecies):  # inherit arithmetics from DerivedSpecies
    name=StringCol(alternateID=True,
                   unique=True,
                   notNone=True)
    species=SQLRelatedJoin( otherClass='Species',
                            joinColumn='reaction_id', 
                            otherColumn='species_id', 
                            intermediateTable='reactant', 
                            orderBy='stoichiometry',
                            createRelatedTable=False)
    calculations=SQLMultipleJoin('Calculation')
    
    def makeFromReactantsProducts(self, reactants, products):
        # give it a list of reactant names and a list of product names
        # assume stoichiometry is one for each
        for reactantname in reactants:
            if reactantname=='ts': 
                continue # discard the 'ts_' at the beginning of the filename
                print "throwing away the 'reactant' called 'ts'"
            try: 
                species=Species.byName(reactantname)
            except SQLObjectNotFound:
                print "species %s doesn't exist yet. Creating now."%reactantname
                species=Species(name=reactantname)
            i=Reactant(stoichiometry=-1,
                       species=species, # this assumes it already exists
                       reaction=self)
            print "adding reactant %s to reaction %d"%(reactantname,self.id)
        for productname in products:
            try: 
                species=Species.byName(productname)
            except SQLObjectNotFound:
                print "species %s doesn't exist yet. Creating now."%productname
                species=Species(name=productname)
            i=Reactant(stoichiometry=1,
                       species=species,
                       reaction=self)
            print "adding product %s to reaction %d"%(productname,self.id)
    
    
class Reactant(SQLObject):
    stoichiometry=IntCol()
    species=ForeignKey('Species')
    reaction=ForeignKey('Reaction')


class DerivedCalculation():
    """Like a calculation, but derived from a calculation(s) and not stored in database"""
    def __init__(self,copySettingsFrom=None, jobType='Mixed', method='CBS-QB3', basisSet='CBS-QB3'):
        self.species=None
        self.reaction=None
        self.jobType=jobType
        self.method=method
        self.basisSet=basisSet
        self.commandLine=''
        self.results=list()
        if copySettingsFrom:
            try:
                self.jobType=copySettingsFrom.jobType
                self.method=copySettingsFrom.method
                self.basisSet=copySettingsFrom.basisSet
            except AttributeError:
               print "can't copy settings from ",copySettingsFrom
               raise
    def __repr__(self):
        return "<DerivedCalculation jobType='%s' method='%s' basisSet='%s' >"%(self.jobType, self.method, self.basisSet)
    def sameAs(self,other):
        if (re.sub('^[RU]','(U/R)',self.method)==re.sub('^[RU]','(U/R)',other.method)
        and self.jobType==other.jobType
        and self.basisSet==other.basisSet):
            return True
        return False
    def __getitem__(self,key):
        for res in self.results:
            if res.name==key:
                return res # return whole DerivedResult object
        raise IndexError
    def __setitem__(self,key,value):
        # if it's already there then change it
        for res in self.results:
            if res.name==key:
                res.value=value
                return
        # not already there so add it
        self.results.append(DerivedResult(name=key,value=value))
    def setResult(self,result):
        self[result.name]=result.value
    def __add__(self,other):
        assert self.sameAs(other), "tried to sum non-equivalent calculations"
        answer=DerivedCalculation(copySettingsFrom=self)
        for r in self.results:
            try:
                answer.setResult(r+other[r.name])
            except IndexError:
                #print "Couldn't find %s in %s"%(r.name, other)
                pass
        return answer
    def __sub__(self,other):
        answer=self+(-other)
        return answer
    def __neg__(self):
        answer=DerivedCalculation(copySettingsFrom=self)
        for r in self.results:
            answer.setResult(-r)
        return answer
    def __mul__(self,other):
        try:    
            factor=float(other)
        except (ValueError,TypeError):
            return NotImplemented
        answer=DerivedCalculation(copySettingsFrom=self)
        for r in self.results:
            answer.setResult(r*factor)
        return answer
    def __rmul__(self,other):
        return self*other

class Calculation(SQLObject,DerivedCalculation):
    species=ForeignKey('Species', default=None)
    reaction=ForeignKey('Reaction', default=None)
    jobType=StringCol(length=20)
    method=StringCol(length=10)
    basisSet=StringCol(length=20)
    commandLine=StringCol()
    #startTime=DateTimeCol()
    results=SQLMultipleJoin('Result')
    # inherit __add__, __sub__, __neg__ from DerivedCalculation, 
    # over-write __getitem__ and __setitem__ here 
    def __getitem__(self,key):
        res=self.results.filter(Result.q.name==key)
        count=res.count()
        if count<1:
            raise IndexError, "no result with name %s"%key
        if count>1:
            raise IndexError, "%s has %d results"%(key,count)
        return res[0]  # return the whole Result object, not just its value
    def __setitem__(self,key,value):
        raise Exception, "probably shouldn't be changing things in the database here!"
        # also, the __setitem__ that would be inherited from DerivedCalculation will 
        # only work here if the item already exists. New items would not be added correctly.
        
        
class DerivedResult():
    def __init__(self,name=None,value=None):
        self.name=name
        self.value=value
    def __repr__(self):
        return "DerivedResult(name='%s',value='%s')"%(self.name,self.value)
    def __add__(self,other):
        answer=DerivedResult(name=self.name)
        try:
            answer.value=float(self.value)+float(other.value)
        except (ValueError,TypeError):
            #print "Couldn't convert '%s' or '%s' to float. Making list"%(self.value,other.value)
            answer.value=[self.value, other.value ]
        return answer
    def __sub__(self,other):
        answer=self+(-other)
        return answer
    def __neg__(self):
        answer=DerivedResult(name=self.name)
        try:
            answer.value=-float(self.value)
        except (ValueError,TypeError):
            # print "Couldn't convert to float. Leaving alone"
            answer.value=self.value
        return answer        
    def __mul__(self,other):
        answer=DerivedResult(name=self.name)
        try:
            answer.value=other*float(self.value)
        except (ValueError,TypeError):
            # print "Couldn't convert to float. Leaving alone"
            answer.value=self.value
        return answer
    def __rmul__(self,other):
        return self*other
        
class Result(SQLObject,DerivedResult): #  inherit __add__,__sub__,__neg__ from DerivedResult
    calculation=ForeignKey('Calculation')  
    name=StringCol(length=50)
    value=StringCol() 
#    item=StringCol(length=50)
#    units=StringCol(length=20)
#    value=FloatCol()
#    uncertainty=FloatCol()
        


def connectToDatabase(reset=False,debug=False):
    db_filename = os.path.abspath('data.db')
    if reset:
#        try:
#            connection=sqlhub.getConnection()
#            connection.close
#            print "trying to close connection"
#        except AttributeError: 
#            print "couldn't close connection (maybe there is none)"
        if os.path.exists(db_filename):
            os.unlink(db_filename)
    connection_string = 'sqlite:' + db_filename
    if debug: 
        connection_string+='?debug=true'
    connection = connectionForURI(connection_string)
    sqlhub.processConnection = connection # The sqlhub.processConnection assignment means that all classes will, by default, use this connection we've just set up.
    
    if reset:
        Species.dropTable(ifExists=True)
        Reaction.dropTable(ifExists=True)
        Reactant.dropTable(ifExists=True)
        Calculation.dropTable(ifExists=True)
        Result.dropTable(ifExists=True)
    
    Species.createTable(ifNotExists=True)
    Reaction.createTable(ifNotExists=True)
    Reactant.createTable(ifNotExists=True)
    Calculation.createTable(ifNotExists=True)
    Result.createTable(ifNotExists=True)

    nSpecies=Species.select().count()
    nReaction=Reaction.select().count()
    nReactant=Reactant.select().count()
    nCalculation=Calculation.select().count()
    nResult=Result.select().count()
    print "I have the following in the database:"
    print " %d Species \n %d Reactions \n %d Reactants \n %d Calculations \n %d Results"%(nSpecies,nReaction,nReactant,nCalculation,nResult)


def dreloadall():
    dreload(Species)
    dreload(Reaction)
    dreload(Reactant)
    dreload(DerivedCalculation)
    dreload(Calculation)
    dreload(DerivedResult)
    dreload(Result)
    dreload(connectToDatabase)
    
    
def main():
    print "nothing"
    ## USAGE:
    rx=Reaction(name='Methane decomposition')
    ch4=Species(name="Methane")
    Reactant(species=ch4,reaction=rx,stoichiometry=1)
    assert rx.species[0].species==ch4
    assert ch4.reactions[0].reaction==rx


if __name__ == "__main__":
    main()
    