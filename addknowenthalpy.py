# Written by Markus Sander
# Modified by Tim Totton (5/03/2009)
# Modified by Tim Totton (11/03/2009)

#####################################################################################################
# User inputs

# Known enthalpy of formation (298K) MySQL table (written in kJ/mol)
knownentSQL = 'Si_knownent_B971'
#####################################################################################################

# Import libraries 
import MySQLdb;
import os.path
import string;
from os import system 

# Import the library           export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib/
import openbabel, pybel
import private  # Database info

def dir_list(dir_name):
    outputList = []
    for root, dirs, files in os.walk(dir_name):
        outputList.append(root)
        for d in dirs:
            outputList.append('/'.join([root, d]))
        for f1 in files:
            outputList.append('/'.join([root, f1]))
    return outputList


conn = MySQLdb.connect (host = private.defaulthost,
                        user = private.defaultuser,
                        passwd = private.defaultpasswd,
                        db = private.defaultdb);
cursor = conn.cursor ()
dirlist=dir_list(".")
units=raw_input("Units kcal/mol or kJ/mol? >")
for filename in dirlist:
    if filename.endswith(".mol") and '\known' in filename or filename.endswith(".mol") and '\knownmp2' in filename:
        mol = pybel.readfile("mol", filename).next()
        totalspin=0
        for atom in mol.atoms:
            valenceunsatisfied=atom.OBAtom.ImplicitHydrogenCount()
            if valenceunsatisfied:
                atom.OBAtom.SetSpinMultiplicity(valenceunsatisfied+1)
                        # multiplicity of 2 if 1 missing H atom
                totalspin+=valenceunsatisfied
                #reset formula
        mol.OBMol.SetFormula( mol.OBMol.GetSpacedFormula(1,'') )
        smile=mol.write('smi').strip()

        name=smile.split("\t")[0]
        sql="select * from %s where name='%s'" %(knownentSQL,name)
        cursor.execute(sql)
        if cursor.rowcount>0:
            print 'Species %s exists already in database' %(filename)
        else:
            print 'H(298) of' +name+ '\nFilename= '+filename        
            
            if units == "kJ/mol":
                enthalpy=float(raw_input("Value >"))
                sql="INSERT into %s (Name, Enthalpy, path) VALUES ('%s','%f','%s')" %(knownentSQL,name,float(enthalpy),filename)
                cursor.execute(sql)
            elif units == "kcal/mol":
                enthalpykcal=float(raw_input("Value >"))
                enthalpykJ=enthalpykcal*4.18400
                sql="INSERT into %s (Name, Enthalpy, path) VALUES ('%s','%f','%s')" %(knownentSQL,name,float(enthalpykJ),filename)
                cursor.execute(sql)
            else:
                print "unknown units"
                
#cursor.execute ("SELECT primeID,cas FROM speciescas ")
#while (1):
# row = cursor.fetchone ()
# if row == None:
#  break
# print "%s, %s" % (row[0], row[1])
# file=row[0]+".struct"
# if (os.path.exists(file)):
#  print "file %s exists" % (file)
# else:
#  get="wget -O"+row[0]+".struct http://webbook.nist.gov/cgi/cbook.cgi?Str2File=C"+row[1]
#  print "%s" % (get)
# file=row[0]+".png"
# if (os.path.exists(file)):
#  print "file %s exists" % (file)
# else:
#  get="wget -O"+row[0]+".png http://webbook.nist.gov/cgi/cbook.cgi?Struct=C"+row[1]
#print "Number of rows returned: %d" % cursor.rowcount


