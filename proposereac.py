
# Markus Sander - Juni 2008 


from math import *
import re
import glob
import os
import MySQLdb;

conn = MySQLdb.connect (host = "127.0.0.1",
                           user = "prime",
                           passwd = "m4sterpl4n",
                           db = "prime");
cursor = conn.cursor ()
periodicTable={ 1: 'H',  2: 'He',  3: 'Li',  4: 'Be',  5: 'B',  6: 'C',  7: 'N',  8: 'O',  9: 'F',  10: 'Ne',  11: 'Na',  12: 'Mg',  13: 'Al',  14: 'Si',  15: 'P',  16: 'S',  17: 'Cl',  18: 'Ar',  19: 'K',  20: 'Ca',  21: 'Sc',  22: 'Ti',  23: 'V',  24: 'Cr',  25: 'Mn',  26: 'Fe',  27: 'Co',  28: 'Ni',  29: 'Cu',  30: 'Zn',  31: 'Ga',  32: 'Ge',  33: 'As',  34: 'Se',  35: 'Br',  36: 'Kr',  37: 'Rb',  38: 'Sr',  39: 'Y',  40: 'Zr',  41: 'Nb',  42: 'Mo',  43: 'Tc',  44: 'Ru',  45: 'Rh',  46: 'Pd',  47: 'Ag',  48: 'Cd',  49: 'In',  50: 'Sn',  51: 'Sb',  52: 'Te',  53: 'I',  54: 'Xe',  55: 'Cs',  56: 'Ba',  57: 'La',  58: 'Ce',  59: 'Pr',  60: 'Nd',  61: 'Pm',  62: 'Sm',  63: 'Eu',  64: 'Gd',  65: 'Tb',  66: 'Dy',  67: 'Ho',  68: 'Er',  69: 'Tm',  70: 'Yb',  71: 'Lu',  72: 'Hf',  73: 'Ta',  74: 'W',  75: 'Re',  76: 'Os',  77: 'Ir',  78: 'Pt',  79: 'Au',  80: 'Hg',  81: 'Tl',  82: 'Pb',  83: 'Bi',  84: 'Po',  85: 'At',  86: 'Rn',  87: 'Fr',  88: 'Ra',  89: 'Ac',  90: 'Th',  91: 'Pa',  92: 'U',  93: 'Np',  94: 'Pu',  95: 'Am',  96: 'Cm',  97: 'Bk',  98: 'Cf',  99: 'Es',  100: 'Fm',  101: 'Md',  102: 'No',  103: 'Lr',  104: 'Rf',  105: 'Db',  106: 'Sg',  107: 'Bh',  108: 'Hs',  109: 'Mt',  110: 'Ds',  111: 'Rg',  112: 'Uub',  113: 'Uut',  114: 'Uuq',  115: 'Uup',  116: 'Uuh',  117: 'Uus',  118: 'Uuo'}


#returns -1 if list2 is not in list1 one, otherwise returns the atoms left in list1 after removing the atoms  in list2
def list2inlist1(list1, list2):
    templist=list(i for i in list2)
   # templist1=list1
    result=[]
    for elem in list1:
        if elem in templist:
            templist.remove(elem)
        else:
            result.append(elem)
    if templist==[]:
        return result
    else:
        return -1
    


periodicTable={ 1: 'H',  2: 'He',  3: 'Li',  4: 'Be',  5: 'B',  6: 'C',  7: 'N',  8: 'O',  9: 'F',  10: 'Ne',  11: 'Na',  12: 'Mg',  13: 'Al',  14: 'Si',  15: 'P',  16: 'S',  17: 'Cl',  18: 'Ar',  19: 'K',  20: 'Ca',  21: 'Sc',  22: 'Ti',  23: 'V',  24: 'Cr',  25: 'Mn',  26: 'Fe',  27: 'Co',  28: 'Ni',  29: 'Cu',  30: 'Zn',  31: 'Ga',  32: 'Ge',  33: 'As',  34: 'Se',  35: 'Br',  36: 'Kr',  37: 'Rb',  38: 'Sr',  39: 'Y',  40: 'Zr',  41: 'Nb',  42: 'Mo',  43: 'Tc',  44: 'Ru',  45: 'Rh',  46: 'Pd',  47: 'Ag',  48: 'Cd',  49: 'In',  50: 'Sn',  51: 'Sb',  52: 'Te',  53: 'I',  54: 'Xe',  55: 'Cs',  56: 'Ba',  57: 'La',  58: 'Ce',  59: 'Pr',  60: 'Nd',  61: 'Pm',  62: 'Sm',  63: 'Eu',  64: 'Gd',  65: 'Tb',  66: 'Dy',  67: 'Ho',  68: 'Er',  69: 'Tm',  70: 'Yb',  71: 'Lu',  72: 'Hf',  73: 'Ta',  74: 'W',  75: 'Re',  76: 'Os',  77: 'Ir',  78: 'Pt',  79: 'Au',  80: 'Hg',  81: 'Tl',  82: 'Pb',  83: 'Bi',  84: 'Po',  85: 'At',  86: 'Rn',  87: 'Fr',  88: 'Ra',  89: 'Ac',  90: 'Th',  91: 'Pa',  92: 'U',  93: 'Np',  94: 'Pu',  95: 'Am',  96: 'Cm',  97: 'Bk',  98: 'Cf',  99: 'Es',  100: 'Fm',  101: 'Md',  102: 'No',  103: 'Lr',  104: 'Rf',  105: 'Db',  106: 'Sg',  107: 'Bh',  108: 'Hs',  109: 'Mt',  110: 'Ds',  111: 'Rg',  112: 'Uub',  113: 'Uut',  114: 'Uuq',  115: 'Uup',  116: 'Uuh',  117: 'Uus',  118: 'Uuo'}

import openbabel, pybel
# please cite:
# Pybel: a Python wrapper for the OpenBabel cheminformatics toolkit
# Noel M O'Boyle, Chris Morley and Geoffrey R Hutchison
# Chemistry Central Journal 2008, 2:5
# doi:10.1186/1752-153X-2-5

sql="select id,molefile from Thermoname where project='silver'"
cursor.execute(sql)
if not os.path.isdir("temp"):
   os.mkdir("temp")
while (1):
    row=cursor.fetchone()
    if row == None:
        break
    tempfile="temp/%s.mol" %row[0]
    f=open(tempfile, 'w')
    f.write(row[1])
 #   print row[0]
#    print row[1]
    f.close()

dirlist=os.listdir("temp/")

speciestable={}
atoms={}
speciesid={}
for filename in dirlist:
 if filename.endswith(".mol"):
    atomlist=[]
    path=filename.rpartition('.')[0]
    filename="temp/"+path + ".mol"
#    print filename
    # read in the first molecule in a mol file
    mol = pybel.readfile("mol", filename).next()


    #correct H
    totalspin=0
    for atom in mol.atoms:
        valenceunsatisfied=atom.OBAtom.ImplicitHydrogenCount()
        if valenceunsatisfied:
            atom.OBAtom.SetSpinMultiplicity(valenceunsatisfied+1)
            # multiplicity of 2 if 1 missing H atom
            totalspin+=valenceunsatisfied

 #   print "  sum of spins: %d" % totalspin

    #reset formula
    mol.OBMol.SetFormula( mol.OBMol.GetSpacedFormula(1,'') )



    
#    mol.OBMol.AddHydrogens()
    speciesid[mol.OBMol.GetFormula()]=path
    speciestable[mol.OBMol.GetFormula()]=mol
   # print os.path.basename(filename),
    smile=mol.write('smi').strip()
    for atom in mol.atoms:
       atomlist.append(periodicTable[atom.OBAtom.GetAtomicNum()])
    atoms[mol.OBMol.GetFormula()]=atomlist
#search for unimolecular reactions
print "\n\nUnimolecular reactions:\n"
enthalpies=[]
for s in speciestable:

     for s2 in speciestable:
         if not s2==s:
             reslist2inlist1=list2inlist1(atoms[s],atoms[s2])
             if (not reslist2inlist1==-1):
                # print reslist2inlist1
                 for s3 in speciestable:
                     if list2inlist1(reslist2inlist1,s3)==[]:
                            
                            sql="select Energy from Thermoname where ID='%s'" %speciesid[s]
                            cursor.execute(sql)
                            row=cursor.fetchone()
                            energys=row[0]
                            sql="select Energy from Thermoname where ID='%s'" %speciesid[s2]
                            cursor.execute(sql)
                            row=cursor.fetchone()
                            energys2=row[0]
                            sql="select Energy from Thermoname where ID='%s'" %speciesid[s3]
                            cursor.execute(sql)
                            row=cursor.fetchone()
                            energys3=row[0]
                            dH=(float(energys)-float(energys2)-float(energys3))*2625.5
                            if dH not in enthalpies and -dH not in enthalpies:
                                
                                enthalpies.append(dH)
                                print  s+"="+s2 +"+"+ s3 +"     dH=%f kj/mol" %dH
                              #  sql="insert into reactions (Reac1, Reac2, Prod1, Prod2, A, n, H) values ('%s', '', '%s', '%s', '10000000000000', '0', '%f')" % (speciesid[s],speciesid[s2],speciesid[s3],dH)
                              #  cursor.execute(sql)
                                
print "\n\nBimolecular reactions:\n"
enthalpies=[]
for s in speciestable:
    for s1 in speciestable:
        # copy s to sumlist
        sumlist=templist=list(i for i in atoms[s])
        # add s1
        sumlist.extend(atoms[s1])

        for s2 in speciestable:
         if not s2==s and not s1==s2:
             reslist2inlist1=list2inlist1(sumlist,atoms[s2])
             if (not reslist2inlist1==-1):
                # print reslist2inlist1
                 for s3 in speciestable:
                     if not s3==s1 and not s3==s:
                         if list2inlist1(s3,reslist2inlist1)==[]:
                             
                          #  print list2inlist1(reslist2inlist1,s3) 
                            
                            sql="select Energy from Thermoname where ID='%s'" %speciesid[s]
                            cursor.execute(sql)
                            row=cursor.fetchone()
                            energys=row[0]
                            sql="select Energy from Thermoname where ID='%s'" %speciesid[s1]
                            cursor.execute(sql)
                            row=cursor.fetchone()
                            energys1=row[0]
                            sql="select Energy from Thermoname where ID='%s'" %speciesid[s2]
                            cursor.execute(sql)
                            row=cursor.fetchone()
                            energys2=row[0]
                            sql="select Energy from Thermoname where ID='%s'" %speciesid[s3]
                            cursor.execute(sql)
                            row=cursor.fetchone()
                            energys3=row[0]
                            dH=(float(energys)+float(energys1)-float(energys2)-float(energys3))*2625.5
                            if dH not in enthalpies and -dH not in enthalpies:
                                print s+"+"+s1+"="+s2 +"+"+ s3+"     dH=%f kj/mol" %dH
                                enthalpies.append(dH)
                                #sql="insert into reactions (Reac1, Reac2, Prod1, Prod2, A, n, H) values ('%s', '%s', '%s', '%s', '10000000000000', '0', '%f')" % (speciesid[s],speciesid[s1],speciesid[s2],speciesid[s3],dH)
                                #cursor.execute(sql)


        

