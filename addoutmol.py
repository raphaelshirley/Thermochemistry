import MySQLdb;
import os.path
import string;
from os import system 

# Import the library           export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib/
import openbabel, pybel

def dir_list(dir_name):
    outputList = []
    for root, dirs, files in os.walk(dir_name):
        outputList.append(root)
        for d in dirs:
            outputList.append('/'.join([root, d]))
        for f1 in files:
            outputList.append('/'.join([root, f1]))
    return outputList






conn = MySQLdb.connect (host = "127.0.0.1",
                           user = "prime",
                           passwd = "m4sterpl4n",
                           db = "prime");
cursor = conn.cursor ()
#dirlist=raw_input("Directory?")
#dirlist=dir_list(".")
dirlist=os.listdir(".")
for filename in dirlist:
   # print filename
    if filename.endswith(".outmol"):
        path=filename.rpartition('.')[0]
        outfile=filename.rpartition('/')[2]
        molefilename=path + ".mol"
        if  not os.path.exists(molefilename):
         print 'Molefile for %s' %(filename)
         molefilename=raw_input(">")

        mol = pybel.readfile("mol", molefilename).next()

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
        print name
        molefile=open(molefilename,'r')
#       print [f for f in files if f.endswith(".outmol")] 
        Efound=0
        thermofound=0
        thermoend=0
        freqfound=0
        calculate=''
        basis=''
        mass=0
        I=[0,0,0]
        spin='Auto'
        thermo=[0, 0, 0, 0, 0, 0]
        func=""
        sql="select * from Thermoname where  path='%s'" %(outfile)
        cursor.execute(sql)
        if cursor.rowcount>0:
            print 'Species %s exists already in database' %(outfile)
            thermofound=1
        #    break   
        else:
         print 'Adding %s to database' %outfile   
         file=open(filename, 'r')
         for line in file:
            if 'Calculate                     ' in line:
#              print line.split(" ")[20]
                for i in line.split(" "):
                    if i.isalpha():
                        calculate=i
                        
            if 'Spin                          ' in line:
 #               print line.split(" ")[20]
                for i in line.split(" "):
                    if i.isdigit():
                        spin = i
            if 'Basis                         ' in line:
#                print line.split(" ")[20]
                for i in line.split(" "):
                    if i.isalpha():
                        basis = i
            if '   Molecular Mass:' in line:
                mass=line.split(" ")[7]
                masssql="Update Thermoname set mass=%f where ID='%i'" %(float(mass),id)
                cursor.execute(masssql)
            if '    Eigenvalues --' in line:
               Ia=line.split(" ")[14]
               Ib=line.split(" ")[19]
               Ic=line.split(" ")[24]

                    
            if freqfound==1:
                      if 'Frequencies (cm-1) and normal modes' in line:
                            freqfound=0
                      else:
                       j=0
                #       print line
                #       print line.split("    ")[14]
                       for i in line.split(" "):
                                
                                if len(i)>0 :
                                  j=j+1
                                  if j==3:
                                      sql="INSERT into frequency (id,Name,frequency) VALUES ('%f','%s','%f')" %(id,name,float(i))
                                      cursor.execute(sql)
                                 # print '\n'
                                 # print j
                                #      print 'new freq \n'
                                  #freq=float(i)
                                 # print freq
         #
            if '  mode     au_amu        cm-1      km/mol' in line:
                freqfound=1
            if 'Functional                    ' in line:
                func= line.split(" ")[20]
                            
            if 'Ef' in line and Efound==0:
            #   print line 
                for i in  line.split("   "):
                    if 'Ha' in i:
                        E=float( i.replace('Ha',''))
                        break
            if 'Entering Properties' in line:
                Efound=1
                print outfile
                print "Please specify the name of the species"
                canname=raw_input(">")
                print "Please specify the project"
                project=raw_input(">")
                print "Any comments?"
                comments=raw_input(">")                
                sql="INSERT into Thermoname (Name,path, Energy, Functional,molefile,basis,spin,calculation,Ia,Ib,Ic,mass,project,canname,comments) VALUES ('%s','%s','%f','%s','%s','%s','%s','%s','%f','%f','%f','%f','%s','%s','%s')" %(name,outfile,E,func,molefile.read(),basis,spin,calculate,float(Ia),float(Ib),float(Ic),float(mass),project,canname,comments)
                cursor.execute (sql)
                
                sql="select id from Thermoname where path='%s'" %(outfile)
                cursor.execute(sql)
                row=cursor.fetchone()
                id=row[0]
            if '       T        Entropy   Heat_Capacity   Enthalpy   Free_Energy' in line:
                thermofound=1
            #   print line
            if thermofound<8 and thermofound >0:
                thermofound=thermofound+1
            elif thermofound>0 and thermoend==0: 
                if len(line)<3:
                    thermoend=1
                elif thermoend==0:
                    j=0
                    for i in  line.split(" "): 
                        if len(i)>0 and j<6:    
                            #print i
                            thermo[j]=float(i)
                            j=j+1
                    sql="INSERT into thermodata (id,name, temp, entropy, heatcapacity, enthalpy, freeenergy) VALUES ('%f','%s','%f','%f','%f','%f','%f')" %(id,name,thermo[1],thermo[2],thermo[3],thermo[4],thermo[5])
                    cursor.execute(sql) 
        if thermofound==0:
            print 'Thermo data not found, setting all values to 0'
            sql="INSERT into thermodata (id, temp, entropy, heatcapacity, enthalpy, freeenergy) VALUES ('%s','298.15','0','0','0','0')"%(id)
            cursor.execute(sql)
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


