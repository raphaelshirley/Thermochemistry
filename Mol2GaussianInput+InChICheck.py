#!/usr/bin/python
####################################################################
#### Mol Files to Gaussian Input Files + Deletion of Dulpicates ####
####################################################################
####          Version 1 (13-10-2008)  Author: T Totton          ####
####################################################################
####             Assumes:  Zero charge on molecules             ####
####           Spin States from SpinStates.txt (ANSI)           ####
####        Format: 'Species Name, Multiplicity (number)'       ####
####################################################################

import os.path
import string
from os import system 

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

dirlist=os.listdir(".")

## Change .mol file to .inchi file

for filename in dirlist:

    if filename.endswith(".mol"):
        path=filename.rpartition('.')[0]
        outfile=filename.rpartition('/')[2]
        name=outfile.rpartition('.')[0]
        inchifilename=name + ".inchi"

        mol = pybel.readfile("mol", filename).next()
        mol.write("inchi", inchifilename)

dirlist=os.listdir(".")

## Delete duplicate InChI strings (taking note of potential error
## when trying to open file which has been already deleted)

for filename in dirlist:

    if filename.endswith(".inchi"):
        k1=filename
        try:
            h1=open(k1,'r')
            j1=h1.read()
            
            for filename in dirlist:
                if filename.endswith(".inchi"):

                    k2=filename
                    try:
                        h2=open(k2,'r')
                        j2=h2.read()
                        if k1!=k2:
                            if j1==j2:
                                h2.close()
                                h1.close()
                                os.remove(k2)
                        h2.close()
                    except IOError:
                        pass
                    h2.close()
            h1.close()
        except IOError:
            pass
        h1.close()

dirlist=os.listdir(".")

## Change .mol file (related to saved .inchi file) to .gau file

for filename in dirlist:

    if filename.endswith(".inchi"):
        path=filename.rpartition('.')[0]
        outfile=filename.rpartition('/')[2]
        name=outfile.rpartition('.')[0]
        molfilename=name + ".gau"
        oldmol=path + ".mol"

        mol = pybel.readfile("mol", oldmol).next()
        mol.write("gau", molfilename)

dirlist=os.listdir(".")

## Add Gaussian Run Settings and Spin State

for filename in dirlist:
    
    if filename.endswith(".gau"):
        path=filename.rpartition('.')[0]
        outfile=filename.rpartition('/')[2]
        myname=outfile.rpartition('.')[0]

        infile=open(filename, 'r')

## Read in Spin State from .txt file
        
        spinstates=open('SpinStates.txt','r')

        for line in spinstates:
            if line.startswith(myname): 
                spinstate=line.rpartition(', ')[2]

## Delete initial opening lines to file given by OpenBabel
                
        f1=infile.readline()
        del f1
        f2=infile.readline()
        del f2
        f3=infile.readline()
        del f3
        f4=infile.readline()
        del f4
        f5=infile.readline()
        del f5
        g=infile.read()
            
        myoutfile=open(filename,'w')

## Input settings at top of file
        
        settings = ('#p ub3lyp/6-311+G(d,p) opt=(Tight, NewEstmFC) freq \n\
#GFInput Population=Regular \n\
#Integral(Grid=UltraFine) Guess=Mix NoSymmetry \n \n\
NewComb ub3lyp/6-311+G(d,p) \n \n\
 0  ')
        myoutfile.write(settings)      
        myoutfile.write(spinstate)
        myoutfile.write(g)
        spinstates.close()
        infile.close()
        myoutfile.close()
            
