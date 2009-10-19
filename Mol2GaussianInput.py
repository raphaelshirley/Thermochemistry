#!/usr/bin/python
##########################################################
####         Mol Files to Gaussian Input Files        ####
##########################################################
####      Version 1 (3-10-2008) Author: T Totton      ####
##########################################################
####        Assumes:  Zero charge on molecules        ####
####      Spin States from SpinStates.txt (ANSI)      ####
####    Format: Species Name, Multiplicity (number)   ####
##########################################################

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

## Change .mol file to .gau file

for filename in dirlist:

    if filename.endswith(".mol"):
        path=filename.rpartition('.')[0]
        outfile=filename.rpartition('/')[2]
        name=outfile.rpartition('.')[0]
        molefilename=name + ".gau"

        mol = pybel.readfile("mol", filename).next()
        mol.write("gau", molefilename)

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
            
