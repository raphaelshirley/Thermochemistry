#!/usr/bin/python
######################################################
####      Gaussian Output Files to Mol Files      ####
######################################################
####    Version 1 (3-10-2008) Author: T Totton    ####
######################################################

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

## Change file from .g03 to .mol file using pybel

for filename in dirlist:

    if filename.endswith(".g03"):
        path=filename.rpartition('.')[0]
        outfile=filename.rpartition('/')[2]
        name=outfile.rpartition('.')[0]
        molefilename=name + ".mol"

        print molefilename
        
        mol = pybel.readfile("g03", filename).next()
        mol.write("mol", molefilename)
