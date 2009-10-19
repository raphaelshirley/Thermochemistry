#!/usr/bin/env python
# encoding: utf-8
"""
gaussian_settings_changer.py
version 0.10

Created by Yaoyao Liu on 2009-02-22.
Copyright (c) 2009 . All rights reserved.

Version history:
0.10 - initial release. Changes gaussian settings and ouputs new job file.
 --- THIS IS STILL INCOMPLETE AND DOES NOT WORK. SORRY.


"""

import sys
import os
import fileinput
import pickle

import openbabel, pybel
from optparse import OptionParser


def readlogfile(filename):
	filedata = fileinput.input(filename)
	return filedata
	
	
# Periodic Table lookup: element <-> number
class PeriodicTable(object):
	element = [None,
		'H', 'He',
		'Li', 'Be',
		'B', 'C', 'N', 'O', 'F', 'Ne',
		'Na', 'Mg',
		'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
		'K', 'Ca',
		'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
		'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
		'Rb', 'Sr',
		'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
		'In', 'Sn', 'Sb', 'Te', 'I', 'Xe',
		'Cs', 'Ba',
		'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
		'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
		'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn',
		'Fr', 'Ra',
		'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No',
		'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Uub']
	number = {}
	for n in range(1, len(element)):
		number[element[n]] = n

# User variables
settingsfile = 'gaussian_settings.txt'

# Constants
#N_a = 6.0221417930E23
#E_h = 4.3597441775E-18

# Parse cmd line optionsf
usage = "%prog [options] filename.gau [filenames.gau]"
parser = OptionParser(usage, version="%prog 0.10 by Yaoyao Liu")




parser.add_option("-w", "--write",
					action="store_true", dest="write", default=False,
					help="write new Gaussian input files")
parser.add_option("-m", "--mol",
					action="store_true", dest="mol", default=False,
					help="use .mol file information")
parser.add_option("-t", "--toggle-spin",
					action="store_true", dest="spin", default=False,
					help="toggle spin multiplcity between singlet/triplet, or doublet/quartet (requires .mol file)")
parser.add_option("--settings-file",
					action="store", dest="settings", metavar="file name", default=settingsfile,
					help="name of file with new Gaussian settings template (default: gaussian_settings.txt)")
parser.add_option("-f", "--name-format",
					action="store", dest="format", metavar="fmt", default="smi",
					help="alternative name string format (default: smi)")
parser.add_option("-q", "--quiet",
					action="store_true", dest="quiet", default=False,
					help="do not print output data")

(options, args) = parser.parse_args()
# Old MySQL stuff. Redundant now.
#conn = MySQLdb.connect (host = "127.0.0.1",
#                          user = "prime",
#                          passwd = "m4sterpl4n",
#                          db = "prime");

print '[Gaussian input settings changer script by Yaoyao Liu. See -h for options]'
print ''

# Not needed with optparse module
#args = sys.argv[1:]

newsettings = open(settingsfile, 'r')

# Begin main loop
for filename in args:

	inputfile = readlogfile(filename)
	print  '----------------------------------'
	print '-=> Processing file', filename, '...'
	
	filenamestub = (filename.rpartition('/')[2]).rpartition('.')[0]
	pathfulenamestub = filename.rpartition('.')[0]

###### START PARSING MOL FILE #####
	if options.mol == True or options.spin == True:

		# Parse file data into OpenBabel
		filenamestub = (filename.rpartition('/')[2]).rpartition('.')[0]
		pathfulenamestub = filename.rpartition('.')[0]
		molfilename = pathfulenamestub + '.mol'
		try:
			obdata = pybel.readfile('mol', molfilename).next()
		except:
			print 'ERROR: Could not open .mol file.'
			exit()

		# Correct for unsatisfied valence
		totalspin = 0
		for atom in obdata.atoms:
		
			# Count number of hydrogens needed to fill the implicit valence of this atom
			overvalence = atom.OBAtom.ImplicitHydrogenCount()
			if overvalence:
				atom.OBAtom.SetSpinMultiplicity(overvalence + 1)
				totalspin = totalspin + overvalence

		# Check that we have satisfied all valences explicitly.
		implicithydrogen = 0
		for atom in obdata.atoms:
			implicithydrogen = implicithydrogen + atom.OBAtom.ImplicitHydrogenCount()
		if implicithydrogen != 0:
			print 'Error in spin multiplicity calculations. Aborting now. DEBUG.'
			exit()
	
		# Get spaced formula in order to count electrons
		formulasp = obdata.OBMol.GetSpacedFormula(0,' ').split()

		# Define periodic table lookup
		pt = PeriodicTable()
	
		# Count total number of electrons in molecule
		i = 0
		electrons = 0
		while i < len(formulasp):
			electrons = electrons + (pt.number[formulasp[i]] * int(formulasp[i+1]))
			i = i + 2
	
#	mult2 = obdata.OBMol.GetTotalSpinMultiplicity()

		if not electrons % 2:
			possiblemult = [1,3]
		else:
			possiblemult = [2,4]
	
		# Get SMILES id from OB
		name = obdata.write(options.format).split()[0]
		# Get additioanl data from OB
		mass = obdata.OBMol.GetExactMass()
		formula = obdata.OBMol.GetSpacedFormula(1,'')

###### END PARSING MOL FILE #####


###### START PARSING GAUSSIAN INPUT FILE #####
	# Reset Gaussian input file variables for each file before parsing
	functionalbasisfound = 0
	functional = ''
	basis = ''
	gaucalc = ''
	(charge, mult) = ('unknown', 'unknown')
	
	# Get data from gaussian input file:
	for line in inputfile:
		
		# Match input settings section only
		while line[0:1] == '#':
			
			# Find the calculation method/functional and basis set
			functionalbasischeck = line.partition('/')
			
			if functionalbasischeck[1] == '/' and functionalbasisfound == 0:
				functional = functionalbasischeck[0].split()[-1]
				basis = functionalbasischeck[2].split()[0]
				
				# Flag that we've already found it and won't look again.
				functionalbasisfound = 1
		
			# Determine if geometry optimisation is performed
			if line.find('opt') >= 0:
				if gaucalc == '':
					gaucalc += 'optimize'
				else:
					gaucalc += ', optimize'
			
			# Determine if frequency analysis is performed
			if line.find('freq') >= 0:
				if gaucalc == '':
					gaucalc += 'freq'
				else:
					gaucalc += ', freq'
		
			# Determine integration grid precision
			if line.find('Integral') >= 0:

				# Common cases where, for example, Grid=Fine or FineGrid is used.
				if line.find('UltraFine') >= 0:
					precision = 'UltraFine'
				elif line.find('Fine') >= 0:
					precision = 'Fine'
				elif line.find('Coarse') >= 0:
					precision = 'Coarse'

				# Catch all case with Grid=xxxxxx
				else:
					precision = line[line.find('Integral'):].partition('(')[2].partition(')')[0].rpartition('=')[2]
				
					# If we still can't determine the integration grid precision... then give up
					if precision == '':
						precision = 'unknown'

			line = inputfile.next()

		# Checking for charge and multiplcity
		chargemultcheck = line.split()
		# Try to see if the line contains only two integers <-->  it's the charge and multiplicity line
		if len(chargemultcheck) == 2:
			try:
				charge = int(chargemultcheck[0])
				mult = int(chargemultcheck[1])
			except:
				pass

###### END PARSING GAUSSIAN INPUT FILE #####

	# Printing output
	if options.quiet == False:
		if options.mol == True:
			print  '--------- Species Summary --------'
			print  'Formula: ', formula
			print  'Name string: ', name
			print  'Molecular mass (amu): ', mass
			print  'Number of electrons: ', electrons
#		print  'Total atomic spin multiplicity: ', mult2
		print  '------- Old Input Parameters ------'
		print  'Functional: ', functional
		print  'Basis set: ', basis
		print  'Calculation: ', gaucalc
		print  'Precision: ', precision
		print  'Spin multiplicity: ', mult
		print  '-----------------------------------'


	# Write new Gaussian input file if specified
	if options.write == True:
	
###### START PARSING GAUSSIAN SETTINGS FILE #####
	# Reset Gaussian input file variables for each file before parsing
		newfunctionalbasisfound = 0
		newfunctional = ''
		newbasis = ''
		newgaucalc = ''
		(newcharge, newmult) = ('unknown', 'unknown')

		# Get data from gaussian input file:
		for line in newsettings:

			# Find the calculation method/functional and basis set
			newfunctionalbasischeck = line.partition('/')

			if newfunctionalbasischeck[1] == '/' and newfunctionalbasisfound == 0:
				newfunctional = newfunctionalbasischeck[0].split()[-1]
				newbasis = newfunctionalbasischeck[2].split()[0]

				# Flag that we've already found it and won't look again.
				newfunctionalbasisfound = 1

			# Determine if geometry optimisation is performed
			if line.find('opt') >= 0:
				if newgaucalc == '':
					newgaucalc += 'optimize'
				else:
					newgaucalc += ', optimize'

			# Determine if frequency analysis is performed
			if line.find('freq') >= 0:
				if newgaucalc == '':
					newgaucalc += 'freq'
				else:
					newgaucalc += ', freq'
					
			# Determine integration grid precision
			if line.find('Integral') >= 0:
				
				# Common cases where, for example, Grid=Fine or FineGrid is used.
				if line.find('UltraFine') >= 0:
					newprecision = 'UltraFine'
				elif line.find('Fine') >= 0:
					newprecision = 'Fine'
				elif line.find('Coarse') >= 0:
					newprecision = 'Coarse'

				# Catch all case with Grid=xxxxxx
				else:
				 	newprecision = line[line.find('Integral'):].partition('(')[2].partition(')')[0].rpartition('=')[2]

					# If we still can't determine the integration grid precision... then give up
					if newprecision == '':
						newprecision = 'unknown'


		# Toggle spin multiplicity if specified
		if options.spin == True:
			if mult == 1 or mult == 2:
				newmult = possiblemult[1]
			if mult == 3 or mult == 4:
				newmult = possiblemult[0]
		else:
			newmult = mult

			# Write new filename systematically
			newfilenamestub =  filenamestub.partition('_')[0] + "_spin=%i_" %(newmult) + newfunctional + '_' + newbasis + '_' + newprecision 
		
		print 'Writing to new Gaussian input file', newfilenamestub + '.gau' , '...' 
		print  '------- New Input Parameters ------'
		print  'Functional: ', newfunctional
		print  'Basis set: ', newbasis
		print  'Calculation: ', newgaucalc
		print  'Precision: ', newprecision
		print  'Spin multiplicity: ', newmult
		print  '-----------------------------------'
	
			
# Cleanup
#if options.write == True:
