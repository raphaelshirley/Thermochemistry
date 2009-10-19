#!/usr/bin/env python
# encoding: utf-8
"""
spin_state_calc.py
version 0.21

Created by Yaoyao Liu on 2009-02-22.
Copyright (c) 2009 . All rights reserved.

Version history:
0.10 - initial release. Calculates spin states and generates SpinStates.txt file.
0.11 - small bugfix and tweaks.
0.20 - added ability and option to generate SpinStates.txt from parsing Gaussian input files. Useful if you want to change the functional/basis set with existing converged geometries.
0.21 - made the output slightly prettier when multiple species are specified (skips a line between species)


"""

import sys
import os
import fileinput
import pickle
import MySQLdb

import openbabel, pybel
from optparse import OptionParser


def readlogfile(filename):
	filedata = fileinput.input(filename)
	return filedata


def serverversion():
	conn = MySQLdb.connect(host = options.mysql_host,
							user = options.mysql_user,
							passwd = options.mysql_passwd,
							db = options.mysql_db);
	cursor = conn.cursor()
	cursor.execute("SELECT VERSION()")
	row = cursor.fetchone()
	print "Server reports MySQL version", row[0]
	cursor.close()
	conn.close()

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
# none

# Constants
#N_a = 6.0221417930E23
#E_h = 4.3597441775E-18

# Parse cmd line optionsf
usage = "%prog [options] filename.mol [filenames.mol]"
parser = OptionParser(usage, version="%prog 0.21 by Yaoyao Liu")




parser.add_option("-w", "--write",
					action="store_true", dest="write", default=False,
					help="write SpinStates.txt and SpinStates-2.txt files from .mol files")
parser.add_option("-g", "--gaussian",
					action="store_true", dest="gaussian", default=False,
					help="write SpinStates.txt from Gaussian input files")
parser.add_option("-f", "--name-format",
					action="store", dest="format", metavar="fmt", default="smi",
					help="alternative name string format (default: smi)")
parser.add_option("-q", "--quiet",
					action="store_true", dest="quiet", default=False,
					help="do not print parsed data")

(options, args) = parser.parse_args()
# Old MySQL stuff. Redundant now.
#conn = MySQLdb.connect (host = "127.0.0.1",
#                          user = "prime",
#                          passwd = "m4sterpl4n",
#                          db = "prime");

print '[Spin state calculation script by Yaoyao Liu. See -h for options]'
print ''



# Not needed with optparse module
#args = sys.argv[1:]
if options.write == True or options.gaussian == True:
	spinstatesfile = "SpinStates.txt"
	spinstates = open(spinstatesfile,'w') 

if options.write == True:
	spinstatesfile2 = "SpinStates-2.txt"
	spinstates2 = open(spinstatesfile2,'w')

# Begin main loop
for filename in args:

	inputfile = readlogfile(filename)
	print '\n-=> Processing file', filename, '...'


	filenamestub = (filename.rpartition('/')[2]).rpartition('.')[0]
	pathfulenamestub = filename.rpartition('.')[0]

	# Generate SpinStates.txt and SpinStates-2.txt from .mol files if specified
	if options.gaussian == False:

		# Parse file data into OpenBabel
		obdata = pybel.readfile('mol', filename).next()

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
	
		# Get SMILES id from OB
		name = obdata.write(options.format).split()[0]
		# Get additioanl data from OB
		mass = obdata.OBMol.GetExactMass()
		formula = obdata.OBMol.GetSpacedFormula(1,'')

	#	mult2 = obdata.OBMol.GetTotalSpinMultiplicity()

		if not electrons % 2:
			mult = [1,3]
		else:
			mult = [2,4]
	
	# Generate SpinStates.txt from Gaussian03 input files if specified
	if options.gaussian == True:
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
		if options.gaussian == False:
			print  '--------- Results Summary --------'
			print  'Formula: ', formula
			print  'Name string: ', name
			print  'Molecular mass (amu): ', mass
			print  'Number of electrons: ', electrons
			print  'Spin multiplicity: ', mult
#		print  'Total atomic spin multiplicity: ', mult2
			print  '----------------------------------'

		if options.gaussian == True:
			print  '------- G03 Input Parameters ------'
			print  'Functional: ', functional
			print  'Basis set: ', basis
			print  'Calculation: ', gaucalc
			print  'Precision: ', precision
			print  'Spin multiplicity: ', mult
			print  '-----------------------------------'

	# Write spinstates files
	if options.gaussian == True:
		spinline = "%s, %i\n" %(filenamestub, mult)
		#print spinline
		spinstates.write(spinline)
		print 'Writing to spin multiplicity file', spinstatesfile, '...' 
		
	if options.write == True:
		# Write spin states file 1
		spinline = "%s, %i\n" %(filenamestub, mult[0])
	#	print spinline
		spinstates.write(spinline)
		print 'Writing to spin multiplicity file', spinstatesfile, '...' 

		spinline2 = "%s, %i\n" %(filenamestub, mult[1])
	#	print spinline2
		spinstates2.write(spinline2)
		print 'Writing to spin multiplicity file', spinstatesfile2, '...'

# Cleanup
if options.write == True:
	spinstates.close()