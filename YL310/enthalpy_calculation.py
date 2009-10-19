#!/usr/bin/env python
# encoding: utf-8
"""
enthalpy_calculation.py
version 0.11

Created by Yaoyao Liu on 2009-02-22.
Copyright (c) 2009 . All rights reserved.

Version history:
0.10 - initial release. Grabs data from MySQL db for a given ID string. Calculates formation enthalpy at 298.15 K and 0 K.
0.11 - made the output slightly prettier when multiple species are specified (skips a line between species)

"""

import sys
import os
import subprocess
import fileinput
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


# User variables
thermonamedb = 'YL310thermoname_B3LYP'
thermodatadb = 'YL310thermodata_B3LYP'
frequencydb = 'YL310frequency_B3LYP'
knownentdb = 'YL310knownent'
comoentdb = 'YL310comoent_B3LYP'
projectname = 'aluminium'
reactemp = '298.15'

# reaction scheme a AlCl3  +  b TiCl4  +  c O2 --> Al_a O_x Cl_y Ti_b  +  (3a+4b-y)/2 Cl_2

# Set up species list
reacnames = [['[Al](Cl)(Cl)Cl'], ['Cl[Ti](Cl)(Cl)Cl'], ['O=O']]
prodnames = [['ClCl']]


# Parse cmd line optionsf
usage = "%prog [options] id number [id numbers]"
parser = OptionParser(usage, version="%prog 0.11 by Yaoyao Liu")




parser.add_option("-k", "--knownent",
					action="store", dest="knownent", metavar="knownent", default='',
					help="upload known species enthalpy by ID number")
parser.add_option("-i", "--ID",
					action="store_true", dest="id", default=True,
					help="calculate enthalpy by species ID number")
parser.add_option("-w", "--upload",
					action="store_true", dest="mysql", default=False,
					help="upload enthalpy data to MySQL database")
parser.add_option("-n", "--host",
					action="store", dest="mysql_host", metavar="host", default="127.0.0.1",
					help="specify alternative MySQL connection host IP address")
parser.add_option("-u", "--user",
					action="store", dest="mysql_user", metavar="user", default="prime",
					help="specify alternative MySQL connection user")
parser.add_option("-p", "--password",
					action="store", dest="mysql_passwd", metavar="passwd", default="m4sterpl4n",
					help="specify alternative MySQL connection password")
parser.add_option("-d", "--database",
					action="store", dest="mysql_db", metavar="database", default="prime",
					help="specify alternative MySQL database")
parser.add_option("-s", "--mysql-version",
					action="store_true", dest="mysql_sv", default=False,
					help="print MySQL server version and exit")
parser.add_option("-q", "--quiet",
					action="store_true", dest="quiet", default=False,
					help="do not print parsed data")

(options, args) = parser.parse_args()

# Old MySQL stuff. Redundant now.
#conn = MySQLdb.connect (host = "127.0.0.1",
#                          user = "prime",
#                          passwd = "m4sterpl4n",
#                          db = "prime");

print '[Enthalpy calculation script by Yaoyao Liu. See -h for options]'
print ''

# Print MySQL server version only
if options.mysql_sv == True:
	serverversion()
	exit()

# Not needed with optparse module
#args = sys.argv[1:]

conn = MySQLdb.connect(host = options.mysql_host,
						user = options.mysql_user,
						passwd = options.mysql_passwd,
						db = options.mysql_db);
cursor = conn.cursor()

# Fetch common thermo data for enthalpy calculations
if options.id == True:

	# Get ID and U_0 total energy data from thermonamedb for reactants
	for n in range(0, len(reacnames)):
#		print reacnames[n][0]
		sql = "SELECT id, Name, Energy from %s where Name='%s'" %(thermonamedb, reacnames[n][0])
		cursor.execute(sql)
		renergydata = cursor.fetchall()
		reacnames[n].append(int(renergydata[0][0]))   # ID
		reacnames[n].append(float(renergydata[0][2])) # Electronic + zero point energies

		# Get H(298.15) - H(0) data from thermonamedb for reactants
#		print reacnames[n][1]
		sql = "SELECT ID, name, temp, enthalpy from %s where (ID='%s' AND temp=%s)" %(thermodatadb, reacnames[n][1], reactemp)
		cursor.execute(sql)
		rentcorrdata = cursor.fetchall()
#		print rentcorrdata
		reacnames[n].append(float(rentcorrdata[0][3])) # Enthalpy correction from 0 K to 298.15 K

		# Get known enthalpy from knownentdb for reactants
		sql = "SELECT id, Name, Enthalpy from %s where Name='%s'" %(knownentdb, reacnames[n][0])
		cursor.execute(sql)
		knownentdata = cursor.fetchall()
		reacnames[n].append(float(knownentdata[0][2])) # Standard enthalpy of formation at 298.15 K
#		print reacnames


	# Get ID and U_0 total energy data from thermonamedb for products
	for n in range(0, len(prodnames)):
#		print prodnames[n][0]
		sql = "SELECT id, Name, Energy from %s where Name='%s'" %(thermonamedb, prodnames[n][0])
		cursor.execute(sql)
		penergydata = cursor.fetchall()
		prodnames[n].append(int(penergydata[0][0]))   # ID
		prodnames[n].append(float(penergydata[0][2])) # Electronic + zero point energies

		# Get H(298.15) - H(0) data from thermonamedb for products
#		print prodnames[n][1]
		sql = "SELECT ID, name, temp, enthalpy from %s where (ID='%s' AND temp=%s)" %(thermodatadb, prodnames[n][1], reactemp)
		cursor.execute(sql)
		pentcorrdata = cursor.fetchall()
#		print pentcorrdata
		prodnames[n].append(float(pentcorrdata[0][3]))  # Enthalpy correction to 298.15 K

		# Get known enthalpy from knownentdb for products
		sql = "SELECT id, Name, Enthalpy from %s where Name='%s'" %(knownentdb, prodnames[n][0])
		cursor.execute(sql)
		knownentdata = cursor.fetchall()
		prodnames[n].append(float(knownentdata[0][2])) # Standard enthalpy of formation at 298.15 K
#		print prodnames
		
#	print reacnames
#	print prodnames



# Upload known enthalpy if specified
if options.knownent != '':
	
	# Stop if no enthalpy specified even though ID is specified
	if args == []:
		print 'You forgot to specify the species ID. Stopping now.'
		exit()

	# If both specified, then proceed
	else:
		knownentid = args[0]
		print '-=> Retrieving species from database', thermonamedb,'...'
		sql = "SELECT Name, ID, canname, mass, spin from %s where ID='%s'" %(thermonamedb, knownentid)  
		cursor.execute(sql)
		speciesdata = cursor.fetchall()

		# Check that species ID exists in database
		if cursor.rowcount == 0:
			print 'WARNING: Species ID', knownentid , 'does not exist in database. Abort.'
			exit()
		else:
			(name, id, formula, mass, mult) = speciesdata[0]
	
		# Print basic output if not quietened down:
		if options.quiet == False:
			print  '--------- Species Details --------'
			print  'Formula: ', formula
			print  'Name string: ', name
			print  'Molecular mass (amu): ', mass
			print  'Spin multiplicity: ', mult
			print  'Standard ∆H_f (kJ/mol): ', options.knownent
			print  '----------------------------------'
		
		# Upload known enthalpies of formation to SQL
		print 'Inserting standard ∆H_f =', options.knownent , 'kJ/mol into database', knownentdb ,'...'

		sql = "SELECT Name, ID from %s where ID='%s'" %(knownentdb, knownentid)  
		cursor.execute(sql)
		speciesdata = cursor.fetchall()

		# Check that species ID does not already exist in knowntentdb
		if cursor.rowcount != 0:
			print 'WARNING: Species ID', knownentid , 'already exists in database ' + knownentdb + '. Abort.'
			exit()
		else:
			sql = "INSERT into " + knownentdb + " (Name, ID, Enthalpy) VALUES ('%s','%s','%s')" %(name, knownentid, options.knownent)
			cursor.execute(sql)


else:
	for i in args:
		
		# Clear variables
		islinear = 0
		mass = 0
		charge = 0
		mult = 0
		freqs = []
		rotsym = 0
		(Ba, Bb, Bc) = (0, 0, 0)
		zvpe = 0
		u_0 = 0
		eigenvals = 0
		moinertia = 0
		(Ia, Ib, Ic) = (0 ,0, 0)
		formula = ''
		name = ''
		
		# Tell the people what we are doing
		print '\n-=> Processing species ID', i,'...'
	
		# Get data from thermoname db
		# Get by project. CURRENTLY BROKEN.
		if options.id == False:
			sql = "SELECT Name, ID, Energy, Functional, path, molefile, basis, spin, calculation, Ia, Ib, Ic, mass, project, canname, comments, Ba, Bb, Bc, symmetrynumber, calclevel from %s where project='%s'" %(thermonamedb, options.project)  
			cursor.execute(sql)
			speciesdata = cursor.fetchall()

		# Get by ID
		else:
			sql = "SELECT Name, ID, Energy, Functional, path, molefile, basis, spin, calculation, Ia, Ib, Ic, mass, project, canname, comments, Ba, Bb, Bc, symmetrynumber, calclevel from %s where ID='%s'" %(thermonamedb, i)  
			cursor.execute(sql)
			speciesdata = cursor.fetchall()

		(name, id, u_0, functional, filenamestub, molfile, basis, mult, gaucalc1, Ia, Ib, Ic, mass, project, formula, comment, Ba, Bb, Bc, rotsym, precision) = speciesdata[0]

		# ------------------
		# Get data from frequency db
		# Get by project. CURRENTLY BROKEN.
		if options.id == False:
			sql = "SELECT  ID, Name, frequency from %s where project='%s'" %(frequencydb, options.project)   
			cursor.execute(sql)
			frequencydata = cursor.fetchall()

		# Get by ID
		else:
			sql = "SELECT  ID, Name, frequency from %s where ID='%s'" %(frequencydb, i)   
			cursor.execute(sql)
			frequencydata = cursor.fetchall()

		# Set up list of freqs
		freqs = []
		f = 0
		while f < len(frequencydata):
			freq = frequencydata[f][2]
			freqs.append(freq)
			f = f + 1

		# Generate better formatted frequency list
		freqlist = ""
		f = 0
		while f < len(freqs):
			freqlist = "%s  %s" %(freqlist, freqs[f])
			f = f + 1

		# ------------------
		# Printing output
		if options.quiet == False:

			print  '--------- Species Details --------'
			print  'Formula: ', formula
			print  'Name string: ', name
			print  'Molecular mass (amu): ', mass
			print  'Spin multiplicity: ', mult
			print  'Harmonic frequencies (cm^-1):', freqlist
			print  'Symmetry number:', rotsym
			print  'Rotational constants (GHz): ', Ba,' ',Bb,' ',Bc
			print  'Moments of inertia (au): ', Ia,' ',Ib,' ',Ic
			print  'Total electronic + ZPVE at 0K (kJ/mol): ', u_0
			print  '------- Gaussian Parameters ------'
			print  'Functional: ',functional
			print  'Basis set: ', basis
			print  'Calculation: ', gaucalc1#+',',gaucalc2
			print  'Precision: ', precision
			print  '----------------------------------'


		# Calculate enthalpy of formation at 0K

		# Set up unknown product's data list
		comoprod = []
		
		# Get ID and U_0 total energy data from thermonamedb for products
		sql = "SELECT Name, id, Energy from %s where ID='%s'" %(thermonamedb, i)
		cursor.execute(sql)
		puenergydata = cursor.fetchall()
		comoprod.extend([puenergydata[0][0], int(puenergydata[0][1]), float(puenergydata[0][2])])

		# Get H(298.15) - H(0) data from thermonamedb for products
		sql = "SELECT ID, name, temp, enthalpy from %s where (ID='%s' AND temp=%s)" %(thermodatadb, i, reactemp)
		cursor.execute(sql)
		puentcorrdata = cursor.fetchall()
		comoprod.append(float(puentcorrdata[0][3]))
		
		
		# Read molfile data into OpenBabel
		obdata = pybel.readstring('mol', molfile)
		
		# Get spaced formula in order to count electrons
		formulasp = obdata.OBMol.GetSpacedFormula(0,' ').split()
		
		# Setup list of atoms in specified product species
		productatoms = [['Al'], ['Cl'], ['O'], ['Ti']]
		
		# Count number of atoms in specified product species
		j = 0
		while j < len(productatoms):
			atom = productatoms[j][0]
			if formulasp.count(atom) > 0:
				k = formulasp.index(atom)
				productatoms[j].append(formulasp[k+1]) 
			else:
				productatoms[j].append('0')
			j = j + 1
			
		# reaction scheme a AlCl3  +  b TiCl4  +  c/2 O2 --> Al_a  Cl_d  O_c Ti_b  +  (3a+4b-d)/2 Cl_2
		# Set up more human-friendly stoichoimetric reaction coefficients
		(r_a, r_b, r_c, r_d) = (int(productatoms[0][1]), int(productatoms[3][1]), int(productatoms[2][1]), int(productatoms[1][1]))
		
#		print productatoms
		# Calculate ∆H_r(298.15 K) using H(products at 298.15 K) - H(reactants at 298.15 K)
		H_products = comoprod[2] + comoprod[3] + 0.5 * (3*r_a + 4* r_b - r_d) * (prodnames[0][2] + prodnames[0][3])
		H_reactants = r_a * (reacnames[0][2] + reacnames[0][3]) + r_b * (reacnames[1][2] + reacnames[1][3]) + 0.5 * r_c * (reacnames[2][2] + reacnames[2][3])
		
		# Elementary chemistry step
		dH_r = H_products - H_reactants
		print 'Calculated ∆H_r at', reactemp, 'K =', dH_r, 'kJ/mol.'
		
		# Equate absolute enthalpy and formation enthalpy calculations of ∆H_r
		dH_f = dH_r + r_a * reacnames[0][4] + r_b * reacnames[1][4] + 0.5*r_c * reacnames[2][4] - 0.5 * (3*r_a + 4* r_b - r_d) * prodnames[0][4]
		print 'Calculated ∆H_f at', reactemp, 'K =', dH_f, 'kJ/mol.'
		
		# Subtract thermal contribution to enthalpy to obtain 0 K standard formation enthalpy
		dH_f0 = dH_f - comoprod[3]
		print 'Calculated ∆H_f at 0 K =', dH_f0, 'kJ/mol.'
		print  '----------------------------------'
		
		
		# Upload data to MySQL comoentdb if specified
		if options.mysql == True:
			print 'Inserting species into database', comoentdb,'...'

			# Upload known enthalpies of formation to SQL
			print 'Inserting enthalpies into database', comoentdb ,'...'

			sql = "SELECT Name, ID from %s where ID='%s'" %(comoentdb, comoprod[1])  
			cursor.execute(sql)
			speciesdata = cursor.fetchall()

			# Check that species ID does not already exist in knowntentdb
			if cursor.rowcount != 0:
				print 'WARNING: Species ID', comoprod[1] , 'already exists in database ' + comoentdb + '. Abort.'
				exit()
			else:
				sql = "INSERT into " + comoentdb + " (Name, ID, Enthalpy, StdEnthalpy) VALUES ('%s','%s','%s','%s')" %(name, comoprod[1], dH_f0, dH_f)
				cursor.execute(sql)


# Clean up
cursor.close()
conn.close()


