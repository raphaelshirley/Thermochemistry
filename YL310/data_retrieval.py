#!/usr/bin/env python
# encoding: utf-8
"""
data_retrieval.py
version 0.10

Created by Yaoyao Liu on 2009-01-30.
Copyright (c) 2009 . All rights reserved.

Version history:
0.10 - initial release. grabs data from MySQL db for a given ID string. Generates CSV file with data.

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
thermonamedb = 'YL310thermoname_mPWPW91'
thermodatadb = 'YL310thermodata_mPWPW91'
frequencydb = 'YL310frequency_mPWPW91'
knownentdb = 'YL310knownent'
comoentdb = 'YL310comoent_mPWPW91'
projectname = 'aluminium'
reactemp = '298.15'

# Parse cmd line optionsf
usage = "%prog [options] id number [id numbers]"
parser = OptionParser(usage, version="%prog 0.10 by Yaoyao Liu")





parser.add_option("--csv",
					action="store", dest="csv", default="",
					help="generate CSV file of thermochemistry data")
parser.add_option("-i", "--ID",
					action="store_true", dest="id", default=True,
					help="fetch data by species ID number")
parser.add_option("--project",
					action="store", dest="project", metavar="name", default=projectname,
					help="fetch data by project name")
parser.add_option("-w", "--upload",
					action="store_true", dest="mysql", default=False,
					help="upload thermo data to MySQL database")
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

print '[Data retrieval script by Yaoyao Liu. See -h for options]'
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


if options.csv != "":
	
	print 'Generating CSV file data_retrieved.csv ...'
	csvfile = open(options.csv + '.txt','w')
	csvfile.write("'ID'\t 'Species'\t Al\t Cl\t O\t Ti\t 'SMI'\t 'multiplcity'\t 'H_f (kJ/mol)'\t 'S (J/mol.K)'\t 'R.C. (GHz)'\t 'R.C. (GHz)'\t 'R.C. (GHz)'\t 'Vib freqs (cm^-1)' \n")

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
		
		
	# Get formation enthalpy at 298K from comoent database	
	sql = "SELECT Name, ID, Enthalpy from %s where ID='%s'" %(comoentdb, i)  
	cursor.execute(sql)
	comoentdata = cursor.fetchall()
	
	
	# Get H(298.15) - H(0) and S(298.15) data from thermonamedb
	sql = "SELECT ID, name, temp, enthalpy, entropy from %s where (ID='%s' AND temp=%s)" %(thermodatadb, i, reactemp)
	cursor.execute(sql)
	thermodata = cursor.fetchall()

	
	

	# Set up list of freqs
	freqs = []
	f = 0
	while f < len(frequencydata):
		freq = int(round(float(frequencydata[f][2]),0))
		
		# Check that frequencies are positive
		if freq > 0:
			freqs.append(freq)
		else:
			print 'WARNING: Discarding negative frequency', freq
		f = f + 1

	# Generate better formatted frequency list
	freqlist = ''
	f = 0
	if len(freqs)>1:
		while f < len(freqs):
			if freqlist == "":
				freqlist = "%s" %(freqs[f])
			else:
				freqlist = "%s,  %s" %(freqlist, freqs[f])
			f = f + 1
			freqlist = freqlist + ''
	else:
		freqlist = freqs[0]

	h_f = comoentdata[0][2] + thermodata[0][3]
	print h_f

	entropy = thermodata[0][4]
	print entropy

	# ------------------
	# Printing output
	if options.quiet == False:

		print  '--------- Results Summary --------'
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
		print  'Functional: ', functional
		print  'Basis set: ', basis
		print  'Calculation: ', gaucalc1#+',',gaucalc2
		print  'Precision: ', precision
		print  '----------------------------------'


	# New systematic filename
#	newfilenamestub =  formula + "_spin=%i_" %(mult) + functional + '_' + basis + '_' + precision 

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
	
	
	
	
	# Write CSV file of thermo data:
	if options.csv != "" :
		

		print 'Adding', name, 'to CSV file data_retrieved.csv ...'
		csvline = "%s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s \n" %(id, formula, int(productatoms[0][1]), int(productatoms[1][1]), int(productatoms[2][1]), int(productatoms[3][1]), name, mult, round(float(h_f),2), round(float(entropy*1000),2) , round(float(Ba),4), round(float(Bb),4), round(float(Bc),4), freqlist)
		print csvline
		csvfile.write(csvline) 
		


if options.csv == True:
	# Clean up CSV file
	csvfile.close
	
# Clean up
cursor.close()
conn.close()
