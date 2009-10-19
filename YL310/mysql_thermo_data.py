#!/usr/bin/env python
# encoding: utf-8
"""
mysql_thermo_data.py
version 0.26

Created by Yaoyao Liu on 2009-01-30.
Copyright (c) 2009 . All rights reserved.

Version history:
0.10 - initial release. grabs data from MySQL db for a given ID string. Generates input file for THERMO.PL and/or .mol file.
0.20 - major new release:
		- fixed bug in THERMO.PL input file for linear molecules where rotational constants may be zero.
		- runs THERMO.PL and captures output.
		- optionally uploads output to thermodatadb in MySQL.
		- optionally generates CSV file of results.
0.21 - added total electronic + ZPVE variables.
0.22 - clears variables for each iteration 'just in case'(TM)
0.23 - discards negative frequencies from calculation
0.26 - added ZPVE. Changed units for MySQL database uploads to kJ/mol. 
       - also made the output slightly prettier when multiple species are specified (skips a line between species)


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
knownentdb = 'YL310knownent_B3LYP'
comoentdb = 'YL310comoent_B3LYP'
projectname = 'aluminium'

# Parse cmd line optionsf
usage = "%prog [options] id number [id numbers]"
parser = OptionParser(usage, version="%prog 0.26 by Yaoyao Liu")





parser.add_option("-t", "--thermo",
					action="store_true", dest="thermo", default=False,
					help="generate data file for THERMO.PL")
parser.add_option("-r", "--run-thermo",
					action="store_true", dest="thermorun", default=False,
					help="run THERMO.PL script")
parser.add_option("--csv",
					action="store_true", dest="csv", default=False,
					help="generate CSV file of thermochemistry data")
parser.add_option("-e",
					action="store_true", dest="enthalpy", default=False,
					help="calculate enthalpy of formation at 0K")
parser.add_option("-m", "--mol",
					action="store_true", dest="mol", default=False,
					help="generate .Mol file from database")
parser.add_option("-i", "--ID",
					action="store_true", dest="id", default=True,
					help="fetch data by species ID number. Note you can use bash shell expansions, e.g., 23{1,2} = 231 232")
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

print '[Data retrieval and thermochemistry script by Yaoyao Liu. Requires THERMO.PL. See -h for options]'


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
		sql = "SELECT Name, ID, Energy, Functional, path, molefile, basis, spin, calculation, Ia, Ib, Ic, mass, project, canname, comments, Ba, Bb, Bc, symmetrynumber, calclevel, zpve from %s where project='%s'" %(thermonamedb, options.project)  
		cursor.execute(sql)
		speciesdata = cursor.fetchall()

	# Get by ID
	else:
		sql = "SELECT Name, ID, Energy, Functional, path, molefile, basis, spin, calculation, Ia, Ib, Ic, mass, project, canname, comments, Ba, Bb, Bc, symmetrynumber, calclevel, zpve from %s where ID='%s'" %(thermonamedb, i)  
		cursor.execute(sql)
		speciesdata = cursor.fetchall()

	(name, id, u_0, functional, filenamestub, molfile, basis, mult, gaucalc1, Ia, Ib, Ic, mass, project, formula, comment, Ba, Bb, Bc, rotsym, precision, zpve) = speciesdata[0]

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
		
		# Check that frequencies are positive
		if freq > 0:
			freqs.append(freq)
		else:
			print 'WARNING: Discarding negative frequency', freq
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
	newfilenamestub =  formula + "_spin=%i_" %(mult) + functional + '_' + basis + '_' + precision 

	# Write THERMO.PL data file
	if options.thermo == True or options.thermorun == True or options.mysql == True or  options.csv == True:

		thermofilename = newfilenamestub+'.dat'
		thermofile = open(thermofilename, 'w')

		line1 = 'MASS\t%f\n'%(mass)

		# Only include non-zero rotational constants
		line2 = line2 = 'GHZ'
		if Ba > 0:
			line2 = line2 + '\t%f' %(Ba)
		if Bb > 0:
			line2 = line2 + '\t%f' %(Bb)
		if Bc > 0:
			line2 = line2 + '\t%f' %(Bc)
		line2 = line2 + '\n'

		line3 = 'VIB\t'
		f = 0
		while f < len(freqs):
			line3 = '%s %f' %(line3, freqs[f])
			f = f + 1
		line3 = '%s\n' %(line3)

		line4 = 'SYMNO\t%i\n' %(rotsym)
		line5 = 'ELEC\t%i\t 0.0' %(mult)
		thermodata = [line1, line2, line3, line4, line5]
	
		print 'Generating',  thermofilename, ' ...'
		t = 0
		while t < len(thermodata):
			thermofile.write(thermodata[t])
			t = t + 1
		thermofile.close()

	# Write mol file
	if options.mol == True:
	
		# Parse molfile data into OpenBabel
		obdata = pybel.readstring('mol', molfile)
		molfilename = newfilenamestub + '.mol'
		print 'Generating', molfilename, '...'
		obdata.write("mol", molfilename)


	# Execute THERMO.PL script and capture stdout
	if options.thermorun == True or options.mysql == True or options.csv == True:
		
		print "Executing 'perl thermo.pl %s' ..." %(thermofilename) 
		
		thermout = subprocess.Popen("perl thermo.pl '%s'" %(thermofilename) , shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		stdout = thermout.stdout
		
#		print stdout
#		# Stop if no thermo.pl script.
#		if thermout.stdout == "":
#			print 'ERROR: Make sure therml.pl script is in the same directory as this script. Aborting.'
#			exit()
			
		# Loop over all thermo.pl stdout
		for line in stdout:
			
			# Capture only thermo output
			if line[0:5] == 'T (K)':

				# Set up list of thermo properties
				thermolist = []
				
				# Be quiet if told to do so
				if options.quiet == False:
					print '------- Printing S, Cp, H against temperature -------'
					print line.strip()

				line = stdout.next()

				# The thermo data does not have any whitespace...
				while line.strip() != "":
					
					try:
						# Be quiet if told to do so
						if options.quiet == False:
							print line.strip()

						thermoline = [map(float, line.split())]
						thermolist.extend(thermoline)

						line = stdout.next()
						
					# Exit gracefully at the end of stdout
					except StopIteration:
						break
				print '-------------------------------------------------------'

			if line[0:11] == '**WARNING**':

				print line.strip()
				print stdout.next()

	# Insert thermodata into into thermodatadb
	if options.mysql == True:

		# Check for existing entries
		sql = "select * from " + thermodatadb + " where id='%s'" %(id)
		cursor.execute(sql)
		if cursor.rowcount > 0:
			print 'WARNING: Species ' + name + ' already exists in database. Stopping now.'
			exit()

		print 'Inserting thermochemistry data into database ' + thermodatadb ,'(all quantities in kJ/mol) ...'

		# Loop over all temperatures
		t = 0
		while t < len(thermolist):

			# G = H - T S
			gibbsenergy = (thermolist[t][3]*1000 - thermolist[t][0] * thermolist[t][1])/1000
			
			sql = "INSERT into " + thermodatadb + " (ID, Name, temp, entropy, heatcapacity, enthalpy, freeenergy) VALUES ('%i','%s','%f','%f','%f','%f','%f')" %(id, name, thermolist[t][0], thermolist[t][1]/1000, thermolist[t][2]/1000, thermolist[t][3], gibbsenergy)
			cursor.execute(sql)
			t = t + 1

	# Write CSV file of thermo data:
	if options.csv == True:
		
		print 'Generating CSV file',newfilenamestub + '.csv ...'
		csvfile = open(newfilenamestub + '.csv','w')
		csvfile.write("'Temp (K)', 'S (J/mol.K)', 'Cp (J/mol.K)', 'H (kJ/mol)', 'G (kJ/mol)'\n")
		
		# Loop over all temperatures
		t = 0
		while t < len(thermolist):

			# G = H - T S
			gibbsenergy = (thermolist[t][3]*1000 - thermolist[t][0] * thermolist[t][1])/1000

			csvline = '%f, %f, %f, %f, %f\n' %(thermolist[t][0],thermolist[t][1],thermolist[t][2],thermolist[t][3],gibbsenergy)
			csvfile.write(csvline)
			t = t + 1
		
		# Clean up
		csvfile.close


# Clean up
cursor.close()
conn.close()
