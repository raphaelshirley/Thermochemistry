#!/usr/bin/env python
# encoding: utf-8
"""
gaussian_data.py
version 0.41

Created by Yaoyao Liu on 2008-12-12.
Copyright (c) 2008 . All rights reserved.

Version history:
0.11 - initial release
0.21 - add basic MySQL upload functionality. Probably includes bugs.
0.25 - add automatic file rename.
0.27 - bug fixes to mol generation - can now generate mol files independently of SQL upload argument.
0.28 - added rotational constants and symmetry number extraction.
0.30 - added ability to generate data input file for THERMO.PL thermochemistry script.
0.31 - fixed bug parsing moments of inertia and rotational constants with linear molecules. Misc cleanups.
0.32 - fixed bug in THERMO.PL input file for linear molecules where rotational constants may be zero.
0.33 - calculates and uploads (total electronic + ZPVE) to SQL.
0.35 - new method for detecting linear molecules such as AlClO using point groups. Also clears out old variables 'just in case'(TM).
0.36 - improved method for detecting linear molecules when Gaussian is confused and shows '*******' for one rotational constant and two other identical constants. Increased precision of moment of inertia parsing to 5dp from 4dp.
0.40 - significantly more robust gaussian input file parsing.
0.41 - made the output slightly prettier when multiple species are specified (skips a line between species)

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


# User variables
thermonamedb = 'YL310thermoname_B971'
thermodatadb = 'YL310thermodata_B971'
frequencydb = 'YL310frequency_B971'
projectname = 'aluminium'


# Constants
N_a = 6.0221417930E23
E_h = 4.3597441775E-18

# Parse cmd line optionsf
usage = "%prog [options] filename [filenames]"
parser = OptionParser(usage, version="%prog 0.41 by Yaoyao Liu")



parser.add_option("-t", "--thermo",
					action="store_true", dest="thermo", default=False,
					help="Generate data file for THERMO.PL")
parser.add_option("-m", "--mol",
					action="store_true", dest="mol", default=False,
					help="Generate .Mol file from Gaussian file")
parser.add_option("-r", "--rename",
					action="store_true", dest="rename", default=False,
					help="Automatically rename files systematically")
parser.add_option("-f", "--name-format",
					action="store", dest="format", metavar="fmt", default="smi",
					help="alternative name string format (default: smi)")
parser.add_option("--project",
					action="store", dest="project", metavar="name", default=projectname,
					help="name of project field in MySQL database (default: aluminium)")
parser.add_option("--comment",
					action="store", dest="comment", metavar="text", default="",
					help="text of comment field in MySQL database")
parser.add_option("-w", "--upload",
					action="store_true", dest="mysql", default=False,
					help="upload Gaussian parameters to MySQL database")
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

print '[Gaussian parsing + MySQL script by Yaoyao Liu. See -h for options]'


# Print MySQL server version only
if options.mysql_sv == True:
	serverversion()
	exit()
	
	
	
if options.mysql == True:
	conn = MySQLdb.connect(host = options.mysql_host,
							user = options.mysql_user,
							passwd = options.mysql_passwd,
							db = options.mysql_db);
	cursor = conn.cursor()

# Not needed with optparse module
#args = sys.argv[1:]

# Begin main loop
for filename in args:

	inputfile = readlogfile(filename)
	print '\n-=> Processing file', filename, '...'

	# Define our variable for reading the file line by line
#	line = inputfile.readline()

#	# Let people know we are dealing with Gaussian files only. BROKEN.
#	for line in inputfile:
#		
#		print line.find("Entering Gaussian System")
#		if line.find("Entering Gaussian System") >= 0:
#			print "Gaussian input file detected. Proceeding..."
#			break
#		else:
#			print "Gaussian input file not detected. Stopping now."
#			exit ()
#			break
	
	# Clear out old variables
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
	
	# Begin file-reading loop
	for line in inputfile:
		
		# Molecular Mass (amu)
		if line[1:16] == 'Molecular mass:':
			mass = float(line.split()[2])
			
		# Charge and multiplicity.
		if line[1:7] == 'Charge' and line.find("Multiplicity") >= 0:
			broken = line.split()
			charge = int(broken[2])
			mult   = int(broken[-1])

		# Harmonic frequencies (cm^-1).
		if line[1:21] == "Harmonic frequencies":
			
			# Set up list of freqs
			freqs = []
		    
			# The freq data does not have any whitespace... hah. Use this to our advantage.
			while line.strip() != "":
				if line[1:15] == "Frequencies --":
						freq = [float(f) for f in line[15:].split()]
						freqs.extend(freq)
				line = inputfile.next()
		
		# Rotational symmetry number
		if line[1:27] == 'Rotational symmetry number':
			rotsym = int(float(line.split()[3]))
		
		# Check if molecule is linear using point groups
		if line[1:17] == 'Full point group':
			fullpointgroup = line[18:].strip()
			
			# Linear <==> D*H or C*V point groups
			if fullpointgroup.partition('*')[1] == '*':
				islinear = True
#				print islinear
				
		# Rotational constants (GHz)
		if line[1:23] == 'Rotational temperature':
			line = inputfile.next() 
			if line[1:20] == 'Rotational constant':
				rotfreqs = line[29:].split()
				
				# Remove useless data if linear molecule
				if islinear == True or (line[38:41]=='***' and rotfreqs[1] == rotfreqs[2]):
					
					if len(rotfreqs) == 3:
						(Ba, Bb, Bc) = map(float,(rotfreqs[2], 0, 0))

					if len(rotfreqs) == 2:
						(Ba, Bb, Bc) = map(float,(rotfreqs[1], 0, 0))

					if len(rotfreqs) == 1:
						(Ba, Bb, Bc) = map(float,(rotfreqs[0], 0, 0))
				
				# Molecule is not linear
				else:
					(Ba, Bb, Bc) = map(float,(rotfreqs[0], rotfreqs[1], rotfreqs[2]))


		# Zero-point vibrational energy (kJ/mol)
		if line[1:30] == 'Zero-point vibrational energy':
			zpve = float(line[31:].split()[0])/1000
			
		if line[1:42] == 'Sum of electronic and zero-point Energies':
			u_0 = float(line[44:].split()[0]) * E_h * N_a / 1000
			

		# Principal moments of inertia (atomic units)
		if line[1:15] == 'Principal axes':
			line = inputfile.next()
			# Blank space in this block
			while line[1:5].strip() == "":
				
				if line[5:19] == "EIGENVALUES --":
					eigenvals =  line[20:].split()
					
					# Check that the line can be spliced into 3 values
					if len(eigenvals) == 3:
						moinertia = map(float, eigenvals)
					
					# If not, the numbers probably run into each other. Splicing manually...
					else:
						moinertia  = map(float, [line[21:31], line[31:41], line[41:51]])
					(Ia, Ib, Ic) = (moinertia[0], moinertia[1], moinertia[2])

				line = inputfile.next()
	
	
	# Parse file data into OpenBabel
	filenamestub = (filename.rpartition('/')[2]).rpartition('.')[0]
	obdata = pybel.readfile('g03', filename).next()
	pathfulenamestub = filename.rpartition('.')[0]

	# Get SMILES id from OB
	name = obdata.write(options.format).split()[0]
	# Get additioanl data from OB
	energy = obdata.OBMol.GetEnergy()
	formula = obdata.OBMol.GetSpacedFormula(1,'')

###### START PARSING GAUSSIAN INPUT FILE #####
	# Open .gau input file for parameters
	gaufilename = pathfulenamestub + '.gau'
	gau = open(gaufilename, 'r')

	# Reset Gaussian input file variables for each file before parsing
	functionalbasisfound = 0
	functional = ''
	basis = ''
	gaucalc = ''
	precision = ''

	for line in gau:

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

			line = gau.next()

		# Checking for charge and multiplcity
		chargemultcheck = line.split()
		# Try to see if the line contains only two integers <-->  it's the charge and multiplicity line
		if len(chargemultcheck) == 2:
			try:
				charge = int(chargemultcheck[0])
				mult = int(chargemultcheck[1])
			except:
				pass
			
	gau.close()
###### END PARSING GAUSSIAN INPUT FILE #####

	# Generate better formatted frequency list
	freqlist = ""
	f = 0
	while f < len(freqs):
		freqlist = "%s  %s" %(freqlist, freqs[f])
		f = f + 1

	# Printing output
	if options.quiet == False:

		print  '--------- Results Summary --------'
		print  'Formula: ', formula
		print  'Name string: ', name
		print  'Molecular mass (amu): ', mass
		print  'Spin multiplicity: ', mult
		print  'Harmonic frequencies (cm^-1):', freqlist
		print  'Symmetry number: ', rotsym
		print  'Rotational constants (GHz): ', Ba,' ',Bb,' ',Bc
		print  'Moments of inertia (au): ', Ia,' ',Ib,' ',Ic
		print  'Zero-point vibrational energy (kJ/mol):', zpve
		print  'Total electronic + ZPVE at 0K (kJ/mol): ', u_0
		print  'OpenBabel energy', energy
		print  '------- Gaussian Parameters ------'
		print  'Functional: ', functional
		print  'Basis set: ', basis
		print  'Calculation: ', gaucalc
		print  'Precision: ', precision
		print  '----------------------------------'

#	# Save output using pickle
#	pickle_data = [mass,mult,freqs,moinertia]
#	pickle_name = os.path.splitext(fileinput.filename())[0]+'.pk1'
#	pickle_file = open(pickle_name,'wb')
#	
#	print 'Saved to', pickle_name
#	pickle.dump(pickle_data,pickle_file)
#	pickle_file.close

	# New systematic filename
	newfilenamestub =  formula + "_spin=%i_" %(mult) + functional + '_' + basis + '_' + precision 
	
	
	# Write mol files
	if options.mysql == True or options.mol == True:
		# Write new automatic filename
		if options.rename == True:
			molfilename = newfilenamestub + '.mol'
		else:
			molfilename = filenamestub + '.mol'
	
		print 'Generating', molfilename, '...'
		obdata.write("mol", molfilename)
		molfile = open(molfilename, 'r')
		
	# Rename Gaussian output file
	if options.rename == True and options.mysql == False:
		newfilename = newfilenamestub + '.g03'
		print 'Renaming',filename,'to',newfilename,'...'
		os.rename(filename,newfilename)
		print 'Renaming', gaufilename, 'to', newfilenamestub + '.gau ...'
		os.rename(gaufilename, newfilenamestub + '.gau')
		filenamestub = newfilenamestub
	
	# Write THERMO.PL output file
	if options.thermo == True:
		thermofile = open(newfilenamestub+'.dat', 'w')
		
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
		
		print 'Generating ' + newfilenamestub+'.dat ...'
		t = 0
		while t < len(thermodata):
			thermofile.write(thermodata[t])
			t = t + 1
		thermofile.close()
		
		
	# Upload data to MySQL database	
	if options.rename == False:
		if options.mysql == True:

			# Check for duplicates in database
			sql = "select * from " + thermonamedb + " where  path='%s'" %(filenamestub)
			cursor.execute(sql)
			if cursor.rowcount>0:
				print 'WARNING: Species ' + filenamestub + ' already exists in database. Stopping now.'
				exit()
		
			# Insert into thermonamedb
			print 'Inserting' , filenamestub , 'to database', thermonamedb, '...'
			sql = "INSERT into " + thermonamedb + " ( Name, Energy, Functional, path, molefile, basis, spin, calculation, Ia, Ib, Ic, mass, project, canname, comments, Ba, Bb, Bc, symmetrynumber, calclevel, zpve ) VALUES ('%s', '%s','%s','%s','%s','%s','%i','%s','%f','%f','%f','%f','%s','%s','%s','%f','%f','%f','%i','%s','%f')" %(name, u_0, functional, filenamestub, molfile.read(), basis, mult, gaucalc, Ia, Ib, Ic, mass, options.project, formula, options.comment, Ba, Bb, Bc, rotsym, precision, zpve)
			cursor.execute(sql)
		
			# Get the new ID number
			cursor.execute("SELECT LAST_INSERT_ID()")
			id = cursor.fetchone()[0]
			print 'Auto-assigned unique id', id
		
			# Insert frequencies into frequencydb
			print 'Inserting frequencies into database ' + frequencydb,'...'
			f = 0
			while f < len(freqs):
				sql = "INSERT into " + frequencydb + " (ID, Name, frequency) VALUES ('%f','%s','%f')" %(id, name, freqs[f])
	#			print freqs[f]
				cursor.execute(sql)
				f = f + 1


#		sql="INSERT into " + thermodatadb + " (id, temp, entropy, heatcapacity, enthalpy, freeenergy) VALUES ('%s','298.15','0','0','0','0')" %(id)
#		cursor.execute(sql)

# Clean up
if options.mysql == True:
	cursor.close()
	conn.close()
	molfile.close()
			