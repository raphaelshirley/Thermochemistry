#!/usr/bin/env python
# encoding: utf-8
"""
mysql_data_copy.py
version 0.10

Created by Yaoyao Liu on 2009-04-30.
Copyright (c) 2009 . All rights reserved.

Version history:
0.10 - initial release. Copies data from specified fields between two tables in MySQL.
        -> Useful for copying symmetry numbers or other data between different tables/calculation methods.

"""

import sys
import os
import subprocess
import fileinput
import MySQLdb

import openbabel, pybel
from optparse import OptionParser
import private

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

projectname = 'aluminium'
reactemp = '298.15'

# Parse cmd line optionsf
usage = "%prog [options] id number [id numbers]"
parser = OptionParser(usage, version="%prog 0.10 by Yaoyao Liu")


parser.add_option("-i", "--ID",
					action="store_true", dest="id", default=True,
					help="fetch data by species ID number")
parser.add_option("--project",
					action="store", dest="project", metavar="name", default=projectname,
					help="fetch data by project name")
parser.add_option("--field",
					action="store", dest="field", metavar="name", default="",
					help="source database to copy from")
parser.add_option("--from",
					action="store", dest="sourcedb", metavar="name", default="",
					help="source database to copy from")
parser.add_option("--to",
					action="store", dest="destdb", metavar="name", default="",
					help="destination database to copy to")
parser.add_option("-n", "--host",
					action="store", dest="mysql_host", metavar="host", default=private.defaulthost,
					help="specify alternative MySQL connection host IP address")
parser.add_option("-u", "--user",
					action="store", dest="mysql_user", metavar="user", default=private.defaultuser,
					help="specify alternative MySQL connection user")
parser.add_option("-p", "--password",
					action="store", dest="mysql_passwd", metavar="passwd", default=private.defaultpasswd,
					help="specify alternative MySQL connection password")
parser.add_option("-d", "--database",
					action="store", dest="mysql_db", metavar="database", default=private.defaultdb,
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

print '[MySQL data copying script by Yaoyao Liu. See -h for options]'


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

	id = ""
	name = ""
	mult = ""
	formula = ""
	functional = ""
	basis = ""
	precision = ""
	copyfield = ""
	destfield = ""
	destid = ""
	
	# Tell the people what we are doing
	print '-=> Processing species ID', i,'...'
	
	# ------------------
	# Get data from data from sourcedb
	sql = "SELECT ID, name, spin, canname, functional, basis, calculation, calclevel, %s from %s where (ID='%s')" %(options.field, options.sourcedb , i)
	cursor.execute(sql)
	sourcedata = cursor.fetchall()
	
	(id, name, mult, formula, functional, basis, gaucalc, precision, fielddata ) = sourcedata[0]
#	print sourcedata[0]

	# ------------------
	# Printing output
	if options.quiet == False:
		print  'SOURCE database:', options.sourcedb
		print  'SOURCE data', options.field , '=', fielddata
		print  '---------- SOURCE SPECIES ---------'
		print  'Formula: ', formula
		print  'Name string: ', name
		print  'Spin multiplicity: ', mult
		print  '------ SOURCE Gaussian Params -----'
		print  'Functional: ', functional
		print  'Basis set: ', basis
		print  'Calculation: ', gaucalc
		print  'Precision: ', precision
		print  '----------------------------------'
		
	# Uplad data to sourcedb in MySQL:
	if options.destdb != "":
		
		print 'Checking for identical species in', options.destdb, '...'
		
		# Get the destination species data
		sql = "SELECT ID, name, spin, canname, functional, basis, calculation, calclevel, %s from %s where (name='%s' AND spin='%s' AND canname='%s')" %(options.field, options.destdb , name, mult, formula)
		cursor.execute(sql)
		destdata = cursor.fetchall()

		# Check that a suitable species exists in the destiantion database
		if cursor.rowcount == 0:
			print 'WARNING: Could not find matching species in ' + options.destdb + '. Nothing changed.'
			pass
		
		if cursor.rowcount > 1:
			print 'WARNING: Multiple matching species in ' + options.destdb + '. Nothing changed.'
			pass

		if cursor.rowcount == 1:
			print '... Found an entry for species ' + formula + ' in ' + options.destdb + '. Hurray!'
			# Set up data vector for destination data
			(destid, destname, destmult, destformula, destfunctional, destbasis, destgaucalc , destprecision, destfielddata ) = destdata[0]
#			print destdata
			if options.quiet == False:
#				print  'Old Destination data', options.field , '=', destfielddata
				print  '--- DESTINATION Gaussian Params ---'
				print  'Functional: ', destfunctional
				print  'Basis set: ', destbasis
				print  'Calculation: ', destgaucalc
				print  'Precision: ', destprecision
				print  '----------------------------------'
			
			# Overwrite the new value into table
			print 'Overwriting old', options.field , 'value of',  destfielddata , 'with new value' , fielddata, '...'
			sql = "UPDATE %s SET %s='%s'  where ID = '%s' " %(options.destdb, options.field , fielddata, destid)
			cursor.execute(sql)
	print  '\n'


# Clean up
cursor.close()
conn.close()
