import MySQLdb;
import os.path
import string;
from os import system 


conn = MySQLdb.connect (host = private.defaulthost,
                        user = private.defaultuser,
                        passwd = private.defaultpasswd,
                        db = private.defaultdb);
cursor = conn.cursor ()
print "Please enter ID to delete"
id=raw_input(">")
sql="select * from Thermoname where ID=%s" %id
cursor.execute(sql)
if cursor.rowcount==0:
    print "ID not found"
else:    
    row=cursor.fetchone()
    print "Do you want to delete:"
    print "Smile: %s " %row[0]
    print "Outfile: %s " %row[4]
    print "Are you sure (y/n)"
    ans=raw_input(">")
    if (ans=='y'):
        sql="delete from frequency where ID=%s" %id
        cursor.execute(sql)
        sql="delete from thermodata where ID=%s" %id
        cursor.execute(sql)
        sql="delete from Thermoname where ID=%s" %id
        cursor.execute(sql)        
        print "Entry deleted"
