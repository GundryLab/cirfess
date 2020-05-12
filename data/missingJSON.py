#! /usr/bin/python3

import sqlite3
import glob

conn = sqlite3.connect('../cirfess.db')
cursor = conn.cursor()

files = glob.glob('/home/jack/work/CirfessBranches/cirfess/www/data/*.json')
print('First file: ' + files[0])

sql = 'SELECT Accession FROM prot;'
cursor.execute(sql)

r = cursor.fetchall()

for acc in r:
  if '/home/jack/work/CirfessBranches/cirfess/www/data/' + acc[0] + '.json' not in files:
    print(acc[0])

conn.close()
