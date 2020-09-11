#! /usr/bin/python3

import argparse

description = '''
This program is meant for updating tmhmm files used 
as proteomic annotation.  Updating from one proteome 
version to another can be difficult. The idea with 
this program is to have one primary tmhmm results
file that gets updated from version to version. Each 
time a new proteome is released, this program can be 
run to update the primary version using the latest 
version. After running update_fasta.py use the fasta
diff file with tmhmm and use the output with this 
program.'''

parser = argparse.ArgumentParser(description=description)
parser.add_argument('primaryfile',  help='The primary tmhmm file that will be updated')
parser.add_argument('updatefile', help='Upudates to the tmhmm file')
args = parser.parse_args()

primaryfile = args.primaryfile
updatefile = args.updatefile

def read_seqs (infile) :
    fo = open(infile, 'r')   ## create a file obj from the specified file

    seqs = {}

    for line in fo:
        ID = line.split('|')[1]
        seqs[ID] = line

    fo.close()
    return seqs

primary = read_seqs(primaryfile)
update = read_seqs(updatefile)

additions = []

for accession in update.keys() :
    primary[accession] = update[accession]
    additions.append(accession)

fo = open(primaryfile, 'w')
for accession in primary.keys() :
    fo.write(primary[accession])
fo.close()

print('Accessions added/update to primary tmhmm file:')
for accession in additions:
    print(accession)
