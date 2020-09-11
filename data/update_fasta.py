#! /usr/bin/python3

import argparse

description = '''
This program is meant for updating fasta files used 
as input to sequence analysis programs that are then
used as proteomic annotation. Updating from one 
proteome version to another can be difficult. The 
idea with this program is to have one primary fasta 
file that gets updated from version to version. Each 
time a new proteome is released, this program can be 
run to update the primary version using the latest 
version. In doing so, a diff file will be created. 
This diff file can be submitted to sequence analysis 
programs such as phobius instead of submitting the 
entire latest realease.'''

parser = argparse.ArgumentParser(description=description)
parser.add_argument('primaryfile',  help='The primary fasta file that will be updated')
parser.add_argument('updatefile', help='The latest version of the fasta file')
parser.add_argument('difffile', help='File created with the sequences different between the primary and update fasta files')
args = parser.parse_args()

primaryfile = args.primaryfile
updatefile = args.updatefile
difffile = args.difffile


def read_seqs (infile) :
    fasta_fileobj = open(infile, 'r')   ## create a file obj from the specified file

    cnt = 0
    seq = ''
    seqs = {}

    for line in fasta_fileobj:
        if line.startswith('>'):
            if cnt > 0:
                seqs[ID] = {}
                seqs[ID]['seq'] = seq
                seqs[ID]['version'] = version
            ID = line.split('|')[1]
            version = line.split()[-1]
            cnt += 1
            seq = line
        else: 
            seq += line
        seqs[ID] = {}
        seqs[ID]['seq'] = seq
        seqs[ID]['version'] = version

    fasta_fileobj.close()
    return seqs

primary = read_seqs(primaryfile)
update = read_seqs(updatefile)

additions = []
updates = []

for accession in update.keys() :
    if accession not in primary.keys() :
        primary[accession] = update[accession]
        additions.append(accession)
    elif update[accession]['version'] != primary[accession]['version'] :
        primary[accession] = update[accession]
        updates.append(accession)

fasta_fileobj = open(difffile, 'w')
print('New Accessions added to fasta file:')
for accession in additions:
    print(accession)
    fasta_fileobj.write(update[accession]['seq'])
print('Accessions that were updated:')
for accession in updates:
    print(accession)
    fasta_fileobj.write(update[accession]['seq'])
fasta_fileobj.close()

fasta_fileobj = open(primaryfile, 'w')  
for accession in primary.keys() :
    fasta_fileobj.write(primary[accession]['seq'])
fasta_fileobj.close()

