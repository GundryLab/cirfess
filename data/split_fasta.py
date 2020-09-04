#! /usr/bin/python3

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-s', '--start', dest='start', default=1, type=int)
parser.add_argument('-e', '--end', dest='end', type=int)
parser.add_argument('-f', '--infile', dest='infile')
args = parser.parse_args()

first = args.start
last = args.end
infile = args.infile

fasta_fileobj = open(infile, 'r')   ## create a file obj from the specified file

cnt = 0
seq = ''
seqs = []

for line in fasta_fileobj:
    if line.startswith('>'):
        if cnt > 0:
            seqs.append(seq)
        cnt += 1
        seq = line
    else: 
        seq += line
seqs.append(seq)

# if they request more sequences than there are, just send what you have.
if last > len(seqs):
    last = len(seqs)
for i in range(first -1, last):
    print(seqs[i], end='')

    
            

