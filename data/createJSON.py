#! /usr/bin/python3

import sqlite3
import json
import re
import urllib.request
from os import path
from os import stat

conn = sqlite3.connect('../cirfess.db')
c = conn.cursor()

def createAccession(accession):
    sql = """SELECT seq FROM prot where Accession = ?;"""
    c.execute( sql, (accession,) )
    r = c.fetchone()
    d={}
    d['Accession'] = accession
    d['sequence'] = r[0]
    d['taxid'] = '9606'
    d['features'] = []
    return(d)

def addMotifs(d):
    motifs = ['NXS', 'NXT', 'NXC', 'NXV', 'C', 'K']
    #E7DFEE
    colors = {'NXS':'#FFA07A', 'NXT':'#4DB07E', 'NXC':'#00BFFF', 'NXV':'#ADF5D1', 'C':'#9987C8', 'K':'#5CEBA3'}
    for motif in motifs:
        sql = 'SELECT motifsLoc' + motif + ' FROM prot where Accession = ?;'
        c.execute( sql, (d['Accession'],) )
        r = c.fetchone()
        if r:
            locs = r[0].split(',')
            for loc in locs:
                m = {}
                m['type'] = motif
                m['category'] = "MOTIF"
                m['begin'] = loc
                m['end'] = loc
                m['color'] = colors[motif]
                d['features'].append(m)
    return(d)

def weighPeptide(s):
    mass_table = {
            "A" : 71.03711,
            "R" : 156.10111,
            "N" : 114.04293,
            "D" : 115.02694,
            "C" : 103.00919 + 57.02146,
            "E" : 129.04259,
            "Q" : 128.05858,
            "G" : 57.02146,
            "H" : 137.05891,
            "I" : 113.08406,
            "L" : 113.08406,
            "K" : 128.09496,
            "M" : 131.04049,
            "F" : 147.06841,
            "P" : 97.05276,
            "S" : 87.03203,
            "T" : 101.04768,
            "W" : 186.07931,
            "Y" : 163.06333,
            "U" : 168.053,
            "V" : 99.06841
            }
    mass = 0

    for aa in s:
        try:
            mass += mass_table[aa]
        except:
            return(0)
    return(mass)

def addPeptides(d):
    sql = 'SELECT OKforMS, range, numMotifsNXS, numMotifsNXT, numMotifsNXC, numMotifsNXV, numMotifsC, numMotifsK, pepSeq from pep WHERE accession = ? AND numMissedCleavages = 0;'
    c.execute( sql, (d['Accession'],) )
    peptides = c.fetchall()
    for peptide in peptides:
        ok = 'No' if peptide[0] == '0' else 'Yes'
        start, end = peptide[1].split('-')
        hasMotif = 0
        if peptide[2] > 0 or peptide[3] > 0 or peptide[4] > 0 or peptide[5] > 0:
            hasMotif = 1
            m = {}
            m['type'] = 'N-GLYCO-CAPTURE'
            m['category'] = "PEPTIDE"
            m['begin'] = start
            m['end'] = end
            m['description'] = 'OK for MS: ' + ok  + ' | |  Mass: ' + str(weighPeptide(peptide[8]))
            if peptide[0] != '0':
                m['color'] = '4DB07E'
            else:
                m['color'] = '#606060'
            d['features'].append(m)
        if peptide[6] > 0:
            hasMotif = 1
            m = {}
            m['type'] = 'C-CAPTURE'
            m['category'] = "PEPTIDE"
            m['begin'] = start
            m['end'] = end
            m['description'] = 'OK for MS: ' + ok  + ' | |  Mass: ' + str(weighPeptide(peptide[8]))
            if peptide[0] != '0':
                m['color'] = '4DB07E'
            else:
                m['color'] = '#606060'
            d['features'].append(m)
        if peptide[7] > 0:
            hasMotif = 1
            m = {}
            m['type'] = 'K-CAPTURE'
            m['category'] = "PEPTIDE"
            m['begin'] = start
            m['end'] = end
            m['description'] = 'OK for MS: ' + ok  + ' | |  Mass: ' + str(weighPeptide(peptide[8]))
            if peptide[0] != '0':
                m['color'] = '4DB07E'
            else:
                m['color'] = '#606060'
            d['features'].append(m)
        if hasMotif == 0:
            m = {}
            m['type'] = 'NOT-CAPTURED'
            m['category'] = "PEPTIDE"
            m['begin'] = start
            m['end'] = end
            m['description'] = 'OK for MS: ' + ok  + ' | |  Mass: ' + str(weighPeptide(peptide[8]))
            if peptide[0] != '0':
                m['color'] = '4DB07E'
            else:
                m['color'] = '#606060'
            d['features'].append(m)
    return(d)

def parsePhobius(pstring, length):
    i1 = list()
    i2 = list()
    o1 = list()
    o2 = list()
    topo = pstring
    if topo.startswith('i'):
        i1.append(1)
    elif topo.startswith('o'):
        o1.append(1)
    else:
        match = re.search(r'\S+/(\S+?)([io])(\S*)', topo)
        if match:
            topo = match.group(2) + match.group(3)
        else:
            print(topo)
        if topo.startswith('i'):
            i1.append(int(match.group(1)))
        elif topo.startswith('o'):
            o1.append(int(match.group(1)))
    for match in re.finditer(r'i(\d+?)-(\d+?)o', topo):
        i2.append(int(match.group(1))-1)
        o1.append(int(match.group(2))+1)
    for match in re.finditer(r'o(\d+?)-(\d+?)i', topo):
        i1.append(int(match.group(2))+1)
        o2.append(int(match.group(1))-1)
    if len(i1) > len(i2):
        i2.append(length)
    elif len(o1) > len(o2):
        o2.append(length)
    l = []
    for i in range(len(o1)):
        r = str(o1[i]) + '-' + str(o2[i]+1)
        l.append(r)
    return l

def addDomains(d):
    sql = 'SELECT StringOutPhobius, Length FROM prot WHERE Accession = ?;'
    c.execute( sql, (d['Accession'],) )
    phobius = c.fetchone()
    domains = parsePhobius(phobius[0],phobius[1])
    for domain in domains:
        start, stop = domain.split('-')
        m = {}
        m['type'] = 'PHOBIUS'
        m['category'] = "EXTRACELLULAR_DOMAIN"
        m['begin'] = start
        m['end'] = stop
        d['features'].append(m)
    sql = 'SELECT StringOutTMHMM, Length FROM prot WHERE Accession = ?;'
    c.execute( sql, (d['Accession'],) )
    tmhmm = c.fetchone()
    if tmhmm[0] != 'n/a':
        domains = parsePhobius(tmhmm[0],tmhmm[1])
        for domain in domains:
            start, stop = domain.split('-')
            m = {}
            m['type'] = 'TMHMM'
            m['category'] = "EXTRACELLULAR_DOMAIN"
            m['begin'] = start
            m['end'] = stop
            d['features'].append(m)
    return(d)

def checkForSignal(d):
    u = 'https://www.ebi.ac.uk/proteins/api/features/' + d['Accession'] + '.json'
    with urllib.request.urlopen(u) as url:
        data = json.loads(url.read().decode())
    for feature in data['features']:
        if feature['type']  == 'SIGNAL':
            feature['category'] = 'SIGNAL_PROCESSING'
            feature['color'] = '#FFA07A'
            d['features'].append(feature)
    return(d)

sql = "SELECT Accession FROM prot;"
c.execute( sql )
accs = c.fetchall()
cnt = 0
for lacc in accs:
    cnt += 1
    if cnt % 250 == 0:
        print( str(cnt) )
    acc = lacc[0]
    if path.exists('../www/data/' + acc + '.json'):
        continue
    else:
        p = createAccession(acc)
        q = addMotifs(p)
        q = addPeptides(p)
        q = addDomains(p)
    try:
        q = checkForSignal(p)
    except urllib.error.HTTPError:
        print(acc)
        continue
    f = open('../www/data/' + acc + '.json', 'w')
    f.write(json.dumps(p))
    f.close()

conn.close()
