#! /usr/bin/python3

import sqlite3
import re
import createCirfess

# extract data from TMHMM file
# get file from
def parse_TMHMM(file_name):
    with open(file_name, 'r') as fo:
        tm_dict = dict()
        for line in fo:
            line = line.rstrip()
            IDstring = line.split('\t')[0]
            ID = IDstring.split('|')[1]
            length = int(line.split('\t')[1].split('=')[1])
            topo_str = line.split('\t')[5].split('=')[1]
            num_TM = int(line.split('\t')[4].split('=')[1])

            if len(topo_str) == 1:
                tm_dict[ID] = { 'length' : length, 'num_TM' : num_TM , 'topo' : 'n/a' }
            else:
                tm_dict[ID] = { 'length' : length, 'num_TM' : num_TM , 'topo' : topo_str }


        for key in tm_dict:
            inside = list()
            outside = list()
            i1 = list()
            i2 = list()
            o1 = list()
            o2 = list()

            topo = tm_dict[key]['topo']

            if topo:

                if topo.startswith('i'):
                   i1.append(1)
                elif topo.startswith('o'):
                   o1.append(1)

                for match in re.finditer(r'i(\d+?)-(\d+?)o', topo):
                    i2.append(int(match.group(1))-1)
                    o1.append(int(match.group(2))+1)
                for match in re.finditer(r'o(\d+?)-(\d+?)i', topo):
                    i1.append(int(match.group(2))+1)
                    o2.append(int(match.group(1))-1)

                if len(i1) > len(i2):
                    i2.append(tm_dict[key]['length'])
                elif len(o1) > len(o2):
                    o2.append(tm_dict[key]['length'])

            inside_residues = list()
            outside_residues = list()

            for i in range(len(i1)):
                inside_residues.extend(list(range(i1[i],i2[i]+1)))

            for i in range(len(o1)):
                outside_residues.extend(list(range(o1[i],o2[i]+1)))

            tm_dict[key]['inside'] = inside_residues
            tm_dict[key]['outside'] = outside_residues


    return tm_dict

def parse_Phobius(file_name):
    import re
    with open(file_name, 'r') as fo:

        tm_dict = dict()

        with open('sources/length.tab' , 'r') as fo2:
            length_dict = dict()
#            header = True
            for line in fo2:
                line = line.rstrip()
                ID = line.split()[0]
                length = line.split('\t')[1]
                length_dict[ID] = int(length)

        for line in fo:
            line = line.rstrip()
            ID, num_TM, SP, topo_str = line.split()
            ID = ID.split('|')[1]

            if bool(length_dict.get(ID)):
                length = length_dict[ID]
            else:
                continue

            if len(topo_str) < 1:
                tm_dict[ID] = { 'length' : length, 'topo' : 'n/a' , \
                                'num_TM' : int(num_TM), 'SP' : SP}
            else:
                tm_dict[ID] = { 'length' : length, 'topo' : topo_str ,  \
                                'num_TM' : int(num_TM), 'SP' : SP}

        for key in tm_dict:
            inside = list()
            outside = list()
            i1 = list()
            i2 = list()
            o1 = list()
            o2 = list()

            topo = tm_dict[key]['topo']

            if topo:

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
                    i2.append(tm_dict[key]['length'])
                elif len(o1) > len(o2):
                    o2.append(tm_dict[key]['length'])

            inside_residues = list()
            outside_residues = list()

            for i in range(len(i1)):
                inside_residues.extend(list(range(i1[i],i2[i]+1)))

            for i in range(len(o1)):
                outside_residues.extend(list(range(o1[i],o2[i]+1)))

            tm_dict[key]['inside'] = inside_residues
            tm_dict[key]['outside'] = outside_residues

    return tm_dict


def parse_signalP(file_name):

        with open(file_name, 'r') as fo:

            SP_dict = dict()

            for line in fo:
                line = line.rstrip()
                IDstring = line.split()[0]
                ID = IDstring.split('_')[1]
                SP = line.split()[1]
                score = line.split()[2]

                if SP.startswith('SP'):
                    SP = 'Y'
                else:
                    SP = 0

                SP_dict[ID] = { 'SP': SP , 'score' : score }

        return SP_dict

def parse_Predisi(file_name):

    with open (file_name, 'r') as fo:

        predisi_dict = dict()

        for line in fo:
            line = line.rstrip()
            IDstring = line.split()[0]
            ID = IDstring.split('|')[1]
            SP = line.split('\t')[-2]
            score  = line.split('\t')[-4]

            predisi_dict[ID] = {'SP' : SP, 'score' : score}


    return predisi_dict

def parse_SPC(seq_dict, file_name):

    with open(file_name, 'r') as fo:

            SPC_dict = dict()
            header = True

            for line in fo:
                    if header:
                        header = False
                    else:
                        line = line.rstrip()
                        ID, SPC, BF, T, dC, DR = line.split(',')

                        SPC_dict[ID] = dict()
                        SPC_dict[ID]['score'] = int(SPC)

                        stringOut = ''

                        if int(BF) == 1 :
                            stringOut += 'BF'
                        if int(T) == 1 :
                            stringOut += ',T'
                        if int(dC) == 1:
                            stringOut += ',dC'
                        if int(DR) == 1:
                            stringOut += ',DR'

                        stringOut = stringOut.lstrip(',')

                        SPC_dict[ID]['stringOut'] = stringOut

            for ID in seq_dict:

                if SPC_dict.get(ID) == None:
                    SPC_dict[ID] = {'score': 0, 'stringOut' : 'n/a'}

    return SPC_dict



def fasta_parser(fasta_filename):


    fasta_fileobj = open(fasta_filename, 'r') ## create a file obj from the specified file
    lfile = open('sources/length.tab', 'w')
    sequence_name = '' ## initialize strings to populate from file object info
    sequence_desc = ''
    sequence_string = ''
    sequence_dict = {}

    for line in fasta_fileobj:   ## iterate through file object with for loop
        line = line.rstrip() ## strip white space on the right side (like a new line character!)

        if line.startswith('>'):

            if len(sequence_string) > 0:
                sequence_dict[sequence_name] = sequence_string
                lfile.write(sequence_name + '\t' + str(len(sequence_string)) + '\n')
                sequence_string = ''   ## reset for the new sequence

            line = line.lstrip('>')   ## remove leading `>` char
            sequence_info = line.split(maxsplit=1)  ## split on only first space
            sequence_name = sequence_info[0].split('|')[1]

            if len(sequence_info) > 1:
                sequence_desc = sequence_info[1]
            else: ## sequence has no description, set to empty
                sequence_desc = ''


            line = line.lstrip('>')   ## remove leading `>` char
            sequence_info = line.split(maxsplit=1)   ## split on only first space

            if len(sequence_info) > 1:
                sequence_desc = sequence_info[1]

            else:
            # sequence has no description, set to empty
                sequence_desc = ''

        else:
            sequence_string += line  # incrementally elongate seq

# When we reach the end of the FASTA file, we drop out of the
# 'for' loop. However, we still have the last sequence record
# stored in memory, which we haven't processed yet, because we
# haven't observed a '>' symbol, so we must copy and paste any
# code that we used to process sequences above to the code block
# below. Check if sequence_string has a non-zero length to
# determine whether to execute the sequence processing code:

    if len(sequence_string) > 0:
        sequence_dict[sequence_name] = sequence_string
        lfile.write(sequence_name + '\t' + str(len(sequence_string)) + '\n')

    lfile.close()
    return sequence_dict


def trypsinize(prot_seq):
    peptides= []
    cut_sites=[0]
    indices = []
    pep = ''

    for i in range(0,len(prot_seq)-1):
        if prot_seq[i] == 'K' and prot_seq[i+1] != 'P':
            cut_sites.append(i+1)
        elif prot_seq[i] == 'R' and prot_seq[i+1] != 'P':
            cut_sites.append(i+1)

    if cut_sites[-1]!=len(prot_seq):
            cut_sites.append(len(prot_seq))

    if len(cut_sites)>2:

        for j in range(0,len(cut_sites)-3):

            pep = prot_seq[cut_sites[j]:cut_sites[j+1]]
            for i in range(cut_sites[j],cut_sites[j+1]):
                indices.append(i+1)
            peptides.append({'seq': pep,'indices': indices, 'missed_cleavages' : 0})
            indices = []

            pep = prot_seq[cut_sites[j]:cut_sites[j+2]]
            for i in range(cut_sites[j],cut_sites[j+2]):
                indices.append(i+1)
            peptides.append({'seq': pep,'indices': indices, 'missed_cleavages' : 1})
            indices = []

            pep = prot_seq[cut_sites[j]:cut_sites[j+3]]
            for i in range(cut_sites[j],cut_sites[j+3]):
                indices.append(i+1)
            peptides.append({'seq': pep,'indices': indices, 'missed_cleavages' : 2})
            indices = []

        pep = prot_seq[cut_sites[-3]:cut_sites[-2]]
        for i in range(cut_sites[-3],cut_sites[-2]):
            indices.append(i+1)
        peptides.append({'seq': pep,'indices': indices, 'missed_cleavages' : 0})
        indices = []

        pep = prot_seq[cut_sites[-3]:cut_sites[-1]]
        for i in range(cut_sites[-3],cut_sites[-1]):
            indices.append(i+1)
        peptides.append({'seq': pep,'indices': indices, 'missed_cleavages' : 1})
        indices = []

        pep = prot_seq[cut_sites[-2]:cut_sites[-1]]
        for i in range(cut_sites[-2],cut_sites[-1]):
            indices.append(i+1)
        peptides.append({'seq': pep,'indices': indices, 'missed_cleavages' : 0})
        indices = []

    else: #there is no trypsin site in the protein sequence
        peptides.append({'seq' : prot_seq, 'indices' : range(1,len(prot_seq)+1), 'missed_cleavages' : 0})

    return peptides

def ion_mim(prot_seq, charge_state):

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
            "V" : 99.06841
            }

    mass = 0

    for aa in mass_table:
        mass += prot_seq.count(aa) * mass_table[aa]


    ion_mass = mass + (charge_state * 1.007276)
    m_z = ion_mass/charge_state

    return m_z


def ok_for_MS(pep_list):

    for pep in pep_list:

        pep['okForMS'] = ''

        if len(pep['seq']) > 5 and (ion_mim(pep['seq'], 2) < 2000):
            pep['okForMS'] += '2'

        if len(pep['seq']) > 5 and (ion_mim(pep['seq'], 3) < 2000):
            pep['okForMS'] += ',3'

        if len(pep['okForMS']) > 0:
            pep['okForMS'] = pep['okForMS'].lstrip(',')
        else:
            pep['okForMS'] = None

    return pep_list

def make_prot_dict(ID1, seq_dict, tm_dict, phob_dict, sp_dict, predisi_dict, SPC_dict):
#def make_prot_dict():
    prot_dict = dict()

    l = [ID1]
    #for ID in seq_dict:
    for ID in l:

      seq = seq_dict[ID]

      if tm_dict.get(ID) and phob_dict.get(ID) and sp_dict.get(ID) and predisi_dict.get(ID):


        pep_list = trypsinize(seq)

        pep_list = ok_for_MS(pep_list)


        glyco_indices_S = list()

        for match in re.finditer(r'N[^P]S', seq ):
            glyco_indices_S.append(match.start()+1)

        glyco_indices_T = list()

        for match in re.finditer(r'N[^P]T', seq ):
            glyco_indices_T.append(match.start()+1)

        glyco_indices_C = list()

        for match in re.finditer(r'N[^P]C', seq ):
            glyco_indices_C.append(match.start()+1)

        glyco_indices_V = list()

        for match in re.finditer(r'N[^P]V', seq ):
            glyco_indices_V.append(match.start()+1)

        K_indices = list()
        C_indices = list()

        for i in range(0, len(seq)):
            if seq[i] == 'K':
                K_indices.append(i+1)
            elif seq[i] == 'C':
                C_indices.append(i+1)




        prot_dict[ID] = {
                            'seq_info' :
                              {
                                  'seq' : seq ,
                                  'seq_len' : len(seq)
                              } ,

                            'topo' :
                              {
                                 'TMHMM' :
                                    { 'inside' : tm_dict[ID]['inside']  ,
                                      'outside' : tm_dict[ID]['outside']  ,
                                      'num_TM' : tm_dict[ID]['num_TM'] ,
                                      'stringOut' : tm_dict[ID]['topo']
                                    } ,

                                 'Phobius' :
                                    { 'inside' : phob_dict[ID]['inside']  ,
                                      'outside' : phob_dict[ID]['outside']  ,
                                      'num_TM' : phob_dict[ID]['num_TM']  ,
                                      'stringOut' : phob_dict[ID]['topo']
                                    }
                              } ,

                            'signal' :
                              {
                                  'Phobius' : {'SP' : phob_dict[ID]['SP'] , 'score' : 0 } ,
                                  'SignalP' : {'SP' : sp_dict[ID]['SP'] , 'score' : sp_dict[ID]['score'] } ,
                                  'PrediSi' : {'SP' : predisi_dict[ID]['SP'] , 'score' : predisi_dict[ID]['score'] }
                              } ,

                            'SPC' : {'score' : int(SPC_dict[ID]['score']) , 'stringOut' : SPC_dict[ID]['stringOut']},

                            'peptides' : pep_list,


                            'motif_sites' :
                               {
                                 'NXS' :
                                   { 'all' : glyco_indices_S , 'extracellular' : dict()  } ,

                                 'NXT' :
                                   { 'all' : glyco_indices_T , 'extracellular' : dict()  } ,

                                 'NXC' :
                                   { 'all' : glyco_indices_C , 'extracellular' : dict()  } ,

                                 'NXV' :
                                   { 'all' : glyco_indices_V , 'extracellular' : dict()  } ,

                                 'C' :
                                   { 'all' : C_indices , 'extracellular' : dict()  } ,

                                 'K' :
                                   { 'all' : K_indices , 'extracellular' : dict()  } ,

                               }

                        }

    return prot_dict

def EC_analysis(prot):

    motif_list = ['NXS','NXT','NXC','NXV','C','K']

    for pep in prot['peptides']:

        for motif in motif_list:

            pep[motif] = dict()
            pep[motif]['extracellular'] = dict()

        for pred in ['TMHMM','Phobius']:

            outside_indices = prot['topo'][pred]['outside']

            for motif in motif_list:

                motif_indices = prot['motif_sites'][motif]['all']
                out_motif = list( set(motif_indices) & set(outside_indices) )

                prot['motif_sites'][motif]['extracellular'][pred] = out_motif

                all = set(motif_indices) & set(pep['indices'])

                if all:
                    pep[motif]['all'] = {'num' : len(all), 'indices' : list(all) }

                    out = set(out_motif) & set(pep['indices'])

                    if out:
                        pep[motif]['extracellular'][pred] = {'num' : len(out), 'indices' : list(out)}
                    else:
                        pep[motif]['extracellular'][pred] = {'num' : '0', 'indices' : 'n/a' }

                else:
                    pep[motif]['all'] = {'num' : '0', 'indices' : 'n/a' }
                    pep[motif]['extracellular'][pred] = {'num' : '0', 'indices' : 'n/a' }


    return prot

def get_pepEntry():
  p = {
    'ID': '',
    'Accession': '',
    'pepSeq': '',
    'range': '',
    'numMissedCleavages': '',
    'OKforMS': '',
    'numMotifsNXS': '',
    'numMotifsNXT': '',
    'numMotifsNXC': '',
    'numMotifsNXV': '',
    'numMotifsC': '',
    'numMotifsK': '',
    'numMotifsPhobiusNXS': '',
    'numMotifsPhobiusNXT': '',
    'numMotifsPhobiusNXC': '',
    'numMotifsPhobiusNXV': '',
    'numMotifsPhobiusC': '',
    'numMotifsPhobiusK': '',
    'numMotifsTMHMMNXS': '',
    'numMotifsTMHMMNXT': '',
    'numMotifsTMHMMNXC': '',
    'numMotifsTMHMMNXV': '',
    'numMotifsTMHMMC': '',
    'numMotifsTMHMMK':    '',
    'motifLocNXS': '',
    'motifLocNXT': '',
    'motifLocNXC': '',
    'motifLocNXV': '',
    'motifLocC': '',
    'motifLocK': '',
    'motifLocPhobiusNXS': '',
    'motifLocPhobiusNXT': '',
    'motifLocPhobiusNXC': '',
    'motifLocPhobiusNXV': '',
    'motifLocPhobiusC': '',
    'motifLocPhobiusK': '',
    'motifLocTMHMMNXS': '',
    'motifLocTMHMMNXT': '',
    'motifLocTMHMMNXC': '',
    'motifLocTMHMMNXV': '',
    'motifLocTMHMMC': '',
    'motifLocTMHMMK': ''}
  return p

def get_protEntry():
    p = {
    'Accession': '',
    'Seq': '',
    'Length': '',
    'StringOutPhobius': '',
    'numTMPhobius': '',
    'numICPhobius': '',
    'numECPhobius': '',
    'StringOutTMHMM': '',
    'numTMTMHMM': '',
    'numICTMHMM': '',
    'numECTMHMM': '',
    'SigPepPhobius': '',
    'ScorePhobius': '',
    'SigPepSignalP': '',
    'ScoreSignalP': '',
    'SigPepPrediSi': '',
    'ScorePrediSi': '',
    'numSPpredictions': '',
    'SPC': '',
    'SPCstring': '',
    'numMotifsNXS': '',
    'motifsLocNXS': '',
    'numMotifsNXT': '',
    'motifsLocNXT': '',
    'numMotifsNXC': '',
    'motifsLocNXC': '',
    'numMotifsNXV': '',
    'motifsLocNXV': '',
    'numMotifsC': '',
    'motifsLocC': '',
    'numMotifsK': '',
    'motifsLocK': '',
    'numMotifsPhobiusNXS': '',
    'motifsLocPhobiusNXS': '',
    'numMotifsPhobiusNXT': '',
    'motifsLocPhobiusNXT': '',
    'numMotifsPhobiusNXC': '',
    'motifsLocPhobiusNXC': '',
    'numMotifsPhobiusNXV': '',
    'motifsLocPhobiusNXV': '',
    'numMotifsPhobiusC': '',
    'motifsLocPhobiusC': '',
    'numMotifsPhobiusK': '',
    'motifsLocPhobiusK': '',
    'numMotifsTMHMMNXS': '',
    'motifsLocTMHMMNXS': '',
    'numMotifsTMHMMNXT': '',
    'motifsLocTMHMMNXT': '',
    'numMotifsTMHMMNXC': '',
    'motifsLocTMHMMNXC': '',
    'numMotifsTMHMMNXV': '',
    'motifsLocTMHMMNXV': '',
    'numMotifsTMHMMC': '',
    'motifsLocTMHMMC': '',
    'numMotifsTMHMMK': '',
    'motifsLocTMHMMK': '',
    'numPepMC0': '',
    'numPepMC1': '',
    'numPepMC2': '',
    'numPepOKMC0': '',
    'numPepOKMC1': '',
    'numPepOKMC2': '',
    'numPepMC0NXS': '',
    'numPepMC0NXT': '',
    'numPepMC0NXC': '',
    'numPepMC0NXV': '',
    'numPepMC0C': '',
    'numPepMC0K': '',
    'numPepMC1NXS': '',
    'numPepMC1NXT': '',
    'numPepMC1NXC': '',
    'numPepMC1NXV': '',
    'numPepMC1C': '',
    'numPepMC1K': '',
    'numPepMC2NXS': '',
    'numPepMC2NXT': '',
    'numPepMC2NXC': '',
    'numPepMC2NXV': '',
    'numPepMC2C': '',
    'numPepMC2K': '',
    'numPepMC0PhobiusNXS': '',
    'numPepMC0PhobiusNXT': '',
    'numPepMC0PhobiusNXC': '',
    'numPepMC0PhobiusNXV': '',
    'numPepMC0PhobiusC': '',
    'numPepMC0PhobiusK': '',
    'numPepMC1PhobiusNXS': '',
    'numPepMC1PhobiusNXT': '',
    'numPepMC1PhobiusNXC': '',
    'numPepMC1PhobiusNXV': '',
    'numPepMC1PhobiusC': '',
    'numPepMC1PhobiusK': '',
    'numPepMC2PhobiusNXS': '',
    'numPepMC2PhobiusNXT': '',
    'numPepMC2PhobiusNXC': '',
    'numPepMC2PhobiusNXV': '',
    'numPepMC2PhobiusC': '',
    'numPepMC2PhobiusK': '',
    'numPepMC0TMHMMNXS': '',
    'numPepMC0TMHMMNXT': '',
    'numPepMC0TMHMMNXC': '',
    'numPepMC0TMHMMNXV': '',
    'numPepMC0TMHMMC': '',
    'numPepMC0TMHMMK': '',
    'numPepMC1TMHMMNXS': '',
    'numPepMC1TMHMMNXT': '',
    'numPepMC1TMHMMNXC': '',
    'numPepMC1TMHMMNXV': '',
    'numPepMC1TMHMMC': '',
    'numPepMC1TMHMMK': '',
    'numPepMC2TMHMMNXS': '',
    'numPepMC2TMHMMNXT': '',
    'numPepMC2TMHMMNXC': '',
    'numPepMC2TMHMMNXV': '',
    'numPepMC2TMHMMC': '',
    'numPepMC2TMHMMK': ''}
    return p





#file_name = 'test.fasta'
file_name = 'proteome.fasta'
path = 'sources/'

seqs_dict = fasta_parser(path + file_name)
print('sequnces:',len(seqs_dict))

sp = parse_signalP(path + 'signalp.txt')
print('SignalP:',len(sp))

TMHMM_dict  = parse_TMHMM(path + 'tmhmm.txt')
print('TMHMM:',len(TMHMM_dict))

Phobius_dict = parse_Phobius(path + 'phobius.txt')
print('Phobius:',len(Phobius_dict))

SPC_dict = parse_SPC(seqs_dict, path + 'SPC_by_Source.csv')
print('SPC:',len(SPC_dict))

predisi_dict = parse_Predisi(path + 'predisi.txt')
print('Predisi:',len(predisi_dict))


###############################################################################
#
# Begin Generating the database
#
###############################################################################

print(' Creating database')

createCirfess.create_cirfess_structure()
conn = sqlite3.connect('../cirfess.db')
cursor = conn.cursor()


###############################################################################
# prot table
###############################################################################
print(' Creating prot table')

counter = 0

motif_list = ['NXS','NXT','NXC','NXV','C','K']

for ID in seqs_dict:
    if ID not in seqs_dict or ID not in TMHMM_dict or ID not in Phobius_dict or ID not in sp or ID not in predisi_dict or ID not in SPC_dict:
        continue
#    print(ID)
    prot_dict = make_prot_dict(ID, seqs_dict, TMHMM_dict, Phobius_dict, sp, predisi_dict, SPC_dict)

    counter += 1
    if counter%1000 == 1 :
        print ( str(counter) )

    protdb = get_protEntry()
    out = list()
    prot = prot_dict[ID]
    protdb['Accession'] = ID

    si = prot['seq_info']

    protdb['Seq'] = si['seq']
    protdb['Length'] = si['seq_len']

    for meth in ['Phobius', 'TMHMM']:
        topo = prot['topo'][meth]

        protdb['StringOut'+meth] = topo['stringOut']
        protdb['numTM'+meth] = str(topo['num_TM'])
        protdb['numIC'+meth] = str(len(topo['inside']))
        protdb['numEC'+meth] = str(len(topo['outside']))
    sppc = 0
    for meth in ['Phobius', 'SignalP', 'PrediSi']:
        sig = prot['signal'][meth]

        if sig['SP'] == 'Y':
            protdb['SigPep'+meth] = 'Y'
            sppc += 1
        else:
            protdb['SigPep'+meth] = 'N'

        if meth == 'Phobius':
            protdb['ScorePhobius'] = 'NA'
        else:
            protdb['Score'+meth] = str(sig['score'])

    protdb['numSPpredictions'] = sppc
    spc = prot['SPC']
    protdb['SPC'] = spc['score']
    protdb['SPCstring'] = spc['stringOut']

    prot = EC_analysis(prot)

    ms = prot['motif_sites']
    for loc in ['all', 'Phobius', 'TMHMM']:
        for motif in motif_list:
            if loc == 'all':
                sites = ms[motif][loc]

            else:
                sites = ms[motif]['extracellular'][loc]

            sites.sort()
            ns = len(sites)
            if loc == 'all' :
                protdb['numMotifs'+motif] = ns
                protdb['motifsLoc'+motif] = ','.join(map(str, sites)) if ns else 'NA'
            else:
                protdb['numMotifs'+loc+motif] = ns
                protdb['motifsLoc'+loc+motif] = ','.join(map(str, sites)) if ns else 'NA'

    peps = prot['peptides']

    pep0 = 0
    pep1 = 0
    pep2 = 0

    ok0 = 0
    ok1 = 0
    ok2 = 0

    counts = {'all' : dict(), 'Phobius': dict(), 'TMHMM' :dict()}
    for key in counts:
        counts[key] = [dict(), dict(), dict()]

        for i in counts[key]:
            for motif in motif_list:
                i[motif] = 0

    for pep in peps:
        mc = pep['missed_cleavages']
        ok = pep['okForMS']

        if mc == 0:
            pep0 += 1
            pep1 += 1
            pep2 += 1
            if ok:
                ok0 += 1
                ok1 += 1
                ok2 += 1

        elif mc == 1:
            pep1 += 1
            pep2 += 1
            if ok:
                ok1 += 1
                ok2 += 1

        else:
            pep2 += 1
            if ok:
                ok2 += 1

        if ok:
            for motif in motif_list:
                for loc in ['all', 'Phobius', 'TMHMM']:
                    if loc == 'all':
                        if int(pep[motif][loc]['num']) > 0:
                            counts[loc][mc][motif] += 1
                    else:
                        if int(pep[motif]['extracellular'][loc]['num']) > 0:
                            counts[loc][mc][motif] += 1

    protdb['numPepMC0'] = pep0
    protdb['numPepMC1'] = pep1
    protdb['numPepMC2'] = pep2
    protdb['numPepOKMC0'] = ok0
    protdb['numPepOKMC1'] = ok1
    protdb['numPepOKMC2'] = ok2

    for loc in ['all', 'Phobius', 'TMHMM']:

        for i in [0,1,2]:

            for motif in motif_list:

                if i == 0:
                    count = counts[loc][i][motif]

                elif i == 1:
                    count = counts[loc][i][motif] + counts[loc][0][motif]

                else:
                    count = counts[loc][i][motif] + counts[loc][1][motif] + counts[loc][0][motif]


                if loc == 'all' :
                    protdb['numPepMC'+str(i)+motif] = str(count)
                else:
                    protdb['numPepMC'+str(i)+loc+motif] = str(count)

    fields = ', '.join(protdb.keys())
    qmarks = ', '.join( list( '?'*len(protdb.keys()) ) )
    values = tuple( protdb.values() )
    insertSQL = 'INSERT INTO prot (' + fields + ') VALUES (' + qmarks + ')'
    try:
        cursor.execute(insertSQL, values)
    except sqlite3.Error as e:
        print(str(e.args[0]))
conn.commit()


###############################################################################
# pep table
###############################################################################
print(' Creating pep table')

counter = 0

topo_list = ['all', 'Phobius', 'TMHMM']
motif_list = ['NXS','NXT','NXC','NXV','C','K']

for ID in seqs_dict:
    if ID not in seqs_dict or ID not in TMHMM_dict or ID not in Phobius_dict or ID not in sp or ID not in predisi_dict or ID not in SPC_dict:
        continue
    prot_dict = make_prot_dict(ID, seqs_dict, TMHMM_dict, Phobius_dict, sp, predisi_dict, SPC_dict)

    counter += 1
    if counter%1000 == 1 :
        print ( str(counter) )

    prot_dict[ID] = EC_analysis(prot_dict[ID])

    peps = prot_dict[ID]['peptides']

    for pep in peps:
        pepdb = get_pepEntry()
        pepdb['Accession'] = ID
        pep_range = str(pep['indices'][0]) + '-' + str(pep['indices'][-1])
        uniq = ID + '_' + pep_range

        pepdb['ID']=uniq
        pepdb['pepSeq'] = pep['seq']
        pepdb['range'] = pep_range
#        pepdb['start'] = pep['indices'][0]
#        pepdb['stop'] = pep['indices'][-1]
        pepdb['numMissedCleavages'] = pep['missed_cleavages']

        if pep['okForMS']:
            pepdb['OKforMS'] = pep['okForMS']
        else:
            pepdb['OKforMS'] = 0

        for topo in topo_list:

            for motif in motif_list:
                if pep[motif]:
                    if pep[motif].get(topo): # if topo is 'all'
                        pepdb['numMotifs'+motif] = str(pep[motif][topo]['num'])
                    elif pep[motif]['extracellular'].get(topo): # else if topo is Phobius or TMMHM
                        pepdb['numMotifs'+topo+motif] = str(pep[motif]['extracellular'][topo]['num'])
                else:  # The motif was not found in the sequences provided
                    out.append('0')
                    if(topo == 'all'):
                        pepdb['numMotifs'+motif] = 0
                    else:
                        pepdb['numMotifs'+topo+motif] = 0

            for motif in motif_list:
                if pep[motif]:
                    if pep[motif].get(topo): # we are in 'all' land here, see above
                        indices = pep[motif][topo]['indices']
                        if indices == 'n/a':
                            pepdb['motifLoc'+motif] = 'NA'
                        else:
                            indices.sort()
                            pepdb['motifLoc'+motif] = ','.join(map(str, indices))
                    elif pep[motif]['extracellular'].get(topo): # we are in phobius and tmhmm land here
                        indices = pep[motif]['extracellular'][topo]['indices']
                        if indices == 'n/a':
                            pepdb['motifLoc'+topo+motif] = 'NA'
                        else:
                            indices.sort()
                            pepdb['motifLoc'+topo+motif] = ','.join(map(str, indices))
                else:
                    print('Got here 6')
                    if(topo == 'all'):
                        pepdb['motifLoc'+motif] = 'NA'
                    else:
                        pepdb['motifLoc'+topo+motif] = 'NA'


        fields = ', '.join(pepdb.keys())
        qmarks = ', '.join( list( '?'*len(pepdb.keys()) ) )
        values = tuple( pepdb.values() )
        insertSQL = 'INSERT INTO pep (' + fields + ') VALUES (' + qmarks + ')'
        try:
            cursor.execute(insertSQL, values)
        except sqlite3.Error as e:
            print(str(e.args[0]))

conn.commit()
conn.close()
